#include "moveageclonal.h"

//#define DEBUG

using namespace std;
namespace weakarg
{

MoveAgeClonal::MoveAgeClonal(Param * p, double a,bool allowclonal,int r,double t)
        : Move(p,a)
{
    description= "Update of the age of a clonal node";
    name= "AgeClonal";
    numtopo=0;
    numtopoaccept=0;
    reps=r;
    tempering=t;
    s1=param->getClonalStepSize();
    s0=param->getClonalStepSize();
    usenormalprop=true;
    backuptree = NULL;
    changeTopology=allowclonal;// we can forbid topology changes with this (no tempering is then needed)
}

inline double fround(double n, double d)
{
    return floor(n * pow(10., d) + .5) / pow(10., d);
}

int MoveAgeClonal::move()
{
    // Propose to change a nodes age
    // reorder the nodes correctly
    // change the recedge times
    // account for factors of 2 in passing other clonal nodes
    // accept/reject step
    // if reject, put things back correctly
    bool usetemper=false;// an indicator as to whether tempering is needed
    if(!changeTopology)s1=0;// probabilities of the return move won't work unless we do this (because 2 reflecting boundaries)


    double temperllorig=0,temperllto=0,temperllbase=0,temperllfinal=0;
    double temperpriororig=0,temperpriorto=0,temperpriorbase=0,temperpriorfinal=0;

    vector<double> store(param->getData()->getL());
    for (int i=0;i<param->getData()->getL();i++)
        store[i]=param->getLLsite(i);
    double ll=param->getLL(), lprior=param->getRecTree()->prior(param);

    if(backuptree == NULL)
    	backuptree = new RecTree(*(param->getRecTree()));
    else
    	backuptree->assign(*(param->getRecTree()));
    RecTree* rectree = backuptree;

    if(rectree==0) throw("Cannot create RecTree in MoveAgeClonal!");
    RecTree * oldrectree=param->getRecTree();
    param->setRecTree(rectree,true);

#if defined DEBUG
	for(int site=0;site<param->getData()->getL();site++) {
	    if(store[site]!=param->getLLsite(site)) 	cout<<"LL site "<<site<<" recomputed: "<<param->getLLsite(site)<<" difference "<<store[site]-param->getLLsite(site)<<endl;
	}
	if(ll!=param->getLL())
	cout<<"LL before "<<ll<<" recomputed: "<<param->getLL()<<" difference "<<ll-param->getLL()<<endl;
#endif

    double oldTT=rectree->getTTotal();
    int which=rectree->getN() + gsl_rng_uniform_int(rng,rectree->getN()-1),whichto;

    // reflecting boundaries at daughter and parent nodes to preserve reversability easily
    double disttodaughter=max(rectree->getNode(which)->getLeft()->getAge(),rectree->getNode(which)->getRight()->getAge())-rectree->getNode(which)->getAge();
    double maxdown= -rectree->getNode(which)->getAge(),maxup=10000.0;
    if(!changeTopology) {
	maxdown= disttodaughter;
	if (which<rectree->getN()*2-2) maxup=rectree->getDist(which);
    }
    double dist,logqratfort;
    if(usenormalprop) {// this is the default now.  Control sigma using s0 and s1 (set in the constructor)
    	dist = gsl_ran_gaussian(rng,s0+s1*rectree->getNode(which)->getAge());// normal distance with variable sigma
	logqratfort=logQrat(rectree->getNode(which)->getAge(),rectree->getNode(which)->getAge()+dist);
    }else {
    	dist = param->getClonalStepSize() * (gsl_rng_uniform(rng) * 2.0 -1.0);// uniform distance
	logqratfort=0;// do this for the uniform proposal
    }
    while(dist< maxdown || dist>maxup)
    {
        if(dist< maxdown)
            dist = 2.0 * maxdown - dist;
        else if(dist > maxup)
            dist = 2.0 * maxup - dist;
    }

// figure out if we have to change topology
    if(dist< disttodaughter || (dist>rectree->getDist(which) && which!=rectree->getN()*2-2)) {
	numtopo++;
	usetemper=true;
    }else {
	numcalls++;
	usetemper=false;
    }
    dlog(1)<<"Proposing move of node "<< which <<" from age "<<rectree->getNode(which)->getAge()<<" to age "<<rectree->getNode(which)->getAge()+dist<<"...";

    // if tempering, set up the moves
    Metropolis *movetemp;
    vector<int> samplespace;
    vector<int> listLR;
    if(usetemper) {
    	temperllbase=param->getLL()/tempering;
	temperpriorbase=rectree->prior(param);
	rectree->chooseNodePath(which, dist, &listLR, &samplespace);
        movetemp= new Metropolis(param);
	movetemp->move(reps,tempering,&samplespace);
    	temperllorig=param->getLL()/tempering;
	temperpriororig=rectree->prior(param);
    }
// Do the move, which returns the samplespace (i.e. affected clonal edges)
    int propfactor=rectree->moveNodeTime(which,&whichto,dist,-1,&listLR,&samplespace);// propfactor:number of recedges moved to younger nodes - number moved to older nodes
    param->computeLikelihood();
    rectree->computeTTotal();
    double logqratio = propfactor*log(2.0)+rectree->numRecEdge()*log(rectree->getTTotal()/oldTT)+logqratfort;
    if(usetemper) {// and complete the tempering
        temperllto=param->getLL()/tempering;
	temperpriorto=rectree->prior(param);
	movetemp->move(reps,tempering,&samplespace);
        if(movetemp!=NULL) delete(movetemp);
	else throw("Error in MoveAgeClonal: movetemp was null!");
    	temperllfinal=param->getLL()/tempering;
	temperpriorfinal=rectree->prior(param);
    }
    double llend=param->getLL();
#if defined DEBUG
    //test its ok so far
    try
    {
        rectree->testNodeAges();
    }
    catch(char * x)
    {
        cout<<x<<endl<<"Moveageclonal: proposed new node label was "<<whichto<<endl;
        exit(1);
    }
    try
    {
        param->testTree();
    }
    catch(char * x)
    {
        cout<<x<<endl<<"Moveageclonal: broke the log liks"<<endl;
        exit(1);
    }
#endif
    double newlprior=rectree->prior(param);
    // accept/reject step
//cout<<"origll="<<ll<<" temperllorig="<<temperllorig<<" temperllto="<<temperllto<<" llend="<<llend<<" temperpriorto="<<temperpriorto<<" temperpriororig="<<temperpriororig<<" lprior="<<lprior<<" newlprior="<<newlprior<<endl;
    if(log(gsl_rng_uniform(rng))>llend-ll+newlprior-lprior + temperllto + temperpriorto - temperllorig - temperpriororig + temperllbase + temperpriorbase - temperllfinal - temperpriorfinal + logqratio)
    {
    	dlog(1)<<" Rejected!"<<endl;
	// put back the rectree as it was
	param->setRecTree(oldrectree,true);
	rectree=oldrectree;
	// put back the loglikelihoods as they were
       for (int i=0;i<param->getData()->getL();i++)
            param->setlocLL(i,store[i]);
        param->setLL(ll);
        // everything should now be back how it was
#if defined DEBUG
        //test its still ok
	int tmp=0;
        try
        {
            rectree->testNodeAges();
        }
        catch(char * x)
        {
            cout<<x<<endl<<"Moveageclonal restore: proposed new node label was "<<whichto<<endl;
            exit(1);
        }
        try
        {
            param->testTree();
        }
        catch(char * x)
        {
            cout<<x<<endl<<"Moveageclonal restore: broke the log liks"<<endl;
            exit(1);
        }
        if(tmp==1)
        {
            cerr<<"Problem replacing after move"<<endl;
            throw("Move not reversed correctly");
        }
#endif
	return(0);
    }
    else {
    	backuptree=oldrectree;	// becomes our new backup tree and we avoid delete
	dlog(1)<<"Accepted!"<<endl;// accept the modified tree
	if(usetemper) numtopoaccept++;
	else numaccept++;
	return(1);
    }
}

double MoveAgeClonal::logQrat(double tfrom,double tto)
{
//	return(log((exp(logdens(tto,tfrom))+exp(logdens(tto,-tfrom)))/(exp(logdens(tfrom,tto))+exp(logdens(tfrom,-tto))))); // less efficient
	return(log((dens(tto,tfrom)+dens(tto,-tfrom))/(dens(tfrom,tto)+dens(tfrom,-tto))));

}


MoveAgeClonal::~MoveAgeClonal()
{}



} // end namespace weakarg
