#include "moveregraftclonal.h"

using namespace std;
namespace weakarg
{

//#define DEBUG

MoveRegraftClonal::MoveRegraftClonal(Param * p,double a,int r,double t)
        : Move(p,a)
{
    description="Update to regraft a clonal node to a new parent";
    name="regraft";
    reps=r;
    tempering=t;
    backuptree=NULL;
}

inline double fround(double n, double d)
{
    return floor(n * pow(10., d) + .5) / pow(10., d);
}


int MoveRegraftClonal::move()
{
    numcalls++;
	dlog(1)<<"Proposing to regraft tree...";

    if(backuptree == NULL)
    	backuptree = new RecTree(*(param->getRecTree()));
    else
    	backuptree->assign(*(param->getRecTree()));
    RecTree* rectree = backuptree;

    RecTree * oldrectree=param->getRecTree();
    param->setRecTree(rectree,true);

    int which=2*rectree->getN()-2;
    if(which<=4){ dlog(1)<<" No regrafting of 3 or fewer nodes!"<<endl;return(-1);}
    while(!validNode(which)) {which=gsl_rng_uniform_int(rng,2*rectree->getN()-2);}
    int origfather=rectree->getNode(which)->getFather()->getId();
    int notmovednode=rectree->otherChild(origfather,which);
    int whichto=getAlive(which);
	if(whichto<0){cerr<<"Error in moveregraftclonal: can't find new place to regraft!"<<endl;throw("Regraft error");};

    vector<double> store(param->getData()->getL());
    for (int i=0;i<param->getData()->getL();i++)
        store[i]=param->getLLsite(i);
    double ll=param->getLL(), lprior=rectree->prior(param);

    dlog(1)<<"Proposing move of node "<< which <<" from parent "<<origfather<<" to parent "<<whichto<<"...";

    Metropolis *movetemp= new Metropolis(param);
    vector <int> samplespace=rectree->getAllChildren(which);

    double temperllbase=param->getLL()/tempering;
    double temperpriorbase=rectree->prior(param);
    movetemp->move(reps,tempering,&samplespace);
    double temperllorig=param->getLL()/tempering;
    double temperpriororig=rectree->prior(param);

        regraft(which,whichto);
    param->computeLikelihood();

    rectree->computeTTotal();
    double temperllto=param->getLL()/tempering;
    double temperpriorto=rectree->prior(param);
    movetemp->move(reps,tempering,&samplespace);
    double llend=param->getLL();
    if(movetemp!=NULL) delete(movetemp);
    else throw("Error in MoveRegraftClonal: movetemp was null!");
    double temperllfinal=param->getLL()/tempering;
    double temperpriorfinal=rectree->prior(param);


#if defined DEBUG
    //test its ok so far
    try
    {
        param->testTree();
    }
    catch(char * x)
    {
        dlog(1)<<x<<endl<<"Moveregraftclonal: broke the log liks"<<endl;
        exit(1);
    }
#endif
    double newlprior=rectree->prior(param);
    // accept/reject step
   if(log(gsl_rng_uniform(rng))>llend-ll+newlprior-lprior + temperllto + temperpriorto - temperllorig - temperpriororig + temperllbase + temperpriorbase - temperllfinal - temperpriorfinal)
    {
    	dlog(1)<<" Rejected!"<<endl;
	param->setRecTree(oldrectree,true);
	rectree=oldrectree;

        for (int i=0;i<param->getData()->getL();i++)
            param->setlocLL(i,store[i]);
        param->setLL(ll);
        // everything should now be back how it was
#if defined DEBUG
        //test its still ok
	int tmp=0;
        for (int i=0;i<param->getData()->getL();i++)
            if(fround(param->getLLsite(i),5)!=fround(store[i],5))
            {
                tmp=1;
                dlog(1)<<"Site "<<i<<" has ll before "<<store[i]<<" and after "<<param->getLLsite(i)<<" ";
                for (int j=0;j<rectree->numRecEdge();j++)
                {
                    if (rectree->getRecEdge(j)->affectsSite(i))
                        dlog(1)<<j<<" ";
                }
                dlog(1)<<endl;
            }
        if(fround(ll,5)!=fround(param->getLL(),5))
        {
            dlog(1)<<"Total ll before "<<ll<<" and after "<<param->getLL()<<endl;
        }
        try
        {
            param->testTree();
        }
        catch(char * x)
        {
            dlog(1)<<x<<endl<<"Moveageclonal restore: broke the log liks"<<endl;
            exit(1);
        }
#endif
	return(0);
    }
    else {
    	backuptree=oldrectree;	// becomes our new backup tree and we avoid delete
	dlog(1)<<"Accepted!"<<endl;// accept the modified tree
	numaccept++;
	return(1);
    }
}

MoveRegraftClonal::~MoveRegraftClonal()
{}

int MoveRegraftClonal::getAlive(int which)
{
	RecTree * rectree=param->getRecTree();
	vector<int> allalive;
	allalive=rectree->alive(rectree->getNode(which)->getFather()->getAge());
	if(allalive.size()<=1) return(-1);
	int rnode=which;
	while(rnode==rectree->getNode(which)->getFather()->getLeft()->getId() || rnode==rectree->getNode(which)->getFather()->getRight()->getId() || rnode==rectree->getNode(which)->getFather()->getId()) rnode=allalive[gsl_rng_uniform_int(rng,allalive.size())];
	return(rnode);
}

void MoveRegraftClonal::regraft(int which,int whichto)
{
	RecTree * rectree=param->getRecTree();
	Node *fanode=rectree->getNode(which)->getFather(), *wnode=rectree->getNode(which), *wtonode=rectree->getNode(whichto);
	int oldfather=fanode->getId();
	int notmovednode=rectree->otherChild(oldfather,which);
	rectree->moveEdges(oldfather,notmovednode);
	rectree->moveEdges(whichto,oldfather,fanode->getAge(),-1.0);
	rectree->swapFather(oldfather,notmovednode);
	rectree->swapFather(oldfather,whichto);
}

} // end namespace weakarg
