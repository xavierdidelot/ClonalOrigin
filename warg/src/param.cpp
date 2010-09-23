#include "param.h"
#include "move.h"
#include "weakarg.h"
#include <fstream>
#include "movetheta.h"
#include "moverho.h"
#include "movedelta.h"
#include "moveremedge.h"
#include "moveaddedge.h"
#include "movesitechange.h"
#include "movetimechange.h"
#include "moveedgechange.h"
#include "moveageclonal.h"
#include "movescaletree.h"
#include "moveregraftclonal.h"
#include "movegreedytree.h"


//#define DEBUG

using namespace std;
namespace weakarg
{

inline double fround(double n, double d)
{
    return floor(n * pow(10., d) + .5) / pow(10., d);
}

Param::Param()
{
    hRho=0.0;
    tempering=1.0;
    rectree=NULL;
    data=NULL;
    theta=1.0;
    rho=1.0;
    delta=500.0;
    ll=0.0;
    moveClonalStepSize=0.02;
    moveScaleTreeSize=0.1;
    esttheta=estthetasq=estrho=estrhosq=estdelta=estdeltasq=estnumrecedge=estnumrecedgesq=estedgeden=estedgedensq=estedgepb=estedgepbsq=estvaredgepb=estvaredgepbsq=0;
    locll=vector<double>();
    computeLikelihood();
}

Param::Param(RecTree * rectree,Data * data)
{
    hRho=0.0;
    tempering=1.0;
    this->rectree=rectree;
    this->data=data;
    theta=1.0;
    rho=1.0;
    delta=500.0;
    ll=0.0;
    moveClonalStepSize=0.1;
    moveScaleTreeSize=0.1;
    esttheta=estthetasq=estrho=estrhosq=estdelta=estdeltasq=estnumrecedge=estnumrecedgesq=estedgeden=estedgedensq=estedgepb=estedgepbsq=estvaredgepb=estvaredgepbsq=0;
    if (data!=NULL)
    {
        locll=vector<double>(data->getL(),0);
        computeLikelihood();
    }
}

Param::~Param()
{
    for (unsigned int i=0;i<mall.size();i++) if(mall.at(i)!=NULL) delete(mall.at(i));
}

void Param::readParamsFromFile(WargXml *infile,std::streampos itstart)
{
	string p1=infile->getParam("theta",itstart);
	string p2=infile->getParam("rho",itstart);
	string p3=infile->getParam("delta",itstart);
	theta=atof(p1.c_str());
	rho=atof(p2.c_str());
	delta=atof(p3.c_str());
	dlog(1)<<"Read theta="<<theta<<" rho="<<rho<<" delta="<<delta<<endl;
}

void Param::readProgramOptions()
{
    if (opt().theta<0.0) setTheta(2.0*data->numPoly()/rectree->getTTotal()); else setTheta(opt().theta);
    if (opt().rho>=0.0) setRho(opt().rho);
    if (opt().delta>=0.0) setDelta(opt().delta);
    if (opt().thetaPerSite) theta*=data->getL();
    if (opt().rhoPerSite) rho*=(delta*data->getB()+data->getL()-data->getB());
    dlog(1)<<"Initiating moves..."<<endl;
    // These should be ordered by movep for efficiency really.
    if (opt().rho==0.0) opt().temperreps=0;
    if (opt().rho<0.0) mall.push_back(new MoveRho(this,opt().movep[0]));
    if (opt().delta<0.0) mall.push_back(new MoveDelta(this,opt().movep[1]));
    if (opt().theta<0.0) mall.push_back(new MoveTheta(this,opt().movep[2]));
    mall.push_back(new MoveRemEdge(this,opt().movep[3]));
    mall.push_back(new MoveAddEdge(this,opt().movep[4]));
    mall.push_back(new MoveSiteChange(this,opt().movep[5]));
    mall.push_back(new MoveTimeChange(this,opt().movep[6]));
    mall.push_back(new MoveEdgeChange(this,opt().movep[7]));
    mall.push_back(new MoveAgeClonal(this,opt().movep[8],opt().allowclonal,opt().temperreps,opt().temperT));
    mall.push_back(new MoveScaleTree(this,opt().movep[9]));
    if(opt().allowclonal) mall.push_back(new MoveRegraftClonal(this,opt().movep[10],opt().temperreps,opt().temperT));
    if(opt().greedyWeight>0) mall.push_back(new MoveGreedyTree(this,opt().greedyWeight));
}

void Param::simulateData(vector<int> blocks)
{
    if (rectree==NULL)
        return;
    if (data!=NULL)
        delete(data);
    int n=rectree->getN();
    int L=blocks.back();
    data=new Data(n,blocks);
    double thetaPerSite=theta/L;
    for (int site=0;site<L;site++)
    {
        rectree->makeLocalTree(site,thetaPerSite);
        vector<char> state(2*n-1,'Z');
        int y,z; // children of this site
        double ey,ez;
        state[2*n-2]=(char)floor(gsl_rng_uniform(rng)*4);
        for (int i=2*n-2;i>=n;i--)
        {
            y=rectree->getSon(i,site,true ,&ey);
            z=rectree->getSon(i,site,false,&ez);
            ey =0.25*(1.0+3.0*ey);
            ez =0.25*(1.0+3.0*ez);
            state[y]=state[i];
            if (gsl_rng_uniform(rng)>ey)
                while (state[y]==state[i])
                    state[y]=(char)floor(gsl_rng_uniform(rng)*4);
            state[z]=state[i];
            if (gsl_rng_uniform(rng)>ez)
                while (state[z]==state[i])
                    state[z]=(char)floor(gsl_rng_uniform(rng)*4);
        }
        for (int i=0;i<n;i++)
            data->set(i,site,state[i]);
    }
    locll=vector<double>(data->getL(),0);
    computeLikelihood();
}

void Param::metropolis(string comment)
{
    double totweight=0;
    for(unsigned int i=0;i<mall.size();i++)
        totweight+=mall[i]->getAlpha();
    if(totweight<=0)
    {
        cerr<<"Error: no move weightings are positive!"<<endl;
        throw;
    }

    mpiofstream output(opt().outfile);
    ostream* csv;
    if(opt().csvoutfile.length()>0)
    	csv = new ofstream(opt().csvoutfile.c_str());
    else
    	csv = new nullstream();

    exportXMLbegin(output(true),comment);
    output.flush();
    startDiagnostics(*csv);

    computeLikelihood();
    long int iterations = opt().burnin+opt().additional;
    for (long int i=0;i<iterations;i++)
    {
    	dlog(1)<<"Iteration "<<i<<endl;

		if (iterations>50 && (i)%((iterations)/50)==0)
		{
			if (100l*i/(iterations)<10)
				cout<<"\b\b\b\b#  "<<100l*(double)i/(iterations)<<"%"<<flush;
			else
				cout<<"\b\b\b\b# "<<100l*(double)i/(iterations)<<"%"<<flush;
		}
		if (i+1==iterations)
			cout<<"\b\b\b\b# 100%"<<endl<<flush;
        for (unsigned int j=0;j<mall.size();j++)
        {
            int j2=chooseMove(&mall,totweight);
            dlog(1)<<"("<<i<<":"<<mall[j2]->getName()<<") "<<flush;
#if defined DEBUG
            dlog(1)<<endl;
#endif
            try{mall[j2]->move();
	    }catch(char * x){ cout<<"Error in Param: "<<x<<endl;}
#if defined DEBUG
            testTree();
#endif
        }
//	if (iterations>50 && (i)%((iterations)/50)==0) greedyTreeMove();
        if(i>=opt().burnin)
        {
            updateDiagnostics(i-opt().burnin,*csv);
            if((i-opt().burnin)%(opt().thinin)==0)
            {
                exportXMLiter(output(false),i);
                output.flush();
            }
        }
    }
    if((iterations-opt().burnin)%(opt().thinin)==0) updateDiagnostics(iterations-opt().burnin,*csv);
    exportXMLiter(output(false),iterations+1);
    exportXMLend(output(true));
    output.flush();
    if(csv!=NULL) delete csv;	// forces a file close, if it was ever open
}

int Param::chooseMove(vector<Move*> *mall,double totweight)
{
    double r=gsl_rng_uniform(rng);
    for(unsigned int i=0;i<mall->size();i++)
    {
        r-=(*mall)[i]->getAlpha()/totweight;
        if(r<0)
            return(i);
    }
    cout<<"Error in choosing move probabilities"<<endl;
    throw;
}

void Param::exportXMLbegin(ostream& out,string comment)
{
    out << "<?xml version = '1.0' encoding = 'UTF-8'?>" << endl;
    out << "<outputFile>" << endl;
    out << "<Blocks>" << endl;
    vector<int> * blocks=data->getBlocks();
    for (unsigned int i=0;i<blocks->size()-1;i++)
        out << blocks->at(i)<<",";
    out << blocks->at(blocks->size()-1)<<endl<<"</Blocks>"<<endl;
    out<< "<comment>"<<comment<<"</comment>"<<endl;
    out<< "<nameMap>";
    data->printNames(out);
    out<<"</nameMap>"<<endl;
    out<<"<regions>";
    vector<int> * regions=data->getRegions();
    for (unsigned int i=0;i<regions->size()-1;i++)
        out << regions->at(i)<<",";
    out << regions->at(regions->size()-1)<<"</regions>"<<endl;
}

void Param::exportXMLiter(ostream& out,long i)
{
    out << "<Iteration>" << endl;
    out << "<Tree>" << endl;
    out << rectree->newick() << endl;
    out << "</Tree>" << endl;
    out << "<number>"<<i<<"</number>"<<endl;
    out << "<ll>"<<ll<<"</ll>"<<endl;
    out << "<prior>"<<rectree->prior(this)<<"</prior>"<<endl;
    out << "<theta>"<<theta<<"</theta>"<<endl;
    out << "<rho>"<<rho<<"</rho>"<<endl;
    out << "<delta>"<<delta<<"</delta>"<<endl;
    out << "<tmrca>"<<rectree->getNode(rectree->getN()*2-2)->getAge()<<"</tmrca>"<<endl;

    out << "<esttheta>"<<esttheta<<"</esttheta>"<<endl;
    out << "<estvartheta>"<<estthetasq-esttheta*esttheta<<"</estvartheta>"<<endl;
    out << "<estrho>"<<estrho<<"</estrho>"<<endl;
    out << "<estvarrho>"<<estrhosq-estrho*estrho<<"</estvarrho>"<<endl;
    out << "<estdelta>"<<estdelta<<"</estdelta>"<<endl;
    out << "<estvardelta>"<<estdeltasq-estdelta*estdelta<<"</estvardelta>"<<endl;
    out << "<estnumrecedge>"<<estnumrecedge<<"</estnumrecedge>"<<endl;
    out << "<estvarnumrecedge>"<<estnumrecedgesq-estnumrecedge*estnumrecedge<<"</estvarnumrecedge>"<<endl;
    out << "<estedgeden>"<<estedgeden<<"</estedgeden>"<<endl;
    out << "<estvaredgeden>"<<estedgedensq-estedgeden*estedgeden<<"</estvaredgeden>"<<endl;
    out << "<estedgepb>"<<estedgepb<<"</estedgepb>"<<endl;
    out << "<estvaredgepb>"<<estedgepbsq-estedgepb*estedgepb<<"</estvaredgepb>"<<endl;
    out << "<estedgevarpb>"<<estvaredgepb<<"</estedgevarpb>"<<endl;
    out << "<estvaredgevarpb>"<<estvaredgepbsq-estvaredgepb*estvaredgepb<<"</estvaredgevarpb>"<<endl;
    int clonal=-1;
    for (unsigned int i=0;i<mall.size();i++)   
    {
        if(mall[i]->getName()=="AgeClonal") clonal=i;
  	out << "<acc"<<mall[i]->getName()<<">"<<((double)mall[i]->getAcc())/max(1,mall[i]->getCounts())<<"</acc"<<mall[i]->getName()<<">"<<endl;
    }
    if(clonal>=0) 	
    {
        out << "<accAgeClonalTopo>"<<((double)mall[clonal]->getTopoAcc())/max(1,mall[clonal]->getTopoCounts())<<"</accAgeClonalTopo>"<<endl;
        out << "<propAgeClonalTopo>"<<((double)mall[clonal]->getTopoCounts())/(mall[clonal]->getTopoCounts()+mall[clonal]->getCounts())<<"</propAgeClonalTopo>"<<endl;
    }

    int rprec=cout.precision();

    for (int i=0;i<rectree->numRecEdge();i++)
    {
        RecEdge*r=rectree->getRecEdge(i);
        out <<"<recedge>";
        out <<"<start>"<<r->getStart   ()<<"</start><end>"<<r->getEnd   ()<<"</end>";
        out <<"<efrom>"<<r->getEdgeFrom()<<"</efrom><eto>"<<r->getEdgeTo()<<"</eto>";
        out <<"<afrom>"<<std::setprecision(16)<<r->getTimeFrom()<<"</afrom><ato>"<<r->getTimeTo()<<"</ato>"<<std::setprecision(rprec);
        out <<"</recedge>"<<endl;
    }
    out << "</Iteration>" << endl;
    out.flush();
}

void Param::exportXMLend(ostream& out)
{
    out << "</outputFile>" << endl;
}

void Param::startDiagnostics(ostream& out)
{
    out<<"theta,delta,rho,numrecedge,edgeden,edgepb,edgevarpb,tmrca"<<endl;
}

void Param::updateDiagnostics(int iter,ostream& out)
{
    // routine diagnostics:
    esttheta=(double)(iter*esttheta + theta)/(iter+1.0);
    estthetasq=(double)(iter*estthetasq + theta*theta)/(iter+1.0);
    estrho=(double)(iter*estrho + rho)/(iter+1.0);
    estrhosq=(double)(iter*estrhosq + rho*rho)/(iter+1.0);
    estdelta=(double)(iter*estdelta + delta)/(iter+1.0);
    estdeltasq=(double)(iter *estdeltasq + delta*delta)/(iter+1.0);
    estnumrecedge=(double)(iter*estnumrecedge +rectree->numRecEdge())/(iter+1.0);
    estnumrecedgesq=(double)(iter *estnumrecedgesq + rectree->numRecEdge()*rectree->numRecEdge())/(iter+1.0);
    double curedgeden=rectree->numRecEdge()/rectree->getTTotal();
    estedgeden=(double)(iter*estedgeden + curedgeden)/(iter+1.0);
    estedgedensq=(double)(iter*estedgedensq + curedgeden*curedgeden)/(iter+1.0);

    double tmpexp=0,tmpexpsq=0;
    int numedges=rectree->getN()*2-1;
    for(int b=0;b<numedges;b++)
    {
        tmpexp+=(double)rectree->numRecEdgeOnBranch(b)/numedges;
        tmpexpsq+=(double)rectree->numRecEdgeOnBranch(b)*rectree->numRecEdgeOnBranch(b)/numedges;
    }
    estedgepb=(double)(iter*estedgepb + tmpexp)/(iter+1.0);
    estedgepbsq=(double)(iter*estedgepbsq + tmpexp*tmpexp)/(iter+1.0);
    estvaredgepb=(double)(iter*estvaredgepb + tmpexpsq-tmpexp*tmpexp)/(iter+1.0);
    estvaredgepbsq=(double)(iter*estvaredgepbsq + (tmpexpsq-tmpexp*tmpexp)*(tmpexpsq-tmpexp*tmpexp))/(iter+1.0);

    out<<theta<<","<<delta<<","<<rho<<","<<rectree->numRecEdge()<<","<<curedgeden<<","<<tmpexp<<","<<tmpexpsq-tmpexp*tmpexp<<","<<rectree->getNode(rectree->getN()*2-2)->getAge()<<endl;
}

void Param::computeSiteLL(unsigned int site,bool makeLT)
{
    //locll[site]=1.0;return;
    if (data==NULL)
        return;
    int n=data->getN();
    if (f.size()==0)
        f=vector<vector<double> >(n*2-1,vector<double>(4,0.0));
    if (makeLT)
        try
        {
            rectree->makeLocalTree(site,theta/data->getL());
        }
        catch(char * x)
        {
            cerr<<"Param::"<<x<<endl;
            throw("computesiteLL::makelocaltreeerror");
        }
    for (int i=n-1;i>=0;i--)
    {
        if (data->get(i,site)>3)
            for (int a=0;a<4;a++)
                f[i][a]=1.0;
        else
        {
            f[i][0]=0.0;
            f[i][1]=0.0;
            f[i][2]=0.0;
            f[i][3]=0.0;
            f[i][data->get(i,site)]=1.0;
        }
    }
    int y,z; // children of this site
    double ey,ez,eyeq,ezeq,eyneq,ezneq;
    for (int i=n;i<2*n-1;i++)
    {
        y=rectree->getSon(i,site,true ,&ey);
        z=rectree->getSon(i,site,false,&ez);
	if(y<0||z<0||y>=(int)f.size()||z>=(int)f.size()){
		// periodically warn the user about these problems, not every time it happens.
		static long badTreeWarning=0;
		if(badTreeWarning <= 100){	
			cerr<<"Error in computeSiteLL: Invalid Local Tree, happened " << badTreeWarning << " times" <<endl;
			if(badTreeWarning==100)  cerr << "Not reporting any more errors, fix this program!!\n";
		}
		badTreeWarning++;
		locll[site]=0;
		return;
/// WARNING!  This does NOT throw an error but ignores the site!
//		throw("computeSiteLL Invalid local tree Error");
	}
        eyeq =0.25*(1.0+3.0*ey);
        ezeq =0.25*(1.0+3.0*ez);
        eyneq=0.25*(1.0-ey);
        ezneq=0.25*(1.0-ez);
        for (unsigned int j=0;j<4;j++)
        {
            int j1=(j+1)%4,j2=(j+2)%4,j3=(j+3)%4;
            double YnotJ=f[y][j1]+f[y][j2]+f[y][j3];
            double ZnotJ=f[z][j1]+f[z][j2]+f[z][j3];
            f[i][j]=(f[z][j]*ezeq+ZnotJ*ezneq)*(f[y][j]*eyeq+YnotJ*eyneq);
        }
    }
    locll[site]=LOG025+log(f[2*n-2][0]+f[2*n-2][1]+f[2*n-2][2]+f[2*n-2][3]);
    //locll[site]=0; //*** NOTE: uncomment this to disconnect the likelihood
    if(isnan(locll[site]))
    {
        cerr<<"ERROR: computeSiteLL:Site "<<site<<"has NaN lik!"<<endl;
        for(int i=n;i<2*n-1;i++)
        {
            cout<<"i="<<i<<endl;
            for (unsigned int j=0;j<4;j++)
                cout<<"j="<<j<<" f[i][j]="<<f[i][j]<<endl;
        }
        cout<<"locll["<<site<<"]=log(0.25)+log("<<f[2*n-2][0]+f[2*n-2][1]+f[2*n-2][2]+f[2*n-2][3]<<")"<<endl;
        throw("NaN lik");
    }

}

void Param::computeLikelihood(unsigned int start, unsigned int end)
{
    if (data==NULL)
        return;
    bool doneOne=false;
    for (unsigned int site=start;site<end;site++)
    {
        ll-=locll[site];
        int prev=site-1;
        if (data->isPoly(site))
            goto noopti;
        //Shortcut. Find previous non-poly site, see if it has same local tree, and if so steal its ll
        while (prev>=0 && data->isPoly(prev))
            prev--;
        if (prev<0)
            goto noopti;//No previous non-poly site found
        for (unsigned int i=prev+1;i<=site;i++)
            if (!rectree->sameLocalTreeAsPrev(i))
                goto noopti;//prev and site don't have the same local tree
        if(isnan(locll[prev]))
        {
            cout<<"Site "<<prev<<" has NaN lik.  Tried to use it for site "<<site<< " isPoly="<<data->isPoly(prev)<<" sameasprev="<<rectree->sameLocalTreeAsPrev(prev)<<endl;
            throw "NaN site lik";
        }
        locll[site]=locll[prev];
        ll+=locll[site];
        if(isnan(locll[site]))
        {
            cerr<<"site "<<site<<" has NaN lik! isPoly="<<data->isPoly(site)<<" sameasprev="<<rectree->sameLocalTreeAsPrev(site)<<endl;
            throw;
        }
        continue;
noopti:
        if (site>start&&doneOne&&rectree->sameLocalTreeAsPrev(site))
	{
            computeSiteLL(site,false);
        }else
        {
            computeSiteLL(site,true);
            doneOne=true;
        }
        ll+=locll[site];
    }
}

void Param::computeLikelihood()
{
    if (data==NULL)
        return;
    computeLikelihood(0,data->getL());
}

double Param::getRlenPrior(int edge)
{
    double lprior=0;
    int start=rectree->getRecEdge(edge)->getStart();
    int end=rectree->getRecEdge(edge)->getEnd();
    int blockin=getData()->inblock(start);
    if(start==getData()->getBlocks()->at(blockin))
    {
        lprior=log((double)getDelta()/(getData()->getB()*(getDelta()-1)+getData()->getL()));
    }
    else
    {
        lprior=log(1.0/(getData()->getB()*(getDelta()-1)+getData()->getL()));
    }
    if(end==getData()->getBlocks()->at(blockin+1))
    {
        lprior+=(end-start)*log(1.0-1.0/getDelta());
    }
    else
    {
        lprior+=(end-start-1.0)*log(1.0-1.0/getDelta())-log(getDelta());
    }// note: end is 1 past the final element of the block
    return lprior;
}

void Param::testTree()
{
    dlog(1)<<"Testing tree in param"<<endl;
    rectree->testTree();
    double ll=getLL();
    dlog(1)<<"Testing local likelihoods in param"<<endl;
    for(int site=0;site<rectree->getL();site++)
    {
        double test=getLLsite(site);
        computeLikelihood(site,site+1);
        if(fround(test,5)!=fround(getLLsite(site),5))
        {
            cout<<"WARNING in param:testTree: ll for site "<<site<<" is not consistent! before="<<test<<" after="<<getLLsite(site)<<endl;
            for(int e =0;e<rectree->numRecEdge();e++)
            {
                cout<<"Edge "<<e<<" "<<rectree->getRecEdge(e)->getStart()<<":"<<rectree->getRecEdge(e)->getEnd()<<endl;
            }
            ;
            throw;
        }
    }
    if(fround(ll,5)!=fround(getLL(),5))
    {
        cout<<"Error in param:testTree: initial LogLL "<<ll<<" not equal to final LogLL "<<getLL()<<endl;
        throw;
    }
}

void Param::greedyTreeMove(){
	cout<<"Performing greedy move"<<endl;
	MoveGreedyTree* mgreedy = new MoveGreedyTree(this,0.0);
	mgreedy->move();
	delete(mgreedy);
}

vector<double> Param::greedyCalcDists(){
	MoveGreedyTree* mgreedy = new MoveGreedyTree(this,0.0);
	vector<double> ret=mgreedy->calcDists();
	delete(mgreedy);
	return(ret);
}

vector<double> Param::greedyCalcDists(vector<double> nmuts,vector<double> nsites){
	MoveGreedyTree* mgreedy = new MoveGreedyTree(this,0.0);
	vector<double> ret=mgreedy->calcDists(nmuts,nsites);
	delete(mgreedy);
	return(ret);
}

void Param::greedyApply(vector<double> dists){
	MoveGreedyTree* mgreedy = new MoveGreedyTree(this,0.0);
	vector<double> newage=mgreedy->calcAges(dists);
	mgreedy->applyChanges(newage);
	delete(mgreedy);
}

vector< vector<double> > Param::greedyDetails(vector< vector<double> > * pairwise){
	int N=rectree->getN();
	vector< vector<double> > ret;
	MoveGreedyTree* mgreedy = new MoveGreedyTree(this,0.0);
	for(int c1=0;c1<8;c1++) ret.push_back(vector<double>(2*N - 1,0.0));

	for(int c1=N;c1<2*N - 1;c1++) {
	vector<int> all2=rectree->getAllSampledSeqs(rectree->getNode(c1)->getLeft()->getId());
	vector<int> all3=rectree->getAllSampledSeqs(rectree->getNode(c1)->getRight()->getId());
	for(unsigned int c2=0;c2<all2.size();c2++) {
		for(unsigned int c3=0;c3<all3.size();c3++) {
			int curindex=mgreedy->getPairIndex(all2[0],all3[0],N);
			ret[0][c1]+=pairwise->at(0)[curindex]/all2.size()/all3.size();
			ret[1][c1]+=pairwise->at(1)[curindex]/all2.size()/all3.size();
			ret[2][c1]+=pairwise->at(2)[curindex]/all2.size()/all3.size();
			ret[3][c1]+=pairwise->at(3)[curindex]/all2.size()/all3.size();
			ret[4][c1]+=pairwise->at(0)[curindex]*pairwise->at(0)[curindex]/all2.size()/all3.size();
			ret[5][c1]+=pairwise->at(1)[curindex]*pairwise->at(1)[curindex]/all2.size()/all3.size();
			ret[6][c1]+=pairwise->at(2)[curindex]*pairwise->at(2)[curindex]/all2.size()/all3.size();
			ret[7][c1]+=pairwise->at(3)[curindex]*pairwise->at(3)[curindex]/all2.size()/all3.size();
	}
	}}
	delete(mgreedy);
	return(ret);
}

vector< vector<double> > * Param::greedyPairwiseDetails(){
	int N=rectree->getN();
	MoveGreedyTree* mgreedy = new MoveGreedyTree(this,0.0);
	vector<vector<bool> >* lcf = mgreedy->localClonalFrame();

	vector<vector<double> > * ret = new vector <vector<double> > (4,vector<double>(lcf->size(),0.0)) ;
	// this is the sum over: nsites, nmuts, nrecsites, nrecmuts
	vector<int> nrecs(lcf->size(),0);

	long pairon=0;
	for(int c1=0;c1<N;c1++){ for(int c2=c1+1;c2<N;c2++){
	  nrecs[pairon]=mgreedy->recCount(c1,c2);
	  for(unsigned int c3=0;c3<lcf->at(pairon).size();c3++) {
	    if(lcf->at(pairon)[c3]){
		ret->at(0)[pairon]+=1.0;
		if(getData()->get(c1,c3)!=getData()->get(c2,c3)) ret->at(1)[pairon]+=1.0;
	    }else{
		ret->at(2)[pairon]+=1.0;
		if(getData()->get(c1,c3)!=getData()->get(c2,c3))ret->at(3)[pairon]+=1.0;
	    }
	  }
	  pairon++;
	}}
	delete(lcf);
	delete(mgreedy);
	return(ret);
}

double Param::empiricalTheta(vector< vector<double> > * mutpairwise){
	return(getTheta());
	//greedyApply(dists);
	double sum=0;
	return(sum);
}


} // end namespace weakarg
