#include "movesitechange.h"

#define MAX_MOVESIZE 10
//The number of sites we'll consider for a "nearly" Gibbs move

using namespace std;
namespace weakarg
{

MoveSiteChange::MoveSiteChange(Param * p,double a)
        : Move(p,a)
{
   description= "Update moving the location of the recedge on the genome";
    name="SiteChange";
}

int MoveSiteChange::move(vector<int> * samplespace)
{
    RecTree * rectree=param->getRecTree();
   int which=rectree->sampleEdge(samplespace);//floor(gsl_rng_uniform(rng)*rectree->numRecEdge());
    if (which<0)
    {
    	dlog(1)<<"No valid edges to change!"<<endl;
        return(-1);
    }
    numcalls++;
    int dsign=gsl_rng_uniform_int(rng,2)*2-1;//plus or minus one
    int movestart=0; // whether the start or end points are moved
    if(gsl_rng_uniform(rng)<0.5)
        movestart=1;
    double ll=param->getLL();

    int start =rectree->getRecEdge(which)->getStart();
    int end   =rectree->getRecEdge(which)->getEnd();
    int blockin=param->getData()->inblock(start);
    vector<double> store1,store2,store3,store4;// the likelihoods at intermediate calculations
    // calculate the probability of moving to a site
    vector<double> v1=calcPost(&store1,&store2,which,dsign,movestart,MAX_MOVESIZE,0);// move but don't change. store2 contains the liks of after the change
    double sv1=sumVec(v1);
    // choose from the posterior distribution of the moves
    double x=gsl_rng_uniform(rng);
    for(int i=0;i<(int)v1.size();i++)
    {
        x-=v1[i]/sv1;
        if(x<0)
        {
            //cout<<"tt1:";param->testTree();
            vector<double> v2=calcPost(&store3,&store4,which,-dsign,movestart,MAX_MOVESIZE-i,0);
            // calculate the probability of the EXTRA moves in the reverse direction
            double sv2=sumVec(v2) -v2[0] + sumVec(v1,i+1);

            // activate the proposed move
            if(movestart)
            {
                rectree->setStart(which,start+i*dsign);
            }
            else
            {
                rectree->setEnd(which,end+i*dsign);
            }
            param->setLL(ll+sumVec(store2,i)-sumVec(store1,i));
            RestoreLik(&store2,dsign,start,end,movestart,0,i);
            //cout<<"tt2:";param->testTree();
            if(gsl_rng_uniform(rng)<sv1/sv2)
            {// rejection step
            	dlog(1)<<"Moved edge "<<which<<" genetic location from "<<start<<":"<<end<<" to "<<start+i*dsign*movestart<< ":"<<end+i*dsign*(1.0-movestart)<<" in block from "<< param->getData()->getBlocks()->at(blockin)<<":"<<param->getData()->getBlocks()->at(blockin+1)<<endl;
		numaccept++;
            	return(1);
            }
            else
            { // put the edge back as it was
            	dlog(1)<<"Rejected moving edge "<<which<<" genetic location from "<<start<<":"<<end<<" to "<<start+i*dsign*movestart<< ":"<<end+i*dsign*(1.0-movestart)<<" in block from "<< param->getData()->getBlocks()->at(blockin)<<":"<<param->getData()->getBlocks()->at(blockin+1)<<endl;
                rectree->setStart(which,start);
                rectree->setEnd(which,end);
                param->setLL(ll);
                RestoreLik(&store1,dsign,start,end,movestart);
            	return(0);
            }
            //cout<<"tt3:";param->testTree();
        }
    }
    if(x>0)
        cerr<<"Error in movesitechange: probability<1 observed!"<<endl;
    throw;
    return(-1);
}

double MoveSiteChange::sumVec(vector<double> vec,int vmax)
{
    double sum=0;
    if(vmax<0)
        vmax = (int)vec.size();
    for(int i=0;i<vmax;i++)
    {
        sum+=vec[i];
    }
    return(sum);
}

vector<double> MoveSiteChange::calcPost(vector<double> *store1,vector<double> *store2,int which, int dsign,int movestart,int maxdist, bool keeplik)
{
    RecTree * rectree=param->getRecTree();
    int dmod=0;
    if(dsign>0)
    {
        dmod=-1;
    }// modification for direction: the ll of the site to be updated isn't always start/end site
    int start=rectree->getRecEdge(which)->getStart();
    int end  =rectree->getRecEdge(which)->getEnd();
    int blockin=param->getData()->inblock(start);
    double oldv=0;// the last value the ll taken
    vector<double> post;//the posteriors and likelihood

    // calculate the initial log posterior
    int imax=maxdist;
    double ll=param->getLL();
    if(dsign>0)
    {
        if(movestart)
        {
            imax=min(imax,end-start-1);
        }
        else
        {
            imax=min(imax,param->getData()->getBlocks()->at(blockin+1)-end);
        }
    }
    else
    {
        if(movestart)
        {
            imax=min(imax, start-param->getData()->getBlocks()->at(blockin));
        }
        else
        {
            imax=min(imax,end-start-1);
        }
    }
    post.push_back(exp(param->getRlenPrior(which)));
    // calculate the log posterior for the moves available
    for (int i=1;i<=imax;i++)
    {
        if(movestart)
        {
            store1->push_back(param->getLLsite(i*dsign+start+dmod));
            oldv-=store1->back();
            rectree->setStart(which,start+i*dsign);
            param->computeLikelihood(i*dsign+start+dmod,i*dsign+start+dmod+1);
            store2->push_back(param->getLLsite(i*dsign+start+dmod));
            oldv+=param->getLLsite(i*dsign+start+dmod);
        }
        else
        {
            store1->push_back(param->getLLsite(i*dsign+end+dmod));
            oldv-=param->getLLsite(i*dsign+end+dmod);
            rectree->setEnd(which,end+i*dsign);
            param->computeLikelihood(i*dsign+end+dmod,i*dsign+end+dmod+1);
            store2->push_back(param->getLLsite(i*dsign+end+dmod));
            oldv+=param->getLLsite(i*dsign+end+dmod);
        }
        double lprior=param->getRlenPrior(which);
        post.push_back(exp(oldv+lprior));
    }
    if(!keeplik)
    {//restore the lik structure
        rectree->setStart(which,start);
        rectree->setEnd(which,end);
        param->setLL(ll);
        RestoreLik(store1,dsign,start,end,movestart);
    }
    return(post);//return the posterior distribution
}

void MoveSiteChange::RestoreLik(vector<double> *store,int dsign,int start,int end,int movestart,int imin,int imax)
{
    int dmod=0;
    if(dsign>0)
    {
        dmod=-1;
    }// modification for direction: the ll of the site to be updated isn't always start/end site
    if(imax<0)
        imax=store->size();
    for(int i=imin;i<imax;i++)
    {
        if(movestart)
            param->setlocLL((i+1)*dsign+start+dmod,store->at(i));
        else
            param->setlocLL((i+1)*dsign+end+dmod,store->at(i));
    }
}

MoveSiteChange::~MoveSiteChange()
{}


} // end namespace weakarg
