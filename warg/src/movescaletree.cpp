#include "movescaletree.h"

//#define DEBUG

using namespace std;
namespace weakarg
{

MoveScaleTree::MoveScaleTree(Param * p,double a)
        : Move(p,a)
{
    description="Update all nodes and edges in the tree to a new rescaled time";
    name="ScaleTree";
}

inline double fround(double n, double d)
{
    return floor(n * pow(10., d) + .5) / pow(10., d);
}

int MoveScaleTree::move()
{
    numcalls++;
    // Propose to change all nodes and recedges by a scaled amount
    RecTree * rectree=param->getRecTree();
    double lprior=rectree->prior(param)+param->logPriorOfRho()+param->logPriorOfTheta();
    // update the tree
    double logscale=param->getScaleTreeSize() *(gsl_rng_uniform(rng) * 2.0 -1.0);
    double scale =exp(logscale);
    dlog(1)<<"Proposing to scale tree TMRCA by "<<scale<<"...";
    scaleTree(scale);
    param->setRho(param->getRho()/scale);
    param->setTheta(param->getTheta()/scale);
    double newlprior=rectree->prior(param) + param->logPriorOfRho() + param->logPriorOfTheta();
    if(log(gsl_rng_uniform(rng))>newlprior-lprior+(2.0*rectree->numRecEdge()+rectree->getN()-3.0)*logscale)
    {
    	dlog(1)<<" Rejected!"<<endl;
        scaleTree(1.0/scale);
    param->setRho(param->getRho()*scale);
    param->setTheta(param->getTheta()*scale);
#if defined DEBUG
        //test its still ok
        int tmp=0;
        for (int i=0;i<param->getData()->getL();i++)
            if(fround(param->getLLsite(i),5)!=fround(store[i],5))
            {
                tmp=1;
                cout<<"Site "<<i<<" has ll before "<<store[i]<<" and after "<<param->getLLsite(i)<<" ";
                for (int j=0;j<rectree->numRecEdge();j++)
                {
                    if (rectree->getRecEdge(j)->affectsSite(i))
                        cout<<j<<" ";
                }
                cout<<endl;
            }
        if(fround(ll,5)!=fround(param->getLL(),5))
        {
            cout<<"Total ll before "<<ll<<" and after "<<param->getLL()<<endl;
        }
        try
        {
            param->testTree();
        }
        catch(char * x)
        {
            cout<<x<<endl<<"Movescaletree restore: broke the log liks"<<endl;
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
    else dlog(1)<<" Accepted!"<<endl;// accept the modified tree
    numaccept++;
    return(1);
}

void MoveScaleTree::scaleTree(double scale)
{
    RecTree * rectree=param->getRecTree();
    for(int c1=0;c1<2*rectree->getN()-1;c1++)
    {
        rectree->getNode(c1)->setAge(rectree->getNode(c1)->getAge()*scale);
    }
    rectree->computeTTotal();
    for(int c1=0;c1<rectree->numRecEdge();c1++)
    {
        rectree->getRecEdge(c1)->setTimeFrom(rectree->getRecEdge(c1)->getTimeFrom()*scale);
        rectree->getRecEdge(c1)->setTimeTo(rectree->getRecEdge(c1)->getTimeTo()*scale);
    }
}

MoveScaleTree::~MoveScaleTree()
{}



} // end namespace weakarg

