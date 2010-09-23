#include "movetheta.h"

using namespace std;
namespace weakarg
{

MoveTheta::MoveTheta(Param * p,double a)
        : Move(p,a)
{
    description= "Update of the mutation rate theta";
    name="Theta";
}

int MoveTheta::move()
{ 
    numcalls++;
    double t=param->getTheta();
    double t2=t+(gsl_rng_uniform(rng)-0.5)*10.0;
    if (t2<0.0)
        t2=-t2;
    vector<double> store(param->getData()->getL());
    for (int i=0;i<param->getData()->getL();i++)
        store[i]=param->getLLsite(i);
    double l=param->getLL();
    param->setTheta(t2);
    param->computeLikelihood();
    double l2=param->getLL();
    dlog(1)<<"Proposing to move theta from "<<t<<" to "<<t2<<"...";
    if (log(gsl_rng_uniform(rng))>l2-l)
    {
    	dlog(1)<<"Rejected!"<<endl;
        param->setTheta(t);
        for (int i=0;i<param->getData()->getL();i++)
            param->setlocLL(i,store[i]);
        param->setLL(l);
	return(0);
        //param->computeLikelihood();
    }
    else dlog(1)<<"Accepted!"<<endl;
    numaccept++;
    return(1);
}

MoveTheta::~MoveTheta()
{}


} // end namespace weakarg

