#include "movedelta.h"

using namespace std;
namespace weakarg
{

MoveDelta::MoveDelta(Param * p,double a)
        : Move(p,a)
{
    description= "Update of the mean recombination length delta";
    name= "Delta";
}

int MoveDelta::move()
{
    double olddelta=param->getDelta();
    numcalls++;
    int X=0,Y=0;
    for (int i=0;i<param->getRecTree()->numRecEdge();i++)
    {
        Y+=param->getRecTree()->getRecEdge(i)->getEnd()-param->getRecTree()->getRecEdge(i)->getStart();
        if (param->getData()->isBegEnd(param->getRecTree()->getRecEdge(i)->getStart()))
            X++;
        if (param->getData()->isBegEnd(param->getRecTree()->getRecEdge(i)->getEnd()))
            Y--;else X--;/*{
            X--;
            Y--;
        }*/
    }
    for (int rep=0;rep<100;rep++)
    {
        double t=param->getDelta();
        double t2=t+(gsl_rng_uniform(rng)-0.5)*10.0;
        if (t2<1.0)
            t2=2.0-t2;
        int b=param->getData()->getB();
        int L=param->getData()->getL();
        int R=param->getRecTree()->numRecEdge();
        double ratio=X*(log(t2)-log(t));
        ratio+=Y*(log(1.0-1.0/t2)-log(1.0-1.0/t));
        ratio-=R*(log(b*t2+L-b)-log(b*t+L-b));
        if (log(gsl_rng_uniform(rng))<ratio)
            param->setDelta(t2);
    }
    dlog(1)<<"Delta was moved from "<<olddelta<<" to "<<param->getDelta()<<"."<<endl;
    numaccept++;
    return(1);
}

MoveDelta::~MoveDelta()
{}


} // end namespace weakarg
