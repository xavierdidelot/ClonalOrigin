#include "moverho.h"

using namespace std;
namespace weakarg
{

MoveRho::MoveRho(Param * p,double a)
        : Move(p,a)
{
    description= "Update of the recombination rate rho";
    name= "Rho";
}

int MoveRho::move()
{
    int r=param->getRecTree()->numRecEdge();
    double T=param->getRecTree()->getTTotal();
    double rho=param->getRho();
    double rho2=gsl_ran_gamma(rng,1.0+r,1.0/(param->hyperPriorOfRho()+T*0.5));
    param->setRho(rho2);
    dlog(1)<<"Gibbs update of rho from "<<rho<<" to "<<rho2<<"..."<<endl;
    numcalls++;numaccept++;
    return(1);
}

MoveRho::~MoveRho()
{}


} // end namespace weakarg
