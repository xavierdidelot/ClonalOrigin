#ifndef MOVEAGECLONAL_H
#define MOVEAGECLONAL_H
//
#include "move.h"
#include "metropolis.h"
#define SQRT2PI 2.506628274631000241612355239340104162693023681640625
//
namespace weakarg
{
/**
    @brief This move updates the age of a clonal node
*/
class MoveAgeClonal : public Move
{
protected:
    bool changeTopology;// we can forbid topology changes with this (no tempering is then needed)
	RecTree* backuptree;
    bool usenormalprop;
    double s0,s1;
    int numtopo;
    int numtopoaccept;
    int reps;
    double tempering;
    double logQrat(double t1,double t2);
    inline double logdens(double tfrom,double tto){
		double delta=tto-tfrom;
	return( -(delta*delta)/2.0/(tfrom*s1+s0)/(tfrom*s1+s0)  -log(tfrom*s1+s0)-log(sqrt(2*M_PI)));
    }
    inline double dens(double tfrom,double tto){
		double delta=tto-tfrom;
	return( exp(-(delta*delta)/2.0/(tfrom*s1+s0)/(tfrom*s1+s0))/(tfrom*s1+s0)/SQRT2PI);
    }
public:
    MoveAgeClonal(Param * p,double a,bool allowclonal=true,int r=20,double t=2.0);
    Move * clone()
    {
        return new MoveAgeClonal(*this);
    };
    int move();
    ~MoveAgeClonal();
    inline int getTopoCounts()
    {
        return(numtopo);
    }
    inline int getTopoAcc()
    {
        return(numtopoaccept);
    }
};

} // end namespace weakarg
#endif
