#ifndef MOVE_H
#define MOVE_H
#include <cmath>
#include "param.h"
#include "mpiutils.h"
#include <gsl/gsl_rng.h>

namespace weakarg
{
extern gsl_rng * rng;

/**
    @brief This virtual class represents a move of the MCMC
*/
class Move
{
protected:
    std::string name;///<The moves name
    std::string description;///<The moves description
    Param*param;///<Parameter set on which the move is performed
    double alpha;///<Weighting for this move (probability is alpha/sum_i(alpha_i) )
    int numcalls;///<Number of attempts at the move
    int numaccept;///<Number of acceptances of the move
public:
    Move(Param*p,double a=1.0);///<Constructor specifying on which Parameter the move is to be performed
    virtual Move * clone()=0;
    virtual ~Move()=0;
    virtual int move()=0;///<Performs the move once
    virtual inline int getTopoCounts(){return(-1);};
    virtual inline int getTopoAcc(){return(-1);};
    virtual inline int move(vector<int> *samplespace){cout<<"Problem in Move!"<<endl;return(move());};///<Performs the move once
    inline void setParam(Param * param)
    {
        this->param=param;
    }///<Change the Parameter on which the move is to be performed
    inline std::string desc()
    {
        return description;
    }///<Returns a description of the move
    inline std::string getName()
    {
        return name;
    }///<Returns a description of the move
    inline void setName(std::string newname)
    {
        name=newname;
    }
    inline void setAlpha(double a)
    {
        alpha=a;
    }
    ;///<Sets alpha for the move
    inline double getAlpha()
    {
        return(alpha);
    }
    inline int getCounts()
    {
        return(numcalls);
    }
    inline int getAcc()
    {
        return(numaccept);
    }

    ;///<Returns alpha for the move
};

} // end namespace weakarg
#endif
