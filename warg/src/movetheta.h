#ifndef MOVETHETA_H
#define MOVETHETA_H

#include "move.h"

namespace weakarg
{

/**
    @brief This move updates theta
*/
class MoveTheta : public Move
{
public:
    MoveTheta(Param * p,double a);
    Move * clone()
    {
        return new MoveTheta(*this);
    };
    int move();
    ~MoveTheta();


};


} // end namespace weakarg
#endif
