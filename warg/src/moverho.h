#ifndef MOVERHO_H
#define MOVERHO_H

#include "move.h"

namespace weakarg
{

/**
    @brief This move updates rho
*/
class MoveRho : public Move
{
public:
    MoveRho(Param * p,double a);
    Move * clone()
    {
        return new MoveRho(*this);
    };
    int move();
    ~MoveRho();
};

} // end namespace weakarg
#endif
