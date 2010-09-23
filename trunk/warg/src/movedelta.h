#ifndef MOVEDELTA_H
#define MOVEDELTA_H
//
#include "move.h"

namespace weakarg
{

/**
    @brief This move updates delta
*/
class MoveDelta : public Move
{
public:
    MoveDelta(Param * p,double a);
    Move * clone()
    {
        return new MoveDelta(*this);
    };
    int move();
    ~MoveDelta();
};

} // end namespace weakarg
#endif
