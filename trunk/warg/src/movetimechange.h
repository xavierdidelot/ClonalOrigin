#ifndef MOVETIMECHANGE_H
#define MOVETIMECHANGE_H

#include "move.h"

namespace weakarg
{

/**
    @brief This move changes the start and end times of a recedge
*/
class MoveTimeChange : public Move
{
public:
    MoveTimeChange(Param * p,double a);
    Move * clone()
    {
        return new MoveTimeChange(*this);
    }
    int move(vector<int> * samplespace=NULL);
    inline int move(){return(move(NULL));}
    ~MoveTimeChange();

};

} // end namespace weakarg
#endif
