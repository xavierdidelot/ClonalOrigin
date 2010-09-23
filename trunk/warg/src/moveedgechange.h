#ifndef MOVEEDGECHANGE_H
#define MOVEEDGECHANGE_H
//
#include "move.h"

namespace weakarg
{

/**
    @brief This move changes the edge of arrival/departunre without changing the time
*/
class MoveEdgeChange : public Move
{
public:
    MoveEdgeChange(Param * p,double a);
    Move * clone()
    {
        return new MoveEdgeChange(*this);
    }
    int move(vector<int> * samplespace=NULL);
    inline int move(){return(move(NULL));}
    ~MoveEdgeChange();
 
};

} // end namespace weakarg
#endif
