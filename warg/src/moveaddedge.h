#ifndef MOVEADDEDGE_H
#define MOVEADDEDGE_H

#include "move.h"

namespace weakarg
{

/**
    @brief This move adds an edge
*/
class MoveAddEdge : public Move
{
public:
    MoveAddEdge(Param * p,double a);
    Move * clone()
    {
        return new MoveAddEdge(*this);
    }
    int move(vector<int> * samplespace=NULL);
    inline int move(){return(move(NULL));}
    ~MoveAddEdge();

};

} // end namespace weakarg
#endif
