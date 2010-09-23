#ifndef MOVEREMEDGE_H
#define MOVEREMEDGE_H

#include "move.h"

namespace weakarg
{

/**
    @brief This move removes an edge
*/
class MoveRemEdge : public Move
{
public:
    MoveRemEdge(Param * p,double a);
    Move * clone()
    {
        return new MoveRemEdge(*this);
    }
    int move(vector<int> * samplespace=NULL);
    inline int move(){return(move(NULL));}
    ~MoveRemEdge();

};

} // end namespace weakarg
#endif
