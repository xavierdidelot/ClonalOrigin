#ifndef MOVESCALETREE_H
#define MOVESCALETREE_H
//
#include "move.h"
//

namespace weakarg
{

/**
    @brief This move scales the whole tree to a new TMRCA
*/
class MoveScaleTree : public Move
{

public:
    MoveScaleTree(Param * p,double a);
    Move * clone()
    {
        return new MoveScaleTree(*this);
    };
    int move();
    ~MoveScaleTree();

    void scaleTree(double scale);///<Scales the tree by a factor scale
};

} // end namespace weakarg
#endif
