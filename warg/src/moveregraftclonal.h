#ifndef MOVEREGRAFTCLONAL_H
#define MOVEREGRAFTCLONAL_H
//
#include "move.h"
#include "metropolis.h"
//
namespace weakarg
{

/**
    @brief This move regrafts a clonal node to a new parent
*/
class MoveRegraftClonal : public Move
{
protected:
    int reps;
    double tempering;
	RecTree* backuptree;
public:
    MoveRegraftClonal(Param * p,double a,int r=20,double t=2.0);
    Move * clone()
    {
        return new MoveRegraftClonal(*this);
    };
    int move();
    ~MoveRegraftClonal();
    inline void setReps(int r){reps=r;};///< Sets the number of reps
    int getAlive(int which);///<Gets a random OTHER alive node at the time of which
    void regraft(int which,int whichto);///< regrafts which to have father whichto, keeping the same time of the event
    inline bool validNode(int which){
        RecTree * rectree=param->getRecTree();
	if(which==rectree->getRoot()->getId() || which == rectree->getRoot()->getLeft()->getId() || which==rectree->getRoot()->getRight()->getId()) return(false);
	return(true);
    }
};

} // end namespace weakarg
#endif
