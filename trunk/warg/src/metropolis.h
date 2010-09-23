#ifndef METROPOLIS_H
#define METROPOLIS_H
//
#include "moveedgechange.h"
#include "movetimechange.h"
#include "movesitechange.h"
#include "moveaddedge.h"
#include "moveremedge.h"
//
namespace weakarg
{
/**
    @brief Performs n, possibly tempered, possibly local moves
*/
class Metropolis 
{

public:
    Metropolis(Param * p);
    void move(int n, double temper=1.0, vector<int> * samplespace=NULL);
    ~Metropolis();
    inline std::string desc()
    {
        return "Performs n local, tempered moves";
    }
protected:
    int chooseMove(vector<Move*> *mall,double totweight);///<Chooses the next move
    Param * param;
};

} // end namespace weakarg
#endif
