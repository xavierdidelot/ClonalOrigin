#ifndef MOVESITECHANGE_H
#define MOVESITECHANGE_H

#include "move.h"

namespace weakarg
{

/**
    @brief This move changes the start and end locations on the genome by a constant factor according to a gibbs move
*/
class MoveSiteChange : public Move
{
public:
    MoveSiteChange(Param * p,double a);
    Move * clone()
    {
        return new MoveSiteChange(*this);
    }
    int move(vector<int> * samplespace=NULL);
    inline int move(){return(move(NULL));}
    ~MoveSiteChange();
    double sumVec(std::vector<double> vec,int vmax=-1);///<Returns the sum of a vector up to vmax or the end of the vector
    void RestoreLik(std::vector<double> *store,int dsign,int start,int end,int movestart,int imin=0,int imax=-1);///< Restores the likelihood from a stored vector between a given range
    std::vector<double> calcPost(std::vector<double> *store1,std::vector<double> *store2,int which, int dsign,int movestart,int maxdist, bool keeplik);///<Calculates the posterior distribution of a move
 
};

} // end namespace weakarg
#endif
