#ifndef MOVEGREEDYTREE_H
#define MOVEGREEDYTREE_H
//
#include "move.h"
//

namespace weakarg
{

/**
    @brief This move does a greedy update of the tree branches.  It IS NOT AN MCMC MOVE!
*/
class MoveGreedyTree : public Move
{

public:
    MoveGreedyTree(Param * p,double a);
    Move * clone()
    {
        return new MoveGreedyTree(*this);
    };
    int move();
    ~MoveGreedyTree();

    inline vector<int> getAllTips(int s,int N,RecTree * t){
	vector<int> tlist=t->getAllChildren(s);
	vector<int> ret;
	for(unsigned int c1=0;c1<tlist.size();c1++) {
		if(tlist[c1]<N) ret.push_back(tlist[c1]);
	}
	return(ret);
    }///< gets all tips beneath a node on the tree
    inline int getPairIndex(int a, int b, int N){
	if(a==b) throw("getPairIndex:Impossible to pair individual with itself!");
	if(a>b){int c=a;a=b;b=c;}// now b>a
	int ret=0;
	for(int c1=0;c1<N;c1++){
	  for(int c2=c1+1;c2<N;c2++){
		if((c1==a) &&(c2==b)) return(ret);
		ret ++;
	  }
	}
	cerr<<a<<" or "<<b<<" is not less than N="<<N<<endl;
	throw("getPairIndex:Not found a possible list index!");
	return(-1);
    }///< gets the pair index used internally to store the set of all tip pairings as a flat list
    int recCount(int s1, int s2);///< Counts the number of recombinations between two sequences
    vector<int> edgesBetween(int s1, int s2);///< returns a list of the edges that are between two sequences
    vector<double> calcDists(vector<double> nmuts,vector<double> nsites);///calculates distances based on provided details
    vector<double> calcDists(double * esttheta=NULL);///< Calculates the distances between each pair of sequences
    vector<double> calcAges(vector<double> dists);///< Calculate the correct ages based on the distances
    void applyChanges(vector<double> newage);///< Updates the ages of a tree to the newages
    vector<vector<bool> >* localClonalFrame();///< returns the matrix of sites in the clonal frame between each pair of sequences
};

} // end namespace weakarg
#endif
