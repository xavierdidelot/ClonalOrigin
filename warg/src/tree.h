#ifndef TREE_H
#define TREE_H
#include "node.h"
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
#include "data.h"

using namespace std;

namespace weakarg
{

extern gsl_rng * rng;

/**
    @brief Genealogy of the isolates under study
*/
class Tree
{
protected:
    Node * root;///<Root node of the genealogy
    std::vector<Node*> nodes;///<Vector of all the nodes in the genealogy
    int n;///<Number of isolates
    double ttotal;///<Sum of branch lengths
public:
    Tree(std::string newick,bool isFilename=true,bool forceages=true);///<Creates a tree from a Newick file
    Tree(int n);///< Creates a Coalescent tree by simulation
    Tree(Tree *intree);///< Creates a copy of a tree via newick (== slow)
    Tree(Data * data);///< Creates a UPGMA tree
    Tree(const Tree& t);///< Copy constructor
    void assign(const Tree& t);///< Copy a tree while allocating as little new storage as possible
    ~Tree();
    std::string newick(int p=6) const; ///<Returns a Newick description of the tree (precision options)
    std::string newickNoInternalLabels(int p=6) const; ///<Returns a Newick description of the tree (precision options)
    void makeFromNewick(std::string newick,bool forceages=false);///< Creates a tree from a newick string
    inline int getN() const
    {
        return n;
    }///<Returns the number of isolates
    inline Node* getNode(int i) const
    {
        return nodes[i];
    }///<Returns a node given its index
    inline Node* getRoot() const
    {
        return root;
    }///<Returns the root node
    inline double getDist(int i) const
    {
        return nodes[i]->getDist();
    }///<Returns the time of a node given its index
    double prior() const; ///<Returns log-prior of the tree
    inline double getTTotal() const
    {
        return ttotal;
    }///<Returns the sum of branch lengths
    void computeTTotal(); ///<Computes the sum of branch lengths
    int getPoint(double * dist, std::vector<int> * samplespace=NULL) const; ///<Returns a point chosen uniformly at random on the tree, or uniformally on a subspace of specified clonal edges
    std::vector<int> getAllChildren(int e);///<Returns a vector of the indices of the children of a given node
    std::vector<int> getAllSampledSeqs(int e);///<Returns a vector of all the observed sequences of the children of a given node

    void orderNodes(double dist=1);///<orders the list by age
    int orderNodes(int which, double dist);///<places a node in a new place in the list
    void swapFather(int a,int b);///<Swaps the father of two nodes
    void swapNode(int a, int b);///<Swaps two nodes in the list
    int getOldestReversedNode();///<Returns the oldest node that has been age reversed

    void testNodeAges() const; ///<Tests node ages
    double tavare() const;

    inline int otherChild(int par,int onechild){
	if(getNode(par)->getLeft()->getId()==onechild) return(getNode(par)->getRight()->getId());
	else if(getNode(par)->getRight()->getId()==onechild) return(getNode(par)->getLeft()->getId());
	else{cerr<<"Error in RecTree::otherChild: parent doesn't have a suitable child!"<<endl;throw("unknown child error");}
    };
    std::vector<int> getMinEdgeList(std::vector<int> seqs);
    int getNextGroup(std::vector<int> *seqs);
    bool isParentToOnly(int e, std::vector<int> seqs,std::vector<int> *which);

};

} // end namespace weakarg
#endif
