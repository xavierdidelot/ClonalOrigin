#ifndef RECTREE_H
#define RECTREE_H
#include "tree.h"
#include "recedge.h"
#include "weakarg.h"
#include "mpiutils.h"
#include "wargxml.h"

#include <gsl/gsl_randist.h>
#define MAXNODES 100000

using namespace std;

namespace weakarg
{

class Param;

/**
    @brief Recombinant tree using Tree class
*/

class RecTree:public Tree
{
protected:
// Variables
    std::vector<std::vector<int> > tabSons;///<Contains the branching order for the current local tree
    std::vector<std::vector<double> > tabSonsDist;///<Contains the distances for the current local tree
    std::vector<int> tabNode;///<Used only in makeLocalTree (made global for speed)
    std::vector<int> tabRec;///<Used only in makeLocalTree (made global for speed)
    std::vector<int> tabFather;///<Used only in makeLocalTree (made global for speed)
    std::vector<double> age;///<Used only in makeLocalTree (made global for speed)
    std::vector<int> tabEdge;///<Used only in makeLocalTree (made global for speed)
    std::vector<bool> sameLTasPrev;///<Indicates for each site whether the local tree is the same as that of the previous site
    unsigned int L;///< Number of sites in the sequences
    std::vector<RecEdge *> edge;///< The recombinant edges
// Functions
    int randomActiveNode(std::vector<int> nodelist);///<Returns a random node id from a vector of indicators for alive (-1 for not alive at current time)
    int getEdgeCoal(std::vector<int> nl,int k,double* time);///<Returns the node the recombinant edge is from and sets the time

    int moveClonalFixFrom(std::vector<int> ornodes=std::vector<int>(0));///< Fixes all from times to be in to valid nodes and reinstates the vector ornodes. returns the number of recedges moved to younger nodes - number moved to older nodes (or -1 if this is not calculated)

    void updateEdgeTimes(int e,double reldist);///moves an edge when its node changes by proportion reldist
    void scaleEdges(int which, double dist);///<Scales the recedges "to" on edge which by an amount dist
    void fixFromTimes(int affedge,double dist);///<Fixes from times when we've moved edge affedge an amount dist
    void fixToTimes(int which, double dist);///<Fixes "to" times when we've changed an node time implicitly

    void swapEdge(int a, int b);///<Swaps two edges in the edge list (ordering only)
    void swapEdgeTo(int a, int b);///<Swaps the edges going to *two specified --nodes--*

    void orderEdges(int which);///<Checks that the edge which is in the correct place in the list (from 0.. which)

    inline int orderNodes(int which, double dist)
    {
	return(orderNodes(which,dist,-1));
    }///<places a node in a new place in the list (updated to account for recedges);
    int orderNodes(int which, double dist, int oldwhich);///<places a node in a new place in the list (updated to account for recedges); a known position is a positive oldwhich...
    void swapNode(int a, int b);///<Swaps two nodes in the list (updated to maintain recedges going to them)
    int lastCommonAncestor(int s1,int s2);///< USES TABNODES AND TABFATHER!  Returns the last common ancestor index of the tabnodes ASSUMING the current local tree.
    void makeLocalTreeKeepingFathers(unsigned int site);///< Makes a version of the local tree that has the fathers but is unoptimised
    double pairwiseDist(int s1,int s2);///< Returns the pairwise distance between two sequences accounting for recombination.  SAME WARNINGS AS lastCommonAncestor!
public:

// Constructors and destructors:
    RecTree(unsigned int numsites,WargXml *infile,bool addedges=true,bool forceages=true,bool isfinal=true,std::streampos itstart=-1);///<Creates a tree from an output file

    RecTree(unsigned int numsites,std::string newick,bool isFilename=true,bool forceages=true);///<Creates a tree from a Newick file
    RecTree(int n,double rho,double delta,std::vector<int> blocks);///<Creates a Coalescent tree with recombinant edges by simulation
    void dropEdges(double rho,double delta,vector<int> blocks);
    RecTree(Data * data,double rho,double delta,vector<int> blocks);///<Creates UPGMA tree
    RecTree(RecTree *intree,string newick,bool isFilename=false,bool forceages=true);/// Copies a rectree
    RecTree(const RecTree& rt);///< Copy constructor
    void assign(const RecTree& rt);
    ~RecTree();
// Part of the simulation initialisation:
    void setBlock(unsigned int* gstart,unsigned int* gend,double delta,std::vector<int> *blocks);///<Sets the start and end positions of a recombinant block based on geometric distribution and independent blocks
    int getEdgeCoal(double* time);///<Returns the node the recombinant edge is from and sets the time
// Part of the likelihood/prior calculation:
    double priorEdge(int e,Param * param) const;///<Returns the LOG prior of a given edge
    double priorEdge(double tFrom,double tTo) const;///<Returns the UNLOGGED prior of a given edge
    double prior(Param * param) const;///<Returns the value of the prior function

    void makeLocalTree(unsigned int site,double thetaPerSite);///<Evaluates the local tree at the given site
// Tests:
    void testEdges() const;///<Tests recedges arrive in between the ages of the nodes they go to
    void testTree();///<Checks ages increase in the tree for all sites.
    inline bool sameLocalTreeAsPrev(unsigned int site) const
    {
        return sameLTasPrev[site];
    }///<Returns true if the local tree of site is the same as that of site-1
    void calcSameLTasPrev(unsigned int site);///<Estimate if the LT is the same as the previous one
// Information functions:
    inline Node * getRoot(){return(root);}
    inline double getEdgeTimeAbsFrom(int e) const
    {
        return edge[e]->getTimeFrom()+nodes[edge[e]->getEdgeFrom()]->getAge();
    }
    ;///<Returns the absolute age of the origin of a recombinant edge
    inline double getEdgeTimeAbsTo  (int e) const
    {
        return edge[e]->getTimeTo  ()+nodes[edge[e]->getEdgeTo  ()]->getAge();
    }
    ;///<Returns the absolute age of the destination of a recombinant edge
    inline int getL() const
    {
        return L;
    }///<Returns the total size of the alignment

    inline RecEdge*getRecEdge(int w) const
    {
        return edge[w];
    }///<Returns the w-th recombinant edge
    inline int getSon (int i, unsigned int site,bool isLeft,double *dist) const
    {
        site=0;
        *dist=tabSonsDist[i][isLeft];
        return tabSons[i][isLeft];
    };
    ///<Returns the left daughter node for a given local site and node (NULL if no left daughter)
    inline int affecting(unsigned int site) const
    {
        int a=0;
        for (unsigned int i=0;i<edge.size();i++)
            if (edge[i]->affectsSite(site))
                a++;
        return a;
    }
    ;///<Returns the number of edges affecting a site
    inline int numRecEdge() const
    {
        return edge.size();
    }///<Returns the number of recombinant edges
    int numRecEdgeOnBranch(int b) const
    {
        int r=0;
        for (unsigned int i=0;i<edge.size();i++)
            if (edge[i]->getEdgeTo()==b)
                r++;
        return r;
    }///<Returns the number of recedges on a given branch
    std::vector<int> alive(double t) const;///<Returns the branches alive at a given time
    bool isThere(int start,int end,int e,double time1,double time2) const;///<Returns yes if there is an edge with departure or arrival on branch e and between time1 and time2, and which affects material that intersects [start;end]
    inline RecEdge* getEdge(int i) const
    {
        return(edge[i]);
    }///<Returns the edge with index i
    double getEdgeTreeTime(int i) const;///<Returns the time an edge spans in the tree
    vector<int> getAffEdges(vector<int> * samplespace);///<Calculates the recedges relevent to a set of clonal edges
    void updateList(int which, int newloc,vector<int> * list);///<Updates a list of edges when which was moved to newloc
    int sampleEdge(vector<int> * samplespace=NULL);///<Samples a random edge that (optionally) interacts with a specific set of clonal edges
// Modification functions:
    void addEdgesFromFile(WargXml *infile,int siteoffset=0);///< adds edges to the rectree from the specified warg XML file (which is set to an iteration). Adds a "siteoffset" to add edges from different runs on different genes
    int addRecEdge(double tfrom,double tto,unsigned int gstart,unsigned int gend,int edgefrom,int edgeto);///<Adds an recombination edge to the tree
    int addRecEdge(unsigned int gstart,unsigned int gend,int edgefrom,int edgeto);///< Adds a recombination edge with random start and end times between given edges
    int addRecEdge(std::string res,int sitesoffset=0);///<Add recedge from our output format
    void remRecEdge(int which);///<Removes a recombinant edge
    void setStart(int edge,int start);///<Sets the starting point of the given recedge
    void setEnd(int edge,int end);///<Sets the ending point of the given recedge
    void changeAge(int which,double dist);///<Changes the "which" node age, maintaining edge times
    inline void changeEdgeNodeFrom(int e,int newnode)
    {
        edge[e]->setTimeFrom(getEdgeTimeAbsFrom(e)-getNode(newnode)->getAge(),newnode);
    }///<changes the node to without changing the absolute time of the edge
    inline void changeEdgeNodeTo(int e,int newnode)
    {
        edge[e]->setTimeTo(getEdgeTimeAbsTo(e)-getNode(newnode)->getAge(),newnode);
    }///<changes the node to without changing the absolute time of the edge
    void chooseNodePath(int which, double dist, vector<int> *listLR, vector<int> *affedges);///< Sets up a specified order for which nodes to swap when they are disordered
    int moveNodeTime(int which, int *whichto,double dist,int oldwhich, vector<int> *listLR,  vector<int> *affedges);///<Updates the node list when the node which is moved a distance dist number of recedges moved to younger nodes - number moved to older nodes. edits affedges to make it the a vector of affected edges

    void moveEdges(int e1,int e2,double t1=0.0,double t2=-1.0);///<Moves all edge events on e1 between absolute times t1..t2 onto edge e2, maintaining their absolute times (t2=-1 means no upper limit)
// Information functions
    std::vector<std::vector<double> > pairwiseDistanceMatrix(long site);///< returns a pairwise distance matrix for a given site

};

} // end namespace weakarg
#endif
