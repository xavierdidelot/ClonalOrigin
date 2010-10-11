#ifndef PARAM_H
#define PARAM_H
#include "data.h"
#include "rectree.h"
#include <vector>
#define LOG025 -1.38629461
#include "weakarg.h"
#include "mpiutils.h"
#include <math.h>
#include "wargxml.h"
#include <iomanip>
/**
    @brief Set of parameters of the model
*/

namespace weakarg
{

class Move;

class Param
{

public:
    Param();
    Param(RecTree* rectree,Data* data);
    virtual ~Param();
    void readParamsFromFile(WargXml * infile,std::streampos itstart=-1);///<Reads the parameter values from an output file
    void readProgramOptions();///<Sets MCMC parameters based on command-line options
    void metropolis(string comment=string(""));///<Performs the Metropolis-Hastings algorithm
    void simulateData(std::vector<int> blocks);///<Simulates data given the weak ARG, using the data structure provided
    inline RecTree * getTree()
    {
        return rectree;
    }///<Returns the tree
    inline void setTree(RecTree * t)
    {
        if(rectree!=NULL) delete(rectree);
        rectree=t;
    }///<Sets the tree
    inline Data * getData()
    {
        return data;
    }///<Returns the Data
    inline void setData(Data*bd)
    {
        if(data!=NULL) delete(data);
        data=bd;
    }///<Sets the data
    void exportXMLbegin(std::ostream& out,string comment=string(""));///<Creates the beginning of an output file
    void exportXMLend(std::ostream& out);///<Creates the end of an output file
    void exportXMLiter(std::ostream& out, long i=-1);///<Creates an iteration in the output file
    void computeLikelihood();///<Calculates the current likelihood
    void computeLikelihood(unsigned int start, unsigned int end);///<Calculates the current likelihood for a range of sites
    double getRlenPrior(int edge);///<Returns the log of the prior due to genetic position for a recedge
    inline double getTheta() const
    {
        return theta;
    }///<Returns the value of theta
    inline void setTheta(double t)
    {
        theta=t;
    }///<Sets the value of theta
    inline double getRho() const
    {
        return rho;
    }///<Returns the value of rho
    inline void setRho(double r)
    {
        rho=r;
    }///<Sets the value of rho
    inline double getDelta() const
    {
        return delta;
    }///<Returns the value of delta
    inline void setDelta(double d)
    {
        delta=d;
    }///<Sets the value of delta
    inline double getLL() const
    {
        return ll/tempering;
    }///<Returns the value of the likelihood
    inline void setLL(double l)
    {
        ll=l*tempering;
    }///<Sets the value of the likelihood
    inline double getLLsite(unsigned int site) const
    {
        return locll[site]/tempering;
    }///<Get the likelihood of a given site
    inline void setlocLL(unsigned int i,double l)
    {
        locll[i]=l*tempering;
    }///<Sets the value of the local likelihood
    inline double getPrior()
    {
	return(rectree->prior(this));
    }///<Returns the value of the prior function

    inline double getClonalStepSize() const
    {
        return moveClonalStepSize;
    }///<Gets the value of the Clonal MCMC step
    inline void setClonalStepSize(double d)
    {
        moveClonalStepSize=d;
    }///<sets the value of the Clonal MCMC step
    inline double getScaleTreeSize() const
    {
        return moveScaleTreeSize;
    }///<Gets the value of the tree scaling MCMC step
    inline void setScaleTreeSize(double d)
    {
        moveScaleTreeSize=d;
    }///<sets the value of the tree scaling MCMC step
    inline RecTree * getRecTree()
    {
        return rectree;
    }
    inline void setRecTree(RecTree *newtree,bool quick=true)
    {
        rectree=newtree;
	if(!quick) computeLikelihood();
    }
    inline void replaceRecTree(RecTree* newrectree,bool quick=true)
    {
	if(rectree!=NULL) delete(rectree);
	rectree=newrectree;
	if(!quick) computeLikelihood();
    }
    inline double logPriorOfRho() {
	if (hRho==0.0) return(0.0);
	return(log(hRho)-hRho*rho);
    }// it is exponential(hRho), ie. uniform if hRho=0
    inline double logPriorOfTheta() {
	return(0);
    }// it is uniform (0,infinity)
    inline double hyperPriorOfRho(){return(hRho);}///<The hyperprior for rho p(rho)=hRho*exp(-rho*hRho), i.e. it is a rate parameter
    int chooseMove(std::vector<Move*> *mall,double totweight);///<Chooses moves proportional to their weights
    void startDiagnostics(std::ostream& out);///<Header for diagnostics
    void updateDiagnostics(int iter,std::ostream& out);///< Updates diagnostic variables
    void testTree();///<Tests the ll vector for consistency and tree structure
    inline void setTempering(double t){tempering=t;}///< Sets the tempering value
    inline double getTempering(){return(tempering);}///< Sets the tempering value
    void greedyTreeMove();///< Performs a non-MCMC "greedy" update of the tree
    vector<double> greedyCalcDists();///< calculates the distances for a greedy move
    vector<double> greedyCalcDists(vector<double> nmuts,vector<double> nsites);///< calculates the distances on the tree given some mutational counts
    void greedyApply(vector<double> dists);///< applies the distances to the tree
    vector< vector<double> > * greedyPairwiseDetails();///< Gets the pairwise distances
    vector< vector<double> > greedyDetails(vector< vector<double> > * pairwise);///< gets all the gory details for r/m calcs (provided by the clonalframe/non-cf details in movegreedytree)
    inline double empiricalDelta(){
	vector<int> * b=data->getBlocks();
	int nsites=0, nedges=0;
	for (int i=0;i<rectree->numRecEdge();i++){
	  unsigned j=0;
	  int ss=1,sf=1;
	  while(j<b->size()){
		if((unsigned)b->at(j)==rectree->getEdge(i)->getStart()) ss=0;
		if((unsigned)b->at(j)==rectree->getEdge(i)->getEnd()) sf=0;
 		if((unsigned)b->at(j)>=rectree->getEdge(i)->getEnd())j=b->size();
		j++;
	  }
	  nedges+=ss+sf;
	  nsites+=rectree->getEdge(i)->getEnd() - rectree->getEdge(i)->getStart();
	}
	return(2.0*nsites/(nedges+1));
    }
    inline double empiricalRho(){
	rectree->computeTTotal();
	return(2.0*rectree->numRecEdge()/rectree->getTTotal());
    }
    double empiricalTheta(vector< vector<double> > * mutpairwise);
protected:
	double tempering;///< Tempering of the likelihood
	std::vector<Move*> mall;///< Vector of moves
    std::vector<std::vector<double> > f;

    double hRho;///<Hyperprior for rho
    double theta;///<Mutation rate
    double rho;///<Recombination rate
    double delta;///<Mean recombination tract length
    double ll;///<Current value of the log-likelihood
    std::vector<double> locll; ///<Current values of the local log-likelihood
    RecTree * rectree;///<Tree on which the parameters apply
    Data * data;///<Data and ancestral states of the ancestral nodes
    void computeSiteLL(unsigned int site,bool makeLT);///<Computes the Log Likelihood of a given site
    // diagnostics:
    double esttheta;///<Running estimate of theta, for diagnostics
    double estthetasq;///<Running estimate of theta squared, for diagnostics
    double estrho;///<Running estimate of rho, for diagnostics
    double estrhosq;///<Running estimate of rho squared, for diagnostics
    double estdelta;///<Running estimate of delta, for diagnostics
    double estdeltasq;///<Running estimate of delta squared, for diagnostics
    double estnumrecedge;///<Running estimate of number of recedges, for diagnostics
    double estnumrecedgesq;///<Running estimate of number of recedges squared, for diagnostics
    double estedgeden;///<Running estimate of edge density, for diagnostics
    double estedgedensq;///<Running estimate of edge density squared, for diagnostics
    double estedgepb;///<Running estimate of edges per branch, for diagnostics
    double estedgepbsq;///<Running estimate of edges per branch squared, for diagnostics
    double estvaredgepb;///<Running estimate of variance in the edges per branch, for diagnostics
    double estvaredgepbsq;///<Running estimate of variance in the edges per branch, squared, for diagnostics

    double moveClonalStepSize;///< The step size for updating the clonal node ages
    double moveScaleTreeSize;///< The step size for scaling the whole tree
};


} // end namespace weakarg
#endif
