#include "moveedgechange.h"

using namespace std;
namespace weakarg
{

// the step size (sigma in a normal distribution) of

MoveEdgeChange::MoveEdgeChange(Param * p,double a)
        : Move(p,a)
{
    description="Update moving the edge of a recedge";
    name="EdgeChange";
}

int MoveEdgeChange::move(vector<int> * samplespace)
{
    RecTree * rectree=param->getRecTree();
    dlog(1)<<"Proposing to change edge sites... ";
    // first choose the move type
    int which=rectree->sampleEdge(samplespace);//floor(gsl_rng_uniform(rng)*rectree->numRecEdge());
    if (which<0)
    {
    	dlog(1)<<"No valid edges to change!"<<endl;
        return(-1);
    }
    numcalls++;
    int movetto=gsl_rng_uniform_int(rng,2); //0 for no or 1 for yes
    // collect the details of the edge
    double l=param->getLL();
    int oldedge;
    double time;
    if (movetto)
    {
        oldedge=rectree->getRecEdge(which)->getEdgeTo  ();
        time=rectree->getEdgeTimeAbsTo  (which);
    }
    else
    {
        oldedge=rectree->getRecEdge(which)->getEdgeFrom();
        time=rectree->getEdgeTimeAbsFrom(which);
    }
    unsigned int start=rectree->getRecEdge(which)->getStart();
    unsigned int end  =rectree->getRecEdge(which)->getEnd();
    dlog(1)<<"edge "<<which<<" ("<<start<<":"<<end<<")"<<flush;
    vector<double> store(end-start);
    for (unsigned int i=start;i<end;i++)
        store[i-start]=param->getLLsite(i);
    // change the edge
    int tmpedge=oldedge;
    vector<int> v=rectree->alive(time);
    if (v.size()==0)
        throw;
    if (v.size()==1)
    {
    	dlog(1)<<"Cannot change chosen edge!"<<endl;
        return(-1);
    }
    while (tmpedge==oldedge)
        tmpedge=v[gsl_rng_uniform_int(rng,v.size())];
    if (movetto)
        rectree->getRecEdge(which)->setTimeTo(time-rectree->getNode(tmpedge)->getAge(),tmpedge);
    else
        rectree->getRecEdge(which)->setTimeFrom(time-rectree->getNode(tmpedge)->getAge(),tmpedge);
    param->computeLikelihood(start,end);
    double l2=param->getLL();
    // acceptance step
    if (log(gsl_rng_uniform_pos(rng))>l2-l)
    {
    	dlog(1)<<"Rejected!"<<endl;
        if (movetto)
            rectree->getRecEdge(which)->setTimeTo(time-rectree->getNode(oldedge)->getAge(),oldedge);
        else
            rectree->getRecEdge(which)->setTimeFrom(time-rectree->getNode(oldedge)->getAge(),oldedge);
        for (unsigned int i=start;i<end;i++)
            param->setlocLL(i,store[i-start]);
        param->setLL(l);
        //param->computeLikelihood(start,end);
        return(0);
    }
    else dlog(1)<<"Accepted!"<<endl;
    numaccept++;
    return(1);
}




MoveEdgeChange::~MoveEdgeChange()
{}

} // end namespace weakarg
