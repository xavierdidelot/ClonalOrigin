#include "moveremedge.h"
//

using namespace std;
namespace weakarg
{

MoveRemEdge::MoveRemEdge(Param * p,double a)
        : Move(p,a)
{
    description= "Update removing a recombinant edge";
    name= "RemEdge";

}


int MoveRemEdge::move(vector<int> * samplespace)
{
    RecTree * rectree=param->getRecTree();
    int which=rectree->sampleEdge(samplespace);//floor(gsl_rng_uniform(rng)*rectree->numRecEdge());
    if (which<0)
    {
    	dlog(1)<<"No valid edges to change!"<<endl;
        return(-1);
    }
    numcalls++;
    double l=param->getLL();
    double tfrom=param->getRecTree()->getRecEdge(which)->getTimeFrom();
    double tto  =param->getRecTree()->getRecEdge(which)->getTimeTo  ();
    int start=param->getRecTree()->getRecEdge(which)->getStart();
    int end  =param->getRecTree()->getRecEdge(which)->getEnd();
    int efrom=param->getRecTree()->getRecEdge(which)->getEdgeFrom();
    int eto  =param->getRecTree()->getRecEdge(which)->getEdgeTo  ();
    rectree->remRecEdge(which);
    vector<double> store(end-start);
    for (int i=start;i<end;i++)
        store[i-start]=param->getLLsite(i);
    param->computeLikelihood(start,end);
    double l2=param->getLL();
    dlog(1)<<"Proposing to remove edge "<<which<<"...";
    if (log(gsl_rng_uniform(rng))>l2-l+log((1.0+rectree->numRecEdge())*2.0/param->getRho()/rectree->getTTotal()))
    {
    	dlog(1)<<"Rejected!"<<endl;
	if(param->getRecTree()->addRecEdge(tfrom,tto,start,end,efrom,eto)<0) throw("MoveRemEdge: Can't restore edge!");
        for (int i=start;i<end;i++)
            param->setlocLL(i,store[i-start]);
        param->setLL(l);
	return(0);
    }
    else dlog(1)<<"Accepted!"<<endl;
    numaccept++;
    return(1);
}

MoveRemEdge::~MoveRemEdge()
{}
//

} // end namespace weakarg
