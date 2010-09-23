#include "movetimechange.h"

#define TIMESTEPSIZE 0.1
// the step size (sigma in a normal distribution) of

using namespace std;
namespace weakarg
{

MoveTimeChange::MoveTimeChange(Param * p,double a)
        : Move(p,a)
{
description="Update moving the time of a recedge";
name= "TimeChange";
}

int MoveTimeChange::move(vector<int> * samplespace)
{
    RecTree * rectree=param->getRecTree();
    int which=rectree->sampleEdge(samplespace);//floor(gsl_rng_uniform(rng)*rectree->numRecEdge());
    if (which<0)
    {
    	dlog(1)<<"No valid edges to change!"<<endl;
        return(-1);
    }
    numcalls++;
    int movetto=gsl_rng_uniform_int(rng,2);//0 for no or 1 for yes
    // collect the details of the edge
    int efrom=rectree->getRecEdge(which)->getEdgeFrom();
    int eto  =rectree->getRecEdge(which)->getEdgeTo  ();
    double tfrom=rectree->getEdgeTimeAbsFrom(which);
    double tto  =rectree->getEdgeTimeAbsTo  (which);
    double tmptime,tfromnew=tfrom,ttonew=tto;
    int start=rectree->getRecEdge(which)->getStart();
    int end  =rectree->getRecEdge(which)->getEnd();
    int tmpedge,etonew=eto,efromnew=efrom;
    double l=param->getLL();
    double lpriorrat=0;// log of (prior ratio * transition ratio)
    double rootage=rectree->getNode(rectree->getN()*2-2)->getAge();
    vector<double> store(end-start);
    for (int i=start;i<end;i++)
        store[i-start]=param->getLLsite(i);
    dlog(1)<<"Proposing to change edge time "<<efrom<<":"<<tfrom<<"->"<<eto<<":"<<tto<<" to"<<flush;

    // change the edge
    if(movetto)
    {// moving the arrival time
        ttonew=tto+gsl_ran_gaussian(rng,TIMESTEPSIZE);
        while(ttonew<0 || ttonew>min(tfrom,rootage))
        {
            if(ttonew<0)
                ttonew=-ttonew;
            if(ttonew>min(tfrom,rootage))
                ttonew=2.0*min(tfrom,rootage)-ttonew;
        }// while loop as can bounce off reflecting boundaries several times
        tmptime=ttonew;
        tmpedge=etonew;
    }
    else
    {// moving the departure time
        tfromnew=tfrom+gsl_ran_gaussian(rng,TIMESTEPSIZE);
        if(tfromnew<tto)
            tfromnew=2.0*tto-tfromnew;
        tmptime=tfromnew;
        tmpedge=efromnew;
    }
    // update the edge index
    while(tmptime<rectree->getNode(tmpedge)->getAge())
    {
        if(gsl_rng_uniform(rng)<0.5)
            tmpedge=rectree->getNode(tmpedge)->getLeft()->getId();
        else
            tmpedge=rectree->getNode(tmpedge)->getRight()->getId();
        lpriorrat+=log(2.0);
    }
    if(tmpedge!=rectree->getN()*2-2)
    {
        while(tmptime>rectree->getNode(tmpedge)->getFather()->getAge())
        {
            tmpedge=rectree->getNode(tmpedge)->getFather()->getId();
            lpriorrat-=log(2.0);
            if(tmpedge==rectree->getN()*2-2)
                break;
        }
    }
    if(movetto)
    {
        ttonew=tmptime;
        etonew=tmpedge;
    }
    else
    {
        tfromnew=tmptime;
        efromnew=tmpedge;
    }
    dlog(1)<<" "<<efromnew<<":"<<tfromnew<<"->"<<etonew<<":"<<ttonew<<"..."<<flush;
    // remove the old edge, add the new, and compute the new likelihood
    lpriorrat-=rectree->priorEdge(which,param);
    rectree->remRecEdge(which);
    bool opti;//Optimization: check that a local tree is changed, otherwise l2=l
    if (tmpedge==rectree->getN()*2-2||(!movetto&&tmpedge!=efrom)||(movetto&&tmpedge!=eto))
        opti=false;
    else if (movetto)
        opti=!rectree->isThere(start,end,tmpedge,tmptime-rectree->getNode(tmpedge)->getAge(),tto  -rectree->getNode(eto  )->getAge());
    else
        opti=false;//opti=!rectree->isThere(start,end,tmpedge,tmptime-rectree->getNode(tmpedge)->getAge(),tfrom-rectree->getNode(efrom)->getAge());
    which=rectree->addRecEdge(tfromnew-rectree->getNode(efromnew)->getAge(),ttonew-rectree->getNode(etonew)->getAge(),start,end,efromnew,etonew);
    if(which<0) throw("Movetimechange: Can't create  edge!");
    double l2=l;
    if (!opti)
    {
        param->computeLikelihood(start,end);
        l2=param->getLL();
    }
    lpriorrat+=rectree->priorEdge(which,param);
    // acceptance step
    if (log(gsl_rng_uniform(rng))>l2-l+lpriorrat)
    {
    	dlog(1)<<"Rejected!"<<endl;
        rectree->remRecEdge(which);
        if(rectree->addRecEdge(tfrom-rectree->getNode(efrom)->getAge(),tto-rectree->getNode(eto)->getAge(),start,end,efrom,eto)<0) throw("Movetimechange: Can't restore edge!");
        for (int i=start;i<end;i++)
            param->setlocLL(i,store[i-start]);
        param->setLL(l);
        //param->computeLikelihood(start,end);
        return(0);
    }
    else dlog(1)<<"Accepted!"<<endl;
    numaccept++;
    return(1);
}




MoveTimeChange::~MoveTimeChange()
{}


} // end namespace weakarg
