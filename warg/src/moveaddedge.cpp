#include "moveaddedge.h"
//
using namespace std;
namespace weakarg
{

MoveAddEdge::MoveAddEdge(Param *p,double a)
        : Move(p,a)
{
       description= "Update adding a recombinant edge";
       name= "AddEdge";
}


int MoveAddEdge::move(vector<int> * samplespace)
{
    RecTree * rectree=param->getRecTree();
    double tfrom,tto;
    unsigned int start,end;
    unsigned int efrom,eto;
    //Draw start and end
    rectree->setBlock(&start,&end,param->getDelta(),param->getData()->getBlocks());
    //Draw eto and tto
    bool insamplespace=false;
    while(!insamplespace){
	eto=rectree->getPoint(&tto);
	//Draw efrom and tfrom
	tfrom=tto+rectree->getNode(eto)->getAge();
	efrom=rectree->getEdgeCoal(&tfrom);
	tfrom-=rectree->getNode(efrom)->getAge();
	if(samplespace==NULL) insamplespace=true;
	else if (samplespace->size()==0) insamplespace=true;
	else
	{
		for(unsigned int i=0;i<samplespace->size();i++) if((int)eto==samplespace->at(i)||(int)efrom==samplespace->at(i))
		{
		insamplespace=true;break;
		}
	}
    }
    double l=param->getLL();
    vector<double> store(end-start);
    for (unsigned int i=start;i<end;i++)
        store[i-start]=param->getLLsite(i);
    int which=rectree->addRecEdge(tfrom,tto,start,end,efrom,eto);
    if(which<0) return(-1);

    numcalls++;
    param->computeLikelihood(start,end);
    double l2=param->getLL();
    dlog(1)<<"Proposing to add edge "<<efrom<<":"<<tfrom<<"->"<<eto<<":"<<tto<<"...";
    if (log(gsl_rng_uniform(rng))>l2-l+log(param->getRho()*rectree->getTTotal()/2.0/rectree->numRecEdge()))
    {
    	dlog(1)<<"Rejected!"<<endl;
        param->getRecTree()->remRecEdge(which);
        for (unsigned int i=start;i<end;i++)
            param->setlocLL(i,store[i-start]);
        param->setLL(l);
	return(0);
        //param->computeLikelihood(start,end);
    }
    else dlog(1)<<"Accepted!"<<endl;
    numaccept++;
    return(1);
}

MoveAddEdge::~MoveAddEdge()
{}
//

} // end namespace weakarg
