#include "metropolis.h"

#define DEBUG

using namespace std;
namespace weakarg
{

Metropolis::Metropolis(Param * p)
{
	param=p;
}

inline double fround(double n, double d)
{
    return floor(n * pow(10., d) + .5) / pow(10., d);
}

void Metropolis::move(int n, double temper, vector<int> * samplespace)
{
	std::vector<Move*> mall;
    mall.push_back(new MoveEdgeChange(param,1.0));
    mall.push_back(new MoveRemEdge(param,1.0));
    mall.push_back(new MoveAddEdge(param,1.0));
    mall.push_back(new MoveSiteChange(param,1.0));
    mall.push_back(new MoveTimeChange(param,1.0));
	double totweight=5.0;



	dlog(1)<<endl;
	double oldtemper=param->getTempering();
	param->setTempering(temper);
	for(int i=0;i<n;i++) {
		int move=chooseMove(&mall,totweight);
		dlog(1)<<"Doing tempered move "<<i<<" of type "<<mall[move]->getName()<<endl;
		try{mall[move]->move(samplespace);
	    	}catch(char * x){ cout<<"Error in Param: "<<x<<endl;}

	}
	param->setTempering(oldtemper);
	for(unsigned int i=0;i<mall.size();i++) if(mall[i]!=NULL) delete(mall[i]);
}

int Metropolis::chooseMove(vector<Move*> *mall,double totweight)
{
    double r=gsl_rng_uniform(rng);
    for(unsigned int i=0;i<mall->size();i++)
    {
        r-=(*mall)[i]->getAlpha()/totweight;
        if(r<0)
            return(i);
    }
    cout<<"Error in choosing move probabilities"<<endl;
    throw;
}

Metropolis::~Metropolis()
{
}



} // end namespace weakarg
