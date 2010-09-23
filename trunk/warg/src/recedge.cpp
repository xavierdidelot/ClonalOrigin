#include "recedge.h"

using namespace std;
namespace weakarg
{

RecEdge::RecEdge(double tfrom,double tto,unsigned int gstart,unsigned int gend,unsigned int edgefrom,unsigned int edgeto)
{
    MoveEdge(tfrom,tto,gstart,gend,edgefrom,edgeto);
}

void RecEdge::MoveEdge(double tfrom,double tto,unsigned int gstart,unsigned int gend,unsigned int edgefrom,unsigned int edgeto)
{
    timeFrom=tfrom;
    timeTo=tto;
    this->gstart=gstart;
    this->gend=gend;
    edgeFrom=edgefrom;
    edgeTo=edgeto;
}

RecEdge::~RecEdge()
{}

} // end namespace weakarg
