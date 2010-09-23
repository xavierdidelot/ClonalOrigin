#include "move.h"

using namespace std;
namespace weakarg
{

Move::Move(Param * p,double a)
{
    param=p;
    alpha=a;
    numcalls=0;
    numaccept=0;
}

Move::~Move()
{}

} // end namespace weakarg
