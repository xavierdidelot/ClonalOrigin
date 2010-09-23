#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

namespace weakarg
{

#define EXP_A (1048576/M_LN2)
#define EXP_C 60801
inline double exponential(double y)
{
    union
    {
        double d;
#ifdef LITTLE_ENDIAN
        struct
        {
            int j, i;
        }
        n;
#else
        struct
        {
            int i, j;
        }
        n;
#endif
    }
    _eco;
    _eco.n.i = (int)(EXP_A*(y)) + (1072693248 - EXP_C);
    _eco.n.j = 0;
    return _eco.d;
}

} // end namespace weakarg
#endif
