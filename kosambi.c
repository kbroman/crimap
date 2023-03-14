#if vms
#include stdio
#else
#include <stdio.h>
#include <math.h>
#endif

#include "defs.h"

double kosambi(t)
     double t;
{
     double log();
     double x;

     if(t == .5) return(100);

     x = 25*log( (1+2*t)/(1-2*t) );

     return(x>100?100:x);

}
