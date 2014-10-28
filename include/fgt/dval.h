/*
  Evaluate weighted gaussian RBF between vectors x,y
 
  Author : Sébastien PARIS : sebastien.paris@lsis.org
  Modified for Python by David PFAU : david.pfau@gmail.com
  -------
*/

#include <math.h>
#ifdef OMP
 #include <omp.h>
#endif

#ifndef MAX_THREADS
#define MAX_THREADS 64
#endif
#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif

void c_dval(double * , double * , double * , double  , double * , int  , int  , int );