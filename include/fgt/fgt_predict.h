/*
  Fast Gauss Transform approximation of the test points y given the model \theta=(xc,A_k)
  Author : S�bastien PARIS : sebastien.paris@lsis.org
  Modified for Python by David PFAU : david.pfau@gmail.com
*/

#include <stdlib.h>
#include <math.h>
#include <time.h>
#define min(a , b) ((a) <= (b) ? (a) : (b))
#define max(a , b) ((a) >= (b) ? (a) : (b))

int invnchoosek(int , int );
void fgt_predict(double * , double * , double * , int , double  , int , double , int , int , 
           double * , 
           double * , double *  , int *);
