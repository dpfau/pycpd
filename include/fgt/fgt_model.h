/*  
  Author : S�bastien PARIS : sebastien.paris@lsis.org
  Adapted to Python by David PFAU : david.pfau@gmail.com
*/

#include <stdlib.h>
#include <math.h>
#include <time.h>
#define min(a , b) ((a) <= (b) ? (a) : (b))
#define max(a , b) ((a) >= (b) ? (a) : (b))

int nchoosek(int  , int );

int idmax(double * , int );

double ddist(double * , double * , int );

void Kcenter(double * , int  , int  , int , 
       double * , int * , int * , int * , 
       double *);

void Compute_C_k(int  , int , double * , int * , int *);

void Compute_A_k(double * , double * , double *, double *  , double  , int , int  , int , int , int , 
         double * , 
         int * , double * , double * , int * );

void fgt_model(double * , double * , double  , int  , int , double ,
         double * , double * , 
         int  , int  ,
         int * , int * , int * , int * , 
         double * , double * , int * , int * , double * , double * , 
         int );
