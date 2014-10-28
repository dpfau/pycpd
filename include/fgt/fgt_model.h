/*
 fgt_model : returns the Fast Gauss Transform Aprroximation Model of a Kernel density

 Usage
 -------
 [xc , Ak]    = fgt_model(x , [w] , [sigma] , [e] , [K] , [p]  );

 Inputs
 -------
 x             Source data (d x Nx)
 w             Weigths (1 x Nx) ( default w = ones(1 , Nx) ) 
 sigma         Noise Standard deviation of the kernel (default sigma = 1)
 e             Ratio of far field (default e = 10)
 K             Number of centers (default K = sqrt(Nx))
 p             Order of truncation (default p = 8)

 Ouputs
 -------
 xc            The K center points of the training set (d x K) 
 Ak            Polynomial coefficient (pd x K), where pd = nchoosek(p + d - 1 , d) = prod(1:p + d - 1)/(prod(1:p - 1)*prod(1:d))

 To compile
 ----------
mex -output fgt_model.dll fgt_model.c
mex  -f mexopts_intelAMD.bat -output fgt_model.dll fgt_model.c

Example 1 
---------
d          = 3;
Nx         = 100;
x          = randn(d , Nx);

[xc , A_k] = fgt_model(x);

Example 2 
---------
d          = 3;
Nx         = 10;
Ny         = 100;
x          = randn(d , Nx);
y          = randn(d , Ny);
w          = rand(1 , Nx);
h          = 2;

e          = 10;
p          = 6;
K          = 5;

v1         = dval(x , y , w , h);

[xc , A_k] = fgt_model(x , w , h , e , K , p);

v2         = fgt_predict(y , xc , A_k , h , e); 

 
Example 3 
---------
d          = 2;
R          = [2 , 0.4 ; 0.4  3];
Nx         = 1000;

h          = 1;
e          = 10;
K          = round(sqrt(Nx));
p          = 6;



vect       = (-5:0.3:5);
Ny         = length(vect);
w          = (1/Nx)*ones(1 , Nx);
  
x          = (chol(R)'*randn(d , Nx));

[X , Y]    = meshgrid(vect);
y          = [X(:) , Y(:)]';

[xc , A_k] = fgt_model(x , w , h , e , K , p);
vy         = fgt_predict(y , xc , A_k , h , e , K , p);

densite    = reshape( vy , Ny , Ny);

figure
set(gcf , 'renderer' , 'opengl');
surfc(X , Y , densite)
shading interp
lighting phong

light
alpha(0.5);
hold on
plot(x(1 , :) , x(2 , :) , 'r+' , xc(1 , :) , xc(2 , :) , 'ko' , 'markersize' , 10);
hold off
view(2)
colorbar


  
Author : S�bastien PARIS : sebastien.paris@lsis.org
-------
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