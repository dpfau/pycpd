/*
 fgt_predict : returns the Fast Gauss Transform approximation of the test points y given the model \theta=(xc,A_k)
 
 Usage
 -------
 v            = fgt_predict(y , xc , A_k , [sigma] , [e] );

 Inputs
 -------
 y             Test point (d x Ny) 
 xc            Kcenter point (d x K)
 A_k           Polynomial coefficient (pd x K), where pd = nchoosek(p + d - 1 , d)
 sigma         Noise Standard deviation of the kernel (default sigma = 1)
 e             Ratio of far field (default e = 10)

 Ouputs
 -------
 v             Density (1 x Ny)

 To compile
 ----------
mex -output fgt_predict.dll fgt_predict.c
mex  -f mexopts_intelAMD.bat -output fgt_predict.dll fgt_predict.c

Example 1 
---------
d          = 3;
Nx         = 300;
Ny         = 1000000;
x          = randn(d , Nx);
w          = ones(1 , Nx);
sigma      = 1;
p          = 7;
K          = round(sqrt(Nx));
e          = 6;
y          = randn(d , Ny);

[xc , A_k] = fgt_model(x , w , sigma , p , K , e);

v          = fgt_predict(y , xc , A_k , sigma , e);

vtrue      = dval(x , y , w);


Example 2 
---------

d       = 3;
Nx      = 10;
Ny      = 100;
x       = randn(d , Nx);
y       = randn(d , Ny);
w       = rand(1 , Nx);
h       = 2;

p       = 6;
K       = 5;
e       = 10;

v1       = dval(x , y , w , h);
v2       = fastgausstransform(x , y , w , h , p , K , e);

 
Example 3 
---------
d          = 2;
R          = [2 , 0.4 ; 0.4  3];
Nx         = 100;
h          = 1;

p          = 6;
K          = 15;
e          = 10;



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

int invnchoosek(int , int );
void fgt_predict(double * , double * , double * , int , double  , int , double , int , int , 
           double * , 
           double * , double *  , int *);
