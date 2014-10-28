#include "fgt_predict.h"
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{	
	double  *y , *xc , *A_k;
	double sigma = 1.0 , e = 10.0;

	double *v;
	int d , pd , K , Ny;
	const int  *dimsxc , *dimsA_k , *dimsy;
		
	int numdimsxc , numdimsA_k , numdimsy ;
	int i , p ;


	double  *dx , *prods;
	int *heads;
		
	/* -------------------------- Parse INPUT  -------------------------------------- */
	if ((nrhs < 3))
	{
		mexErrMsgTxt("Usage : v = fgt_predict(y , xc , A_k, [sigma] , [e] );"); 	
	}


	/* ----- Input 1 ----- */
	y           = mxGetPr(prhs[0]);
  numdimsy    = mxGetNumberOfDimensions(prhs[0]);
	dimsy       = mxGetDimensions(prhs[0]);
	d           = dimsy[0];
	Ny          = dimsy[1];


	/* ----- Input 2 ----- */
	xc           = mxGetPr(prhs[1]);
	numdimsxc    = mxGetNumberOfDimensions(prhs[1]);
	dimsxc       = mxGetDimensions(prhs[1]);
	
	K            = dimsxc[1];

	if (dimsxc[0] != d )
	{
		mexErrMsgTxt("xc must be (d x K)"); 	
	}
	
	/* ----- Input 3 ----- */
	A_k           = mxGetPr(prhs[2]);  
  numdimsA_k    = mxGetNumberOfDimensions(prhs[2]);
	dimsA_k       = mxGetDimensions(prhs[2]);
	pd            = dimsA_k[0];

	if (dimsA_k[1] != K )
	{
		mexErrMsgTxt("A_k must be (pd , K) where pd = nchoosek(p + d - 1 , d)"); 	
	}

	/* ----- Input 4 ----- */	
	if (nrhs > 3)
	{
		sigma = (double)mxGetScalar(prhs[3]);	
	}

	/* ----- Input 5 ----- */
	if (nrhs > 4)
	{
		e = (double)mxGetScalar(prhs[4]);	
	}

	/* -------------------------- Parse OUTPUT  ------------------------------------- */
	
	/* ----- output 1 ----- */
	plhs[0]        = mxCreateDoubleMatrix(1 , Ny , mxREAL); 
	v              = mxGetPr(plhs[0]);

  dx             = (double *)mxMalloc(d*sizeof(double));

  prods          = (double *)mxMalloc(pd*sizeof(double));

  heads          = (int *)mxMalloc((d + 1)*sizeof(int));

	/* ----------------------- MAIN CALL  -------------------------------------------- */	

	fgt_predict(y , xc , A_k , Ny , sigma  , K , e , d , pd , 
		        v , 
			    dx , prods , heads);
	
	/* ------------ END of Mex File ---------------- */

	mxFree(heads);
	mxFree(dx);
	mxFree(prods);
}

void fgt_predict(double *y , double *xc , double *A_k  , int Ny, double sigma , int K , double e , int d , int pd ,  
			     double *v , 
			     double *dy , double *prods , int *heads)
{
	int p , i , j , m , k , t , tail , mbase , kn , xbase , head , ind;
	double sum2 , ctesigma = 1.0/(sigma) , temp , temp1;
	p              = invnchoosek(d , pd);
	for (m=0 ; m < Ny ; m++)
	{	
		temp    = 0.0;
		mbase   = m*d;
		for (kn = 0 ; kn < K ; kn++)
		{
			xbase = kn*d;
			ind   = kn*pd;
			sum2  = 0.0;
			for (i = 0 ; i < d ; i++)
			{
				dy[i]    = (y[i + mbase] - xc[i + xbase])*ctesigma;
				sum2    += dy[i] * dy[i];
				heads[i] = 0;
			}
			if (sum2 > e) continue; /* skip to next kn */
			prods[0] = exp(-sum2);		
			for (k=1, t=1, tail=1 ; k < p ; k++ , tail=t)
			{
				for (i = 0 ; i < d; i++)
				{
					head     = heads[i];
					heads[i] = t;
					temp1    = dy[i];
					for (j = head ; j < tail ; j++ , t++)
					{
						prods[t] = temp1 * prods[j];
					}
				} 
			}
			
			for (i = 0 ; i < pd ; i++)
			{
				temp += A_k[i + ind]*prods[i];
			}
		}
		v[m] = temp;
	}
}	

int invnchoosek(int d , int cnk)
{
	int i , j , cted=1 , ctep , cte , p  ;	
	for(i = 2 ; i <= d ; i++)
	{
		cted *=i;
	}
	
	cte  = cnk*cted;
	p    = 2;
	ctep = p;
	
	for (i = p + 1 ; i < p + d ; i++)
	{
		ctep *=i ; 	
	}
	
	while(ctep != cte)
	{
		ctep = ((p+d)*ctep)/p;	
		p++;		
	}
	
	return p;
}




