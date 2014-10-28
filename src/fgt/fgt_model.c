#include "fgt_model.h"
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
	double *x , *w;
	double sigma = 1.0 , e = 10.0;
	int p = 8 , K = 30;
	double *xc , *A_k;
	const int  *dimsx ;
	int numdimsx  ;
	int i , d , Nx ;
	int pd;
	double *dist_C , *C_k , *dx , *prods;
	int *indxc , *indx , *xhead , *xboxsz , *heads , *cinds;
	
	/* -------------------------- Parse INPUT  -------------------------------------- */
	
	if ((nrhs < 1))
	{
		mexErrMsgTxt("Usage : v = fgt_model(x  , [w] , [sigma] , [p] , [K] , [e] );"); 
	}
	
	/* ----- Input 1 ----- */
	x           = mxGetPr(prhs[0]);
	numdimsx    = mxGetNumberOfDimensions(prhs[0]);
	dimsx       = mxGetDimensions(prhs[0]);
	
	d           = dimsx[0];
	Nx          = dimsx[1];
	
  K           = min(Nx , K);

	/* ----- Input 2 ----- */
	if ((nrhs < 2) || mxIsEmpty(prhs[1]))
	{
		w            = (double *)mxMalloc(Nx*sizeof(double));
		for (i = 0 ; i < Nx ; i++)
		{
			w[i] = 1.0;
		}
	}
	else
	{
		w           = mxGetPr(prhs[1]);
	}
	
	/* ----- Input 3 ----- */
	if (nrhs > 2)
	{
		sigma       = (double)mxGetScalar(prhs[2]);
	}
	
	/* ----- Input 4 ----- */
	if (nrhs > 3)
	{
		e       = (double)mxGetScalar(prhs[3]);
	}
	
	/* ----- Input 5 ----- */
	if (nrhs > 4)
	{
		K         = (int)mxGetScalar(prhs[4]);
		if(K > Nx)
		{
			mexErrMsgTxt("K must be <= Nx"); 
		}
	}
	else
	{
		K = (int) sqrt(Nx);
	}
	
	/* ----- Input 6 ----- */
	if (nrhs > 5)
	{
		p          = (int)mxGetScalar(prhs[5]);
	}
	
	/* -------------------------- Parse OUTPUT  ------------------------------------- */

	pd             = nchoosek(p + d - 1 , d); 

	/* ----- output 1 ----- */
	plhs[0]        = mxCreateDoubleMatrix(d , K , mxREAL); 
	xc             = mxGetPr(plhs[0]);
	
	plhs[1]        = mxCreateDoubleMatrix(pd , K , mxREAL); 
	A_k            = mxGetPr(plhs[1]);
	
  C_k            = (double *)mxMalloc(pd*sizeof(double));
  dist_C         = (double *)mxMalloc(Nx*sizeof(double));
  dx             = (double *)mxMalloc(d*sizeof(double));
  prods          = (double *)mxMalloc(pd*sizeof(double));
	
  indxc          = (int *)mxMalloc(K*sizeof(int));
  indx           = (int *)mxMalloc(Nx*sizeof(int));
  xhead          = (int *)mxMalloc(K*sizeof(int));
  xboxsz         = (int *)mxMalloc(K*sizeof(int));
  heads          = (int *)mxMalloc((d + 1)*sizeof(int));
  cinds          = (int *)mxMalloc(pd*sizeof(int));
	
	/* ----------------------- MAIN CALL  -------------------------------------------- */

	fgt_model(x ,  w , sigma , p , K , e , 
		      xc , A_k , 
			  d , Nx ,
			  indxc , indx , xhead , xboxsz , 
			  dist_C , C_k , heads , cinds , dx , prods ,  
			  pd);

	/* ------------ END of Mex File ---------------- */
	
	mxFree(indxc);
	mxFree(indx);
	mxFree(xhead);
	mxFree(xboxsz);
	mxFree(dist_C);
	mxFree(C_k);
	mxFree(heads);
	mxFree(cinds);
	mxFree(dx);
	mxFree(prods);
	if ((nrhs < 3) || mxIsEmpty(prhs[2]) )
	{
		mxFree(w);
	}
}

void fgt_model(double *x , double *w , double sigma , int p , int K , double e ,
			   double *xc , double *A_k ,
			   int d , int Nx ,
			   int *indxc , int *indx , int *xhead , int *xboxsz , 
			   double *dist_C , double *C_k , int *heads , int *cinds , double *dx , double *prods  ,
			   int pd)
{
	Kcenter(x , d , Nx , K , 
		    xc , indxc , indx , xboxsz , 
		    dist_C);
	Compute_C_k(d , p , 
		        C_k , 
		        heads , cinds);
  Compute_A_k(x , w , xc , C_k , sigma , d , Nx , p , K , pd , 
		        A_k , 
		        indx  , dx , prods , heads );
}

void Kcenter(double *x , int d , int Nx , int K , 
			 double *xc , int *indxc , int *indx , int *xboxsz , 
			 double *dist_C)
{
	double *x_ind , *x_j;
	register double temp ;
	int i , j , ind , nd , ibase;
	
	/* randomly pick one node as the first center. */
	/*	srand( (unsigned)time( NULL ) ); */
	/*	ind      = rand() % Nx; */
	ind      = 1;
	*indxc++ = ind;
	x_j      = x;
	x_ind    = x + ind*d;
	
	for (j = 0 ; j < Nx ; x_j += d , j++)
	{
		dist_C[j] = (j==ind) ? 0.0 : ddist(x_j , x_ind , d);
		indx[j]   = 0;
	}
	
	for(i = 1 ; i < K ; i++)
	{	 
		ind      = idmax(dist_C , Nx);
		*indxc++ = ind; 
		x_j      = x;
		x_ind    = x + ind*d;
		for (j = 0 ; j < Nx ; x_j += d, j++)
		{
			temp = (j==ind) ? 0.0 : ddist(x_j , x_ind , d);
			if (temp < dist_C[j])
			{
				dist_C[j] = temp;
				indx[j]   = i;
			}
		}
	}
	
	for (i = 0 ; i < K ; i++)
	{
		xboxsz[i] = 0;
	}
	
	for (i = 0; i < d*K; i++)
	{
		xc[i] = 0.0;
	}
	
	for (i = 0 , nd = 0 ; i < Nx ; i++ , nd += d)
	{
		xboxsz[indx[i]]++;
		ibase = indx[i]*d;
		for (j = 0 ; j < d; j++)
		{
			xc[j + ibase ] += x[j + nd];
		}
	}
	
	for (i = 0 , ibase = 0 ; i < K ; i++ , ibase += d)
	{
		temp = 1.0/xboxsz[i];
		for (j = 0; j < d; j++)
		{
			xc[j + ibase] *= temp;			
		}
	}	
}

void Compute_C_k(int d , int p , 
				 double *C_k , 
				 int *heads , int *cinds)
{
  int i , k , t , j , tail , head;
	for (i = 0; i < d; i++)
	{
		heads[i] = 0;
	}
	heads[d] = INT_MAX;
	cinds[0] = 0;
	C_k[0]   = 1.0;
	for (k=1 , t=1, tail = 1 ; k < p ; k++ , tail=t)
	{
		for (i = 0; i < d; i++)
		{
			head     = heads[i];
			heads[i] = t;
			for ( j = head ; j < tail ; j++ , t++)
			{
				cinds[t] = (j < heads[i+1]) ? cinds[j] + 1 : 1;
				C_k[t]   = 2.0 * C_k[j];
				C_k[t]  /= (double) cinds[t];
			}
		}
	}
}

void Compute_A_k(double *x , double *w , double *xc, double *C_k , double sigma , int d , int Nx , int p , int K , int pd , 
				 double *A_k , 
				 int *indx  , double *dx , double *prods , int *heads )
{
	int n , i , k , t , tail , j , head , ind;
	int nbase , ix2c , ix2cbase;
	register double sum , ctesigma = 1.0/(sigma) , temp , temp1;

	for (i = 0; i < pd*K; i++)
	{
		A_k[i] = 0.0;
	}
	
	for (n = 0 ; n < Nx ; n++)
	{
		nbase    = n*d;
		ix2c     = indx[n];
		ix2cbase = ix2c*d;
		ind      = ix2c*pd;
    temp     = w[n];
		sum      = 0.0;
		
		for (i = 0 ; i < d ; i++)
		{
			dx[i]    = (x[i + nbase] - xc[i + ix2cbase])*ctesigma;
			sum     += dx[i] * dx[i];
			heads[i] = 0;
		}		
		prods[0] = exp(-sum);
		
		for (k = 1 , t = 1 , tail = 1 ; k < p ; k++ , tail = t)
		{
			for (i = 0 ; i < d; i++)
			{
				head     = heads[i];
				heads[i] = t;
				temp1    = dx[i];
				for ( j = head; j < tail ; j++, t++)
				{
					prods[t] = temp1 * prods[j];
				}
			} 
		}
		
		for (i = 0 ; i < pd ; i++)
		{
			A_k[i + ind] += temp*prods[i];
		}
	}
	
	for (k = 0 ; k < K ; k++)
	{
		ind  = k*pd;
		for (i = 0 ; i < pd ; i++)
		{
			A_k[i + ind] *= C_k[i];
		}
	}
}

int nchoosek(int n , int k)
{
	int i , n_k = n - k , nchsk = 1;
	if (k < n_k)
	{
		k   = n_k;
		n_k = n - k;
	}
	for ( i = 1 ; i <= n_k ; i++)
	{
		nchsk *= (++k);
		nchsk /= i;
	}
	return nchsk;
}

double ddist(double *x , double *y , int d)
{
	int i;
	register double t , s = 0.0;
	for (i = 0 ; i < d ; i++)
	{
		t  = (x[i] - y[i]);
		s += (t * t);
	}
	return s;
}

int idmax(double *x , int N)
{
	int i , k = 0;
	double t = -1.0;
	for (i = 0 ; i < N ; i++ )
	{
		if( t < x[i] )
		{
			t = x[i];
			k = i;
		}
	}
	return k;
}