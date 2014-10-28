#include "fgt/fgt_model.h"
#include <limits.h>

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