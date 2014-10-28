#include "fgt/fgt_predict.h"

void fgt_predict(double *y , double *xc , double *A_k  , int Ny, double sigma , int K , double e , int d , int pd ,  
			     double *v , 
			     double *dy , double *prods , int *heads)
{
	int p , i , j , m , k , t , tail , mbase , kn , xbase , head , ind;
	double sum2 , ctesigma = 1.0/(sigma) , temp , temp1;
	p = invnchoosek(d , pd);
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




