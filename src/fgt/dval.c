/*
  Evaluate weighted gaussian RBF between vectors x,y
  Author : Sébastien PARIS : sebastien.paris@lsis.org
  Modified for Python by David PFAU : david.pfau@gmail.com
*/

void dval(double *x , double *y , double *w , double sigma , double *v , int d , int Nx , int Ny )
{
  int i , j , l , id , jd;  
  register double temp , res ;
  double tempv;
  double cte = -1.0/(sigma*sigma);
  int num_threads;
  
#ifdef OMP 
    num_threads          = min(MAX_THREADS,omp_get_num_procs());
    omp_set_num_threads(num_threads);
#endif

  for (i = 0; i < Ny ; i++) 
  {
    id    = i*d;  
    tempv = 0.0;
#ifdef OMP 
/* #pragma omp parallel for default(none) private(j,l,temp,jd,res) shared(x,y,w,d,id,Nx,cte,tempv) */
#pragma omp parallel for default(none) firstprivate(l,temp,res) lastprivate(j,jd,tempv) shared(x,y,w,d,id,Nx,cte)
#endif
    for (j = 0 ; j < Nx ; j++) 
    {
      jd  = j*d;
      res = 0.0;
      for (l = 0 ; l < d ; l++)
      {
        temp  = (y[l + id] - x[l + jd]);  
        res  +=temp*temp;
      }
      tempv += w[j]*exp(cte*res); 
    }
    v[i] = tempv;
  }
}
