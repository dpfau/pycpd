void comp_p(
    double* x,
    double* y, 
    double* sigma2,
    double* outlier,
    double* P1,
    double* Pt1,
    double* Px,
    double* E,
    int N,
    int M,
    int D)
{
  int   n, m, d;
  double  ksig, diff, razn, outlier_tmp, sp;
  double  *P, *temp_x;
  
  P = (double*) calloc(M, sizeof(double));
  temp_x = (double*) calloc(D, sizeof(double));
  
  ksig = -2.0 * *sigma2;
  outlier_tmp=(*outlier*M*pow (-ksig*3.14159265358979,0.5*D))/((1-*outlier)*N); 
 /* printf ("ksig = %lf\n", *sigma2);*/
  /* outlier_tmp=*outlier*N/(1- *outlier)/M*(-ksig*3.14159265358979); */
  
  
  for (n=0; n < N; n++) {
      
      sp=0;
      for (m=0; m < M; m++) {
          razn=0;
          for (d=0; d < D; d++) {
             diff=*(x+n+d*N)-*(y+m+d*M);  diff=diff*diff;
             razn+=diff;
          }
          
          *(P+m)=exp(razn/ksig);
          sp+=*(P+m);
      }
      
      sp+=outlier_tmp;
      *(Pt1+n)=1-outlier_tmp/ sp;
      
      for (d=0; d < D; d++) {
       *(temp_x+d)=*(x+n+d*N)/ sp;
      }
         
      for (m=0; m < M; m++) {
         
          *(P1+m)+=*(P+m)/ sp;
          
          for (d=0; d < D; d++) {
          *(Px+m+d*M)+= *(temp_x+d)**(P+m);
          }
          
      }
      
   *E +=  -log(sp);     
  }
  *E +=D*N*log(*sigma2)/2;
    
  
  free((void*)P);
  free((void*)temp_x);

  return;
}

void comp_truncate(
    double* x,
    double* y, 
    double* sigma2,
    double* outlier,
    double* P1,
    double* Pt1,
    double* Px,
    double* E,
    int N,
    int M,
    int D,
    double* truncate)
{
  int   n, m, d;
  double  ksig, diff, razn, outlier_tmp, sp;
  double  *P;
  
  P = (double*) calloc(M, sizeof(double));

  ksig = -2.0 * *sigma2;
  *truncate=log(*truncate);
  if (*outlier==0) *outlier=1e-8;
  
  outlier_tmp=(*outlier*M*pow (-ksig*3.14159265358979,0.5*D))/((1-*outlier)*N); 
 /* printf ("ksig = %lf\n", *sigma2);*/
  /* outlier_tmp=*outlier*N/(1- *outlier)/M*(-ksig*3.14159265358979); */
  
  
  for (n=0; n < N; n++) {
      sp=0;
      
      for (m=0; m < M; m++) {
          razn=0;
          for (d=0; d < D; d++) {
             diff=*(x+n+d*N)-*(y+m+d*M);  diff=diff*diff;
             razn+=diff;
          };
          
          razn=razn/ksig;
          
          if (razn< *truncate)
          {*(P+m)=0;}
          else
          {
          *(P+m)=exp(razn);
          sp+=*(P+m);
          };
      }
      
      sp+=outlier_tmp;
      *(Pt1+n)=1-outlier_tmp/ sp;
      
            
      for (m=0; m < M; m++) {
   
          if (*(P+m)==0){}
          else
          {
          *(P1+m)+=*(P+m)/sp;
          for (d=0; d < D; d++) {
          *(Px+m+d*M)+= *(x+n+d*N)/ sp* *(P+m);
          }
          }
          
      }
      
   *E +=  -log(sp);     
  }
  *E +=D*N*log(*sigma2)/2;
    

  free((void*)P);
  

  return;
}

void comp_correspondence(
    double* x,
    double* y, 
    double* sigma2,
    double* outlier,
    double* Pc,
    int N,
    int M,
    int D)
{
  int   n, m, d;
  double  ksig, diff, razn, outlier_tmp,temp_x,sp;
  double  *P, *P1;
  
  
  P = (double*) calloc(M, sizeof(double));
  P1 = (double*) calloc(M, sizeof(double));
  

  ksig = -2.0 * (*sigma2+1e-3);
  outlier_tmp=(*outlier*M*pow (-ksig*3.14159265358979,0.5*D))/((1-*outlier)*N); 
  if (outlier_tmp==0) outlier_tmp=1e-10;
  
  
 /* printf ("ksig = %lf\n", *sigma2);*/
  
  
  for (n=0; n < N; n++) {
      sp=0;
      for (m=0; m < M; m++) {
          razn=0;
          for (d=0; d < D; d++) {
              diff=*(x+n+d*N)-*(y+m+d*M);  diff=diff*diff;
              razn+=diff;
          }
          
          *(P+m)=exp(razn/ksig);
          sp+=*(P+m);
          
      }
      
      sp+=outlier_tmp;
      
      
      for (m=0; m < M; m++) {
          
          temp_x=*(P+m)/ sp;
          
          if (n==0)
          {*(P1+m)= *(P+m)/ sp;
          *(Pc+m)=n+1;};
          
          if (temp_x > *(P1+m))
          {
              *(P1+m)= *(P+m)/ sp;
              *(Pc+m)=n+1;
          }
    
              
      }
      
  }

   
  free((void*)P);
  free((void*)P1);

  return;
}