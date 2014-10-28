cdef extern from "fgt/fgt_predict.h":
  int invnchoosek(int, int)
  void fgt_predict(double *, double *, double *, int, double, int, double, int, int,
           double *, 
           double *, double *, int *)

cdef extern from "fgt/fgt_model.h":
  int nchoosek(int, int)
  int idmax(double *, int)
  double ddist(double *, double *, int)
  void Kcenter(double *, int, int, int, 
       double *, int *, int *, int *, 
       double *)
  void Compute_C_k(int, int, double *, int *, int *)
  void Compute_A_k(double *, double *, double *, double *, double, int, int, int, int, int, 
         double *, 
         int *, double *, double *, int *)
  void fgt_model(double *, double *, double, int, int, double,
         double *, double *,
         int, int,
         int *, int *, int *, int *,
         double *, double *, int *, int *, double *, double *, 
         int)

  cdef extern from "fgt/dval.h":
    void dval(double *, double *, double *, double, double *, int, int, int)