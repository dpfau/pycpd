cdef extern from "fgt/fgt_predict.h":
  int invnchoosek(int , int)
  void fgt_predict(double *, double *, double *, int, double, int, double, int, int,
           double *, 
           double *, double *, int *)