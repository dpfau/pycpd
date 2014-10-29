cdef extern from "comp.h":
  void comp_p(double*, double*, double*, double*, double*, double*, double*, double*, int, int, int)
  void comp_truncate(double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, double*)
  void comp_correspondence(double*, double*, double*, double*, double*, int, int, int)
