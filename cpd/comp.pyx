#distutils: language = c
#distutils: sources = src/comp.c

from comp cimport *
import numpy as np
cimport numpy as np

cpdef p(np.ndarray[np.double_t, ndim=2] x,
        np.ndarray[np.double_t, ndim=2] y,
        np.ndarray[np.double_t, ndim=2] sigma2,
        np.ndarray[np.double_t, ndim=2] outlier):

  cdef int N, M, D
  N = x.shape[0]
  M = y.shape[0]
  D = x.shape[1]
  assert D == y.shape[1]

  cdef np.ndarray[np.double_t, ndim=2] P1, Pt1, Px, E
  P1  = np.ndarray([M,1], dtype=np.double);
  Pt1 = np.ndarray([N,1], dtype=np.double);
  Px  = np.ndarray([M,D], dtype=np.double);
  E   = np.ndarray([1,1], dtype=np.double);

  comp_p(<double*> x.data, <double*> y.data, <double*> sigma2.data, <double*> outlier.data,
         <double*> P1.data, <double*> Pt1.data, <double*> Px.data, <double*> E.data, N, M, D)

  return P1, Pt1, Px, E

cpdef truncate():
  pass

cpdef correspondence():
  pass
