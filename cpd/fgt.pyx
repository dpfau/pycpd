import numpy as np
cimport numpy as np

cdef predict(np.ndarray y, np.ndarray xc, np.ndarray Ak, double sigma=1, double e=10):
  return v

cpdef model(np.ndarray[np.double_t, ndim=2] x,
            np.ndarray[np.double_t, ndim=2] w=None,
            double sigma=1,
            double e=10,
            int K=None,
            int p=8):
  cdef int d, Nx
  cdef int pd
  cdef np.ndarray[np.double_t, ndim=2] xc, Ak
  cdef np.ndarray[np.double_t, ndim=1] C_k, dist_C, dx, prods
  cdef np.ndarray[np.int_t, ndim=1] indxc, indx, xhead, xboxsz, heads, cinds

  d, Nx = x.shape

  if w == None:
    cdef np.ndarray w = np.ones([1,Nx], dtype=DTYPE)
  else:
    cdef int wx, wy
    wx, wy = w.shape
    assert wx == 1 and wy == Nx

  if K == None:
    K = sqrt(Nx)
  else:
    assert K <= Nx

  pd = nchoosek(p + d - 1 , d);
  xc = np.zeros([d,K], dtype=np.double)
  Ak = np.zeros([pd,K], dtype=np.double)

  C_k     = np.zeros([pd], dtype=np.double)
  dist_C  = np.zeros([Nx], dtype=np.double)
  dx      = np.zeros([d], dtype=np.double)
  prods   = np.zeros([pd], dtype=np.double)

  indxc  = np.zeros([K], dtype=np.int)
  indx   = np.zeros([Nx], dtype=np.int)
  xhead  = np.zeros([K], dtype=np.int)
  xboxsz = np.zeros([K], dtype=np.int)
  heads  = np.zeros([d + 1], dtype=np.int)
  cinds  = np.zeros([pd], dtype=np.int)

  fgt_model(<double*> x.data, <double*> w.data, sigma, p, K, e, 
            <double*> xc.data, <double*> Ak.data, 
            d, Nx,
            <int*> indxc.data, <int*> indx.data, <int*> xhead.data, <int*> xboxsz.data, 
            <double*> dist_C.data, <double*> C_k.data, <int*> heads.data, <int*> cinds.data,
            <double*> dx.data, <double*> prods.data,  
            pd)

  return xc, Ak