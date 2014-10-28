import numpy as np
cimport numpy as np

cpdef predict(np.ndarray[np.double_t, ndim=2] y,
              np.ndarray[np.double_t, ndim=2] xc,
              np.ndarray[np.double_t, ndim=2] Ak,
              double sigma=1, double e=10):
  """predict: returns the Fast Gauss Transform approximation of the test points y 
              given the model \theta=(xc,A_k)

  Parameters
  ----------
  y:      test point (d x Ny) numpy array
  xc:     Kcenter point (d x K) numpy array
  Ak:     Polynomial coefficient (pd x K), where pd = nchoosek(p+d-1, d)
  sigma:  Noise Standard deviation of the kernel (default sigma = 1)
  e:      Ratio of far field (default e = 10)

  Returns
  -------
  v:      density (1 x Ny) numpy array
  """
  cdef int d, pd, K, Ny
  cdef np.ndarray[np.double_t, ndim=1] dx, prods
  cdef np.ndarray[np.int_t, ndim=1] heads

  d  = y.shape[0]
  Ny = y.shape[1]
  K = xc.shape[1]
  assert d == xc.shape[0]

  pd = Ak.shape[0]
  assert K == Ak.shape[1]

  cdef np.ndarray[np.double_t, ndim=2] v = np.ndarray([1,Ny], dtype=np.double)

  dx    = np.ndarray([d], dtype=np.double)
  prods = np.ndarray([pd], dtype=np.double)
  heads = np.ndarray([d+1], dtype=np.int)

  fgt_predict(<double*> y.data, <double*> xc.data, <double*> Ak.data, Ny, sigma, K, e, d, pd, 
              <double*> v.data, 
              <double*> dx.data, <double*> prods.data, <int*> heads.data)

  return v

cpdef model(np.ndarray[np.double_t, ndim=2] x,
            np.ndarray w=None,
            double sigma=1, double e=10, int K=None, int p=8):
  """model: returns the Fast Gauss Transform Aprroximation Model of a Kernel density

  Parameters
  ----------
  x:            Source data (d x Nx) numpy array
  w:            Weigths (1 x Nx) numpy array ( default w = ones(1 , Nx) ) 
  sigma:        Noise Standard deviation of the kernel (default sigma = 1)
  e:            Ratio of far field (default e = 10)
  K:            Number of centers (default K = sqrt(Nx))
  p:            Order of truncation (default p = 8)

  Returns
  -------
  xc:           The K center points of the training set (d x K) numpy array
  Ak:           Polynomial coefficient (pd x K) numpy array, 
                where pd = nchoosek(p + d - 1 , d) = prod(1:p + d - 1)/(prod(1:p - 1)*prod(1:d))
  """
  cdef int d, Nx, wx, pd
  cdef np.ndarray[np.double_t, ndim=2] xc, Ak
  cdef np.ndarray[np.double_t, ndim=1] C_k, dist_C, dx, prods
  cdef np.ndarray[np.int_t, ndim=1] indxc, indx, xhead, xboxsz, heads, cinds

  d  = x.shape[0]
  Nx = x.shape[1]

  if w == None:
    w = np.ones([1,Nx], dtype=np.double)
  else:
    assert w.dtype == np.double and w.ndim == 2
    wx = w.shape[0]
    assert wx == 1 and w.shape[1] == Nx

  if K == None:
    K = np.sqrt(Nx)
  else:
    assert K <= Nx

  pd = nchoosek(p+d-1, d);
  xc = np.ndarray([d,K], dtype=np.double)
  Ak = np.ndarray([pd,K], dtype=np.double)

  C_k     = np.ndarray([pd], dtype=np.double)
  dist_C  = np.ndarray([Nx], dtype=np.double)
  dx      = np.ndarray([d], dtype=np.double)
  prods   = np.ndarray([pd], dtype=np.double)

  indxc  = np.ndarray([K], dtype=np.int)
  indx   = np.ndarray([Nx], dtype=np.int)
  xhead  = np.ndarray([K], dtype=np.int)
  xboxsz = np.ndarray([K], dtype=np.int)
  heads  = np.ndarray([d + 1], dtype=np.int)
  cinds  = np.ndarray([pd], dtype=np.int)

  fgt_model(<double*> x.data, <double*> w.data, sigma, p, K, e, 
            <double*> xc.data, <double*> Ak.data, 
            d, Nx,
            <int*> indxc.data, <int*> indx.data, <int*> xhead.data, <int*> xboxsz.data, 
            <double*> dist_C.data, <double*> C_k.data, <int*> heads.data, <int*> cinds.data,
            <double*> dx.data, <double*> prods.data,  
            pd)

  return xc, Ak