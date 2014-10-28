import numpy as np

cdef predict(y, xc, A_k, sigma=1, e=10):
  return v

cdef model(x, w=None, sigma=1, e=10, K=None, p=8):
  d, Nx = x.shape

  if w == None:
    w = np.ones((1,Nx))
  else:
    wx, wy = w.shape
    if wx != 1 or wy != Nx:
      raise Error('w must be (1,Nx) array')

  if K == None:
    K = sqrt(Nx)

  return xc, A_k