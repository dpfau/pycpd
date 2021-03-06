from cpd.fast_gaussian_transform import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def test_fgt_run():
  d   = 3
  Nx  = 100
  x   = np.random.randn(d, Nx)

  xc, Ak = model(x)

def test_fgt_value():
  d  = 3
  Nx = 10
  Ny = 100
  x  = np.random.randn(d, Nx)
  y  = np.random.randn(d, Ny)
  w  = np.random.rand(1 , Nx)
  h  = 2

  e  = 10
  p  = 6
  K  = 5

  xc, Ak = model(x, w, h, e, K, p)
  v = predict(y, xc, Ak, h, e)

 
def test_fgt_plot():
  d    = 2
  R    = np.array([[2 , 0.4], [0.4, 3]])
  Nx   = 1000

  h    = 1
  e    = 10
  K    = np.round(np.sqrt(Nx))
  p    = 6

  vect = np.arange(-5,5,0.3)
  Ny   = len(vect)
  w    = (1.0/Nx)*np.ones((1, Nx))
    
  x    = (np.linalg.cholesky(R).transpose().dot(np.random.randn(d, Nx)))

  X, Y = np.meshgrid(vect, vect)
  y    = np.zeros((d,Ny*Ny))

  y[0,:] = X.flatten()
  y[1,:] = Y.flatten()

  xc, Ak = model(x, w, h, e, K, p);
  vy     = predict(y, xc, Ak, h, e);

  density = vy.reshape((Ny, Ny));

  fig = plt.figure()
  ax = fig.gca(projection='3d')
  surf = ax.plot_surface(X, Y, density, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False, alpha=0.5)
  ax.plot(x[0, :], x[1, :], 'r+', xc[0, :], xc[1, :], 'ko', ms=10)

  fig.colorbar(surf, shrink=0.5, aspect=5)
  plt.show()

def test_dval_1():
  d       = 10
  Nx      = 100
  Ny      = 10000
  x       = np.random.randn(d, Nx)
  y       = np.random.randn(d, Ny)
  v       = dval(x, y)

def test_dval_2():
  d       = 10
  Nx      = 100
  Ny      = 10000
  x       = np.random.randn(d, Nx)
  y       = np.random.randn(d, Ny)
  u       = np.random.rand(1, Nx)
  sigma   = 2
  v       = dval(x, y, u, sigma)
