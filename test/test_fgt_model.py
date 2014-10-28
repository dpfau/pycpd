from cpd.fast_gaussian_transform import model, predict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def test_fgt_model_1():
  d   = 3
  Nx  = 100
  x   = np.random.randn((d, Nx))

  xc, Ak = model(x)

def test_fgt_model_2():
  d  = 3
  Nx = 10
  Ny = 100
  x  = np.random.randn((d, Nx))
  y  = np.random.randn((d, Ny))
  w  = np.random.rand((1 , Nx))
  h  = 2

  e  = 10
  p  = 6
  K  = 5

  xc, Ak = model(x, w, h, e, K, p)
  v = predict(y, xc, Ak, h, e)

 
def test_fgt_model_3():
  d    = 2
  R    = np.array([[2 , 0.4], [0.4, 3]])
  Nx   = 1000

  h    = 1
  e    = 10
  K    = np.round(np.sqrt(Nx))
  p    = 6

  vect = xrange(-5,5,0.3)
  Ny   = len(vect)
  w    = (1.0/Nx)*np.ones((1, Nx))
    
  x    = (np.linalg.chol(R).transpose().dot(np.random.randn((d, Nx))))

  X, Y = np.meshgrid(vect)
  y    = np.concatenate([X.flatten(), Y.flatten()]).transpose()

  xc, Ak = model(x, w, h, e, K, p);
  vy     = predict(y, xc, Ak, h, e, K, p);

  density = vy.reshape((Ny, Ny));

  fig = plt.figure()
  ax = fig.gca(projection='3d')
  surf = ax.plot_surface(X, Y, density, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False, alpha=0.5)
  ax.plot(x[1, :], x[2, :], 'r+', xc[1, :], xc[2, :], 'ko', ms=10)

  fig.colorbar(surf, shrink=0.5, aspect=5)
  plt.show()