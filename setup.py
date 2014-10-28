from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from numpy import get_include

extensions = [
  Extension("cpd.fast_gaussian_transform", ["cpd/fast_gaussian_transform.pyx"],
    include_dirs = ["include", get_include()])
]

setup(
  name = "PyCPD",
  ext_modules = cythonize(extensions)
)
