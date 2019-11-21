from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from glob import glob

# Use python setup.py build_ext --inplace
# to compile


setup(
  name = "ekats2", version='0.1.3', description='wrapper for eKATS2',
  ext_modules = cythonize(["ekats2.pyx"])
)
