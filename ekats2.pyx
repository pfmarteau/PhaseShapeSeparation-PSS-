# distutils: language = c++
# distutils: sources = eKATS2.cpp

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string
import numpy as np

cdef extern from "string" namespace "std":
    cdef cppclass string:
        char* c_str()

# c++ interface to cython
cdef extern from "eKATS2.h" namespace "ekats2":

  cdef cppclass eKats2:
    eKats() except +
    double SqEuclideanDistance(vector[vector[double]] x, vector[vector[double]] y) 
    vector[vector[double]] interpolate(vector[vector[double]])
    double kdtw(vector[vector[double]], vector[vector[double]], double nu)
    vector[vector[double]] iTEKA(vector[vector[double]], vector[vector[vector[double]]], double nu)
    vector[vector[vector[double]]] iTEKA_stdev(vector[vector[double]], vector[vector[vector[double]]], double nu)

# creating a cython wrapper class
cdef class PyeKats2:
    cdef eKats2 *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new eKats2()
    def SqEuclideanDistance(self, sv1, sv2):
        return self.thisptr.SqEuclideanDistance(sv1, sv2)
    def interpolate(self, sv):
        return np.array(self.thisptr.interpolate(sv))
    def kdtw(self, sv1, sv2, nu):
        return self.thisptr.kdtw(sv1, sv2, nu)
    def iTEKA(self, sv1, ds, nu):
        return np.array(self.thisptr.iTEKA(sv1, ds, nu), dtype=object)
    def iTEKA_stdev(self, sv1, ds, nu):
        return np.array(self.thisptr.iTEKA_stdev(sv1, ds, nu), dtype=object)

