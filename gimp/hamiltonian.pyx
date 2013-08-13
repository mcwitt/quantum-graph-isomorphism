import numpy as np
import scipy.sparse
from gimp import *

from libc.stdint cimport uint32_t
from libc.math cimport sqrt
cimport numpy as np

ctypedef uint32_t UINT

def hamiltonian(np.ndarray[int, ndim=2] a,
                np.ndarray[double] h,
                double s):
    
    N = len(a)
    d = 2**N
    nnz = d*N*(N+3)/2   # number of nonzero entries
    
    cdef np.ndarray[int]    rows = np.zeros(nnz, dtype=np.intc)
    cdef np.ndarray[int]    cols = np.zeros(nnz, dtype=np.intc)
    cdef np.ndarray[double] vals = np.zeros(nnz, dtype=np.double)
    cdef double c
    cdef int j, k
    cdef UINT i, m, inz = 0

    # transverse field
    for i in range(d):
        m = 1
        for j in range(N):
            rows[inz] = i
            cols[inz] = i^m
            vals[inz] = 1. - s
            inz += 1
            m <<= 1
            
    # longitudinal field
    for i in range(d):
        for k in range(N):
            rows[inz] = i
            cols[inz] = i
            vals[inz] = s * h[k] * spin(i, k)
            inz += 1
                    
    #  antiferromagnetic exchange
    for i in range(d):
        for j in range(N):
            c = s * spin(i, j)
            for k in range(j):
                if a[j, k] == 1:
                    rows[inz] = i
                    cols[inz] = i
                    vals[inz] = c * spin(i, k)
                    inz += 1

    return scipy.sparse.coo_matrix((vals, (rows, cols)), shape=(d, d))
