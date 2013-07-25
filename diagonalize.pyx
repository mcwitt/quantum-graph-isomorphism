import numpy as np
import scipy.sparse
from libc.stdint cimport uint32_t
cimport numpy as np

cdef inline int spin(uint32_t i, int j):
    return (((1 << j) & i) >> j)*2 - 1

cdef inline uint32_t neighbor(uint32_t i, int j):
    return (1 << j) ^ i

def hamiltonian5(np.ndarray[long, ndim=2] a,
                 np.ndarray[double] h,
                 double s):
    
    N = len(a)
    d = 2**N
    nnz = d*N*(N+3)/2 # number of nonzero entries
    
    cdef np.ndarray[int]    rows = np.zeros(nnz, dtype=np.intc)
    cdef np.ndarray[int]    cols = np.zeros(nnz, dtype=np.intc)
    cdef np.ndarray[double] vals = np.zeros(nnz, dtype=np.double)
    cdef int m, n, s_n
    cdef uint32_t i, inz = 0

    # transverse field
    for i in range(d):
        for m in range(N):
            rows[inz] = i
            cols[inz] = neighbor(i, m)
            vals[inz] = 1. - s
            inz += 1
            
    # longitudinal field
    for i in range(d):
        for n in range(N):
            rows[inz] = i
            cols[inz] = i
            vals[inz] = s * h[n] * spin(i, n)
            inz += 1
                    
    #  antiferromagnetic exchange
    for i in range(d):
        for n in range(N):
            s_n = spin(i, n)
            for m in range(n):
                if a[n, m] == 1:
                    rows[inz] = i
                    cols[inz] = i
                    vals[inz] = s * s_n * spin(i, m)
                    inz += 1

    return scipy.sparse.coo_matrix((vals, (rows, cols)), shape=(d, d))
