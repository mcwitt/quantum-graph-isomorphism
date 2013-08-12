import numpy as np
import scipy.sparse
from libc.stdint cimport uint32_t
from libc.math cimport sqrt
cimport numpy as np

def load_graph(filename):
    f = open(filename, 'r')
    a = []
    for line in f: a.append(list(line.strip()))
    a = np.array(a, dtype=int)
    assert(a.shape[0] == a.shape[1])
    return a

cdef inline int spin(uint32_t i, int j):
    return ((int)((i >> j) & 1) << 1) - 1

cdef double sigma_z(np.ndarray[double] psi, int j):
    cdef double result = 0.
    cdef uint32_t i
    for i in range(len(psi)): result += psi[i]**2 * spin(i, j)
    return result

cdef double sigma2_z(np.ndarray[double] psi, int j, int k):
    cdef double result = 0.
    cdef uint32_t i
    for i in range(len(psi)): result += psi[i]**2 * spin(i, j) * spin(i, k)
    return result

cdef double sigma_x(np.ndarray[double] psi, int j):
    cdef double result = 0.
    cdef uint32_t i, m = 1UL << j
    for i in range(len(psi)): result += psi[i] * psi[i^m]
    return result

def hamiltonian(np.ndarray[long, ndim=2] a,
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
    cdef uint32_t i, m, inz = 0

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

def mag_z(int N, np.ndarray[double] psi):
    result = 0.
    for j in range(N): result += sigma_z(psi, j)
    return result / N

def mag_x(int N, np.ndarray[double] psi):
    result = 0.
    for j in range(N): result += sigma_x(psi, j)
    return result / N

def overlap(int N, np.ndarray[double] psi):
    result = 0.
    for j in range(1, N):
        for k in range(j):
            result += sigma2_z(psi, j, k)**2

    return sqrt(2. * result / N / (N-1))
