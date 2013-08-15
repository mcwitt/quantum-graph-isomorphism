#!python
#cython: embedsignature=True

import numpy as np
import scipy.sparse

cimport numpy as np
from libc.stdint cimport uint32_t
from libc.math cimport sqrt

ctypedef uint32_t UINT

cdef extern from "global.h":
    int n
    UINT ndim

cdef extern from "qaa.h":
    void qaa_compute_diagonals(
        int a[],
        double h[],
        double d[]
        )
    double qaa_minimize_energy(
        double  s,
        double  d[],
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  *edrvr,
        double  psi[],
        double  *psi2,
        double  delta[],
        double  r[],
        double  *r2
        )

N = n
D = ndim

def load_graph(filename):
    f = open(filename, 'r')
    a = []
    for line in f: a.append(list(line.strip()))
    return np.array(a, dtype=np.int)

def tobitstr(a):
    N = a.shape[0]
    b = [a[i, j] for i in range(1, N) for j in range(i)]
    return np.array(b, dtype=int)

def compute_diagonals(
        np.ndarray[np.int_t, ndim=1] b,
        np.ndarray[np.double_t, ndim=1] h,
        np.ndarray[np.double_t, ndim=1] d=None,
        inplace=False
        ):

    if not inplace:
        d = np.empty(D, dtype=np.double)

    qaa_compute_diagonals(
            <int*> b.data,
            <double*> h.data,
            <double*> d.data
            )

    if not inplace: return d

def minimize_energy(
        s,
        np.ndarray[np.double_t, ndim=1] d,
        eps,
        max_iter,
        np.ndarray[np.double_t, ndim=1] psi,
        np.ndarray[np.double_t, ndim=1] delta,
        np.ndarray[np.double_t, ndim=1] resid
        ):

    'returns: energy, edrvr, psi2, r2, num_iter'

    cdef double edrvr, psi2, r2
    cdef int num_iter

    energy = qaa_minimize_energy(
            s,
            <double*> d.data,
            eps,
            max_iter,
            &num_iter,
            &edrvr,
            <double*> psi.data,
            &psi2,
            <double*> delta.data,
            <double*> resid.data,
            &r2
            )

    return energy, edrvr, psi2, r2, num_iter

def hamiltonian(
        np.ndarray[np.int_t, ndim=2] a,
        np.ndarray[np.double_t, ndim=1] h,
        double s
        ):
    
    d = 2**N
    nnz = d*N*(N+3)/2   # number of nonzero entries
    
    cdef np.ndarray[np.int_t, ndim=1]    rows = np.zeros(nnz, dtype=np.int)
    cdef np.ndarray[np.int_t, ndim=1]    cols = np.zeros(nnz, dtype=np.int)
    cdef np.ndarray[np.double_t, ndim=1] vals = np.zeros(nnz, dtype=np.double)
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

def mag_z(np.ndarray[np.double_t, ndim=1] psi):
    result = 0.
    for j in range(N): result += sigma_z(psi, j)
    return result / N

def mag_x(np.ndarray[np.double_t, ndim=1] psi):
    result = 0.
    for j in range(N): result += sigma_x(psi, j)
    return result / N

def overlap(np.ndarray[np.double_t, ndim=1] psi):
    result = 0.
    for j in range(1, N):
        for k in range(j):
            result += sigma2_z(psi, j, k)**2

    return sqrt(2. * result / N / (N-1))

cdef int spin(UINT i, int j):
    return ((int)((i >> j) & 1) << 1) - 1

cdef double sigma_z(np.ndarray[np.double_t, ndim=1] psi, int j):
    cdef double result = 0.
    cdef UINT i
    for i in range(len(psi)): result += psi[i]**2 * spin(i, j)
    return result

cdef double sigma2_z(np.ndarray[np.double_t, ndim=1] psi, int j, int k):
    cdef double result = 0.
    cdef UINT i
    for i in range(len(psi)): result += psi[i]**2 * spin(i, j) * spin(i, k)
    return result

cdef double sigma_x(np.ndarray[np.double_t, ndim=1] psi, int j):
    cdef double result = 0.
    cdef UINT i, m = 1UL << j
    for i in range(len(psi)): result += psi[i] * psi[i^m]
    return result

