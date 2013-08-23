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

    void qaa_compute_diagonals(int a[], double d[])
    void qaa_update_diagonals(double dh, double d[])
    void qaa_update_diagonals_1(int j, double dh, double d[])

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

    double qaa_mag_z(double psi[])
    double qaa_mag_x(double psi[])
    double qaa_overlap(double psi[])


N = n
D = ndim

cdef int spin(UINT i, int j):
    return ((int)((i >> j) & 1) << 1) - 1

def load_graph(filename):
    f = open(filename, 'r')
    a = []
    for line in f: a.append(list(line.strip()))
    f.close()
    return np.array(a, dtype=np.int)

def save_graph(filename, a):
    f = open(filename, 'w')
    for row in a:
        for entry in row:
            f.write('{:d}'.format(entry))
        f.write('\n')
    f.close()

def tobitstr(a):
    N = a.shape[0]
    b = [a[i, j] for i in range(1, N) for j in range(i)]
    return np.array(b, dtype=np.int32)

def compute_diagonals(b, np.ndarray[np.double_t, ndim=1] d=None):
    cdef np.ndarray[np.int32_t, ndim=1] buf = b.astype(np.int32)
    if d is None: d = np.empty(D, dtype=np.double)
    qaa_compute_diagonals(<int*> buf.data, <double*> d.data)
    return d

def update_diagonals(
        double dh,
        np.ndarray[np.double_t, ndim=1] d):

    qaa_update_diagonals(dh, <double*> d.data)

def update_diagonals_1(
        int j,
        double dh,
        np.ndarray[np.double_t, ndim=1] d):

    qaa_update_diagonals_1(j, dh, <double*> d.data)

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
    
    cdef np.ndarray[np.int_t, ndim=1] rows = np.zeros(nnz, dtype=np.int)
    cdef np.ndarray[np.int_t, ndim=1] cols = np.zeros(nnz, dtype=np.int)
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
    return qaa_mag_z(<double*> psi.data)

def mag_x(np.ndarray[np.double_t, ndim=1] psi):
    return qaa_mag_x(<double*> psi.data)

def overlap(np.ndarray[np.double_t, ndim=1] psi):
    return qaa_overlap(<double*> psi.data)

def sigma_z(np.ndarray[np.double_t, ndim=1] psi):
    cdef double c2
    cdef int j
    cdef np.ndarray[np.double_t, ndim=1] result = np.zeros(N, dtype=np.double)
    cdef UINT i

    for i in range(D):
        c2 = psi[i]**2
        for j in range(N):
            result[j] += c2 * spin(i, j)

    return result

def sigma2_z(np.ndarray[np.double_t, ndim=1] psi):
    cdef double c2
    cdef int j, k, s_j
    cdef np.ndarray[np.double_t, ndim=2] result = np.zeros((N, N), dtype=np.double)
    cdef UINT i

    for i in range(D):
        c2 = psi[i]**2
        for j in range(N):
            s_j = spin(i, j)
            for k in range(j):
                result[j, k] += c2 * s_j * spin(i, k)

    result = result + result.T - np.diag(result.diagonal)   # symmetrize result
    return result

def sigma_x(np.ndarray[np.double_t, ndim=1] psi):
    cdef int j
    cdef np.ndarray[np.double_t, ndim=1] result = np.zeros(N, dtype=np.double)
    cdef UINT i, m

    for i in range(D):
        for j in range(N):
            m = 1UL << j
            result[j] += psi[i] * psi[i^m]

    return result

