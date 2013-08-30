#cython: embedsignature=True

import numpy as np
import scipy.sparse

cimport numpy as np
from libc.stdint cimport uint32_t

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

def load_amatrix(filename):
    f = open(filename, 'r')
    a = []
    for line in f: a.append(list(line.strip()))
    f.close()
    return np.matrix(a, dtype=np.dtype('i'))

def save_amatrix(filename, a):
    f = open(filename, 'w')
    for row in a:
        for entry in row:
            f.write('{:d}'.format(entry))
        f.write('\n')
    f.close()

def tobitstr(a):
    N = a.shape[0]
    b = [a[i, j] for i in range(1, N) for j in range(i)]
    return np.array(b, dtype=np.dtype('i'))

def compute_diagonals(int[:] b, double[:] d=None):
    if d is None: d = np.empty(D, dtype=np.dtype('d'))
    qaa_compute_diagonals(&b[0], &d[0])
    return np.asarray(d)

def update_diagonals(double dh, double[:] d):
    qaa_update_diagonals(dh, &d[0])

def update_diagonals_1(int j, double dh, double[:] d):
    qaa_update_diagonals_1(j, dh, &d[0])

def minimize_energy(
        s,
        double[:] d,
        eps = 1e-10,
        max_iter = 300,
        double[:] psi = None,
        double[:] delta = None,
        double[:] resid = None,
        normalize = True,
        full_output = False
        ):

    '''
    returns: energy, psi
    full output: energy, psi, resid, delta, edrvr, num_iter
    '''

    cdef double edrvr, psi2, r2
    cdef int num_iter

    if psi is None: psi = np.random.normal(scale=1./np.sqrt(D), size=D)
    if delta is None: delta = np.empty(D, dtype=np.dtype('d'))
    if resid is None: resid = np.empty(D, dtype=np.dtype('d'))

    energy = qaa_minimize_energy(s, &d[0], eps, max_iter, &num_iter, &edrvr,
            &psi[0], &psi2, &delta[0], &resid[0], &r2)

    if normalize:
        delta /= np.sqrt(np.dot(delta, delta))
        norm = np.sqrt(np.dot(psi, psi))
        psi /= norm
        resid /= norm

    if full_output:
        return energy, np.asarray(psi), np.asarray(resid), np.asarray(delta), edrvr, num_iter

    return energy, np.asarray(psi)

def hamiltonian(int[:, :] a, double[:] h, double s):
    
    cdef int[:] rows, cols
    cdef double[:] vals
    cdef double c, oms = 1. - s
    cdef int j, k
    cdef UINT i, ia, m

    n = (N+1)*D # number of nonzero entries

    rows = np.empty(n, dtype=np.dtype('i'))
    cols = np.empty(n, dtype=np.dtype('i'))
    vals = np.empty(n, dtype=np.dtype('d'))

    # problem hamiltonian
    for i in range(D):
        rows[i] = cols[i] = i
        vals[i] = 0.
        for j in range(N):
            c = s * spin(i, j)
            vals[i] -= c * h[j]
            for k in range(j):
                if a[j, k] == 1:
                    vals[i] += c * spin(i, k)

    ia = D

    # transverse field
    for i in range(D):
        m = 1
        for j in range(N):
            rows[ia] = i
            cols[ia] = i^m
            vals[ia] = oms
            ia += 1
            m <<= 1

    return scipy.sparse.coo_matrix((vals, (rows, cols)), shape=(D, D)).tocsr()

def mag_z(double[:] psi):
    return qaa_mag_z(&psi[0])

def mag_x(double[:] psi):
    return qaa_mag_x(&psi[0])

def overlap(double[:] psi):
    return qaa_overlap(&psi[0])

def sigma_z(double[:] psi):
    cdef double[:] result
    cdef double c2
    cdef int j
    cdef UINT i

    result = np.zeros(N, dtype=np.dtype('d'))

    for i in range(D):
        c2 = psi[i]**2
        for j in range(N):
            result[j] += c2 * spin(i, j)

    return np.asarray(result)

def sigma2_z(double[:] psi):
    cdef double[:,:] result
    cdef double c2
    cdef int j, k, s_j
    cdef UINT i

    result = np.zeros((N, N), dtype=np.dtype('d'))

    for i in range(D):
        c2 = psi[i]**2
        for j in range(N):
            s_j = spin(i, j)
            for k in range(j):
                result[j, k] += c2 * s_j * spin(i, k)

    a = np.asarray(result)
    a = a + a.T - np.diag(a.diagonal()) # symmetrize result
    return a

def sigma_x(double[:] psi):
    cdef double[:] result
    cdef int j
    cdef UINT i, m

    result = np.zeros(N, dtype=np.dtype('d'))

    for i in range(D):
        for j in range(N):
            m = 1UL << j
            result[j] += psi[i] * psi[i^m]

    return np.asarray(result)

