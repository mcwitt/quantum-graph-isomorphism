from libc.stdint cimport uint32_t
from libc.math cimport sqrt
cimport numpy as np

ctypedef uint32_t UINT

cdef int spin(UINT i, int j):
    return ((int)((i >> j) & 1) << 1) - 1

cdef double sigma_z(np.ndarray[double] psi, int j):
    cdef double result = 0.
    cdef UINT i
    for i in range(len(psi)): result += psi[i]**2 * spin(i, j)
    return result

cdef double sigma2_z(np.ndarray[double] psi, int j, int k):
    cdef double result = 0.
    cdef UINT i
    for i in range(len(psi)): result += psi[i]**2 * spin(i, j) * spin(i, k)
    return result

cdef double sigma_x(np.ndarray[double] psi, int j):
    cdef double result = 0.
    cdef UINT i, m = 1UL << j
    for i in range(len(psi)): result += psi[i] * psi[i^m]
    return result

def mag_z(int n, np.ndarray[double] psi):
    result = 0.
    for j in range(n): result += sigma_z(psi, j)
    return result / N

def mag_x(int n, np.ndarray[double] psi):
    result = 0.
    for j in range(n): result += sigma_x(psi, j)
    return result / n

def overlap(int n, np.ndarray[double] psi):
    result = 0.
    for j in range(1, n):
        for k in range(j):
            result += sigma2_z(psi, j, k)**2

    return sqrt(2. * result / n / (n-1))
