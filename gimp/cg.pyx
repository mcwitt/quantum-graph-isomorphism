from gimp import *

ctypedef uint32_t UINT

cdef np.ndarray[double] problem_hamiltonian_diagonals(
        int n,
        UINT ndim,
        double s,
        np.ndarray[int, ndim=2] a,
        np.ndarray[double] h):

    cdef int j, k
    cdef UINT i

    d = np.zeros(ndim)

    for i in range(ndim):
        for j in range(n):
            cdef int s_j = spin(i, j)
            d[i] -= h[j] * s_j
            for k in range(j):
                if a[j][k] == 1:
                    d[i] += s_j * spin(i, k)

    return d

cdef double driver_matrix_element(
        int n,
        UINT ndim,
        np.ndarray[double] u,
        np.ndarray[double] v):
{
    cdef double result = 0.
    cdef int j
    cdef UINT i, m

    for i in range(ndim):
        m = 1
        for j in range(n):
            result += u[i] * v[i^m]
            m <<= 1

    return result
}

cdef double problem_matrix_element(
        int n,
        UINT ndim,
        np.ndarray[double] d,
        np.ndarray[double] u,
        np.ndarray[double] v):
{
    cdef double prod, udotv = 0., result = 0.
    cdef UINT i;

    for i in range(ndim):
        prod = u[i] * v[i];
        udotv += prod;
        result += prod * d[i];

    return result, udotv
}

cdef double matrix_element(
        int n,
        UINT ndim,
        double s,
        np.ndarray[double] d,
        np.ndarray[double] u,
        np.ndarray[double] v):

    pme, udotv = problem_matrix_element(n, ndim, d, u, v)

    return (1. - s) * driver_matrix_element(u, v) + s * pme(d, u, v, udotv),
           udotv

cdef double energy_grad(
        int n,
        UINT ndim,
        double s,
        np.ndarray[double] d,
        np.ndarray[double] psi,
        np.ndarray[double] grad):

    cdef double energy
    cdef UINT i

    # compute the energy
    edrvr = driver_matrix_element(n, ndim, psi, psi)
    energy, psi2 = problem_matrix_element(n, ndim, d, psi, psi)
    energy = (1. - s) * edrvr + s * energy
    energy /= psi2

    # compute the gradient
    for i in range(ndim):
        cdef int j
        cdef double sum = 0.
        cdef UINT m = 1

        for j in range(n):
            sum += psi[i^m]
            m <<= 1

        grad[i] = 2. * (psi[i] * (s * d[i] - energy) + (1. - s) * sum) / psi2

    return energy, edrvr, psi2

cdef double line_min(
        int n,
        UINT ndim,
        double s,
        np.ndarray[double] d,
        np.ndarray[double] psi,
        np.ndarray[double] delta):

    cdef double psi2, psi_dot_delta, delta2,
           psi_H_psi, psi_H_delta, delta_H_delta,
           a, b, c, coef, sqr, alpha

    psi_H_psi,     psi2          = matrix_element(n, ndim, s, d, psi,   psi)
    psi_H_delta,   psi_dot_delta = matrix_element(n, ndim, s, d, psi,   delta)
    delta_H_delta, delta2        = matrix_element(n, ndim, s, d, delta, delta)

    a = psi_dot_delta * delta_H_delta - delta2 * psi_H_delta
    b = psi2 * delta_H_delta - delta2 * psi_H_psi
    c = psi2 * psi_H_delta - psi_dot_delta * psi_H_psi

    coef = -0.5 / a
    sqr = sqrt(b*b - 4.*a*c)
    alpha = coef * (b - sqr)

    # if critical point is a maximum, use other solution
    if 2*a*alpha + b < 0: alpha = coef * (b + sqr)

    return alpha
