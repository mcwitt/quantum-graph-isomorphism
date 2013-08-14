import numpy as np

cimport numpy as np
from libc.stdint cimport uint32_t

ctypedef uint32_t UINT

cdef extern from "global.h":
    int n
    UINT ndim

cdef extern from "qaa.h":
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

def minimize_energy(
        s,
        np.ndarray[np.double_t,ndim=1] d,
        eps,
        max_iter,
        np.ndarray[np.double_t,ndim=1] psi,
        np.ndarray[np.double_t,ndim=1] delta,
        np.ndarray[np.double_t,ndim=1] resid,
        full_output=False
        ):

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

