/**
 * @file    qaa.h
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

#include "nlcg.h"

typedef struct
{
    nlcg_t cg;
    double s;
    double d[D];    /* diagonal elements of problem Hamiltonian */
    double edrvr;
} qaa_t;

/**
 * Compute problem hamiltonian and initialize minimization algorithm. Initially
 * all \f$h_j\f$ are set to zero. Return initial energy.
 * @param[out]  p   qaa_t instance.
 * @param[in]   b   Independent adjacency matrix entries (A_21, A_31, A_32, A_41, etc.)
 */
double qaa_init(qaa_t *p, const int *b, const double *psi);

/**
 * Update problem hamiltonian when \f$h_j \rightarrow h_j + \delta\f$.
 * @param[in,out]   p   qaa_t instance.
 * @param[in]       dh  Change in field.
 */
void qaa_shift_field(qaa_t *p, double dh);

/**
 * Update problem hamiltonian for a change in field at a single site.
 * @param[in,out]   p   qaa_t instance.
 * @param[in]       j   Index of modified field.
 * @param[in]       dh  Change in field.
 */
void qaa_shift_field_1(qaa_t *p, int j, double dh);

/**
 * Reset the energy minimization algorithm. Return energy.
 * @param   p       qaa_t instance.
 * @param   psi     Wavefunction.
 */
double qaa_reset(qaa_t *p, const double *psi);

/**
 * Do one iteration of the energy-minimization algorithm. Return energy.
 * @param[in,out]   p   qaa_t instance.
 */
double qaa_iterate(qaa_t *p, double *psi);

/**
 * Return the minimum energy.
 * @param[in,out]   p           qaa_t instance.
 * @param[in,out]   psi         Wavefunction.
 * @param[in]       tol         Stop when \f$|r|/|\psi|\f$ is less than *tol*.
 * @param[in]       max_iter    Maximum number of iterations.
 * @param[out]      num_iter    Number of iterations used.
 */
double qaa_minimize(qaa_t *p, double *psi, double tol, int max_iter, int *num_iter);

/**
 * Compute the energy and energy gradient.
 * @param[in]   s       Adiabatic parameter.
 * @param[in]   d       Diagonal elements of \f$H_d\f$.
 * @param[in]   psi     Wavefunction.
 * @param[out]  grad    Energy gradient.
 * @param[out]  psi2    Squared norm of wavefunction.
 * @param[out]  edrvr   Driver (off-diagonal) part of energy.
 */
double qaa_energy_grad(
        double s,
        const double *d,
        const double *psi,
        double *grad,
        double *psi2,
        double *edrvr);

/**
 * Minimize the energy along a specified direction.
 * @param[in]   s       Adiabatic parameter.
 * @param[in]   d       Diagonal elements of \f$H_d\f$.
 * @param[in]   psi     Wavefunction.
 * @param[in]   delta   Search direction.
 */
double qaa_line_min(
        double s,
        const double *d,
        const double *psi,
        const double *delta);

/**
 * Compute matrix element of the driver Hamiltonian, \f$H_d\f$.
 * @param[in]   u   Bra vector.
 * @param[in]   v   Ket vector.
 */
double qaa_me_driver(const double *u, const double *v);

/**
 * Compute matrix element of the problem Hamiltonian, \f$H_p\f$.
 * @param[in]   d       Diagonal elements of \f$H_d\f$.
 * @param[in]   u       Bra vector.
 * @param[in]   v       Ket vector.
 * @param[out]  udotv   \f$\langle u \mid v \rangle\f$.
 */
double qaa_me_problem(
        const double *d,
        const double *u,
        const double *v,
        double *udotv);

double qaa_sigma_z(const double *psi, int j);

double qaa_sigma_x(const double *psi, int j);

double qaa_sigma2_z(const double *psi, int j, int k);

double qaa_mag_z(const double *psi);

double qaa_overlap(const double *psi);
