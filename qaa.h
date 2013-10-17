/**
 * @file    qaa.h
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

#include "nlcg.h"

typedef struct
{
    double s, *edrvr;
    const double *d;
} qaa_args_t;

/**
 * Compute the diagonal elements of the Hamiltonian for h_j = 0.
 * @param[in]   b   independent adjacency matrix entries (A_21, A_31, A_32, A_41, etc.)
 * @param[out]  d   diagonal elements of \f$H_d\f$
 */
void qaa_compute_diagonals(const int *b, double *d);

/**
 * Update diagonal elements when \f$h_j \rightarrow h_j + \delta\f$.
 * @param[in]       dh  change in field
 * @param[in,out]   d   diagonal elements of \f$H_d\f$
 */
void qaa_update_diagonals(double dh, double *d);

/**
 * Update diagonal elements for a change in field at a single site.
 * @param[in]       j   index of modified field
 * @param[in]       dh  change in field
 * @param[in,out]   d   diagonal elements of \f$H_d\f$
 */
void qaa_update_diagonals_1(int j, double dh, double *d);

/**
 * Initialize conjugate gradient algorithm.
 * @param[in]       s       Adiabatic parameter.
 * @param[in]       d       Diagonal elements of \f$H_d\f$.
 * @param[in,out]   psi     Initial guess for wavefunction.
 * @param[out]      edrvr   Driver energy.
 * @param[out]      nlcg    nlcg_t instance.
 */
double qaa_nlcg_init(
        double s,
        const double *d,
        double *psi,
        double *edrvr,
        qaa_args_t *args,
        nlcg_t *nlcg
        );

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
