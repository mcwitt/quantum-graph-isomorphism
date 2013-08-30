/**
 * @file    qaa.h
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

/**
 * Compute the diagonal elements of the Hamiltonian for h_j = 0.
 * @param[in]   b   independent adjacency matrix entries (A_21, A_31, A_32, A_41, etc.)
 * @param[out]  d   diagonal elements of f$H_df$
 */
void qaa_compute_diagonals(const int *b, double *d);

/**
 * Update diagonal elements when \f$h_j \rightarrow h_j + \delta\f$.
 * @param[in]       dh  change in field
 * @param[in,out]   d   diagonal elements of f$H_df$
 */
void qaa_update_diagonals(double dh, double *d);

/**
 * Update diagonal elements for a change in field at a single site.
 * @param[in]       j   index of modified field
 * @param[in]       dh  change in field
 * @param[in,out]   d   diagonal elements of f$H_df$
 */
void qaa_update_diagonals_1(int j, double dh, double *d);

/**
 * Find the ground state using the conjugate gradient algorithm.
 * @param[in]   s           adiabatic parameter
 * @param[in]   d           diagonal elements of f$H_df$
 * @param[in]   eps         maximum norm of residual
 * @param[in]   max_iter    maximum number of iterations
 * @param[out]  num_iter    number of iterations
 * @param[out]  edrvr       driver energy
 * @param[out]  psi         ground state vector
 * @param[out]  psi2        squared norm of wavefunction
 * @param[out]  r           residual vector
 * @param[out]  r2          squared norm of residual
 * @param[out]  delta       CG search direction at last iteration
 */
double qaa_minimize_energy(
        double s,
        const double *d,
        double eps,
        int max_iter,
        int *num_iter,
        double *edrvr,
        double *psi,
        double *psi2,
        double *r,
        double *r2,
        double *delta
        );

/**
 * Compute the energy and energy gradient.
 * @param[in]   s       adiabatic parameter
 * @param[in]   d       diagonal elements of f$H_df$
 * @param[in]   psi     wavefunction
 * @param[out]  grad    energy gradient
 * @param[out]  psi2    squared norm of wavefunction
 * @param[out]  edrvr   driver (off-diagonal) part of energy
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
 * @param[in]   s       adiabatic parameter
 * @param[in]   d       diagonal elements of f$H_df$
 * @param[in]   psi     wavefunction
 * @param[in]   delta   search direction
 */
double qaa_line_min(
        double s,
        const double *d,
        const double *psi,
        const double *delta);

/**
 * Compute matrix element of the driver Hamiltonian, \f$H_d\f$.
 * @param[in]   u   bra vector
 * @param[in]   v   ket vector
 */
double qaa_me_driver(const double *u, const double *v);

/**
 * Compute matrix element of the problem Hamiltonian, \f$H_p\f$.
 * @param[in]   d       diagonal elements of \f$H_d\f$
 * @param[in]   u       bra vector
 * @param[in]   v       ket vector
 * @param[out]  udotv   \f$u \cdot v\f$
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
