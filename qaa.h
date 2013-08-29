/**
 * @file    qaa.h
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

/**
 * Compute the diagonal elements of the Hamiltonian for h_j = 0.
 * @param[in]   b   independent adjacency matrix entries (A_21, A_31, A_32, A_41, etc.)
 * @param[out]  d   diagonal elements of problem Hamiltonian
 */
void qaa_compute_diagonals(const int *b, double *d);

/**
 * Update diagonal elements when \f$h_j \rightarrow h_j + \delta\f$.
 * @param[in]       dh  change in field
 * @param[in,out]   d   diagonal elements of problem Hamiltonian
 */
void qaa_update_diagonals(double dh, double *d);

/**
 * Update diagonal elements for a change in field at a single site.
 * @param[in]       j   index of modified field
 * @param[in]       dh  change in field
 * @param[in,out]   d   diagonal elements of problem Hamiltonian
 */
void qaa_update_diagonals_1(int j, double dh, double *d);

double qaa_minimize_energy(
        double s,
        const double *d,
        double eps,
        int max_iter,
        int *num_iter,
        double *edrvr,
        double *psi,
        double *psi2,
        double *delta,
        double *r,
        double *r2
        );

/**
 * Compute the energy and energy gradient.
 * @param[in]   s       adiabatic parameter
 * @param[in]   d       diagonal elements of problem Hamiltonian
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
 * @param[in]   d       diagonal elements of problem Hamiltonian
 * @param[in]   psi     wavefunction
 * @param[in]   delta   search direction
 */
double qaa_line_min(
        double s,
        const double *d,
        const double *psi,
        const double *delta);

double qaa_sigma_z(const double *psi, int j);

double qaa_sigma_x(const double *psi, int j);

double qaa_sigma2_z(const double *psi, int j, int k);

double qaa_mag_z(const double *psi);

double qaa_mag_x(const double *psi);

double qaa_overlap(const double *psi);
