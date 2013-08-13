/**
 * @file    qaa.h
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

#include "global.h"

/**
 * Compute the diagonal elements of the Hamiltonian.
 * @param[in]   a   adjacency matrix of graph
 * @param[in]   h   field at each site
 * @param[out]  d   diagonal elements of problem Hamiltonian
 */
void qaa_compute_diagonals(int a[N][N], double h[N], double d[D]);

/**
 * Compute the energy and energy gradient.
 * @param[in]   s       adiabatic parameter
 * @param[in]   d       diagonal elements of problem Hamiltonian
 * @param[in]   psi     wavefunction
 * @param[out]  grad    energy gradient
 * @param[out]  psi2    squared norm of wavefunction
 * @param[out]  edrvr   driver (off-diagonal) part of energy
 */
double qaa_energy_grad(double s, double d[D], double psi[D],
        double grad[D], double *psi2, double *edrvr);

/**
 * Minimize the energy along a specified direction.
 * @param[in]   s       adiabatic parameter
 * @param[in]   d       diagonal elements of problem Hamiltonian
 * @param[in]   psi     wavefunction
 * @param[in]   delta   search direction
 */
double qaa_line_min(double s, double d[D], double psi[D], double delta[D]);

double qaa_sigma_z(double psi[D], int j);

double qaa_sigma_x(double psi[D], int j);

double qaa_sigma2_z(double psi[D], int j, int k);

double qaa_mag_z(double psi[D]);

double qaa_mag_x(double psi[D]);

double qaa_overlap(double psi[D]);
