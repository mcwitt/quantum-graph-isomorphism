/**
 * @file   ising.h
 * @author Matt Wittmann
 * @brief  Algorithms related to the quantum Ising model.
 */

#include "global.h"

/**
 * Compute the diagonal elements of the Hamiltonian.
 * @param[in]   a   adjacency matrix of graph
 * @param[in]   h   field at each site
 * @param[out]  d   diagonal elements of Hamiltonian
 */
void ising_diagonals(int a[N][N], double h[N], double d[D]);

/**
 * Compute the energy and energy gradient.
 * @param[in]   d       diagonal elements of Hamiltonian
 * @param[in]   psi     wavefunction
 * @param[out]  grad    energy gradient
 * @param[out]  psi2    squared norm of wavefunction
 */
double ising_energy_grad(double d[D], double psi[D],
        double grad[D], double *psi2);

double ising_sigma_z(double psi[D], int j);

double ising_sigma_x(double psi[D], int j);

double ising_sigma2_z(double psi[D], int j, int k);

double ising_mag_z(double psi[D]);

double ising_mag_x(double psi[D]);

double ising_overlap(double psi[D]);
