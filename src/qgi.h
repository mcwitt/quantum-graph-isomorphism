#include "params.h"

/**
 * Compute the diagonal elements of the problem Hamiltonian.
 * @param[in] a Adjacency matrix of graph.
 * @param[in] h Field at each site.
 * @param[out] d Diagonal elements of problem Hamiltonian.
 */
void qgi_compute_problem_hamiltonian(int a[N][N], double h[N], double d[D]);

/**
 * Compute the energy and energy gradient.
 * @param[in] s Adiabatic parameter.
 * @param[in] d Diagonal elements of problem Hamiltonian.
 * @param[in] psi Wavefunction.
 * @param[out] grad Energy gradient.
 * @param[out] psi2 Squared norm of wavefunction.
 * @param[out] edrvr Driver (off-diagonal) part of energy.
 */
double qgi_energy_grad(double s, double d[D], double psi[D],
        double grad[D], double *psi2, double *edrvr);

/**
 * Minimize the energy along a specified direction.
 * @param[in] s Adiabatic parameter.
 * @param[in] d Diagonal elements of problem Hamiltonian.
 * @param[in] psi Wavefunction.
 * @param[in] delta Search direction.
 */
double qgi_line_min(double s, double d[D], double psi[D], double delta[D]);

int qgi_minimize_energy(double s, double d[D], int max_iter, double eps,
        double *energy, double psi[D]);
