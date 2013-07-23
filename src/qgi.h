/**
 * @file   qgi.h
 * @author Matt Wittmann
 * @brief  Implementations of quantum graph isomorphism algorithms.
 */

#include "params.h"
#include <stdint.h>

typedef uint64_t index_t;

/**
 * A structure to store parameters and intermediate values for the CG energy
 * minimization
 * */
typedef struct
{
    double d[D];        /**< diagonal elements of problem hamiltonian */
    double psi[D];      /**< wavefunction */

    double delta[D];    /**< CG search direction */
    double r[D];        /**< residual */
} qgi_t;

/**
 * Initialize data structures.
 * @param[out]      qgi pointer to qgi_t data structure
 * @param[in]       a   adjacency matrix of graph
 * @param[in]       h   fields at each site (vertex)
 */
void qgi_init(qgi_t *qgi, int a[N][N], double h[N]);

/**
 * Compute the diagonal elements of the problem Hamiltonian.
 * @param[in]   a   adjacency matrix of graph
 * @param[in]   h   field at each site
 * @param[out]  d   diagonal elements of problem Hamiltonian
 */
void qgi_compute_problem_hamiltonian(int a[N][N], double h[N], double d[D]);

/**
 * Compute the energy and energy gradient.
 * @param[in]   s       adiabatic parameter
 * @param[in]   d       diagonal elements of problem Hamiltonian
 * @param[in]   psi     wavefunction
 * @param[out]  grad    energy gradient
 * @param[out]  psi2    squared norm of wavefunction
 * @param[out]  edrvr   driver (off-diagonal) part of energy
 */
double qgi_energy_grad(double s, double d[D], double psi[D],
        double grad[D], double *psi2, double *edrvr);

/**
 * Minimize the energy along a specified direction.
 * @param[in]   s       adiabatic parameter
 * @param[in]   d       diagonal elements of problem Hamiltonian
 * @param[in]   psi     wavefunction
 * @param[in]   delta   search direction
 */
double qgi_line_min(double s, double d[D], double psi[D], double delta[D]);

/**
 * Find the ground state using the nonlinear conjugate gradient algorithm.
 * @param[in,out]   qgi pointer to qgi_t data structure
 * @param[in]       s           adiabatic parameter
 * @param[in]       eps         stop when residual is reduced by this factor
 * @param[out]      num_iter    number of CG iterations used 
 */
double qgi_minimize_energy(qgi_t *qgi, double s, double eps, int *num_iter);

double qgi_sigma_z(double psi[D], int j);

double qgi_sigma_x(double psi[D], int j);

double qgi_sigma2_z(double psi[D], int j, int k);

double qgi_mag_z(double psi[D]);

double qgi_mag_x(double psi[D]);

double qgi_overlap(double psi[D]);
