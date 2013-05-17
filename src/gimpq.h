#include "params.h"

/** Compute the diagonal elements of the problem Hamiltonian. */
void gq_compute_problem_hamiltonian(int a[][N], double h[N], double d[D]);

double gq_energy_grad(double s, double d[D], double psi[D], double grad[D]);

double gq_line_min(double s, double d[D], double psi[D], double delta[D]);

