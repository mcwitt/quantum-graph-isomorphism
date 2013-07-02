#include "params.h"

/** Compute the diagonal elements of the problem Hamiltonian. */
void qgi_compute_problem_hamiltonian(int a[N][N], double h[N], double d[D]);

double qgi_energy_grad(double s, double d[D], double psi[D],
        double grad[D], double *eod);

double qgi_line_min(double s, double d[D], double psi[D], double delta[D]);

