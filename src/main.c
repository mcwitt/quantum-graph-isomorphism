#include <stdlib.h>
#include <stdint.h>
#include "amatrix.h"
#include "gimpq.h"
#include "nlcg.h"
#include "params.h"

double d[D];
double s;

/* wrappers for gradient and line minimization functions
 * to be passed to NLCG routine */
void gradient(double *x, double *grad) { gq_energy_grad(s, d, x, grad); }
double line_min(double *x, double *delta) { return gq_line_min(s, d, x, delta); }

int main(int argc, char *argv[])
{
    char *filename;
    double h[N];
    double x[D], grad[D];
    double energy;
    ULONG i;
    int a[N][N], j;

    filename = "test-graph-16.txt";

    if (! amatrix_load(filename, 16, &a[0][0]))
    {
        fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                argv[0], filename);

        return EXIT_FAILURE;
    }

    for (j = 0; j < N; j++) h[j] = 1.;

    gq_compute_problem_hamiltonian(a, h, d);

    /*
    for (j = 0; j < D; j++)
        printf("%g\n", d[j]);
    exit(0);
    */

    for (i = 0; i < D; i++) x[i] = 0.;
    s = 0.99;

    nlcg_minimize(gradient, line_min, x);

    energy = gq_energy_grad(s, d, x, grad);
    printf("energy = %g\n", energy);

    return EXIT_SUCCESS;
}
