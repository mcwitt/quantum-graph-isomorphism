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
    int a[N][N], j, iter, max_iter = 100;

    filename = "test-graph-8.txt";

    if (! amatrix_load(filename, 8, &a[0][0]))
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

    for (i = 0; i < D; i++) x[i] = 1.;
    s = 0.9;

    iter = nlcg_minimize(gradient, line_min, 100, 1e-18, x);
    printf("stopped after %d iterations. (%d max)\n", iter, max_iter);

    energy = gq_energy_grad(s, d, x, grad);
    printf("energy = %.10g\n", energy);

    return EXIT_SUCCESS;
}
