#include <stdlib.h>
#include <stdint.h>
#include "amatrix.h"
#include "qgi.h"
#include "nlcg.h"
#include "params.h"

#define S           0.9
#define H_0         1.
#define X_0_i       1.
#define MAX_ITER    100
#define EPS         1e-12

double d[D];

/* wrappers for gradient and line minimization functions
 * to be passed to NLCG routine */
void gradient(double *x, double *grad) { qgi_energy_grad(S, d, x, grad); }
double line_min(double *x, double *delta) { return qgi_line_min(S, d, x, delta); }

int main(int argc, char *argv[])
{
    char **file_names;
    double h[N];
    double x[D], grad[D];
    double energy;
    index_t i;
    int a[N][N], ifile, j, iter;

    if (argc < 2)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    file_names = &argv[1];

    printf("#%11s %12s %12s\n", "file", "E_0", "iterations");

    for (ifile = 0; ifile < argc - 1; ifile++)
    {
        if (! amatrix_load(file_names[ifile], N, &a[0][0]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], file_names[ifile]);

            return EXIT_FAILURE;
        }

        for (j = 0; j < N; j++) h[j] = H_0;
        qgi_compute_problem_hamiltonian(a, h, d);
        for (i = 0; i < D; i++) x[i] = X_0_i;
        iter = nlcg_minimize(gradient, line_min, MAX_ITER, EPS, x);
        energy = qgi_energy_grad(S, d, x, grad);
        printf("%12s %12.9g %12d\n", file_names[ifile], energy, iter);
    }

    return EXIT_SUCCESS;
}
