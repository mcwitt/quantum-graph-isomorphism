#include <stdlib.h>
#include <stdint.h>
#include "amatrix.h"
#include "qgi.h"
#include "params.h"

#define DS          0.01
#define SMIN        DS
#define SMAX        1.

#define MAX_ITER    1000
#define EPS         1e-16

double d[D];    /* diagonal elements of H_p */
double psi[D];  /* wavefunction */

int main(int argc, char *argv[])
{
    char **file_names;
    double h[N];    /* fields */
    double s, energy, mx, mz, qz;
    index_t i;
    int a[N][N], ifile, j, iter;

    if (argc < 2)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    file_names = &argv[1];

    printf("%16s %6s %12s %12s %12s %12s %12s\n",
            "file", "s", "iterations", "E0", "mz2", "mx2", "q");

    for (ifile = 0; ifile < argc - 1; ifile++)
    {
        if (! amatrix_load(file_names[ifile], N, &a[0][0]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], file_names[ifile]);

            return EXIT_FAILURE;
        }

        for (s = SMIN; s < SMAX; s += DS)
        {
            for (j = 0; j < N; j++) h[j] = 1.;
            qgi_compute_problem_hamiltonian(a, h, d);
            for (i = 0; i < D; i++) psi[i] = 1.;
            iter = qgi_minimize_energy(s, d, MAX_ITER, EPS, &energy, psi);
            mz = qgi_mag_z(psi);
            mx = qgi_mag_x(psi);
            qz = qgi_overlap(psi);

            printf("%16s %6.3f %12d %12g %12g %12g %12g\n",
                    file_names[ifile], s, iter, energy, mz*mz, mx*mx, qz);
        }
    }

    return EXIT_SUCCESS;
}
