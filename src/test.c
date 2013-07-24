#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "amatrix.h"
#include "nlcg.h"
#include "qgi.h"
#include "params.h"

#define DS          0.005
#define SMIN        DS
#define SMAX        1.

#define MAX_ITER    1000
#define EPS         1e-16

double d[D];        /* diagonal elements of problem hamiltonian */
double psi[D];      /* wavefunction */
double delta[D];    /* CG search direction */
double r[D];        /* residual */
double s;           /* adiabatic parameter */
double psi2;        /* squared norm of wavefunction */
double edrvr;       /* driver part of energy */

double obj_grad(double psi[D], double grad[D])
{
    return qgi_energy_grad(s, d, psi, grad, &psi2, &edrvr);
}

double line_min(double psi[D], double delta[D])
{
    return qgi_line_min(s, d, psi, delta);
}

int main(int argc, char *argv[])
{
    char **file_names;
    double h[N];    /* fields */
    double energy, mx, mz, q2;
    index_t i;
    int a[N][N], ifile, ih, iter, j;

    int nh = 5;
    double h0[] = {0.25, 0.5, 1., 2., 4.};

    if (argc < 2)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    file_names = &argv[1];

    printf("%16s %6s %6s %12s %12s %12s %12s %12s\n",
            "file", "h0", "s", "iterations", "energy", "mz2", "mx2", "q2");

    for (ifile = 0; ifile < argc - 1; ifile++)
    {
        if (! amatrix_load(file_names[ifile], N, &a[0][0]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], file_names[ifile]);

            return EXIT_FAILURE;
        }

        for (i = 0; i < D; i++) psi[i] = 1.;

        for (ih = 0; ih < nh; ih++)
        {
            for (j = 0; j < N; j++) h[j] = h0[ih];
            qgi_compute_problem_hamiltonian(a, h, d);
                
            for (s = SMIN; s < SMAX; s += DS)
            {
                energy = nlcg_minimize(obj_grad, line_min, EPS, MAX_ITER,
                        &iter, psi, delta, r);

                /* normalize wavefunction */
                psi2 = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= psi2;

                mz = qgi_mag_z(psi);
                mx = qgi_mag_x(psi);
                q2 = qgi_overlap(psi);

                printf("%16s %6.3f %6.3f %12d %12g %12g %12g %12g\n",
                        file_names[ifile], h0[ih], s, iter, energy, mz*mz, mx*mx, q2);
            }
        }
    }

    return EXIT_SUCCESS;
}
