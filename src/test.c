#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "amatrix.h"
#include "nlcg.h"
#include "qgi.h"
#include "params.h"

#define SMIN    0.01
#define SMAX    0.99
#define DS      0.01

#define HMIN    0.01
#define NDEC    4
#define NH_DEC  25

#define MAX_ITER    300
#define EPS         1e-16

#define NS      (int) (1 + (SMAX - SMIN) / DS)
#define NH      1 + NH_DEC * NDEC

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
    double h0, h0mult, energy, mx, mz, q2, r2;
    index_t i;
    int a[N][N], ifile, ih, is, iter, j;

    if (argc < 2)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    file_names = &argv[1];

    printf("%16s %9s %9s %12s %12s %12s %12s %12s %12s\n",
            "file", "h0", "s", "iterations", "res2", "energy", "mz", "mx", "q2");

    for (ifile = 0; ifile < argc - 1; ifile++)
    {
        if (! amatrix_load(file_names[ifile], N, &a[0][0]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], file_names[ifile]);

            return EXIT_FAILURE;
        }

        for (i = 0; i < D; i++) psi[i] = 1.;

        h0 = HMIN;
        h0mult = pow(10., 1./NH_DEC);

        for (ih = 0; ih < NH; ih++)
        {
            for (j = 0; j < N; j++) h[j] = h0;
            qgi_compute_problem_hamiltonian(a, h, d);

            s = SMIN;
                
            for (is = 0; is < NS; is++)
            {
                energy = nlcg_minimize(obj_grad, line_min, EPS, MAX_ITER,
                        &iter, psi, delta, r, &r2);

                /* normalize wavefunction */
                psi2 = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= psi2;

                mz = qgi_mag_z(psi);
                mx = qgi_mag_x(psi);
                q2 = qgi_overlap(psi);

                printf("%16s %9g %9g %12d %12g %12g %12g %12g %12g\n",
                        file_names[ifile], h0, s, iter, r2, energy, mz, mx, q2);

                s += DS;
            }

            h0 *= h0mult;
        }
    }

    return EXIT_SUCCESS;
}
