#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <libgen.h>
#include "amatrix.h"
#include "global.h"
#include "nlcg.h"
#include "params.h"
#include "qaa.h"

double d[D];        /* diagonal elements of problem hamiltonian */
double psi[D];      /* wavefunction */
double delta[D];    /* CG search direction */
double r[D];        /* residual */

int main(int argc, char *argv[])
{
    const gsl_rng_type *T;
    gsl_rng *rng;
    params_t p;
    double h[N];    /* fields */
    double edrvr, energy, h0, mx, mz, psi2, q2, r2, s;
    int a[N][N], ifile, ip, is, iter, j;
    UINT i;

    params_from_cmd(&p, argc, argv);

    if (p.num_files == 0)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    printf("%16s %12s %12s %12s %12s %12s %12s %12s %12s\n",
            "graph", "h0", "s", "iterations", "res2", "energy", "mz", "mx", "q2");

    for (ifile = 0; ifile < p.num_files; ifile++)
    {
        if (! amatrix_load(p.files[ifile], N, &a[0][0]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], p.files[ifile]);

            return EXIT_FAILURE;
        }

        /* generate random initial wavefunction */
        gsl_rng_env_setup();
        T = gsl_rng_default;
        rng = gsl_rng_alloc(T);
        for (i = 0; i < D; i++) psi[i] = gsl_ran_gaussian(rng, 1.) /sqrt(D);

        for (ip = 0; ip < p.np; ip++)
        {
            h0 = pow(10., p.emin + (p.emax - p.emin)*ip/(p.np - 1.)) ;
            for (j = 0; j < N; j++) h[j] = h0;
            qaa_compute_diagonals(a, h, d);

            for (is = 0; is < p.ns; is++)
            {
                s = p.smin;
                if (p.ns > 1) s += (p.smax - p.smin)*is/(p.ns - 1.);

                energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                        &edrvr, psi, &psi2, delta, r, &r2);

                /* normalize wavefunction */
                psi2 = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= psi2;

                mz = qaa_mag_z(psi);
                mx = qaa_mag_x(psi);
                q2 = qaa_overlap(psi);

                printf("%16s %12g %12g %12d %12g %12g %12g %12g %12g\n",
                        basename(p.files[ifile]), h0, s, iter, r2, energy, mz, mx, q2);
            }
        }
    }

    return EXIT_SUCCESS;
}
