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
double s;           /* adiabatic parameter */
double edrvr;       /* driver part of energy */

double obj_x2_grad(double psi[D], double *x2, double grad[D])
{
    return qaa_energy_grad(s, d, psi, grad, x2, &edrvr);
}

double line_min(double psi[D], double delta[D])
{
    return qaa_line_min(s, d, psi, delta);
}

int main(int argc, char *argv[])
{
    const gsl_rng_type *T;
    gsl_rng *rng;
    params_t p;
    double h[N];    /* fields */
    double energy, h0, mx, mz, psi2, q2, r2;
    UINT i;
    int a[N][N], ifile, ip, is, iter, j;

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

                energy = nlcg_minimize_norm_ind(obj_x2_grad, line_min, p.eps,
                        p.itermax, &iter, psi, &psi2, delta, r, &r2);

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
