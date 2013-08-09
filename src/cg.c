#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "amatrix.h"
#include "defs.h"
#include "nlcg.h"
#include "params.h"
#include "qgi.h"

double d[D];        /* diagonal elements of problem hamiltonian */
double psi[D];      /* wavefunction */
double delta[D];    /* CG search direction */
double r[D];        /* residual */
double s;           /* adiabatic parameter */
double edrvr;       /* driver part of energy */

double obj_x2_grad(double psi[D], double *x2, double grad[D])
{
    return qgi_energy_grad(s, d, psi, grad, x2, &edrvr);
}

double line_min(double psi[D], double delta[D])
{
    return qgi_line_min(s, d, psi, delta);
}

int main(int argc, char *argv[])
{
    const gsl_rng_type *T;
    gsl_rng *rng;
    params_t p;
    double h[N];    /* fields */
    double energy, h0, h0mult, mx, mz, psi2, q2, r2;
    index_t i;
    int a[N][N], ifile, ih, is, iter, j, nh, ns;

    params_from_cmd(&p, argc, argv);

    if (p.num_files == 0)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    ns = (int) (1 + (p.smax - p.smin) / p.ds);
    nh = 1 + p.nh_dec * p.ndec;

    printf("%16s %9s %9s %12s %12s %12s %12s %12s %12s\n",
            "file", "h0", "s", "iterations", "res2", "energy", "mz", "mx", "q2");

    for (ifile = 0; ifile < argc - 1; ifile++)
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

        h0 = p.hmin;
        h0mult = pow(10., 1./p.nh_dec);

        for (ih = 0; ih < nh; ih++)
        {
            for (j = 0; j < N; j++) h[j] = h0;
            qgi_compute_problem_hamiltonian(a, h, d);

            s = p.smin;
                
            for (is = 0; is < ns; is++)
            {
                energy = nlcg_minimize_norm_ind(obj_x2_grad, line_min, p.eps,
                        p.itermax, &iter, psi, &psi2, delta, r, &r2);

                /* normalize wavefunction */
                psi2 = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= psi2;

                mz = qgi_mag_z(psi);
                mx = qgi_mag_x(psi);
                q2 = qgi_overlap(psi);

                printf("%16s %9g %9g %12d %12g %12g %12g %12g %12g\n",
                        p.files[ifile], h0, s, iter, r2, energy, mz, mx, q2);

                s += p.ds;
            }

            h0 *= h0mult;
        }
    }

    return EXIT_SUCCESS;
}
