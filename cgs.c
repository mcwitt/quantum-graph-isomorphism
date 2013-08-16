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

#define DH(h0)  1e-10

double d[D];        /* diagonal elements of problem hamiltonian */
double psi[D];      /* wavefunction */
double psi1[D];
double delta[D];    /* CG search direction */
double r[D];        /* residual */

int main(int argc, char *argv[])
{
    const double sqrt_D = sqrt(D);
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *rng;
    params_t p;
    double h[N], fj[N];
    double edrvr, energy, f, h0, mx, mz, norm, psi2, q2, r2, s;
    double ds = 0., mh0 = 0.;
    int a[N*(N-1)/2], ifile, ih, is, iter, j, k;
    UINT i;

    params_from_cmd(&p, argc, argv);

    if (p.num_files == 0)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (p.ns > 1) ds  = (p.smax - p.smin)/(p.ns - 1.);
    if (p.nh > 1) mh0 = pow(10., (p.emax - p.emin)/(p.nh - 1.)) ;

    rng = gsl_rng_alloc(T);

    printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
            "ver", "graph", "h0", "s",
            "iterations", "res2", "energy", "mz", "mx", "q2", "q2p");

    for (ifile = 0; ifile < p.num_files; ifile++)
    {
        if (! amatrix_load(p.files[ifile], a))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], p.files[ifile]);

            return EXIT_FAILURE;
        }

        /* generate random initial wavefunction */
        gsl_rng_env_setup();
        for (i = 0; i < D; i++) psi[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;
        h0 = pow(10., p.emin);

        for (ih = 0; ih < p.nh; ih++, h0 *= mh0)
        {
            for (j = 0; j < N; j++) h[j] = h0;
            qaa_compute_diagonals(a, h, d);
            s = p.smin;

            for (is = 0; is < p.ns; is++, s += ds)
            {
                energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                        &edrvr, psi, &psi2, delta, r, &r2);

                norm = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= norm;

                mz = qaa_mag_z(psi);
                mx = qaa_mag_x(psi);
                q2 = qaa_overlap(psi);

                f = 0.;

                /* vary the field at each site to approximate susceptibilities */
                for (j = 0; j < N; j++)
                {
                    qaa_update_diagonals(j, 0.5 * DH(h0), d);
                    for (i = 0; i < D; i++) psi1[i] = psi[i];
                    energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                            &edrvr, psi1, &psi2, delta, r, &r2);
                    norm = sqrt(psi2);
                    for (i = 0; i < D; i++) psi1[i] /= norm;
                    for (k = 0; k < N; k++) fj[k] = qaa_sigma_z(psi1, k);

                    qaa_update_diagonals(j, -DH(h0), d);
                    for (i = 0; i < D; i++) psi1[i] = psi[i];
                    energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                            &edrvr, psi1, &psi2, delta, r, &r2);
                    norm = sqrt(psi2);
                    for (i = 0; i < D; i++) psi1[i] /= norm;
                    for (k = 0; k < N; k++)
                        fj[k] = (fj[k] - qaa_sigma_z(psi1, k)) / DH(h0);

                    for (k = 0; k < N; k++) f += fj[k] * fj[k];
                    qaa_update_diagonals(j, 0.5 * DH(h0), d);
                }

                f = sqrt(f) / N;

                printf("%12s %12s %12g %12g " \
                       "%12d %12g %12g %12g %12g %12g %12g\n",
                        VERSION, basename(p.files[ifile]), h0, s,
                        iter, r2, energy, mz, mx, q2, f);
            }
        }
    }

    gsl_rng_free(rng);
    return EXIT_SUCCESS;
}