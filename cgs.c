#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdlib.h>
#include <libgen.h>
#include "global.h"
#include "graph.h"
#include "params.h"
#include "qaa.h"

double d[D];        /* diagonal elements of problem hamiltonian */
double psi_0[D];    /* ground state wavefunction for h_j = h_0 */
double psi[D];      /* for h_j = h_0 +/- dh/2 */
double delta[D];    /* CG search direction */
double r[D];        /* residual */

int main(int argc, char *argv[])
{
    const double sqrt_D = sqrt(D);
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *rng;
    params_t p;
    double fj[N];
    double edrvr, energy, h, mx, mz, norm, psi2, q2, q2p, r2, r2_0;
    int c, ih, is, iter, iter_0, j, k;
    UINT i;

    c = params_from_cmd(&p, argc, argv);

    if (c > PARAMS_ERR)
    {
        fprintf(stderr, "%s: %s\n", argv[0], params_errmsg[c]);
        exit(EXIT_FAILURE);
    }

    if (c > PARAMS_SUC)
    {
        fprintf(stderr, "%s\n", params_usage);
        exit(EXIT_SUCCESS);
    }
    rng = gsl_rng_alloc(T);

    printf("%16s %16s %16s %16s %16s %16s " \
           "%16s %16s %16s %16s %16s %16s %16s\n",
            "ver", "N", "graph", "dh", "h", "s",
            "iterations", "resid", "energy", "mz", "mx", "q2", "q2p");

    /* encode graph into diagonal elements of hamiltonian */
    qaa_compute_diagonals(p.graph.b, d);
    h = 0.;

    /* generate random initial wavefunction */
    gsl_rng_env_setup();
    for (i = 0; i < D; i++) psi_0[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;

    for (ih = 0; ih < p.nh; ih++)
    {
        qaa_update_diagonals(p.h[ih] - h, d);
        h = p.h[ih];

        for (is = 0; is < p.ns; is++)
        {
            energy = qaa_minimize_energy(p.s[is], d, p.eps, p.itermax, &iter_0,
                    &edrvr, psi_0, &psi2, r, &r2_0, delta);

            r2_0 /= psi2;
            norm = sqrt(psi2);
            for (i = 0; i < D; i++) psi_0[i] /= norm;

            mz = qaa_mag_z(psi_0);
            mx = 2. * qaa_me_driver(psi_0, psi_0);
            q2 = qaa_overlap(psi_0);

            q2p = 0.;

            /* vary the field at each site to approximate susceptibilities */
            for (j = 0; j < N; j++)
            {
                /* h_j = h_0 + dh/2 */
                /* use solution at midpoint as initial guess */
                for (i = 0; i < D; i++) psi[i] = psi_0[i];
                qaa_update_diagonals_1(j, 0.5 * p.dh, d);
                qaa_minimize_energy(p.s[is], d, p.eps, p.itermax, &iter,
                        &edrvr, psi, &psi2, r, &r2, delta);
                norm = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= norm;
                for (k = 0; k < N; k++) fj[k] = qaa_sigma_z(psi, k);

                /* h_j = h_0 - dh/2 */
                for (i = 0; i < D; i++) psi[i] = psi_0[i];
                qaa_update_diagonals_1(j, -p.dh, d);
                qaa_minimize_energy(p.s[is], d, p.eps, p.itermax, &iter,
                        &edrvr, psi, &psi2, r, &r2, delta);
                norm = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= norm;
                for (k = 0; k < N; k++)
                    fj[k] = (fj[k] - qaa_sigma_z(psi, k)) / p.dh;

                for (k = 0; k < N; k++) q2p += fj[k] * fj[k];
                qaa_update_diagonals_1(j, 0.5 * p.dh, d);
            }

            q2p = sqrt(q2p) / N;

            printf("%16s %16d %16s %16g %16g %16g " \
                   "%16d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                    VERSION, N, basename(p.fname), p.dh, p.h[ih], p.s[is],
                    iter_0, sqrt(r2_0), energy, mz, mx, q2, q2p);
        }
    }

    gsl_rng_free(rng);
    params_free(&p);
    return EXIT_SUCCESS;
}
