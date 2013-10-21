#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdlib.h>
#include <libgen.h>
#include "global.h"
#include "graph.h"
#include "params.h"
#include "qaa.h"

params_t p;
qaa_t qaa;
double psi_0[D];    /* ground state wavefunction for h_j = h_0 */
double psi[D];      /* for h_j = h_0 +/- dh/2 */
UINT i;

int main(int argc, char *argv[])
{
    const double sqrt_D = sqrt(D);
    const gsl_rng_type *T = gsl_rng_default;
    double fj[N];
    double energy, h, mx, mz, norm, q2, q2p, r2_0;
    gsl_rng *rng;
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
        fprintf(stderr, params_usage, argv[0]);
        exit(EXIT_SUCCESS);
    }

    rng = gsl_rng_alloc(T);

    /* generate random initial wavefunction */
    gsl_rng_env_setup();
    for (i = 0; i < D; i++) psi_0[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;

    qaa_init(&qaa, p.graph.b, psi_0);
    h = 0.;

    printf("%8s %6s %8s %16s %16s %6s %6s "\
           "%16s %16s %16s %16s %16s %16s\n",
            "ver", "N", "tol", "graph", "h", "s", "iter",
            "resid", "energy", "mz", "mx", "q2", "q2p");

    for (ih = 0; ih < p.nh; ih++)
    {
        qaa_shift_field(&qaa, p.h[ih] - h);
        h = p.h[ih];

        for (qaa.s=p.s[is=0]; is < p.ns; qaa.s=p.s[++is], qaa_reset(&qaa, psi))
        {
            energy = qaa_minimize(&qaa, psi_0, p.tol, p.itermax, &iter_0);

            r2_0 = qaa.cg.r2 / qaa.cg.x2;
            norm = sqrt(qaa.cg.x2);
            for (i = 0; i < D; i++) psi_0[i] /= norm;

            mz = qaa_mag_z(psi_0);
            mx = 2. / N * qaa_me_driver(psi_0, psi_0);
            q2 = qaa_overlap(psi_0);

            q2p = 0.;

            /* vary the field at each site to approximate susceptibilities */
            for (j = 0; j < N; j++)
            {
                /* h_j = h_0 + dh/2 */
                /* use solution at midpoint as initial guess */
                for (i = 0; i < D; i++) psi[i] = psi_0[i];
                qaa_shift_field_1(&qaa, j, 0.5 * p.dh);
                qaa_reset(&qaa, psi);
                energy = qaa_minimize(&qaa, psi, p.tol, p.itermax, &iter);
                norm = sqrt(qaa.cg.x2);
                for (i = 0; i < D; i++) psi[i] /= norm;
                for (k = 0; k < N; k++) fj[k] = qaa_sigma_z(psi, k);

                /* h_j = h_0 - dh/2 */
                for (i = 0; i < D; i++) psi[i] = psi_0[i];
                qaa_shift_field_1(&qaa, j, -p.dh);
                qaa_reset(&qaa, psi);
                energy = qaa_minimize(&qaa, psi, p.tol, p.itermax, &iter);
                norm = sqrt(qaa.cg.x2);
                for (i = 0; i < D; i++) psi[i] /= norm;
                for (k = 0; k < N; k++)
                    fj[k] = (fj[k] - qaa_sigma_z(psi, k)) / p.dh;

                for (k = 0; k < N; k++) q2p += fj[k] * fj[k];
                qaa_shift_field_1(&qaa, j, 0.5 * p.dh);
            }

            q2p = sqrt(q2p) / N;

            printf("%8s %6d %8g %16s %16g %6g %6d " \
                   "%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                    VERSION, N, p.tol, basename(p.file), p.h[ih], p.s[is], iter_0,
                    sqrt(r2_0), energy, mz, mx, q2, q2p);
        }
    }

    gsl_rng_free(rng);
    params_free(&p);
    return EXIT_SUCCESS;
}
