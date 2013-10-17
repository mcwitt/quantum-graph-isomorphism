#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdlib.h>
#include <libgen.h>
#include "global.h"
#include "graph.h"
#include "nlcg.h"
#include "params.h"
#include "qaa.h"

double d[D];        /* diagonal elements of problem hamiltonian */
double psi_0[D];    /* ground state wavefunction for h_j = h_0 */
double psi[D];      /* for h_j = h_0 +/- dh/2 */
nlcg_t nlcg;
qaa_args_t args;    /* QAA-specific parameters */

int main(int argc, char *argv[])
{
    const double sqrt_D = sqrt(D);
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *rng;
    params_t p;
    double fj[N];
    double edrvr, energy, h, mx, mz, norm, q2, q2p, r2_0;
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

    printf("%8s %6s %8s %16s %16s %6s %6s "\
           "%16s %16s %16s %16s %16s %16s\n",
            "ver", "N", "tol", "graph", "h", "s", "iter",
            "resid", "energy", "mz", "mx", "q2", "q2p");

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
            qaa_nlcg_init(p.s[is], d, psi_0, &edrvr, &args, &nlcg);
            energy = nlcg_minimize(&nlcg, p.tol, p.itermax, &iter_0);

            r2_0 = nlcg.r2 / nlcg.x2;
            norm = sqrt(nlcg.x2);
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
                qaa_update_diagonals_1(j, 0.5 * p.dh, d);
                qaa_nlcg_init(p.s[is], d, psi, &edrvr, &args, &nlcg);
                energy = nlcg_minimize(&nlcg, p.tol, p.itermax, &iter);
                norm = sqrt(nlcg.x2);
                for (i = 0; i < D; i++) psi[i] /= norm;
                for (k = 0; k < N; k++) fj[k] = qaa_sigma_z(psi, k);

                /* h_j = h_0 - dh/2 */
                for (i = 0; i < D; i++) psi[i] = psi_0[i];
                qaa_update_diagonals_1(j, -p.dh, d);
                qaa_nlcg_init(p.s[is], d, psi, &edrvr, &args, &nlcg);
                energy = nlcg_minimize(&nlcg, p.tol, p.itermax, &iter);
                norm = sqrt(nlcg.x2);
                for (i = 0; i < D; i++) psi[i] /= norm;
                for (k = 0; k < N; k++)
                    fj[k] = (fj[k] - qaa_sigma_z(psi, k)) / p.dh;

                for (k = 0; k < N; k++) q2p += fj[k] * fj[k];
                qaa_update_diagonals_1(j, 0.5 * p.dh, d);
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
