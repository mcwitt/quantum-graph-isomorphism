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
double psi[D];      /* wavefunction */
double delta[D];    /* CG search direction */
double r[D];        /* residual */

int main(int argc, char *argv[])
{
    const double sqrt_D = sqrt(D);
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *rng;
    params_t p;
    double edrvr, energy, h, mx, mz, norm, psi2, q2, r2;
    int c, ih, is, iter;
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

    printf("%16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n",
            "ver", "N", "graph", "h", "s",
            "iterations", "resid", "energy", "mz", "mx", "q2");
        
    /* encode graph into diagonal elements of hamiltonian */
    qaa_compute_diagonals(p.graph.b, d);
    h = 0.;

    /* generate random initial wavefunction */
    gsl_rng_env_setup();
    for (i = 0; i < D; i++) psi[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;

    for (ih = 0; ih < p.nh; ih++)
    {
        qaa_update_diagonals(p.h[ih] - h, d);
        h = p.h[ih];

        for (is = 0; is < p.ns; is++)
        {
            energy = qaa_minimize_energy(p.s[is], d, p.eps, p.itermax, &iter,
                    &edrvr, psi, &psi2, r, &r2, delta);

            /* normalize wavefunction */
            r2 /= psi2;
            norm = sqrt(psi2);
            for (i = 0; i < D; i++) psi[i] /= norm;

            mz = qaa_mag_z(psi);
            mx = 2. * qaa_me_driver(psi, psi);
            q2 = qaa_overlap(psi);

            printf("%16s %16d %16s %16g %16g " \
                   "%16d %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                    VERSION, N, basename(p.fname), p.h[ih], p.s[is],
                    iter, sqrt(r2), energy, mz, mx, q2);
        }
    }

    gsl_rng_free(rng);
    params_free(&p);
    return EXIT_SUCCESS;
}
