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
double psi[D];
double h, energy;
int iter;
UINT i;

void write_output()
{
    double norm, mx, mz, q2;

    /* normalize wavefunction */
    norm = sqrt(qaa.cg.x2);
    for (i = 0; i < D; i++) psi[i] /= norm;

    mz = qaa_mag_z(psi);
    mx = 2. / N * qaa_me_driver(psi, psi);
    q2 = qaa_overlap(psi);

    printf("%8s %6d %8g %16s %16g %6g %6d " \
           "%16.9e %16.9e %16.9e %16.9e %16.9e\n",
            VERSION, N, p.tol, basename(p.file), h, qaa.s, iter,
            sqrt(qaa.cg.r2) / norm, energy, mz, mx, q2);
}

int main(int argc, char *argv[])
{
    const double sqrt_D = sqrt(D);
    const gsl_rng_type *T = gsl_rng_default;
    double tol2;
    gsl_rng *rng;
    int c, ih, is;

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

    tol2 = p.tol * p.tol;
    rng = gsl_rng_alloc(T);

    /* generate random initial wavefunction */
    gsl_rng_env_setup();
    for (i = 0; i < D; i++) psi[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;

    qaa_init(&qaa, p.graph.b, psi);
    h = 0.;

    printf("%8s %6s %8s %16s %16s %6s %6s "\
           "%16s %16s %16s %16s %16s\n",
            "ver", "N", "tol", "graph", "h", "s", "iter",
            "resid", "energy", "mz", "mx", "q2");

    for (ih = 0; ih < p.nh; ih++)
    {
        qaa_shift_field(&qaa, p.h[ih] - h);
        h = p.h[ih];

        for (qaa.s=p.s[is=0]; is < p.ns; qaa.s=p.s[++is], qaa_reset(&qaa, psi))
        {
            if (p.fullout)
            {
                for (iter = 0; iter < p.itermax; iter++)
                {
                    if ((qaa.cg.r2 / qaa.cg.x2) < tol2) break;
                    energy = qaa_iterate(&qaa, psi);
                    write_output();
                }
            }
            else
            {
                for (iter = 0; iter < p.itermax; iter++)
                {
                    if ((qaa.cg.r2 / qaa.cg.x2) < tol2) break;
                    energy = qaa_iterate(&qaa, psi);
                }

                write_output();
            }
        }
    }

    gsl_rng_free(rng);
    params_free(&p);

    return EXIT_SUCCESS;
}
