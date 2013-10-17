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

double d[D];        /* diagonal elements of problem Hamiltonian */
double psi[D];      /* wavefunction */
nlcg_t nlcg;
qaa_args_t args;    /* QAA-specific parameters */
params_t p;
double h, s, edrvr, energy;
int iter;
UINT i;

void write_output()
{
    double norm, mx, mz, q2;

    /* normalize wavefunction */
    norm = sqrt(nlcg.x2);
    for (i = 0; i < D; i++) psi[i] /= norm;

    mz = qaa_mag_z(psi);
    mx = 2. / N * qaa_me_driver(psi, psi);
    q2 = qaa_overlap(psi);

    printf("%8s %6d %8g %16s %16g %6g %6d " \
           "%16.9e %16.9e %16.9e %16.9e %16.9e\n",
            VERSION, N, p.tol, basename(p.file), h, s, iter,
            sqrt(nlcg.r2) / norm, energy, mz, mx, q2);
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

    printf("%8s %6s %8s %16s %16s %6s %6s "\
           "%16s %16s %16s %16s %16s\n",
            "ver", "N", "tol", "graph", "h", "s", "iter",
            "resid", "energy", "mz", "mx", "q2");
        
    /* encode graph into diagonal elements of hamiltonian */
    qaa_compute_diagonals(p.graph.b, d);
    h = 0.;

    /* generate random initial wavefunction */
    gsl_rng_env_setup();
    for (i = 0; i < D; i++) psi[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;

    for (h=p.h[ih=0]; ih < p.nh; h=p.h[++ih])
    {
        qaa_update_diagonals(p.h[ih] - h, d);
        h = p.h[ih];

        for (s=p.s[is=0]; is < p.ns; s=p.s[++is])
        {
            qaa_nlcg_init(p.s[is], d, psi, &edrvr, &args, &nlcg);

            if (p.fullout)
            {
                for (iter = 0; iter < p.itermax; iter++)
                {
                    if ((nlcg.r2 / nlcg.x2) < tol2) break;
                    energy = nlcg_iterate(&nlcg);
                    write_output();
                }
            }
            else
            {
                for (iter = 0; iter < p.itermax; iter++)
                {
                    if ((nlcg.r2 / nlcg.x2) < tol2) break;
                    energy = nlcg_iterate(&nlcg);
                }

                write_output();
            }
        }
    }

    gsl_rng_free(rng);
    params_free(&p);

    return EXIT_SUCCESS;
}
