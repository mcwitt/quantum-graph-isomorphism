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
double psi[D];      /* wavefunction */
double delta[D];    /* CG search direction */
double r[D];        /* residual */

int main(int argc, char *argv[])
{
    const double sqrt_D = sqrt(D);
    const gsl_rng_type *T = gsl_rng_default;
    graph_t g;
    gsl_rng *rng;
    params_t p;
    double edrvr, energy, h, mx, mz, norm, psi2, q2, r2;
    int ifile, ih, is, iter;
    UINT i;

    params_from_cmd(&p, argc, argv);

    if (p.num_files == 0)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    rng = gsl_rng_alloc(T);

    printf("%16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n",
            "ver", "graph", "h", "s",
            "iterations", "res2", "energy", "mz", "mx", "q2");

    for (ifile = 0; ifile < p.num_files; ifile++)
    {
        if (! graph_read_amatrix(&g, p.files[ifile]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], p.files[ifile]);

            return EXIT_FAILURE;
        }
        
        /* encode graph into diagonal elements of hamiltonian */
        qaa_compute_diagonals(g.b, d);
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
                        &edrvr, psi, &psi2, delta, r, &r2);

                /* normalize wavefunction */
                r2 /= psi2;
                norm = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= norm;

                mz = qaa_mag_z(psi);
                mx = qaa_mag_x(psi);
                q2 = qaa_overlap(psi);

                printf("%16s %16s %16g %16g " \
                       "%16d %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                        VERSION, basename(p.files[ifile]), p.h[ih], p.s[is],
                        iter, r2, energy, mz, mx, q2);
            }
        }
    }

    gsl_rng_free(rng);
    params_free(&p);
    return EXIT_SUCCESS;
}
