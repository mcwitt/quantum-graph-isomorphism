#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
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
    double edrvr, energy, h, mx, mz, norm, psi2, q2, r2, s;
    double ds = 0., mh = 0., hprev = 0.;
    int ifile, ih, is, iter;
    UINT i;

    params_from_cmd(&p, argc, argv);

    if (p.num_files == 0)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (p.ns > 1) ds = (p.smax - p.smin)/(p.ns - 1.);
    if (p.nh > 1) mh = pow(10., (p.emax - p.emin)/(p.nh - 1.)) ;

    rng = gsl_rng_alloc(T);

    printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
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

        /* generate random initial wavefunction */
        gsl_rng_env_setup();
        for (i = 0; i < D; i++) psi[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;
        h = pow(10., p.emin);
        s = p.smin;

        for (ih = 0; ih < p.nh; ih++, hprev = h, h *= mh)
        {
            qaa_update_diagonals(h - hprev, d);

            for (is = 0; is < p.ns; is++, s += ds)
            {
                energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                        &edrvr, psi, &psi2, delta, r, &r2);

                /* normalize wavefunction */
                norm = sqrt(psi2);
                for (i = 0; i < D; i++) psi[i] /= norm;

                mz = qaa_mag_z(psi);
                mx = qaa_mag_x(psi);
                q2 = qaa_overlap(psi);

                printf("%12s %12s %12g %12g " \
                       "%12d %12g %12g %12g %12g %12g\n",
                        VERSION, basename(p.files[ifile]), h, s,
                        iter, r2, energy, mz, mx, q2);
            }

            ds = -ds; s += ds;  /* reverse s scan direction for next row*/
        }
    }

    gsl_rng_free(rng);
    return EXIT_SUCCESS;
}
