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

#define DH(h)   1e-3

double d[D];        /* diagonal elements of problem hamiltonian */
double psi0[D];     /* ground state wavefunction for h_j = h_0 */
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
    double edrvr, energy, h, mx, mz, norm, psi2, q2, q2p, r2, s;
    double ds = 0., mh = 0., hprev = 0.;
    int a[N*(N-1)/2], ifile, ih, is, iter, j, k;
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

    printf("%12s %12s %12s %12s %12s " \
           "%12s %12s %12s %12s %12s %12s %12s\n",
            "ver", "graph", "dh", "h", "s",
            "iterations", "res2", "energy", "mz", "mx", "q2", "q2p");

    for (ifile = 0; ifile < p.num_files; ifile++)
    {
        if (! amatrix_load(p.files[ifile], a))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], p.files[ifile]);

            return EXIT_FAILURE;
        }
        
        /* encode graph into diagonal elements of hamiltonian */
        qaa_compute_diagonals(a, d);

        /* generate random initial wavefunction */
        gsl_rng_env_setup();
        for (i = 0; i < D; i++) psi0[i] = gsl_ran_gaussian(rng, 1.) / sqrt_D;
        h = pow(10., p.emin);

        for (ih = 0; ih < p.nh; ih++, hprev = h, h *= mh)
        {
            qaa_update_diagonals(h - hprev, d);
            s = p.smin;

            for (is = 0; is < p.ns; is++, s += ds)
            {
                energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                        &edrvr, psi0, &psi2, delta, r, &r2);

                norm = sqrt(psi2);
                for (i = 0; i < D; i++) psi0[i] /= norm;

                mz = qaa_mag_z(psi0);
                mx = qaa_mag_x(psi0);
                q2 = qaa_overlap(psi0);

                q2p = 0.;

                /* vary the field at each site to approximate susceptibilities */
                for (j = 0; j < N; j++)
                {
                    /* h_j = h_0 + dh/2 */
                    qaa_update_diagonals_1(j, 0.5 * DH(h), d);
                    for (i = 0; i < D; i++) psi[i] = psi0[i];
                    energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                            &edrvr, psi, &psi2, delta, r, &r2);
                    norm = sqrt(psi2);
                    for (i = 0; i < D; i++) psi[i] /= norm;
                    for (k = 0; k < N; k++) fj[k] = qaa_sigma_z(psi, k);

                    /* h_j = h_0 - dh/2 */
                    qaa_update_diagonals_1(j, -DH(h), d);
                    for (i = 0; i < D; i++) psi[i] = psi0[i];
                    energy = qaa_minimize_energy(s, d, p.eps, p.itermax, &iter,
                            &edrvr, psi, &psi2, delta, r, &r2);
                    norm = sqrt(psi2);
                    for (i = 0; i < D; i++) psi[i] /= norm;
                    for (k = 0; k < N; k++)
                        fj[k] = (fj[k] - qaa_sigma_z(psi, k)) / DH(h);

                    for (k = 0; k < N; k++) q2p += fj[k] * fj[k];
                    qaa_update_diagonals_1(j, 0.5 * DH(h), d);
                }

                q2p = sqrt(q2p) / N;

                printf("%12s %12s %12g %12g %12g " \
                       "%12d %12g %12g %12g %12g %12g %12g\n",
                        VERSION, basename(p.files[ifile]), DH(h), h, s,
                        iter, r2, energy, mz, mx, q2, q2p);
            }
        }
    }

    gsl_rng_free(rng);
    return EXIT_SUCCESS;
}
