#include <stdlib.h>
#include <stdint.h>
#include "amatrix.h"
#include "qgi.h"
#include "params.h"

#define DS          0.01
#define SMIN        DS
#define SMAX        1.

#define MAX_ITER    1000
#define EPS         1e-16

qgi_t qgi;

int main(int argc, char *argv[])
{
    char **file_names;
    double h[N];    /* fields */
    double s, energy, mx, mz, q2;
    index_t i;
    int a[N][N], ifile, ih, iter, j;

    int nh = 5;
    double h0[] = {0.25, 0.5, 1., 2., 4.};

    if (argc < 2)
    {
        fprintf(stderr, "%s: no input files\n", argv[0]);
        return EXIT_FAILURE;
    }

    file_names = &argv[1];

    printf("%16s %6s %6s %12s %12s %12s %12s %12s\n",
            "file", "h0", "s", "iterations", "energy", "mz2", "mx2", "q2");

    for (ifile = 0; ifile < argc - 1; ifile++)
    {
        if (! amatrix_load(file_names[ifile], N, &a[0][0]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], file_names[ifile]);

            return EXIT_FAILURE;
        }

        for (i = 0; i < D; i++) qgi.psi[i] = 1.;

        for (ih = 0; ih < nh; ih++)
        {
            for (j = 0; j < N; j++) h[j] = h0[ih];
            qgi_init(&qgi, a, h);
                
            for (s = SMIN; s < SMAX; s += DS)
            {
                energy = qgi_minimize_energy(&qgi, s, EPS, &iter);
                mz = qgi_mag_z(qgi.psi);
                mx = qgi_mag_x(qgi.psi);
                q2 = qgi_overlap(qgi.psi);

                printf("%16s %6.3f %6.3f %12d %12g %12g %12g %12g\n",
                        file_names[ifile], h0[ih], s, iter, energy, mz*mz, mx*mx, q2);
            }
        }
    }

    return EXIT_SUCCESS;
}
