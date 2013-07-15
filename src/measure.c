#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include "amatrix.h"
#include "qgi.h"
#include "nlcg.h"
#include "params.h"
#include "time.h"

#define DEFAULT_INPUT_SVALS     "s_vals.in"
#define DEFAULT_OUTPUT_SUSC     "susc.out"

#define DEFAULT_FIELD       1.
#define DEFAULT_MAX_ITER    1000
#define DEFAULT_TOL         1e-16

#define MAX_S_VALS  1000

double d[D];
double s[MAX_S_VALS];
int s_index;
double psi2, tmp;

/* wrappers for gradient and line minimization functions
 * to be passed to NLCG routine */

void gradient(double *psi, double *grad) {
    qgi_energy_grad(s[s_index], d, psi, grad, &psi2, &tmp);
}

double line_min(double *psi, double *delta) {
    return qgi_line_min(s[s_index], d, psi, delta);
}

void process_graph(a, h, s)
{
    qgi_compute_problem_hamiltonian(a, h, d);

    for (s_index = 0; s_index < num_s_vals; s_index++)
    {
        for (i = 0; i < D; i++) psi[i] = rand() / (RAND_MAX + 1.) * 2. - 1.;
        iter = nlcg_minimize(gradient, line_min, max_iter, etol, psi);
        energy = qgi_energy_grad(s[s_index], d, psi, grad, &eod);

        printf("%12s %6g %12d %16.9g %16.9g\n",
                file_names[ifile], s[s_index], iter, energy, eod);
    }
}

int main(int argc, char *argv[])
{
    extern char *optarg;
    extern int optind;
    int c, err = 0;

    FILE *fp;

    char *input_svals, *output_susc, **file_names;
    double h0, etol;
    int max_iter;

    double h[N];
    double psi[D], grad[D];
    double energy, eod;
    index_t i;
    int a[N][N], ifile, s_index, iter, j, num_s_vals, num_files;

    input_svals = DEFAULT_INPUT_SVALS;
    output_susc = DEFAULT_OUTPUT_SUSC;
    h0          = DEFAULT_FIELD;
    max_iter    = DEFAULT_MAX_ITER;
    etol        = DEFAULT_TOL;

    while ((c = getopt(argc, argv, "e:h:i:o:s:")) != -1)
    {
        switch (c)
        {
            case 'e':
                etol = atof(optarg);
                break;
            case 'h':
                h0 = atof(optarg);
                break;
            case 'i':
                max_iter = atoi(optarg);
                break;
            case 'o':
                output_susc = optarg;
                break;
            case 's':
                input_svals = optarg;
                break;
            case '?':
                err = 1;
                break;
        }
    }

    if (err || (optind == argc))
    {
        fprintf(stderr,
                "Usage: %s [-ehios value] file [file ...]\n" \
                "  -e   error tolerance\n" \
                "  -h   external field\n" \
                "  -i   max CG iterations\n" \
                "  -o   output file for susceptibilities\n" \
                "  -s   file with s values\n",
                argv[0]);

        return EXIT_FAILURE;
    }

    srand(time(NULL));
    /*srand(123);*/

    file_names = &argv[optind];
    num_files = argc - optind;

    if ((fp = fopen(input_svals, "r")) == NULL)
    {
        fprintf(stderr, "%s: error opening file %s\n", argv[0], input_svals);
        return EXIT_FAILURE;
    }

    /* read s-values from file */
    for (num_s_vals = 0; num_s_vals < MAX_S_VALS; num_s_vals++)
        if (fscanf(fp, "%lf", &s[num_s_vals]) == EOF) break;

    printf("#%11s %6s %12s %16s %16s\n", "file", "s", "iterations", "E_0", "E_0(offd)");

    for (ifile = 0; ifile < num_files; ifile++)
    {
        if (! amatrix_load(file_names[ifile], N, &a[0][0]))
        {
            fprintf(stderr, "%s: error loading adjacency matrix from file \"%s\"\n",
                    argv[0], file_names[ifile]);

            return EXIT_FAILURE;
        }

        for (j = 0; j < N; j++) h[j] = h0;

    }

    return EXIT_SUCCESS;
}
