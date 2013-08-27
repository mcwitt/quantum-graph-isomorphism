#include "params.h"
#include "global.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static void print_usage(int argc, char *argv[])
{
    fprintf(stderr,
            "Conjugate gradient energy minimization " \
            "for graphs with N=%d vertices.\n" \
            "usage: %s [-option arg ...] file [file ...]\n" \
            "options:\n" \
            "  -s, --smin    : minimum value of adiabatic parameter s\n" \
            "  -S, --smax    : max value of s\n" \
            "  -n, --ns      : number of s values\n" \
            "  -e, --emin    : log_10 of minimum magnetic field\n" \
            "  -E, --emax    : log_10 of maximum magnetic field\n" \
            "  -m, --nh      : number of magnetic field values\n" \
            "  -d, --dh      : delta for finite differences\n" \
            "  -i, --itermax : maxiumum number of CG iterations\n" \
            "  -t, --tol     : error tolerance\n",
            N, argv[0]);
}

static double *linspace(double left, double right, int n)
{
    double *a, delta = (right - left) / (n - 1);
    int i;

    a = (double*) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) a[i] = left + i*delta;
    return a;
}

static double *logspace(double left, double right, int n)
{
    double *a;
    int i;

    a = linspace(left, right, n);
    for (i = 0; i < n; i++) a[i] = pow(10., a[i]);
    return a;
}

void params_defaults(params_t *p)
{
    p->s = NULL;
    p->h = NULL;
    p->smin = 0.02;
    p->smax = 0.98;
    p->emin = -2;
    p->emax = 1;
    p->ns = 49;
    p->nh = 51;
    p->dh = 1e-3;
    p->itermax = 300;
    p->eps = 1e-12;
    p->num_files = 0;
    p->files = NULL;
}

void params_from_cmd(params_t *p, int argc, char *argv[])
{
    extern char *optarg;
    extern int optind;
    int c, long_index;

    static struct option long_options[] = {
        {"smin",    required_argument,  0,  's' },
        {"smax",    required_argument,  0,  'S' },
        {"ns",      required_argument,  0,  'n' },
        {"emin",    required_argument,  0,  'e' },
        {"emax",    required_argument,  0,  'E' },
        {"nh",      required_argument,  0,  'm' },
        {"dh",      required_argument,  0,  'd' },
        {"itermax", required_argument,  0,  'i' },
        {"tol",     required_argument,  0,  't' },
        {"help",    no_argument,        0,  'h' },
        {NULL, 0, 0, 0}
    };

    params_defaults(p);

    while ((c = getopt_long(argc, argv, "s:S:n:e:E:m:d:i:t:h",
                    long_options, &long_index)) != -1)
    {
        switch (c)
        {
             case 's' : p->smin    = atof(optarg); break;
             case 'S' : p->smax    = atof(optarg); break;
             case 'n' : p->ns      = atoi(optarg); break;
             case 'e' : p->emin    = atof(optarg); break;
             case 'E' : p->emax    = atof(optarg); break;
             case 'm' : p->nh      = atoi(optarg); break;
             case 'd' : p->dh      = atof(optarg); break;
             case 'i' : p->itermax = atoi(optarg); break;
             case 't' : p->eps     = atof(optarg); break;
             default: print_usage(argc, argv); exit(EXIT_FAILURE);
        }
    }

    p->s = linspace(p->smin, p->smax, p->ns);
    p->h = logspace(p->emin, p->emax, p->nh);
    p->num_files = argc - optind;
    p->files = argv + optind;
}

void params_free(params_t *p)
{
    free(p->s);
    free(p->h);
}
