#include "params.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>

static void print_usage(int argc, char *argv[])
{
    fprintf(stderr,
            "usage: %s [-option arg ...] file [file ...]\n" \
            "options:\n" \
            "  -s, --smin    : minimum value of adiabatic parameter s\n" \
            "  -S, --smax    : max value of s\n" \
            "  -n, --ns      : number of s values\n" \
            "  -e, --emin    : log_10 of minimum magnetic field\n" \
            "  -E, --emax    : log_10 of maximum magnetic field\n" \
            "  -m, --np      : number of magnetic field values\n" \
            "  -i, --itermax : maxiumum number of CG iterations\n" \
            "  -t, --tol     : error tolerance\n",
            argv[0]);
}

void params_defaults(params_t *p)
{
    p->smin = 0.01;
    p->smax = 0.99;
    p->ns   = 99;

    p->emin = -1.;
    p->emax = 1.;
    p->np   = 51;

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
        {"np",      required_argument,  0,  'm' },
        {"itermax", required_argument,  0,  'i' },
        {"tol",     required_argument,  0,  't' },
        {NULL, 0, 0, 0}
    };

    params_defaults(p);

    while ((c = getopt_long(argc, argv, "s:S:n:e:E:m:i:t:",
                    long_options, &long_index)) != -1)
    {
        switch (c)
        {
             case 's' : p->smin    = atof(optarg); break;
             case 'S' : p->smax    = atof(optarg); break;
             case 'n' : p->ns      = atoi(optarg); break;
             case 'e' : p->emin    = atof(optarg); break;
             case 'E' : p->emax    = atof(optarg); break;
             case 'm' : p->np      = atoi(optarg); break;
             case 'i' : p->itermax = atoi(optarg); break;
             case 't' : p->eps     = atof(optarg); break;
             default: print_usage(argc, argv); exit(EXIT_FAILURE);
        }
    }

    p->num_files = argc - optind;
    p->files = argv + optind;
}
