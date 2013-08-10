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
            "  -p, --pmin    : log_10 of minimum magnetic field\n" \
            "  -P, --pmax    : log_10 of maximum magnetic field\n" \
            "  -d, --np      : number of magnetic field values\n" \
            "  -i, --itermax : maxiumum number of CG iterations\n" \
            "  -e, --eps     : error tolerance\n",
            argv[0]);
}

void params_defaults(params_t *p)
{
    p->smin = 0.01;
    p->smax = 0.99;
    p->ns   = 99;

    p->pmin = -1.;
    p->pmax = 1.;
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
        {"pmin",    required_argument,  0,  'p' },
        {"pmax",    required_argument,  0,  'P' },
        {"np",      required_argument,  0,  'd' },
        {"itermax", required_argument,  0,  'i' },
        {"eps",     required_argument,  0,  'e' },
        {NULL, 0, 0, 0}
    };

    params_defaults(p);

    while ((c = getopt_long(argc, argv, "s:S:n:p:P:d:i:e:",
                    long_options, &long_index)) != -1)
    {
        switch (c)
        {
             case 's' : p->smin    = atof(optarg); break;
             case 'S' : p->smax    = atof(optarg); break;
             case 'n' : p->ns      = atoi(optarg); break;
             case 'p' : p->pmin    = atof(optarg); break;
             case 'P' : p->pmax    = atof(optarg); break;
             case 'd' : p->np      = atoi(optarg); break;
             case 'i' : p->itermax = atoi(optarg); break;
             case 'e' : p->eps     = atof(optarg); break;
             default: print_usage(argc, argv); exit(EXIT_FAILURE);
        }
    }

    p->num_files = argc - optind;
    p->files = argv + optind;
}
