#include "params.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>

static void print_usage(int argc, char *argv[])
{
    fprintf(stderr,
            "usage: %s [-option arg ...] file [file ...]\n" \
            "options:\n" \
            "  -s, --smin    : minimum value of adiabatic parameter\n" \
            "  -S, --smax    : max value\n" \
            "  -d, --ds      : step size\n" \
            "  -h, --hmin    : minimum value of magnetic field\n" \
            "  -n, --nh      : number of magnetic field values per decade\n" \
            "  -D, --ndec    : number of decades in magnetic field\n" \
            "  -i, --itermax : maxiumum number of CG iterations\n" \
            "  -e, --eps     : error tolerance\n",
            argv[0]);
}

void params_defaults(params_t *p)
{
    p->smin = 0.01;
    p->smax = 0.99;
    p->ds   = 0.01;

    p->hmin   = 0.1;
    p->nh_dec = 50;
    p->ndec   = 2;

    p->itermax = 300;
    p->eps = 1e-16;

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
        {"ds",      required_argument,  0,  'd' },
        {"hmin",    required_argument,  0,  'h' },
        {"nh",      required_argument,  0,  'n' },
        {"ndec",    required_argument,  0,  'D' },
        {"itermax", required_argument,  0,  'i' },
        {"eps",     required_argument,  0,  'e' }
    };

    params_defaults(p);

    while ((c = getopt_long(argc, argv, "s:S:d:h:n:D:i:e:",
                    long_options, &long_index)) != -1)
    {
        switch (c)
        {
             case 's' : p->smin    = atof(optarg); break;
             case 'S' : p->smax    = atof(optarg); break;
             case 'd' : p->ds      = atof(optarg); break;
             case 'h' : p->hmin    = atof(optarg); break;
             case 'n' : p->nh_dec  = atoi(optarg); break;
             case 'D' : p->ndec    = atoi(optarg); break;
             case 'i' : p->itermax = atoi(optarg); break;
             case 'e' : p->eps     = atof(optarg); break;
             default: print_usage(argc, argv); exit(EXIT_FAILURE);
        }
    }

    p->num_files = argc - optind;
    p->files = argv + optind;
}
