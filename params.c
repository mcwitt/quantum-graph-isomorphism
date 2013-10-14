#include "params.h"
#include "global.h"
#include "graph.h"
#include "range.h"
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char *params_errmsg[] = {
    [PARAMS_ERR_LOAD]   = "couldn't read input file",
    [PARAMS_ERR_USAGE]  = "incorrect usage (use -h for info)"
};

char *params_usage = \
        "Conjugate gradient energy minimization.\n" \
        "Usage: %s [-option arg ...] " \
        "[(-f amatrix_file | -b bit_string | -x hex_string)]\n" \
        "Reads graph in bit-string format from stdin if " \
        "-f, -b or -x not specified.\n" \
        "Options:\n" \
        "  -s, --smin    : minimum value of adiabatic parameter s\n" \
        "  -S, --smax    : max value of s\n" \
        "  -n, --ns      : number of s values\n" \
        "  -e, --emin    : log_10 of minimum magnetic field\n" \
        "  -E, --emax    : log_10 of maximum magnetic field\n" \
        "  -m, --nh      : number of magnetic field values\n" \
        "  -d, --dh      : delta for finite differences\n" \
        "  -i, --itermax : maxiumum number of CG iterations\n" \
        "  -t, --tol     : error tolerance\n";

void params_defaults(params_t *p)
{
    p->s = NULL;
    p->h = NULL;
    p->smin = 0.02;
    p->smax = 0.98;
    p->emin = -2.;
    p->emax = 1.;
    p->ns = 49;
    p->nh = 52;
    p->dh = 1e-3;
    p->itermax = 300;
    p->eps = 1e-12;
    p->file = "";
}

enum read_mode { READ_FILE, READ_BITS, READ_HEXS };

static int read_file(graph_t *g, char *file)
{
    FILE *fp;

    fp = fopen(file, "r");
    if (fp == NULL || graph_read_amatrix(g, fp) != 0) return 1;
    fclose(fp);
    return 0;
}

static int read_stdin(graph_t *g)
{
    char bits[GRAPH_BITS_LEN+1];

    fgets(bits, GRAPH_BITS_LEN+1, stdin);
    return graph_read_bits(g, bits);
}

int params_from_cmd(params_t *p, int argc, char *argv[])
{
    extern char *optarg;
    extern int optind;
    enum read_mode mode = READ_BITS;
    char *arg;
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

    while ((c = getopt_long(argc, argv, "fxs:S:n:e:E:m:d:i:t:h",
                    long_options, &long_index)) != -1)
    {
        switch (c)
        {
            case 'f' : mode       = READ_FILE; break;
            case 'x' : mode       = READ_HEXS; break;
            case 's' : p->smin    = atof(optarg); break;
            case 'S' : p->smax    = atof(optarg); break;
            case 'n' : p->ns      = atoi(optarg); break;
            case 'e' : p->emin    = atof(optarg); break;
            case 'E' : p->emax    = atof(optarg); break;
            case 'm' : p->nh      = atoi(optarg); break;
            case 'd' : p->dh      = atof(optarg); break;
            case 'i' : p->itermax = atoi(optarg); break;
            case 't' : p->eps     = atof(optarg); break;
            case 'h' : return PARAMS_SUC_USAGE;   break;
            default  : return PARAMS_ERR_USAGE;
        }
    }

    p->s = linspace(p->smin, p->smax, p->ns);
    p->h = logspace(p->emin, p->emax, p->nh);

    if (argc > optind)  /* graph specified in command */
    {
        arg = argv[optind];

        switch(mode)
        {
            case READ_BITS:
                if (graph_read_bits(&p->graph, arg) != 0)
                    return PARAMS_ERR_LOAD;
                graph_to_hexs(&p->graph, p->hexs);
                break;
            case READ_HEXS:
                if (graph_read_hexs(&p->graph, arg) != 0)
                    return PARAMS_ERR_LOAD;
                strcpy(p->hexs, arg);
                break;
            case READ_FILE:
                if (read_file(&p->graph, arg) != 0) 
                   return PARAMS_ERR_LOAD;
                graph_to_hexs(&p->graph, p->hexs);
                p->file = arg;
        }
    }
    else if (read_stdin(&p->graph) != 0) return PARAMS_ERR_LOAD;

    return PARAMS_SUC;
}

void params_free(params_t *p)
{
    free(p->s);
    free(p->h);
}
