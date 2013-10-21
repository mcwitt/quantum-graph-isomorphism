#include "params.h"
#include "global.h"
#include "graph.h"
#include "range.h"
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char *params_errmsg[] = {
    [PARAMS_ERR_LOAD]   = "couldn't read input file",
    [PARAMS_ERR_USAGE]  = "incorrect usage (use -h for info)"
};

char *params_usage = \
        "Conjugate gradient energy minimization.\n" \
        "Usage: %s [-option arg ...]\n" \
        "Reads graph in bit-string format from stdin if " \
        "-f, -b or -x not specified.\n" \
        "Flags:\n" \
        "  -v, --fullout    show full output including all iterations\n"
        "Parameters:\n" \
        "  -f, --file       read adjacency matrix from file\n" \
        "  -b, --bits       read graph in bit-string format \n" \
        "  -x, --hexs       read graph in hex-string format \n" \
        "  -s, --smin       set minimum value of adiabatic parameter s\n" \
        "  -S, --smax       set max value of s\n" \
        "  -n, --ns         set number of s values\n" \
        "  -e, --emin       set log_10 of minimum magnetic field\n" \
        "  -E, --emax       set log_10 of maximum magnetic field\n" \
        "  -m, --nh         set number of magnetic field values\n" \
        "  -d, --dh         set delta for finite differences\n" \
        "  -c, --cutoff     set maxiumum number of CG iterations\n" \
        "  -t, --tol        set error tolerance\n";

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
    p->tol = 1e-12;
    p->fullout = 0;
    p->file = NULL;
}

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
    int c, long_index, is_graph_loaded = 0;

    static struct option long_options[] = {
        {"file",    required_argument, 0, 'f'},
        {"bits",    required_argument, 0, 'b'},
        {"hexs",    required_argument, 0, 'x'},
        {"smin",    required_argument, 0, 's'},
        {"smax",    required_argument, 0, 'S'},
        {"ns",      required_argument, 0, 'n'},
        {"emin",    required_argument, 0, 'e'},
        {"emax",    required_argument, 0, 'E'},
        {"nh",      required_argument, 0, 'm'},
        {"dh",      required_argument, 0, 'd'},
        {"cutoff",  required_argument, 0, 'c'},
        {"tol",     required_argument, 0, 't'},
        {"fullout", no_argument,       0, 'v'},
        {"help",    no_argument,       0, 'h'},
        {NULL, 0, 0, 0}
    };

    params_defaults(p);

    while ((c = getopt_long(argc, argv, "f:b:x:s:S:n:e:E:m:d:c:t:vh",
                    long_options, &long_index)) != -1)
    {
        switch (c)
        {
            case 'f': 
                if (read_file(&p->graph, optarg) != 0)
                    return PARAMS_ERR_LOAD;
                graph_to_hexs(&p->graph, p->hexs);
                p->file = optarg;
                is_graph_loaded = 1;
                break;
            case 'b':
                if (graph_read_bits(&p->graph, optarg) != 0)
                    return PARAMS_ERR_LOAD;
                graph_to_hexs(&p->graph, p->hexs);
                is_graph_loaded = 1;
                break;
            case 'x':
                if (graph_read_hexs(&p->graph, optarg) != 0)
                    return PARAMS_ERR_LOAD;
                strcpy(p->hexs, optarg);
                is_graph_loaded = 1;
                break;
            case 's' : p->smin    = atof(optarg); break;
            case 'S' : p->smax    = atof(optarg); break;
            case 'n' : p->ns      = atoi(optarg); break;
            case 'e' : p->emin    = atof(optarg); break;
            case 'E' : p->emax    = atof(optarg); break;
            case 'm' : p->nh      = atoi(optarg); break;
            case 'd' : p->dh      = atof(optarg); break;
            case 'c' : p->itermax = atoi(optarg); break;
            case 't' : p->tol     = atof(optarg); break;
            case 'v' : p->fullout = 1;            break;
            case 'h' : return PARAMS_SUC_USAGE;   break;
            default  : return PARAMS_ERR_USAGE;
        }
    }

    if ((! is_graph_loaded) && (read_stdin(&p->graph) != 0))
        return PARAMS_ERR_LOAD;

    p->s = linspace(p->smin, p->smax, p->ns);
    p->h = logspace(p->emin, p->emax, p->nh);

    return PARAMS_SUC;
}

void params_free(params_t *p)
{
    free(p->s);
    free(p->h);
}
