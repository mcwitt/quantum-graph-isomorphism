#include "graph.h"

enum
{
    PARAMS_SUC = 0,
    PARAMS_SUC_USAGE,
    PARAMS_ERR,
    PARAMS_ERR_LOAD,
    PARAMS_ERR_USAGE,
    PARAMS_NUM_MSG 
};

typedef struct
{
    double *s;
    double *h;
    double smin, smax;
    double emin, emax;
    int ns, nh;
    double dh;
    int itermax;
    double eps;
    char hexs[GRAPH_BITS_LEN / 4];
    char *file;
    graph_t graph;
} params_t;

extern char *params_usage;
extern char *params_errmsg[PARAMS_NUM_MSG];

void params_defaults(params_t *p);
int params_from_cmd(params_t *p, int argc, char *argv[]);
void params_free(params_t *p);
