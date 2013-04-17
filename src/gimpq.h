#include "params.h"

typedef struct
{
    double d[D];    /* diagonal elements of H_P */
    double s;       /* adiabatic parameter */
} GQProblem;

void gq_init(GQProblem *p, int g[N][N], double h[N]);

void gq_grad_energy(gimpq_t *g, double c[D], double grad[D]);
