#include "params.h"

typedef struct
{
    double d[D];    /* diagonal elements of H_P */
    double s;       /* adiabatic parameter */
} GQProblem;

void gq_init(GQProblem *p, int g[N][N], double h[N]);

void gq_grad_energy(GQProblem *p, double s, double c[D], double grad[D]);

double gq_line_min(GQProblem *p, double s, double c[D], double delta[D]);

