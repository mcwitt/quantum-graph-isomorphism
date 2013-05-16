#include <stdlib.h>
#include "amatrix.h"
#include "gimpq.h"
#include "nlcg.h"
#include "params.h"

GQProblem g;

/* wrappers for gradient and line minimization functions
 * to be passed to NLCG routine */
void gradient(double *x, double *grad) { gq_grad_energy(&g, x, grad); }
void line_min(double *x, double *d) { gq_line_min(&g, x, d); }

int main()
{
    double h[N];
    int a[N][N], j;

    amatrix_load("test-graph-16.txt", 16, a);

    g.s = 1.;
    for (j = 0; j < N; j++) h[j] = 1.;

    gq_init(&g, a, h);

    nlcg_minimize(gradient, line_min, x);

    return EXIT_SUCCESS;
}
