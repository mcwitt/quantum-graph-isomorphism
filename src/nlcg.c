#include "nlcg.h"
#include <stdio.h>

int nlcg_minimize(
        void   (*gradient)(double*, double*),
        double (*line_min)(double*, double*),
        int max_iter,
        double eps,
        double x[D]
        )
{
    index_t i;
    double d[D], r0[D], r1[D];
    double a, b, r0norm2, r1norm2, rtol;
    int iter;

    gradient(x, r0);
    r0norm2 = 0.;

    for (i = 0; i < D; i++)
    {
        d[i] = r0[i];   /* set initial search direction */
        r0norm2 += r0[i] * r0[i];
    }

    rtol = eps * r0norm2;   /* stop when we've reached some small fraction of
                               the initial residual */

    for (iter = 0; iter < max_iter; iter++)
    {
        a = line_min(x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        gradient(x, r1);
        r1norm2 = 0.; for (i = 0; i < D; i++) r1norm2 += r1[i] * r1[i];
        if (r1norm2 < rtol) break;   /* are we done? */
        b = r1norm2 / r0norm2;  /* Fletcher-Reeves method */
        for (i = 0; i < D; i++) d[i] = r1[i] + b*d[i];  /* update search direction */
        for (i = 0; i < D; i++) r0[i] = r1[i];  /* save residual for next iteration */
        r0norm2 = r1norm2;
    }

    return iter;
}
