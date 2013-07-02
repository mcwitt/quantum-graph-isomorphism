#include "nlcg.h"
#include <stdio.h>

void nlcg_init(
        void    (*gradient)(double*, double*),
        double  x[D],
        double  d[D],
        double  r[D],
        double  *r2)
{
    index_t i;

    gradient(x, r);
    *r2 = 0.;

    for (i = 0; i < D; i++)
    {
        d[i] = r[i];   /* set initial search direction */
        *r2 += r[i] * r[i];
    }
}

void nlcg_iterate(
        void    (*gradient)(double*, double*),
        double  (*line_min)(double*, double*),
        double  x[D],
        double  d[D],
        double  rprev[D],
        double  *r2prev
        )
{
    double r[D], a, b, r2;
    index_t i;

    a = line_min(x, d);
    for (i = 0; i < D; i++) x[i] += a * d[i];
    gradient(x, r);
    r2 = 0.; for (i = 0; i < D; i++) r2 += r[i] * r[i];
#ifdef USE_FLETCHER_REEVES
    /* compute b (\beta) using the simpler Fletcher-Reeves method */
    b = r2 / *r2prev;
#else
    /* use Polak-Ribiere (generally converges in fewer iterations) */
    b = 0.;
    for (i = 0; i < D; i++) b += r[i] * (r[i] - rprev[i]);
    b /= *r2prev;
    if (b < 0.) b = 0.;
#endif
    for (i = 0; i < D; i++) d[i] = r[i] + b*d[i];  /* update search direction */
    for (i = 0; i < D; i++) rprev[i] = r[i];  /* save residual for next iteration */
    *r2prev = r2;
}

int nlcg_minimize(
        void    (*gradient)(double*, double*),
        double  (*line_min)(double*, double*),
        int     max_iter,
        double  eps,
        double  x[D]
        )
{
    double d[D], r[D], r2, r2max;
    int iter;

    nlcg_init(gradient, x, d, r, &r2);
    r2max = eps * r2;   /* stop when we've reached some small fraction of
                           the initial residual */

    for (iter = 0; iter < max_iter; iter++)
    {
        nlcg_iterate(gradient, line_min, x, d, r, &r2);
        if (r2 < r2max) break;   /* are we done? */
    }

    return iter;
}
