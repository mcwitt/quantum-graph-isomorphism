#include "nlcg.h"

typedef uint64_t index_t;

double nlcg_minimize(
        double  (*obj_grad)(double*, double*),
        double  (*line_min)(double*, double*),
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D]
        )
{
    double a, b, obj, r2, r2prev, r2stop;
    index_t i;

    obj = obj_grad(x, r);
    r2prev = 0.; for (i = 0; i < D; i++) r2prev += r[i] * r[i];
    /* r2stop = eps * r2prev; */
    r2stop = eps;
    for (i = 0; i < D; i++) d[i] = r[i];

    for (*num_iter = 0; *num_iter < max_iter; *num_iter += 1)
    {
        a = line_min(x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        obj = obj_grad(x, r);
        r2 = 0.; for (i = 0; i < D; i++) r2 += r[i] * r[i];
        if (r2 < r2stop) break;   /* are we done? */
        /* else update search direction and continue... */
        b = r2 / r2prev;    /* Fletcher-Reeves */
        for (i = 0; i < D; i++) d[i] = r[i] + b*d[i];
        r2prev = r2;
    }

    return obj;
}

double nlcg_minimize_pr(
        double  (*obj_grad)(double*, double*),
        double  (*line_min)(double*, double*),
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D],
        double  rprev[D]
        )
{
    double a, b, obj, r2, r2prev, r2stop;
    index_t i;

    obj = obj_grad(x, rprev);
    r2prev = 0.; for (i = 0; i < D; i++) r2prev += rprev[i] * rprev[i];
    r2stop = eps * r2prev;
    for (i = 0; i < D; i++) d[i] = rprev[i];

    for (*num_iter = 0; *num_iter < max_iter; *num_iter += 1)
    {
        a = line_min(x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        obj = obj_grad(x, r);
        r2 = 0.; for (i = 0; i < D; i++) r2 += r[i] * r[i];
        if (r2 < r2stop) break;   /* are we done? */
        /* else update search direction and continue... */

        /* Polak-Ribiere update (generally converges faster) */
        b = 0.; for (i = 0; i < D; i++) b += r[i] * rprev[i];
        b = (r2 - b) / r2prev;
        if (b < 0.) b = 0.;

        for (i = 0; i < D; i++) d[i] = r[i] + b*d[i];
        for (i = 0; i < D; i++) rprev[i] = r[i];
        r2prev = r2;
    }

    return obj;
}
