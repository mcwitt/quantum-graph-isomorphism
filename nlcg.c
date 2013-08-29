#include "nlcg.h"

double nlcg_minimize(
        double  (*obj_grad)(void*, const double*, double*),
        double  (*line_min)(void*, const double*, double*),
        void    *arg,
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D],
        double  *r2
        )
{
    double a, b, obj, r2prev, r2stop;
    UINT i;

    obj = obj_grad(arg, x, d);
    r2prev = 0.; for (i = 0; i < D; i++) r2prev += d[i] * d[i];
    r2stop = eps * eps * r2prev;

    for (*num_iter = 0; *num_iter < max_iter; *num_iter += 1)
    {
        a = line_min(arg, x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        obj = obj_grad(arg, x, r);
        *r2 = 0.; for (i = 0; i < D; i++) *r2 += r[i] * r[i];
        if (*r2 < r2stop) break;    /* are we done? */
        /* else update search direction and continue... */
        b = *r2 / r2prev;   /* Fletcher-Reeves */
        for (i = 0; i < D; i++) d[i] = r[i] + b*d[i];
        r2prev = *r2;
    }

    return obj;
}

double nlcg_minimize_norm_ind(
        double  (*obj_x2_grad)(void*, const double*, double*, double*),
        double  (*line_min)(void*, const double*, double*),
        void    *arg,
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  *x2,
        double  d[D],
        double  r[D],
        double  *r2
        )
{
    double a, b, eps2, obj, r2prev;
    UINT i;

    eps2 = eps * eps;
    obj = obj_x2_grad(arg, x, x2, d);
    r2prev = 0.; for (i = 0; i < D; i++) r2prev += d[i] * d[i];

    for (*num_iter = 0; *num_iter < max_iter; *num_iter += 1)
    {
        a = line_min(arg, x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        obj = obj_x2_grad(arg, x, x2, r);
        *r2 = 0.; for (i = 0; i < D; i++) *r2 += r[i] * r[i];
        if (*r2 < eps2 * (*x2)) break;  /* are we done? */
        /* else update search direction and continue... */
        b = *r2 / r2prev;   /* Fletcher-Reeves */
        for (i = 0; i < D; i++) d[i] = r[i] + b*d[i];
        r2prev = *r2;
    }

    return obj;
}

double nlcg_minimize_pr(
        double  (*obj_grad)(void*, const double*, double*),
        double  (*line_min)(void*, const double*, double*),
        void    *arg,
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D],
        double  *r2,
        double  rprev[D]
        )
{
    double a, b, obj, r2prev, r2stop;
    UINT i;

    obj = obj_grad(arg, x, rprev);
    r2prev = 0.; for (i = 0; i < D; i++) r2prev += rprev[i] * rprev[i];
    r2stop = eps * eps * r2prev;
    for (i = 0; i < D; i++) d[i] = rprev[i];

    for (*num_iter = 0; *num_iter < max_iter; *num_iter += 1)
    {
        a = line_min(arg, x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        obj = obj_grad(arg, x, r);
        *r2 = 0.; for (i = 0; i < D; i++) *r2 += r[i] * r[i];
        if (*r2 < r2stop) break;   /* are we done? */
        /* else update search direction and continue... */

        /* Polak-Ribiere update (generally converges faster) */
        b = 0.; for (i = 0; i < D; i++) b += r[i] * rprev[i];
        b = (*r2 - b) / r2prev;
        if (b < 0.) b = 0.;

        for (i = 0; i < D; i++) d[i] = r[i] + b*d[i];
        for (i = 0; i < D; i++) rprev[i] = r[i];
        r2prev = *r2;
    }

    return obj;
}
