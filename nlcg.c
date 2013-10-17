#include "nlcg.h"
#include "global.h"

#define DOT(i, u, v, uv) { \
    (uv) = 0.; for (i = 0; i < D; i++) (uv) += (u)[i] * (v)[i]; \
}

double nlcg_init(
        nlcg_t  *p,
        double  (*obj_x2_grad)(void*, const double*, double*, double*),
        double  (*line_min)(void*, const double*, double*),
        void    *arg,
        double  *x
        )
{
    double obj;
    UINT i;

    p->obj_x2_grad = obj_x2_grad;
    p->line_min = line_min;
    p->arg = arg;
    p->x = x;

    /* compute initial search direction (residual) */
    obj = p->obj_x2_grad(p->arg, p->x, &p->x2, p->d);
    DOT(i, p->d, p->d, p->r2);
    p->r2p = p->r2;

    return obj;
}

double nlcg_iterate(nlcg_t *p)
{
    double a, b, obj;
    UINT i;

    a = p->line_min(p->arg, p->x, p->d);
    for (i = 0; i < D; i++) p->x[i] += a * p->d[i];
    obj = p->obj_x2_grad(p->arg, p->x, &p->x2, p->r);
    DOT(i, p->r, p->r, p->r2);
    /* update search direction */
    b = p->r2 / p->r2p; /* Fletcher-Reeves */
    for (i = 0; i < D; i++) p->d[i] = p->r[i] + b * p->d[i];
    p->r2p = p->r2;

    return obj;
}

double nlcg_minimize(nlcg_t *p, double tol, int max_iter, int *num_iter)
{
    double obj = 0., tol2 = tol * tol;
    
    for (*num_iter = 0; *num_iter < max_iter; *num_iter += 1)
    {
        if ((p->r2)/(p->x2) < tol2) break;
        obj = nlcg_iterate(p);
    }

    return obj;
}

