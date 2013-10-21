#include "nlcg.h"
#include "global.h"

double nlcg_init(
        nlcg_t          *p,
        double          (*obj_x2_grad)(void*, const double*, double*, double*),
        double          (*line_min)(void*, const double*, double*),
        void            *arg,
        const double    *x
        )
{
    p->obj_x2_grad = obj_x2_grad;
    p->line_min = line_min;
    p->arg = arg;
    return nlcg_reset(p, x);
}

double nlcg_reset(nlcg_t *p, const double *x)
{
    double obj;
    UINT i;

    obj = p->obj_x2_grad(p->arg, x, &p->x2, p->d);
    p->r2 = 0.; for (i = 0; i < D; i++) p->r2 += p->d[i] * p->d[i];
    p->r2p = p->r2;
    return obj;
}

double nlcg_iterate(nlcg_t *p, double *x)
{
    double a, b, obj;
    UINT i;

    a = p->line_min(p->arg, x, p->d);
    for (i = 0; i < D; i++) x[i] += a * p->d[i];
    obj = p->obj_x2_grad(p->arg, x, &p->x2, p->r);
    p->r2 = 0.; for (i = 0; i < D; i++) p->r2 += p->r[i] * p->r[i];
    /* update search direction */
    b = p->r2 / p->r2p; /* Fletcher-Reeves */
    for (i = 0; i < D; i++) p->d[i] = p->r[i] + b * p->d[i];
    p->r2p = p->r2;
    return obj;
}

double nlcg_minimize(nlcg_t *p, double *x, double tol, int max_iter, int *num_iter)
{
    double obj = 0., tol2 = tol * tol;
    
    for (*num_iter = 0; *num_iter < max_iter; *num_iter += 1)
    {
        if ((p->r2)/(p->x2) < tol2) break;
        obj = nlcg_iterate(p, x);
    }

    return obj;
}

