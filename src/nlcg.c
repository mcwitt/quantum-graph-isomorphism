#include "nlcg.h"

int nlcg_minimize(
        void   (*gradient)(double*, double*),
        double (*line_min)(double*, double*),
        int    max_iter,
        double eps,
        double x[D]
        )
{

    double d[D], r[D], rprev[D];
    double a, b, r2, r2prev, r2stop;
    index_t i;
    int iter;

    gradient(x, rprev);
    r2prev = 0.; for (i = 0; i < D; i++) r2prev += rprev[i] * rprev[i];
    r2stop = eps * r2prev;
    for (i = 0; i < D; i++) d[i] = rprev[i];

    for (iter = 0; iter < max_iter; iter++)
    {
        a = line_min(x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        gradient(x, r);
        r2 = 0.; for (i = 0; i < D; i++) r2 += r[i] * r[i];
        if (r2 < r2stop) break;   /* are we done? */
        /* else update search direction and continue... */
#ifdef USE_POLAK_RIBIERE
        /* compute b (\beta) using the simpler Fletcher-Reeves method */
        b = r2 / r2prev;
#else
        /* use Polak-Ribiere (generally converges faster) */
        b = 0.; for (i = 0; i < D; i++) b += r[i] * rprev[i];
        b = (r2 - b) / r2prev;
        if (b < 0.) b = 0.;
#endif
        for (i = 0; i < D; i++) d[i] = r[i] + b*d[i];
        for (i = 0; i < D; i++) rprev[i] = r[i];
        r2prev = r2;
    }

    return iter;
}
