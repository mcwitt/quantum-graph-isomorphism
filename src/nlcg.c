#include "nlcg.h"

#define MAX_ITER 100
#define EPS 1e-6

void nlcg_minimize(
        void   (*gradient)(double*, double*),
        double (*line_min)(double*, double*),
        double x[D])
{
    ULONG i;
    double d[D], r0[D], r1[D];
    double a, b, r0norm2, r1norm2, tol;
    int iter;

    gradient(x, r0);
    r0norm2 = 0.;
    for (i = 0; i < D; i++) { d[i] = r0[i]; r0norm2 += r0[i] * r0[i]; }
    tol = EPS * r0norm2;

    for (iter = 0; iter < MAX_ITER; iter++)
    {
        a = line_min(x, d);
        for (i = 0; i < D; i++) x[i] += a * d[i];
        gradient(x, r1);
        r1norm2 = 0.; for (i = 0; i < D; i++) r1norm2 += r1[i] * r1[i];
        if (r1norm2 < tol) break;
        b = r1norm2 / r0norm2;  /* Fletcher-Reeves method */
        for (i = 0; i < D; i++) d[i] = r1[i] + b*d[i];
        for (i = 0; i < D; i++) r0[i] = r1[i];
        r0norm2 = r1norm2;
    }
}
