#include "nlcg.h"

#define MAX_ITER 100
#define EPS 1e-6

void nlcg_minimize(
        void   (*gradient)(double*, double*),
        double (*line_min)(double*, double*),
        double x[D])
{
    double d[D], r0[D], r1[D];
    double a, b, r0norm2, r1norm2, tol;
    int i, j;

    gradient(x, r0);
    r0norm2 = 0.;
    for (j = 0; j < D; j++) { d[j] = r0[j]; r0norm2 += r0[j] * r0[j]; }
    tol = EPS * r0norm2;

    for (i = 0; i < MAX_ITER; i++)
    {
        a = line_min(x, d);
        for (j = 0; j < D; j++) x[j] += a * d[j];
        gradient(x, r1);
        r1norm2 = 0.; for (j = 0; j < D; j++) r1norm2 += r1[j] * r1[j];
        if (r1norm2 < tol) break;
        b = r1norm2 / r0norm2;  /* Fletcher-Reeves method */
        for (j = 0; j < D; j++) d[j] = r1[j] + b*d[j];
        for (j = 0; j < D; j++) r0[j] = r1[j];
        r0norm2 = r1norm2;
    }
}
