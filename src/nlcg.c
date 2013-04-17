#define MAX_ITER 100
#define EPS 1e-6

void nlcg_minimize(
        void   (*gradient)(double*, double*),
        double (*line_min)(double*, double*),
        double x[CG_DIM])
{
    double d[CG_DIM], r0[CG_DIM], r1[CG_DIM];
    double a, b, r0norm2, r1norm2, goal;

    gradient(x, r0);
    for (i = 0; i < CG_DIM; i++) d[i] = r0[i];
    r0norm2 = 0.; for (i = 0; i < CG_DIM; i++) r0norm2 += r0[i] * r0[i];
    goal = EPS * r0norm2;

    for (iter = 0; iter < MAX_ITER; iter++)
    {
        a = line_min(x, r0);
        for (i = 0; i < CG_DIM; i++) x[i] += a * r0[i];
        gradient(x, r1);
        r1norm2 = 0.; for (i = 0; i < CG_DIM; i++) r1norm2 += r1[i] * r1[i];
        if (r1norm2 < goal) break;
        b = r1norm2 / r0norm2;  /* Fletcher-Reeves method */
        for (i = 0; i < CG_DIM; i++) d[i] = r1[i] + b*d[i];
    }
}
