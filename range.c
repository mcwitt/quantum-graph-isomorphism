#include "range.h"

double *linspace(double left, double right, int n)
{
    double *a, delta;
    int i;

    delta = (n > 1) ? (right - left) / (n - 1) : 0.;
    a = (double*) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) a[i] = left + i*delta;
    return a;
}

double *logspace(double left, double right, int n)
{
    double *a;
    int i;

    a = linspace(left, right, n);
    for (i = 0; i < n; i++) a[i] = pow(10., a[i]);
    return a;
}

