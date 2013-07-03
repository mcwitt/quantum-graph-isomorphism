/*
 * N        number of vertices (spins)
 * D        dimension of Hilbert space
 * a        adjacency matrix describing graph G
 * psi      state vector in the \sigma^z basis
 * d        diagonal elements of H_P
 */

#include "qgi.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/*
 * NEIGHBOR(i, j)   returns the jth neighbor state of |i>
 * SPIN(i, j)       returns the eigenvalue of sigma^z_j for state |i>
 */

#define PSI2_MAX    10e6

#define NEIGHBOR(i, j)  ((1UL << (j)) ^ (i))
#define SPIN(i, j)      (((int) (((1UL << (j)) & (i)) >> (j)))*2 - 1)

typedef uint64_t index_t;

void qgi_compute_problem_hamiltonian(int a[N][N], double h[N], double d[D])
{
    index_t i;
    int m, n;

    for (i = 0; i < D; i++)
    {
        d[i] = 0.;

        for (n = 0; n < N; n++)
        {
            int s_n = SPIN(i, n);

            d[i] += h[n] * s_n;

            for (m = 0; m < n; m++)
                if (a[n][m] == 1)
                    d[i] += s_n * SPIN(i, m);
        }
    }
}

double qgi_driver_matrix_element(double u[D], double v[D])
{
    double result = 0.;
    index_t i;
    int j;

    for (i = 0; i < D; i++)
        for (j = 0; j < N; j++)
            result += u[i] * v[NEIGHBOR(i, j)];

    return result;
}

double qgi_problem_matrix_element(double d[D], double u[D], double v[D], double *udotv)
{
    double udotv_i, result = 0.;
    index_t i;

    *udotv = 0.;

    for (i = 0; i < D; i++)
    {
        udotv_i = u[i] * v[i];
        *udotv += udotv_i;
        result += d[i] * udotv_i;
    }

    return result;
}

double qgi_matrix_element(double s, double d[D], double u[D], double v[D], double *udotv)
{
    return (1. - s) * qgi_driver_matrix_element(u, v)
               + s  * qgi_problem_matrix_element(d, u, v, udotv);
}

double qgi_energy_grad(double s, double d[D], double psi[D],
        double grad[D], double *psi2, double *edrvr)
{
    double energy;
    int j, k;

    /* compute the energy */
    *edrvr = qgi_driver_matrix_element(psi, psi);
    energy = qgi_problem_matrix_element(d, psi, psi, psi2);
    energy = (1. - s) * (*edrvr) + s * energy;
    energy /= *psi2;

    /* compute the gradient */
    for (k = 0; k < D; k++)
    {
        double sum = 0.;

        for (j = 0; j < N; j++) sum += psi[NEIGHBOR(k, j)];
        grad[k] = 2. * (psi[k] * (s * d[k] - energy) + (1. - s) * sum) / *psi2;
    }

    return energy;
}

double qgi_line_min(double s, double d[D], double psi[D], double delta[D])
{
    double psi2, psi_dot_delta, delta2,
           psi_H_psi, psi_H_delta, delta_H_delta,
           a, b, c, coef, sqr, alpha;

    psi_H_psi     = qgi_matrix_element(s, d, psi,   psi,   &psi2);
    psi_H_delta   = qgi_matrix_element(s, d, psi,   delta, &psi_dot_delta);
    delta_H_delta = qgi_matrix_element(s, d, delta, delta, &delta2);

    a = psi_dot_delta * delta_H_delta - delta2 * psi_H_delta;
    b = psi2 * delta_H_delta - delta2 * psi_H_psi;
    c = psi2 * psi_H_delta - psi_dot_delta * psi_H_psi;

    coef = -0.5 / a;
    sqr = sqrt(b*b - 4.*a*c);
    alpha = coef * (b - sqr);

    /* if critical point is a maximum, use other solution */
    if (2*a*alpha + b < 0.) alpha = coef * (b + sqr);

    return alpha;
}

int qgi_minimize_energy(double s, double d[D], int max_iter, double eps,
        double *energy, double psi[D])
{
    double delta[D], r[D], rprev[D];
    double a, b, edrvr, psi2, r2, r2prev, r2stop;
    index_t i;
    int iter;

    *energy = qgi_energy_grad(s, d, psi, rprev, &psi2, &edrvr);
    r2prev = 0.; for (i = 0; i < D; i++) r2prev += rprev[i] * rprev[i];
    r2stop = eps * r2prev;
    for (i = 0; i < D; i++) delta[i] = rprev[i];

    for (iter = 0; iter < max_iter; iter++)
    {
        a = qgi_line_min(s, d, psi, delta);
        for (i = 0; i < D; i++) psi[i] += a * delta[i];
        *energy = qgi_energy_grad(s, d, psi, r, &psi2, &edrvr);
        r2 = 0.; for (i = 0; i < D; i++) r2 += r[i] * r[i];
        if (r2 < r2stop) break;   /* are we done? */
        /* else update search direction and continue... */
#ifdef USE_FLETCHER_REEVES
        /* compute b (\beta) using the simpler Fletcher-Reeves method */
        b = r2 / r2prev;
#else
        /* use Polak-Ribiere (generally converges faster) */
        b = 0.;
        for (i = 0; i < D; i++) b += r[i] * (r[i] - rprev[i]);
        b /= r2prev;
        if (b < 0.) b = 0.;
#endif
        for (i = 0; i < D; i++) delta[i] = r[i] + b * delta[i];
        for (i = 0; i < D; i++) rprev[i] = r[i];
        r2prev = r2;

        /* renormalize wavefunction if necessary */
        if (psi2 > PSI2_MAX)
        {
            psi2 = sqrt(psi2);
            for (i = 0; i < D; i++) psi[i] /= psi2;
        }
    }

    return iter;
}
