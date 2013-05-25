/*
 * N        number of vertices (spins)
 * D        dimension of Hilbert space
 * a        adjacency matrix describing graph G
 * psi      state vector in the \sigma^z basis
 * d        diagonal elements of H_P
 */

#include "gimpq.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/* DEBUG */
#include <assert.h>
#include <stdio.h>

/*
 * NEIGHBOR(i, j)   returns the jth neighbor state of |i>
 * SPIN(i, j)       returns the eigenvalue of sigma^z_j for state |i>
 */

#define NEIGHBOR(i, j)  ((1UL << (j)) ^ (i))
#define SPIN(i, j)      (((int) (((1UL << (j)) & (i)) >> (j)))*2 - 1)

typedef uint64_t index_t;

void gq_compute_problem_hamiltonian(int a[N][N], double h[N], double d[D])
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

double gq_driver_matrix_element(double u[D], double v[D])
{
    double result = 0.;
    index_t i;
    int j;

    for (i = 0; i < D; i++)
        for (j = 0; j < N; j++)
            result += u[i] * v[NEIGHBOR(i, j)];

    return result;
}

double gq_problem_matrix_element(double d[D], double u[D], double v[D], double *udotv)
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

double gq_matrix_element(double s, double d[D], double u[D], double v[D], double *udotv)
{
    return (1. - s) * gq_driver_matrix_element(u, v)
               + s  * gq_problem_matrix_element(d, u, v, udotv);
}

double gq_energy_grad(double s, double d[D], double psi[D], double grad[D])
{
    double psi2, energy;
    int j, k;

    /* compute the energy */
    energy = gq_matrix_element(s, d, psi, psi, &psi2);
    energy /= psi2;

    /* compute the gradient */
    for (k = 0; k < D; k++)
    {
        double sum = 0.;

        for (j = 0; j < N; j++) sum += psi[NEIGHBOR(k, j)];
        grad[k] = 2. * (psi[k] * (s * d[k] - energy) + (1. - s) * sum) / psi2;
    }

    return energy;
}

double gq_line_min(double s, double d[D], double psi[D], double delta[D])
{
    double psi2, psi_dot_delta, delta2,
           psi_H_psi, psi_H_delta, delta_H_delta,
           a, b, c, coef, sqr, alpha;

    psi_H_psi     = gq_matrix_element(s, d, psi,   psi,   &psi2);
    psi_H_delta   = gq_matrix_element(s, d, psi,   delta, &psi_dot_delta);
    delta_H_delta = gq_matrix_element(s, d, delta, delta, &delta2);

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
