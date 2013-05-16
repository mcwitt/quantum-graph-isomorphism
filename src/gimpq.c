/*
 * N        number of vertices (spins)
 * D        dimension of Hilbert space
 * a        off-diagonal part of Hamiltonian
 * g        adjacency matrix describing graph G
 * c        state vector in the \sigma^z basis
 * d        diagonal elements of H_P
 */

#include "gimpq.h"

#include <stdint.h>
#include <stdlib.h>

/*
 * NEIGHBOR(i, j)   returns the jth neighbor state of |i>
 * SPIN(i, j)       returns the eigenvalue of sigma^z_j for state |i>
 */

#define NEIGHBOR(i, j)  ((i) ^ (1UL << (j)))
#define SPIN(i, j)      (2 * (((1UL << (j)) & (i)) >> (j)) - 1)

typedef uint64_t ULONG;

void gq_init(GQProblem *p, int g[N][N], double h[N])
{
    ULONG i;
    int n, m;

    s = 0.;

    /* compute diagonal elements of H_P */

    for (i = 0; i < D; i++)
    {
        p->d[i] = 0.;

        for (n = 0; n < N; n++)
        {
            int s_n = SPIN(i, n);

            p->d[i] += h[n] * s_n;

            for (m = 0; m < N; m++)
            {
                if (g[n][m] == 0) continue;
                p->d[i] += s_n * SPIN(i, m);
            }
        }
    }
}


double gq_driver_matrix_element(double u[D], double v[D])
{
    double result = 0.;
    ULONG i;
    int j;

    for (i = 0; i < D; i++)
        for (j = 0; j < N; j++)
            result += u[i] * v[NEIGHBOR(i, j)];
}

void gq_problem_matrix_element(double d[D], double u[D], double v[D], double *udotv)
{
    double udotv_i, result = 0.;
    ULONG i;

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

void gq_grad_energy(double s, double d[D], double psi[D], double grad[D])
{
    int k;
    double psi2, energy;

    /* compute the energy */
    energy = gq_matrix_element(s, d, psi, psi, &psi2);
    energy /= psi2;

    /* compute the gradient */
    for (k = 0; k < D; k++)
    {
        double apsi_i = 0.;   /* component of A psi */

        for (j = 0; j < N; j++) apsi_i += psi[NEIGHBOR(k, j)];
        grad[k] = 2. * (psi[k] * (s * p->d[k] - energy) + (1. - s) * apsi_i) / psi2;
    }
}

void gq_grad_energy(GQProblem *p, double s, double c[D], double grad[D])
{
    ULONG i;
    int j;
    double c2, energy, dsum, asum;

    c2 = 0.;
    dsum = 0.;
    asum = 0.;

    /* compute the energy */
    for (i = 0; i < D; i++)
    {
        double ci2 = c[i] * c[i];

        c2 += ci2;
        dsum += ci2 * p->d[i];

        for (j = 0; j < N; j++)
            asum += c[i] * c[NEIGHBOR(i, j)];
    }

    energy = ((1. - s) * asum + s * dsum) / c2;

    /* compute the gradient */
    for (i = 0; i < D; i++)
    {
        double ac_i = 0.;   /* component of A c */

        for (j = 0; j < N; j++) ac_i += c[NEIGHBOR(i, j)];
        grad[i] = 2. * (c[i] * (s * p->d[i] - energy) + (1. - s) * ac_i) / c2;
    }
}

double gq_line_min(GQProblem *p, double s, double psi[D], double delta[D])
{
    double psi2, psi_dot_delta, delta2, psi_H_psi, psi_H_delta, delta_H_delta,
           alpha, beta, gamma;

    psi_H_psi = gq_matrix_element(s, d, psi, psi, &psi2);
    psi_H_delta = gq_matrix_element(s, d, psi, delta, &psi_dot_delta);
    delta_H_delta = gq_matrix_element(s, d, delta, delta, &delta2);

    alpha = psi_dot_delta * delta_H_delta - delta2 * psi_H_delta;
    beta  = psi2 * delta_H_delta - delta2 * psi_H_psi;
    gamma = psi2 * delta_H_psi - psi_dot_delta * psi_H_psi;

    disc = beta*beta - 4.*alpha*gamma;
    return 0.5*(sqrt(disc) - beta) / alpha;
}
