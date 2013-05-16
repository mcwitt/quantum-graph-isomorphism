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
    double psi_dot_delta = 0.,
           psi2 = 0.,
           delta2 = 0.,
           psi_H_delta = 0.,
           psi_H_psi = 0.,
           delta_H_delta = 0.;
    int i;

    for (i = 0; i < D; i++)
    {
        psi2 += 
        psi_dot_delta += psi[i] * delta[i];
    }




    

}
