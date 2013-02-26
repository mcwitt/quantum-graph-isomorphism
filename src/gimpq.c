/*
 * N        number of vertices (spins)
 * D        dimension of Hilbert space
 * a        off-diagonal part of Hamiltonian
 * g        adjacency matrix describing graph G
 * c        state vector in the \sigma^z basis
 * d        diagonal elements of H_P
 */

#include <stdint.h>
#include <stdlib.h>

#include "amatrix.h"

/*
 * NEIGHBOR(i, j)   returns the jth neighbor state of |i>
 * SPIN(i, j)       returns the eigenvalue of sigma^z_i for state |i>
 */

#define NEIGHBOR(i, j)  ((i) ^ (1UL << (j)))
#define SPIN(i, j)      (2 * (((1UL << (j)) & (i)) >> (j)) - 1)

typedef uint64_t ULONG;

void compute_gradient(double c[D], double d[D], double s, double grad[D])
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
        dsum += ci2 * d[i];

        for (j = 0; j < N; j++)
            asum += c[i] * c[NEIGHBOR(i, j)];
    }

    energy = ((1. - s) * asum + s * dsum) / c2;

    /* compute the gradient */
    for (i = 0; i < D; i++)
    {
        double ac_i = 0.;   /* component of A acted on c */

        for (j = 0; j < N; j++) ac_i += c[NEIGHBOR(i, j)];
        grad[i] = 2. * (c[i] * (s * d[i] - energy) + (1. - s) * ac_i) / c2;
    }
}

static void compute_diagonal_elements(int g[N][N], double h[N], double d[D])
{
    ULONG i;
    int n, m;

    for (i = 0; i < D; i++)
    {
        d[i] = 0.;

        for (n = 0; n < N; n++)
        {
            int s_n = SPIN(i, n);

            d[i] += h[n] * s_n;

            for (m = 0; m < N; m++)
            {
                if (g[n][m] == 0) continue;
                d[i] += s_n * SPIN(i, m);
            }
        }
    }
}
