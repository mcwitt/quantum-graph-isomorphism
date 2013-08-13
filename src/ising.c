/**
 * @file    ising.c
 * @author  Matt Wittmann
 * @brief  Algorithms related to the quantum Ising model.
 */

#include "ising.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* SPIN(i, j) returns the eigenvalue of sigma^z_j for state |i> */
#define SPIN(i, j)      ((int) ((((i) >> (j)) & 1) << 1) - 1)

void ising_diagonals(int a[N][N], double h[N], double d[D]);
{
    int j, k;
    UINT i;

    for (i = 0; i < D; i++)
    {
        d[i] = 0.;

        for (j = 0; j < N; j++)
        {
            int s_j = SPIN(i, j);

            d[i] -= h[j] * s_j;

            for (k = 0; k < j; k++)
                if (a[j][k] == 1)
                    d[i] += s_j * SPIN(i, k);
        }
    }
}

double ising_matrix_element(double d[D], double u[D], double v[D],
        double *udotv)
{
    double prod, result = 0.;
    UINT i;

    *udotv = 0.;

    for (i = 0; i < D; i++)
    {
        prod = u[i] * v[i];
        *udotv += prod;
        result += prod * d[i];
    }

    return result;
}

double ising_energy_grad(double d[D], double psi[D],
        double grad[D], double *psi2)
{
    double energy;
    UINT i, m;

    /* compute the energy */
    *edrvr = qgi_driver_matrix_element(psi, psi);
    energy = qgi_problem_matrix_element(d, psi, psi, psi2);
    energy = (1. - s) * (*edrvr) + s * energy;
    energy /= *psi2;

    /* compute the gradient */
    for (i = 0; i < D; i++)
    {
        double sum = 0.;

        for (m = 1; m < D; m <<= 1) sum += psi[i^m];
        grad[i] = 2. * (psi[i] * (s * d[i] - energy) + (1. - s) * sum) / *psi2;
    }

    return energy;
}

double ising_sigma_z(double psi[D], int j)
{
    double result = 0.;
    UINT i;

    for (i = 0; i < D; i++)
        result += psi[i] * psi[i] * SPIN(i, j);

    return result;
}

double ising_sigma2_z(double psi[D], int j, int k)
{
    double result = 0.;
    UINT i;

    for (i = 0; i < D; i++)
        result += psi[i] * psi[i] * SPIN(i, j) * SPIN(i, k);

    return result;
}

double ising_sigma_x(double psi[D], int j)
{
    double result = 0.;
    UINT i, m = 1UL << j;

    for (i = 0; i < D; i++) result += psi[i] * psi[i^m];

    return result;
}

double ising_mag_z(double psi[D])
{
    double result = 0.;
    int j;

    for (j = 0; j < N; j++) result += ising_sigma_z(psi, j);

    return result / N;
}

double ising_mag_x(double psi[D])
{
    double result = 0.;
    int j;

    for (j = 0; j < N; j++) result += ising_sigma_x(psi, j);

    return result / N;
}

double ising_overlap(double psi[D])
{
    double m, result = 0.;
    int j, k;

    for (j = 1; j < N; j++)
    {
        for (k = 0; k < j; k++)
        {
            m = ising_sigma2_z(psi, j, k);
            result += m*m;
        }
    }

    return sqrt(2. * result / N / (N-1));
}
