/**
 * @file    qaa.c
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

#include "qaa.h"
#include "nlcg.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* SPIN(i, j) returns the eigenvalue of sigma^z_j for state |i> */
#define SPIN(i, j)      ((int) ((((i) >> (j)) & 1) << 1) - 1)

void qaa_compute_diagonals(int a[N][N], double h[N], double d[D])
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

typedef struct { double s, *edrvr, *d; } arg_t;

static double obj_x2_grad(void *arg, double psi[D], double *x2, double grad[D])
{

    arg_t args = *((arg_t*) arg);
    return qaa_energy_grad(args.s, args.d, psi, grad, x2, args.edrvr);
}

static double line_min(void *arg, double psi[D], double delta[D])
{
    arg_t args = *((arg_t*) arg);
    return qaa_line_min(args.s, args.d, psi, delta);
}

double qaa_minimize_energy(
        double  s,
        double  d[D],
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  *edrvr,
        double  psi[D],
        double  *psi2,
        double  delta[D],
        double  r[D],
        double  *r2
        )
{
    arg_t args = {s, edrvr, d};

    return nlcg_minimize_norm_ind(obj_x2_grad, line_min, &args, eps,
            max_iter, num_iter, psi, psi2, delta, r, r2);
}

double qaa_driver_matrix_element(double u[D], double v[D])
{
    double result = 0.;
    UINT i, m;

    for (i = 0; i < D; i++)
        for (m = 1; m < D; m <<= 1)
            result += u[i] * v[i^m];

    return result;
}

double qaa_problem_matrix_element(double d[D], double u[D], double v[D], double *udotv)
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

double qaa_matrix_element(double s, double d[D], double u[D], double v[D], double *udotv)
{
    return (1. - s) * qaa_driver_matrix_element(u, v)
               + s  * qaa_problem_matrix_element(d, u, v, udotv);
}

double qaa_energy_grad(double s, double d[D], double psi[D],
        double grad[D], double *psi2, double *edrvr)
{
    double energy;
    UINT i, m;

    /* compute the energy */
    *edrvr = qaa_driver_matrix_element(psi, psi);
    energy = qaa_problem_matrix_element(d, psi, psi, psi2);
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

double qaa_line_min(double s, double d[D], double psi[D], double delta[D])
{
    double psi2, psi_dot_delta, delta2,
           psi_H_psi, psi_H_delta, delta_H_delta,
           a, b, c, coef, sqrd, x;

    psi_H_psi     = qaa_matrix_element(s, d, psi,   psi,   &psi2);
    psi_H_delta   = qaa_matrix_element(s, d, psi,   delta, &psi_dot_delta);
    delta_H_delta = qaa_matrix_element(s, d, delta, delta, &delta2);

    a = psi_dot_delta * delta_H_delta - delta2 * psi_H_delta;
    b = psi2 * delta_H_delta - delta2 * psi_H_psi;
    c = psi2 * psi_H_delta - psi_dot_delta * psi_H_psi;

    coef = -0.5 / a;
    sqrd = sqrt(b*b - 4.*a*c);
    x = coef * (b - sqrd);

    /* if critical point is a maximum, use other solution */
    if (2*a*x + b < 0.) x = coef * (b + sqrd);

    return x;
}


double qaa_sigma_z(double psi[D], int j)
{
    double result = 0.;
    UINT i;

    for (i = 0; i < D; i++)
        result += psi[i] * psi[i] * SPIN(i, j);

    return result;
}

double qaa_sigma2_z(double psi[D], int j, int k)
{
    double result = 0.;
    UINT i;

    for (i = 0; i < D; i++)
        result += psi[i] * psi[i] * SPIN(i, j) * SPIN(i, k);

    return result;
}

double qaa_sigma_x(double psi[D], int j)
{
    double result = 0.;
    UINT i, m = 1UL << j;

    for (i = 0; i < D; i++) result += psi[i] * psi[i^m];

    return result;
}

double qaa_mag_z(double psi[D])
{
    double result = 0.;
    int j;

    for (j = 0; j < N; j++) result += qaa_sigma_z(psi, j);

    return result / N;
}

double qaa_mag_x(double psi[D])
{
    double result = 0.;
    int j;

    for (j = 0; j < N; j++) result += qaa_sigma_x(psi, j);

    return result / N;
}

double qaa_overlap(double psi[D])
{
    double m, result = 0.;
    int j, k;

    for (j = 1; j < N; j++)
    {
        for (k = 0; k < j; k++)
        {
            m = qaa_sigma2_z(psi, j, k);
            result += m*m;
        }
    }

    return sqrt(2. * result / N / (N-1));
}
