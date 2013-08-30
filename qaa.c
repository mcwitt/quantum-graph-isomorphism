/**
 * @file    qaa.c
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

#include "qaa.h"
#include "global.h"
#include "nlcg.h"
#include <math.h>
#include <stdlib.h>

/* SPIN(i, j) returns the eigenvalue of sigma^z_j for state |i> */
#define SPIN(i, j)  ((int) ((((i) >> (j)) & 1) << 1) - 1)

typedef struct { double s, *edrvr; const double *d; } arg_t;

static double obj_x2_grad(void *arg, const double *psi, double *x2, double *grad);
static double line_min(void *arg, const double *psi, double *delta);

void qaa_compute_diagonals(const int *b, double *d)
{
    int ib, j, k, sj;
    UINT i;

    for (i = 0; i < D; i++)
    {
        d[i] = 0.;
        ib = 0;

        for (j = 1; j < N; j++)
        {
            sj = SPIN(i, j);

            for (k = 0; k < j; k++)
                if (b[ib++] == 1)
                    d[i] += sj * SPIN(i, k);
        }
    }
}

void qaa_update_diagonals(double dh, double *d)
{
    int j;
    UINT i;

    for (i = 0; i < D; i++)
        for (j = 0; j < N; j++)
            d[i] -= dh * SPIN(i, j);
}

void qaa_update_diagonals_1(int j, double dh, double *d)
{
    UINT i;
    for (i = 0; i < D; i++) d[i] -= dh * SPIN(i, j);
}

double qaa_minimize_energy(
        double s,
        const double *d,
        double eps,
        int max_iter,
        int *num_iter,
        double *edrvr,
        double *psi,
        double *psi2,
        double *r,
        double *r2,
        double *delta
        )
{
    arg_t args = {s, edrvr, d};

    return nlcg_minimize_norm_ind(obj_x2_grad, line_min, &args, eps,
            max_iter, num_iter, psi, psi2, r, r2, delta);
}

double qaa_me_driver(const double *u, const double *v)
{
    double result = 0.;
    UINT i, m;

    for (i = 0; i < D; i++)
        for (m = 1; m < D; m <<= 1)
            result += u[i] * v[i^m];

    return 0.5 * result;
}

double qaa_me_problem(
        const double *d,
        const double *u,
        const double *v,
        double *udotv)
{
    double p, result = 0.;
    UINT i;

    *udotv = 0.;

    for (i = 0; i < D; i++)
    {
        p = u[i] * v[i];
        *udotv += p;
        result += p * d[i];
    }

    return result;
}

double qaa_energy_grad(
        double s,
        const double *d,
        const double *psi,
        double *grad,
        double *psi2,
        double *edrvr)
{
    double energy;
    UINT i, m;

    /* compute the energy */
    *edrvr = qaa_me_driver(psi, psi);
    energy = qaa_me_problem(d, psi, psi, psi2);
    energy = (1. - s) * (*edrvr) + s * energy;
    energy /= *psi2;

    /* compute the gradient */
    for (i = 0; i < D; i++)
    {
        double sum = 0.;

        for (m = 1; m < D; m <<= 1) sum += psi[i^m];
        grad[i] = (2. * psi[i] * (s * d[i] - energy) + (1. - s) * sum) / *psi2;
    }

    return energy;
}

double qaa_line_min(
        double s,
        const double *d,
        const double *psi,
        const double *delta)
{
    double psi2, psi_dot_delta, delta2,
           psi_H_psi, psi_H_delta, delta_H_delta,
           a, b, c, coef, sqrd, x, oms = 1. - s;

    psi_H_psi     = s * qaa_me_problem(d, psi,   psi,   &psi2);
    psi_H_delta   = s * qaa_me_problem(d, psi,   delta, &psi_dot_delta);
    delta_H_delta = s * qaa_me_problem(d, delta, delta, &delta2);

    psi_H_psi     += oms * qaa_me_driver(psi, psi);
    psi_H_delta   += oms * qaa_me_driver(psi, delta);
    delta_H_delta += oms * qaa_me_driver(delta, delta);

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

double qaa_sigma_z(const double *psi, int j)
{
    double result = 0.;
    UINT i;

    for (i = 0; i < D; i++)
        result += psi[i] * psi[i] * SPIN(i, j);

    return result;
}

double qaa_sigma2_z(const double *psi, int j, int k)
{
    double result = 0.;
    UINT i;

    for (i = 0; i < D; i++)
        result += psi[i] * psi[i] * SPIN(i, j) * SPIN(i, k);

    return result;
}

double qaa_sigma_x(const double *psi, int j)
{
    double result = 0.;
    UINT i, m = 1UL << j;

    for (i = 0; i < D; i++) result += psi[i] * psi[i^m];

    return result;
}

double qaa_mag_z(const double *psi)
{
    double result = 0.;
    int j;

    for (j = 0; j < N; j++) result += qaa_sigma_z(psi, j);

    return result / N;
}

double qaa_overlap(const double *psi)
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

static double obj_x2_grad(void *arg, const double *psi, double *x2, double *grad)
{
    arg_t args = *((arg_t*) arg);
    return qaa_energy_grad(args.s, args.d, psi, grad, x2, args.edrvr);
}

static double line_min(void *arg, const double *psi, double *delta)
{
    arg_t args = *((arg_t*) arg);
    return qaa_line_min(args.s, args.d, psi, delta);
}

