/**
 * @file    qaa.c
 * @author  Matt Wittmann
 * @brief   Subroutines related to the QAA hamiltonian for the GIMP.
 */

#include "qaa.h"
#include "global.h"
#include <math.h>
#include <stdlib.h>

/* SPIN(i, j) returns the eigenvalue of sigma^z_j for state |i> */
#define SPIN(i, j)  ((int) ((((i) >> (j)) & 1) << 1) - 1)

static double obj_x2_grad(void *arg, const double *psi, double *psi2, double *grad);
static double line_min(void *arg, const double *psi, double *delta);

double qaa_init(qaa_t *p, const int *b, const double *psi)
{
    int ib, j, k, sj;
    UINT i;

    for (i = 0; i < D; i++)
    {
        p->d[i] = 0.;
        ib = 0;

        for (j = 1; j < N; j++)
        {
            sj = SPIN(i, j);

            for (k = 0; k < j; k++)
                if (b[ib++] == 1)
                    p->d[i] += sj * SPIN(i, k);
        }
    }

    return nlcg_init(&p->cg, obj_x2_grad, line_min, p, psi);
}

void qaa_shift_field(qaa_t *p, double dh)
{
    int j;
    UINT i;

    for (i = 0; i < D; i++)
        for (j = 0; j < N; j++)
            p->d[i] -= dh * SPIN(i, j);
}

void qaa_shift_field_1(qaa_t *p, int j, double dh)
{
    UINT i;
    for (i = 0; i < D; i++) p->d[i] -= dh * SPIN(i, j);
}

double qaa_reset(qaa_t *p, const double *psi)
{
    return nlcg_reset(&p->cg, psi);
}

double qaa_iterate(qaa_t *p, double *psi)
{
    return nlcg_iterate(&p->cg, psi);
}

double qaa_minimize(qaa_t *p, double *psi, double tol, int max_iter, int *num_iter)
{
    return nlcg_minimize(&p->cg, psi, tol, max_iter, num_iter);
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

    for (j = 0; j < N; j++)
    {
        for (k = 0; k < N; k++)
        {
            m = qaa_sigma2_z(psi, j, k);
            m -= qaa_sigma_z(psi, j) * qaa_sigma_z(psi, k);
            result += m * m;
        }
    }

    return sqrt(result) / N;
}

static double obj_x2_grad(void *arg, const double *psi, double *psi2, double *grad)
{
    qaa_t *p = arg;
    return qaa_energy_grad(p->s, p->d, psi, grad, psi2, &p->edrvr);
}

static double line_min(void *arg, const double *psi, double *delta)
{
    qaa_t *p = arg;
    return qaa_line_min(p->s, p->d, psi, delta);
}

