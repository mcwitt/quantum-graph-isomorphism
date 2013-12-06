/**
 * @file   nlcg.h
 * @author Matt Wittmann
 * @brief  Simple implementations of the nonlinear conjugate gradient algorithm.
 */

#ifndef NLCG_H
#define NLCG_H

#include "global.h"

typedef struct
{
    /** Pointer to objective function */
    double  (*obj_x2_grad)(void*, const double*, double*, double*);

    /** Pointer to line minimization function */
    double  (*line_min)(void*, const double*, double*);

    void        *arg;   /**< Optional first argument to *obj_grad* and *line_min*. */
    double      r[D];   /**< Residual. */
    double      d[D];   /**< Search direction */
    double      x2;     /**< Squared norm of solution vector. */
    double      r2;     /**< Squared norm of the residual. */
    double      r2p;    /**< Squared norm of previous residual */
} nlcg_t;

/**
 * Initialize the algorithm.
 *
 * @param[out]  p
 * Pointer to nlcg_t instance.
 *
 * @param[in]   obj_x2_grad
 * Pointer to the objective function. The four arguments are (1) an optional
 * argument specified as *arg*, (2) the current solution vector, and the
 * memory locations to store the computed (3) norm and the (4) gradient.
 *
 * @param[in]   line_min
 * Pointer to the line minimization function. This function should return
 * \f$\alpha\f$ to minimize \f$f(x + \alpha d)\f$. The three arguments are an
 * optional argument specified as *arg*, \f$x\f$, and \f$d\f$.
 *
 * @param[in]   arg
 * Optional first argument to *obj_grad* and *line_min*.
 *
 * @param[in]   x
 * Initial guess for solution vector.
 */
double nlcg_init(
        nlcg_t *p,
        double (*obj_x2_grad)(void*, const double*, double*, double*),
        double (*line_min)(void*, const double*, double*),
        void *arg,
        const double *x
        );
{
    p->obj_x2_grad = obj_x2_grad;
    p->line_min = line_min;
    p->arg = arg;
    return nlcg_reset(p, x);
}

/**
 * Reset algorithm following a change in the objective function.
 * @param
 */
double nlcg_reset(nlcg_t *p, const double *x);

/**
 * Do one iteration of the algorithm.
 * @param[in,out]   p   Pointer to nlcg_t instance.
 * @param[in,out]   x   Approximate solution vector.
 */
double nlcg_iterate(nlcg_t *p, double *x);

/**
 * Minimize a function using the nonlinear conjugate gradient method with the
 * Fletcher-Reeves update. Returns the minimum value.
 * @param[in,out]   p           Pointer to nlcg_t instance.
 * @param[in,out]   x           Solution vector.
 * @param[in]       tol         Stop when \f$|r|/|x|\f$ is less than *tol*.
 * @param[in]       max_iter    Maximum number of iterations.
 * @param[out]      num_iter    Number of iterations used.
 */
double nlcg_minimize(nlcg_t *p, double *x, double tol, int max_iter, int *num_iter);

#endif
