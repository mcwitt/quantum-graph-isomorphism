/**
 * @file   nlcg.h
 * @author Matt Wittmann
 * @brief  Simple implementations of the nonlinear conjugate gradient algorithm.
 */

#include "global.h"

/**
 * Minimize a function using the nonlinear conjugate gradient method with the
 * Fletcher-Reeves update. Returns the minimum value.
 * @param[in]   obj_grad    Function which returns the objective and gradient.
 *                          First argument is specified as *arg*; second and
 *                          third are the state vector and gradient.
 * @param[in]   line_min    Line minimization routine. Return \f$\alpha\f$ to
 *                          minimize \f$f(x+\alpha d)\f$. Arguments are \f$x\f$
 *                          and \f$d\f$.
 * @param[in]   arg         Passed unchanged as first argument to the functions
 *                          *obj_grad* and *line_min*; eliminates the need for
 *                          global variables.
 * @param[in]   eps         Stop when residual is *eps* times initial residual.
 * @param[in]   max_iter    Maximum number of iterations.
 * @param[out]  num_iter    Number of iterations used to reach specified accuracy.
 * @param[out]  x           Solution vector.
 * @param[out]  d           CG search direction.
 * @param[out]  r           Residual.
 * @param[out]  r2          Squared norm of the residual.
 */
double nlcg_minimize(
        double  (*obj_grad)(void*, double*, double*),
        double  (*line_min)(void*, double*, double*),
        void    *arg,
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D],
        double  *r2
        );

/**
 * Same as nlcg_minimize but specialized for minimizing objective functions
 * that are independent of the norm of the vector. Uses a norm-independent
 * stopping condition.
 * @param[in]   obj_x2_grad Function which returns the objective, squared norm,
 *                          and gradient. First argument is specified as *arg*;
 *                          the remaining arguments are the state vector,
 *                          squared norm, and gradient.
 * @param[in]   line_min    Line minimization routine. Return \f$\alpha\f$ to
 *                          minimize \f$f(x+\alpha d)\f$. Arguments are \f$x\f$
 *                          and \f$d\f$.
 * @param[in]   arg         Passed unchanged as first argument to the functions
 *                          *obj_grad* and *line_min*; eliminates the need for
 *                          global variables.
 * @param[in]   eps         Stop when residual is *eps* times initial residual.
 * @param[in]   max_iter    Maximum number of iterations.
 * @param[out]  num_iter    Number of iterations used to reach specified accuracy.
 * @param[out]  x           Solution vector.
 * @param[out]  d           CG search direction.
 * @param[out]  r           Residual.
 * @param[out]  r2          Squared norm of the residual.
 */
double nlcg_minimize_norm_ind(
        double  (*obj_x2_grad)(void*, double*, double*, double*),
        double  (*line_min)(void*, double*, double*),
        void    *arg,
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  *x2,
        double  d[D],
        double  r[D],
        double  *r2
        );

/**
 * Same as nlcg_minimize but uses the Polak-Ribiere update, which is more
 * complex and requires more memory but generally converges faster.
 * @param[in]   obj_grad    Function which returns the objective and gradient.
 *                          First argument is specified as *arg*; second and
 *                          third are the state vector and gradient.
 * @param[in]   line_min    Line minimization routine. Return \f$\alpha\f$ to
 *                          minimize \f$f(x+\alpha d)\f$. Arguments are \f$x\f$
 *                          and \f$d\f$.
 * @param[in]   arg         Passed unchanged as first argument to the functions
 *                          *obj_grad* and *line_min*; eliminates the need for
 *                          global variables.
 * @param[in]   eps         Stop when residual is *eps* times initial residual.
 * @param[in]   max_iter    Maximum number of iterations.
 * @param[out]  num_iter    Number of iterations used to reach specified accuracy.
 * @param[out]  x           Solution vector.
 * @param[out]  d           CG search direction.
 * @param[out]  r           Residual.
 * @param[out]  r2          Squared norm of the residual.
 * @param[out]  rprev       Residual from previous iteration.
 */
double nlcg_minimize_pr(
        double  (*obj_grad)(void*, double*, double*),
        double  (*line_min)(void*, double*, double*),
        void    *arg,
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D],
        double  *r2,
        double  rprev[D]
        );
