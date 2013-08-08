#include "params.h"
#include <stdint.h>

/**
 * Minimize a function using the nonlinear conjugate gradient method with the
 * Fletcher-Reeves update. Returns the minimum value.
 * @param[in]   obj_grad    function which returns the objective and gradient
 * First argument is state vector; gradient is stored in second argument
 * @param[in]   line_min    line minimization routine
 * Return $\alpha$ to minimize $f(x+\alpha d)$. Arguments are $x$ and $d$.
 * @param[in]   eps         stop when residual is eps times initial residual
 * @param[in]   max_iter    maximum number of iterations
 * @param[out]  num_iter    number of iterations used to reach specified accuracy
 * @param[out]  x           solution vector
 * @param[out]  d           CG search direction
 * @param[out]  r           residual
 * @param[out]  r2          squared norm of the residual
 */
double nlcg_minimize(
        double  (*obj_grad)(double*, double*),
        double  (*line_min)(double*, double*),
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D],
        double  *r2
        );

/**
 * Same as nlcg_minimize but uses the Polak-Ribiere update, which is more
 * complex and requires more memory but generally converges faster.
 * @param[in]   obj_grad    function which returns the objective and gradient
 * First argument is state vector; gradient is stored in second argument
 * @param[in]   line_min    line minimization routine
 * Return $\alpha$ to minimize $f(x+\alpha d)$. Arguments are $x$ and $d$.
 * @param[in]   eps         stop when residual is eps times initial residual
 * @param[in]   max_iter    maximum number of iterations
 * @param[out]  num_iter    number of iterations used to reach specified accuracy
 * @param[out]  x           solution vector
 * @param[out]  d           CG search direction
 * @param[out]  r           residual
 * @param[out]  rprev       residual from previous iteration
 */
double nlcg_minimize_pr(
        double  (*obj_grad)(double*, double*),
        double  (*line_min)(double*, double*),
        double  eps,
        int     max_iter,
        int     *num_iter,
        double  x[D],
        double  d[D],
        double  r[D],
        double  rprev[D]
        );
