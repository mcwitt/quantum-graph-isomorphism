#include "params.h"
#include <stdint.h>

typedef uint64_t index_t;

void nlcg_init(
        void    (*gradient)(double*, double*),
        double  x[D],
        double  d[D],
        double  r[D],
        double  *r2);

void nlcg_iterate(
        void    (*gradient)(double*, double*),
        double  (*line_min)(double*, double*),
        double  x[D],
        double  d[D],
        double  rprev[D],
        double  *r2prev
        );

/**
 * Minimize a function using the nonlinear conjugate gradient method. Returns
 * the number of iterations required to reach specified accuracy.
 * @param[in] gradient Function to compute the gradient.
 * First argument is position vector; result is stored in second argument
 * @param[in] line_min Line minimization routine.
 * Return $\alpha$ to minimize $f(x+\alpha d)$. Arguments are $x$ and $d$.
 * @param[in] max_iter Maximum number of iterations.
 * @param[in] eps Stop when residual is eps times initial residual.
 * @param[out] x Solution vector.
 */
int nlcg_minimize(
        void    (*gradient)(double*, double*),
        double  (*line_min)(double*, double*),
        int     max_iter,
        double  eps,
        double  x[D]
        );
