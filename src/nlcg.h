/**
 * Minimize a function using the nonlinear conjugate gradient method.
 * @param[in] gradient Function to compute the gradient.
 * First argument is location to compute gradient, result is stored in second argument
 * @param[in] line_min Line minimization routine.
 * Return $\alpha$ to minimize $f(x+\alpha d)$. Arguments are $x$ and $d$.
 * @param[out] x Solution vector.
 */

#include "params.h"

void nlcg_minimize(
        void   (*gradient)(double*, double*),
        double (*line_min)(double*, double*),
        double x[D]);
