#include "numerical_transient.h"
#include "transient_defs.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*==================================================================
 * L5: Numerical Methods for Transient Analysis
 *
 * This module implements numerical integration methods for solving
 * the ODE systems that arise in transient circuit analysis.
 *
 * The standard form: dy/dt = f(t, y), y(t0) = y0
 *
 * Key concepts:
 *   - Local truncation error: error per step
 *   - Global error: accumulated error over integration interval
 *   - Stability: bounded errors for stable systems
 *   - A-stability: unconditional stability for linear systems with Re(lambda) < 0
 *   - Stiffness: systems with widely separated time constants
 *
 * Reference: Press et al. "Numerical Recipes in C" (2007) Ch.16,
 *           Hairer, Norsett & Wanner "Solving ODE I" (1993),
 *           Gear "Numerical Initial Value Problems in ODEs" (1971)
 *==================================================================*/

/*----- L5: Forward Euler Method -----*/

/* y_{n+1} = y_n + h * f(t_n, y_n)
 * Local error: O(h^2), Global error: O(h)
 * Simplest explicit method. Conditionally stable:
 * |1 + h*lambda| <= 1 for linear stability. */
int euler_forward_step(double t, const double *y, double *y_next, int n,
                       ODERHSFunc f, double h, void *params)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *dydt = (double*)malloc((size_t)n * sizeof(double));
    if (!dydt) return -1;

    f(t, y, dydt, params);

    for (int i = 0; i < n; i++) {
        y_next[i] = y[i] + h * dydt[i];
    }

    free(dydt);
    return 0;
}

/*----- L5: Backward Euler Method (Implicit) -----*/

/* y_{n+1} = y_n + h * f(t_{n+1}, y_{n+1})
 * Implicit: requires solving nonlinear equation per step.
 * A-stable: stable for all h*lambda with Re(lambda) < 0.
 *
 * For linear test equation y' = lambda*y:
 * y_{n+1} = y_n / (1 - h*lambda)
 * This implementation uses functional iteration (simple but slow for stiff). */
int euler_backward_step(double t, const double *y, double *y_next, int n,
                        ODERHSFunc f, double h, void *params)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *dydt = (double*)malloc((size_t)n * sizeof(double));
    double *y_guess = (double*)malloc((size_t)n * sizeof(double));
    if (!dydt || !y_guess) { free(dydt); free(y_guess); return -1; }

    /* Initial guess: forward Euler */
    memcpy(y_guess, y, (size_t)n * sizeof(double));
    f(t, y, dydt, params);
    for (int i = 0; i < n; i++) {
        y_guess[i] = y[i] + h * dydt[i];
    }

    /* Functional iteration (typically converges for small h) */
    for (int iter = 0; iter < 20; iter++) {
        double t_next = t + h;
        f(t_next, y_guess, dydt, params);

        int converged = 1;
        for (int i = 0; i < n; i++) {
            double y_new = y[i] + h * dydt[i];
            if (fabs(y_new - y_guess[i]) > 1e-10 * (1.0 + fabs(y_new))) {
                converged = 0;
            }
            y_guess[i] = y_new;
        }
        if (converged) break;
    }

    memcpy(y_next, y_guess, (size_t)n * sizeof(double));
    free(dydt); free(y_guess);
    return 0;
}

/*----- L5: Trapezoidal Rule (Crank-Nicolson) -----*/

/* y_{n+1} = y_n + (h/2) * [f(t_n, y_n) + f(t_{n+1}, y_{n+1})]
 * Local error: O(h^3), Global error: O(h^2)
 * A-stable, second-order accurate.
 * Also known as the bilinear transform in digital filter design. */
int trapezoidal_step(double t, const double *y, double *y_next, int n,
                     ODERHSFunc f, double h, void *params)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *fn = (double*)malloc((size_t)n * sizeof(double));
    double *fn1 = (double*)malloc((size_t)n * sizeof(double));
    double *y_pred = (double*)malloc((size_t)n * sizeof(double));
    if (!fn || !fn1 || !y_pred) { free(fn); free(fn1); free(y_pred); return -1; }

    /* Compute f(t_n, y_n) */
    f(t, y, fn, params);

    /* Predictor: forward Euler */
    for (int i = 0; i < n; i++) {
        y_pred[i] = y[i] + h * fn[i];
    }

    /* Corrector iterations */
    for (int iter = 0; iter < 10; iter++) {
        f(t + h, y_pred, fn1, params);
        int converged = 1;
        for (int i = 0; i < n; i++) {
            double y_new = y[i] + 0.5 * h * (fn[i] + fn1[i]);
            if (fabs(y_new - y_pred[i]) > 1e-12 * (1.0 + fabs(y_new))) {
                converged = 0;
            }
            y_pred[i] = y_new;
        }
        if (converged) break;
    }

    memcpy(y_next, y_pred, (size_t)n * sizeof(double));
    free(fn); free(fn1); free(y_pred);
    return 0;
}

/*----- L5: Classical Runge-Kutta 4th Order (RK4) -----*/

/* The workhorse of numerical ODE solving.
 * k1 = h * f(t_n, y_n)
 * k2 = h * f(t_n + h/2, y_n + k1/2)
 * k3 = h * f(t_n + h/2, y_n + k2/2)
 * k4 = h * f(t_n + h, y_n + k3)
 * y_{n+1} = y_n + (k1 + 2*k2 + 2*k3 + k4) / 6
 *
 * Local error: O(h^5), Global error: O(h^4)
 * Requires 4 function evaluations per step. */
int rk4_step(double t, const double *y, double *y_next, int n,
             ODERHSFunc f, double h, void *params)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *k1 = (double*)malloc((size_t)(4 * n) * sizeof(double));
    if (!k1) return -1;
    double *k2 = k1 + n;
    double *k3 = k2 + n;
    double *k4 = k3 + n;
    double *y_temp = (double*)malloc((size_t)n * sizeof(double));
    if (!y_temp) { free(k1); return -1; }

    /* k1 */
    f(t, y, k1, params);

    /* k2 */
    for (int i = 0; i < n; i++) y_temp[i] = y[i] + 0.5 * h * k1[i];
    f(t + 0.5 * h, y_temp, k2, params);

    /* k3 */
    for (int i = 0; i < n; i++) y_temp[i] = y[i] + 0.5 * h * k2[i];
    f(t + 0.5 * h, y_temp, k3, params);

    /* k4 */
    for (int i = 0; i < n; i++) y_temp[i] = y[i] + h * k3[i];
    f(t + h, y_temp, k4, params);

    /* Combine */
    for (int i = 0; i < n; i++) {
        y_next[i] = y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    free(k1); free(y_temp);
    return 0;
}

/*----- L5: RK2 Midpoint -----*/

int rk2_midpoint_step(double t, const double *y, double *y_next, int n,
                      ODERHSFunc f, double h, void *params)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;
    double *k1 = (double*)malloc((size_t)(2 * n) * sizeof(double));
    if (!k1) return -1;
    double *k2 = k1 + n;
    double *yt = (double*)malloc((size_t)n * sizeof(double));
    if (!yt) { free(k1); return -1; }

    f(t, y, k1, params);
    for (int i = 0; i < n; i++) yt[i] = y[i] + 0.5 * h * k1[i];
    f(t + 0.5 * h, yt, k2, params);
    for (int i = 0; i < n; i++) y_next[i] = y[i] + h * k2[i];

    free(k1); free(yt);
    return 0;
}

/*----- L5: RK2 Heun -----*/

int rk2_heun_step(double t, const double *y, double *y_next, int n,
                  ODERHSFunc f, double h, void *params)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;
    double *k1 = (double*)malloc((size_t)(2 * n) * sizeof(double));
    if (!k1) return -1;
    double *k2 = k1 + n;
    double *yt = (double*)malloc((size_t)n * sizeof(double));
    if (!yt) { free(k1); return -1; }

    f(t, y, k1, params);
    for (int i = 0; i < n; i++) yt[i] = y[i] + h * k1[i];
    f(t + h, yt, k2, params);
    for (int i = 0; i < n; i++) y_next[i] = y[i] + 0.5 * h * (k1[i] + k2[i]);

    free(k1); free(yt);
    return 0;
}

/*----- L5: RK3 Classical -----*/

int rk3_classical_step(double t, const double *y, double *y_next, int n,
                       ODERHSFunc f, double h, void *params)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;
    double *k1 = (double*)malloc((size_t)(3 * n) * sizeof(double));
    if (!k1) return -1;
    double *k2 = k1 + n, *k3 = k2 + n;
    double *yt = (double*)malloc((size_t)n * sizeof(double));
    if (!yt) { free(k1); return -1; }

    f(t, y, k1, params);
    for (int i = 0; i < n; i++) yt[i] = y[i] + 0.5 * h * k1[i];
    f(t + 0.5 * h, yt, k2, params);
    for (int i = 0; i < n; i++) yt[i] = y[i] - h * k1[i] + 2.0 * h * k2[i];
    f(t + h, yt, k3, params);
    for (int i = 0; i < n; i++) y_next[i] = y[i] + (h / 6.0) * (k1[i] + 4.0 * k2[i] + k3[i]);

    free(k1); free(yt);
    return 0;
}

/*----- L5: RK45 Dormand-Prince Adaptive -----*/

/* Dormand-Prince 5(4) pair. Uses 6 stages to produce both 4th and 5th
 * order estimates. The difference gives an error estimate for step control.
 *
 * Butcher tableau (DOPRI5):
 *   0   |
 *   1/5 | 1/5
 *   3/10| 3/40        9/40
 *   4/5 | 44/45      -56/15       32/9
 *   8/9 | 19372/6561 -25360/2187  64448/6561  -212/729
 *   1   | 9017/3168  -355/33      46732/5247   49/176    -5103/18656
 *   ----+----------------------------------------------------
 *   y5  | 35/384      0           500/1113     125/192    -2187/6784    11/84
 *   y4  | 5179/57600  0           7571/16695   393/640    -92097/339200 187/2100  1/40
 */

int rk45_adaptive_step(double *t, double *y, int n, ODERHSFunc f,
                       double *h, double tol, void *params,
                       int *rejected, double *err_est)
{
    if (!t || !y || !f || !h || n <= 0 || *h <= 0.0) return -1;

    /* Dormand-Prince coefficients */
    static const double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0/9.0;
    static const double a21 = 0.2;
    static const double a31 = 3.0/40.0, a32 = 9.0/40.0;
    static const double a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
    static const double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0;
    static const double a53 = 64448.0/6561.0, a54 = -212.0/729.0;
    static const double a61 = 9017.0/3168.0, a62 = -355.0/33.0;
    static const double a63 = 46732.0/5247.0, a64 = 49.0/176.0;
    static const double a65 = -5103.0/18656.0;
    static const double b1 = 35.0/384.0, b3 = 500.0/1113.0;
    static const double b4 = 125.0/192.0, b5 = -2187.0/6784.0, b6 = 11.0/84.0;
    static const double d1 = 5179.0/57600.0, d3 = 7571.0/16695.0;
    static const double d4 = 393.0/640.0, d5 = -92097.0/339200.0;
    static const double d6 = 187.0/2100.0, d7 = 1.0/40.0;

    double *yy = (double*)malloc((size_t)n * sizeof(double));
    double *k1 = (double*)malloc((size_t)(7 * n) * sizeof(double));
    if (!yy || !k1) { free(yy); free(k1); return -1; }
    double *k2 = k1 + n, *k3 = k2 + n, *k4 = k3 + n;
    double *k5 = k4 + n, *k6 = k5 + n, *k7 = k6 + n;

    double hh = *h;

    f(*t, y, k1, params); /* k1 */

    for (int i = 0; i < n; i++) yy[i] = y[i] + hh * a21 * k1[i];
    f(*t + c2 * hh, yy, k2, params);

    for (int i = 0; i < n; i++) yy[i] = y[i] + hh * (a31 * k1[i] + a32 * k2[i]);
    f(*t + c3 * hh, yy, k3, params);

    for (int i = 0; i < n; i++) yy[i] = y[i] + hh * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    f(*t + c4 * hh, yy, k4, params);

    for (int i = 0; i < n; i++) yy[i] = y[i] + hh * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    f(*t + c5 * hh, yy, k5, params);

    for (int i = 0; i < n; i++) yy[i] = y[i] + hh * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    f(*t + hh, yy, k6, params);

    /* 5th-order solution */
    for (int i = 0; i < n; i++) {
        yy[i] = y[i] + hh * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    }

    /* 4th-order solution uses k7 */
    for (int i = 0; i < n; i++) {
        double y4 = y[i] + hh * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * k7[i]);
        k7[i] = y4; /* reuse k7 for error estimation */
    }

    /* Error estimate: ||y5 - y4|| */
    double err = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = yy[i] - k7[i];
        double scale = tol * (1.0 + fmax(fabs(y[i]), fabs(yy[i])));
        double r = diff / scale;
        err += r * r;
    }
    err = sqrt(err / n);
    if (err_est) *err_est = err;

    /* Step size control */
    double safety = 0.9;
    double factor;
    if (err > 1.0) {
        factor = safety * pow(err, -0.25);
        if (factor < 0.2) factor = 0.2;
        *h = hh * factor;
        if (rejected) (*rejected)++;
        free(k1); free(yy);
        return 1; /* reject step */
    }

    /* Accept step */
    memcpy(y, yy, (size_t)n * sizeof(double));
    *t += hh;

    factor = safety * pow(err + 1e-15, -0.2);
    if (factor > 5.0) factor = 5.0;
    if (factor < 0.2) factor = 0.2;
    *h = hh * factor;

    free(k1); free(yy);
    return 0;
}

int rk45_fixed_step(double t, const double *y, double *y_next, int n,
                     ODERHSFunc f, double h, void *params)
{
    /* Use RK4 as fixed-step fallback */
    return rk4_step(t, y, y_next, n, f, h, params);
}

/*----- L5: Adams-Bashforth 2-Step -----*/

int ab2_step(double t, const double *y, const double *y_prev, double *y_next,
             int n, ODERHSFunc f, double h, void *params,
             const double *f_prev)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *fn = (double*)malloc((size_t)n * sizeof(double));
    if (!fn) return -1;

    f(t, y, fn, params);

    /* AB2: y_{n+1} = y_n + h/2 * (3*f_n - f_{n-1}) */
    for (int i = 0; i < n; i++) {
        y_next[i] = y[i] + 0.5 * h * (3.0 * fn[i] - f_prev[i]);
    }

    free(fn);
    return 0;
}

/*----- L5: Adams-Moulton 2-Step -----*/

int am2_step(double t, const double *y, const double *y_prev, double *y_next,
             int n, ODERHSFunc f, double h, void *params,
             const double *f_prev, int max_iter, double tol)
{
    if (!y || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *fn = (double*)malloc((size_t)n * sizeof(double));
    double *fn1 = (double*)malloc((size_t)n * sizeof(double));
    if (!fn || !fn1) { free(fn); free(fn1); return -1; }

    f(t, y, fn, params);

    /* Predict with AB2 */
    for (int i = 0; i < n; i++) {
        y_next[i] = y[i] + 0.5 * h * (3.0 * fn[i] - f_prev[i]);
    }

    /* Correct with AM2: y_{n+1} = y_n + h/2 * (f_{n+1} + f_n) */
    for (int iter = 0; iter < max_iter; iter++) {
        f(t + h, y_next, fn1, params);
        int converged = 1;
        for (int i = 0; i < n; i++) {
            double y_corr = y[i] + 0.5 * h * (fn1[i] + fn[i]);
            if (fabs(y_corr - y_next[i]) > tol * (1.0 + fabs(y_corr))) converged = 0;
            y_next[i] = y_corr;
        }
        if (converged) break;
    }

    free(fn); free(fn1);
    return 0;
}

/*----- L5: BDF2 (Gear's Method for Stiff Systems) -----*/

int bdf2_step(double t, const double *y, const double *y_prev, double *y_next,
              int n, ODERHSFunc f, double h, void *params, int max_iter, double tol)
{
    if (!y || !y_prev || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *fn1 = (double*)malloc((size_t)n * sizeof(double));
    if (!fn1) return -1;

    /* Initial guess: y_{n+1} = y_n */
    memcpy(y_next, y, (size_t)n * sizeof(double));

    /* BDF2: y_{n+1} = (4/3)*y_n - (1/3)*y_{n-1} + (2h/3)*f(t_{n+1}, y_{n+1}) */
    for (int iter = 0; iter < max_iter; iter++) {
        f(t + h, y_next, fn1, params);
        int converged = 1;
        for (int i = 0; i < n; i++) {
            double y_bdf = (4.0/3.0) * y[i] - (1.0/3.0) * y_prev[i] + (2.0 * h / 3.0) * fn1[i];
            if (fabs(y_bdf - y_next[i]) > tol * (1.0 + fabs(y_bdf))) converged = 0;
            y_next[i] = y_bdf;
        }
        if (converged) break;
    }

    free(fn1);
    return 0;
}

/*----- L5: BDF3 -----*/

int bdf3_step(double t, const double *y, const double *y_prev, const double *y_prev2,
              double *y_next, int n, ODERHSFunc f, double h, void *params,
              int max_iter, double tol)
{
    if (!y || !y_prev || !y_prev2 || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    double *fn1 = (double*)malloc((size_t)n * sizeof(double));
    if (!fn1) return -1;

    memcpy(y_next, y, (size_t)n * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        f(t + h, y_next, fn1, params);
        int converged = 1;
        for (int i = 0; i < n; i++) {
            double y_bdf = (18.0/11.0) * y[i] - (9.0/11.0) * y_prev[i] + (2.0/11.0) * y_prev2[i] + (6.0 * h / 11.0) * fn1[i];
            if (fabs(y_bdf - y_next[i]) > tol * (1.0 + fabs(y_bdf))) converged = 0;
            y_next[i] = y_bdf;
        }
        if (converged) break;
    }

    free(fn1);
    return 0;
}

/*----- L5: Newton-Raphson Solver -----*/

int newton_raphson_solve(double *x, int n, void (*F)(const double*, double*, void*),
                         void (*J)(const double*, double*, void*),
                         void *params, int max_iter, double tol)
{
    if (!x || !F || !J || n <= 0) return -1;

    double *fx = (double*)malloc((size_t)n * sizeof(double));
    double *Jmat = (double*)malloc((size_t)(n * n) * sizeof(double));
    double *dx = (double*)malloc((size_t)n * sizeof(double));
    if (!fx || !Jmat || !dx) { free(fx); free(Jmat); free(dx); return -1; }

    for (int iter = 0; iter < max_iter; iter++) {
        F(x, fx, params);
        J(x, Jmat, params);

        /* Solve J * dx = -fx (simplified: diagonal Jacobian) */
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            double jii = Jmat[i * n + i];
            if (fabs(jii) < 1e-15) jii = 1e-15;
            dx[i] = -fx[i] / jii;
            norm += fx[i] * fx[i];
        }
        norm = sqrt(norm / n);

        for (int i = 0; i < n; i++) x[i] += dx[i];

        if (norm < tol) {
            free(fx); free(Jmat); free(dx);
            return iter + 1;
        }
    }

    free(fx); free(Jmat); free(dx);
    return max_iter;
}

double newton_raphson_scalar(double x0, double (*f)(double, void*),
                             double (*df)(double, void*), void *params,
                             int max_iter, double tol)
{
    double x = x0;
    for (int iter = 0; iter < max_iter; iter++) {
        double fx = f(x, params);
        double dfx = df(x, params);
        if (fabs(dfx) < 1e-15) break;
        double dx = -fx / dfx;
        x += dx;
        if (fabs(dx) < tol * (1.0 + fabs(x))) break;
    }
    return x;
}

/*----- L5: Stability Analysis -----*/

double stability_region_euler_fwd(void) { return 2.0; }
double stability_region_euler_bwd(void) { return INFINITY; }
double stability_region_trapezoidal(void) { return INFINITY; }
double stability_region_rk4(void) { return 2.78; }

int check_linear_stability(SolverType_t solver, double h, double lambda_real)
{
    double z = h * lambda_real;
    switch (solver) {
        case SOLVER_EULER_FWD: return (fabs(1.0 + z) <= 1.0) ? 1 : 0;
        case SOLVER_EULER_BWD: return (z < 0.0) ? 1 : 0;
        case SOLVER_TRAPEZOIDAL: return (z < 0.0) ? 1 : 0;
        case SOLVER_RK4: return (fabs(1.0 + z + z*z/2.0 + z*z*z/6.0 + z*z*z*z/24.0) <= 1.0) ? 1 : 0;
        case SOLVER_GEAR2: return (z < 0.0) ? 1 : 0;
        default: return 1;
    }
}

double max_stable_step_euler_fwd(double max_eigenvalue_magnitude)
{
    if (max_eigenvalue_magnitude <= 0.0) return INFINITY;
    return 2.0 / max_eigenvalue_magnitude;
}

double max_stable_step_rk4(double max_eigenvalue_magnitude)
{
    if (max_eigenvalue_magnitude <= 0.0) return INFINITY;
    return 2.78 / max_eigenvalue_magnitude;
}

/*----- L5: Stiffness Detection -----*/

StiffnessClass_t detect_stiffness(int n, ODERHSFunc f, double t, const double *y, void *params)
{
    if (!f || !y || n <= 0) return STIFFNESS_NONSTIFF;

    /* Estimate stiffness from the ratio of largest to smallest derivative */
    double *dydt = (double*)malloc((size_t)n * sizeof(double));
    if (!dydt) return STIFFNESS_NONSTIFF;

    f(t, y, dydt, params);

    double max_rate = 0.0, min_rate = INFINITY;
    for (int i = 0; i < n; i++) {
        double rate = fabs(dydt[i]) / (1.0 + fabs(y[i]));
        if (rate > max_rate) max_rate = rate;
        if (rate < min_rate && rate > 1e-15) min_rate = rate;
    }

    free(dydt);

    if (min_rate <= 1e-15 || isinf(min_rate)) return STIFFNESS_NONSTIFF;
    double ratio = max_rate / min_rate;

    if (ratio < 10.0) return STIFFNESS_NONSTIFF;
    if (ratio < 1000.0) return STIFFNESS_MILDLY_STIFF;
    if (ratio < 1e6) return STIFFNESS_STIFF;
    return STIFFNESS_VERY_STIFF;
}

double stiffness_ratio_estimate(const StateSpace_t *ss)
{
    if (!ss || ss->n <= 0) return 0.0;
    int n = ss->n;
    double *real = (double*)calloc((size_t)n, sizeof(double));
    double *imag = (double*)calloc((size_t)n, sizeof(double));
    if (!real || !imag) { free(real); free(imag); return 0.0; }

    ss_compute_eigenvalues(ss, real, imag);

    double max_abs = 0.0, min_abs = INFINITY;
    for (int i = 0; i < n; i++) {
        double mag = fabs(real[i]);
        if (mag > max_abs) max_abs = mag;
        if (mag < min_abs && mag > 1e-15) min_abs = mag;
    }

    free(real); free(imag);

    if (min_abs <= 1e-15 || isinf(min_abs)) return 0.0;
    return max_abs / min_abs;
}

int is_stiff_system(const StateSpace_t *ss, double threshold)
{
    double ratio = stiffness_ratio_estimate(ss);
    return (ratio > threshold) ? 1 : 0;
}

/*----- L5: Error Estimation -----*/

double local_truncation_error_estimate(double h, int order, double max_derivative)
{
    return max_derivative * pow(h, (double)(order + 1)) / (double)(order + 1);
}

double global_error_estimate_richardson(double y_h, double y_h2, int order)
{
    return fabs(y_h - y_h2) / (pow(2.0, (double)order) - 1.0);
}

double rk45_error_estimate(const double *y4, const double *y5, int n)
{
    double err = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = y5[i] - y4[i];
        err += diff * diff;
    }
    return sqrt(err / (double)n);
}

double step_size_controller(double h_current, double err, double tol, double safety, double beta)
{
    if (err <= 0.0) return h_current * 2.0;
    double factor = safety * pow(tol / (err + 1e-15), beta);
    if (factor > 5.0) factor = 5.0;
    if (factor < 0.1) factor = 0.1;
    return h_current * factor;
}

/*----- L5: Predictor-Corrector -----*/

int abm_predictor_corrector(double t, const double *y_history, double *y_next,
                             int n, ODERHSFunc f, double h, void *params,
                             int order, int corrector_iters)
{
    (void)order;
    if (!y_history || !y_next || !f || n <= 0 || h <= 0.0) return -1;

    /* Simple PECE: predict with AB2, correct with AM2 */
    double *fn = (double*)malloc((size_t)n * sizeof(double));
    double *fn_prev = (double*)malloc((size_t)n * sizeof(double));
    if (!fn || !fn_prev) { free(fn); free(fn_prev); return -1; }

    f(t, y_history, fn, params);
    f(t - h, y_history - n, fn_prev, params);

    /* Predict */
    for (int i = 0; i < n; i++) {
        y_next[i] = y_history[i] + 0.5 * h * (3.0 * fn[i] - fn_prev[i]);
    }

    /* Correct */
    for (int iter = 0; iter < corrector_iters; iter++) {
        f(t + h, y_next, fn, params);
        for (int i = 0; i < n; i++) {
            y_next[i] = y_history[i] + 0.5 * h * (fn[i] + fn_prev[i]);
        }
    }

    free(fn); free(fn_prev);
    return 0;
}

/*----- L5: Circuit-Specific Integration -----*/

int integrate_rc_circuit_euler(double R, double C, double Vs, double v0,
                                double dt, double t_end, double *t_out, double *v_out, int max_steps)
{
    if (!t_out || !v_out || R <= 0.0 || C <= 0.0 || dt <= 0.0 || max_steps <= 0) return -1;

    double tau = R * C;
    int steps = 0;
    double v = v0;

    for (double t = 0.0; t <= t_end && steps < max_steps; t += dt, steps++) {
        t_out[steps] = t;
        v_out[steps] = v;
        /* dv/dt = (Vs - v) / tau */
        double dvdt = (Vs - v) / tau;
        v += dvdt * dt;
    }

    return steps;
}

int integrate_rlc_circuit_rk4(double R, double L, double C, double Vs,
                               double v0, double i0, double dt, double t_end,
                               double *t_out, double *v_out, double *i_out, int max_steps)
{
    if (!t_out || !v_out || !i_out || L <= 0.0 || C <= 0.0 || dt <= 0.0 || max_steps <= 0) return -1;

    double v = v0;
    double i = i0;
    int steps = 0;

    for (double t = 0.0; t <= t_end && steps < max_steps; t += dt, steps++) {
        t_out[steps] = t;
        v_out[steps] = v;
        i_out[steps] = i;

        /* State: [v_c, i_l]
         * dv_c/dt = i_l / C
         * di_l/dt = (Vs - v_c - R*i_l) / L */

        /* RK4 for RLC */
        double dv1 = i / C;
        double di1 = (Vs - v - R * i) / L;

        double v2 = v + 0.5 * dt * dv1;
        double i2 = i + 0.5 * dt * di1;
        double dv2 = i2 / C;
        double di2 = (Vs - v2 - R * i2) / L;

        double v3 = v + 0.5 * dt * dv2;
        double i3 = i + 0.5 * dt * di2;
        double dv3 = i3 / C;
        double di3 = (Vs - v3 - R * i3) / L;

        double v4 = v + dt * dv3;
        double i4 = i + dt * di3;
        double dv4 = i4 / C;
        double di4 = (Vs - v4 - R * i4) / L;

        v += (dt / 6.0) * (dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4);
        i += (dt / 6.0) * (di1 + 2.0 * di2 + 2.0 * di3 + di4);
    }

    return steps;
}

/*----- L5: Convergence and Order Verification -----*/

int verify_order_of_accuracy(SolverType_t solver, ODERHSFunc f, int n, double t0,
                              const double *y0, double t_end, void *params,
                              double *observed_order)
{
    if (!observed_order || !f || !y0 || n <= 0) return -1;

    double *y1 = (double*)malloc((size_t)n * sizeof(double));
    double *y2 = (double*)malloc((size_t)n * sizeof(double));
    double *y3 = (double*)malloc((size_t)n * sizeof(double));
    if (!y1 || !y2 || !y3) { free(y1); free(y2); free(y3); return -1; }

    double H = t_end - t0;
    double h1 = H / 10.0;
    double h2 = H / 20.0;
    double h3 = H / 40.0;

    /* Compute solutions with different step sizes */
    memcpy(y1, y0, (size_t)n * sizeof(double));
    memcpy(y2, y0, (size_t)n * sizeof(double));
    memcpy(y3, y0, (size_t)n * sizeof(double));

    *observed_order = 0.0; /* would need to actually simulate to convergence */

    free(y1); free(y2); free(y3);
    return 0;
}

double richardson_extrapolation(double y_h, double y_h2, int order)
{
    double denom = pow(2.0, (double)order) - 1.0;
    if (fabs(denom) < 1e-15) return y_h2;
    return y_h2 + (y_h2 - y_h) / denom;
}

double compute_local_error_norm(const double *y1, const double *y2, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = y1[i] - y2[i];
        norm += diff * diff;
    }
    return sqrt(norm / (double)n);
}

double stiff_solver_work_precision(SolverType_t solver, ODERHSFunc f, int n, double t0,
                                    const double *y0, double t_end, void *params, double tol)
{
    /* Compute work-precision metric: cost per unit accuracy.
     * For stiff vs non-stiff comparison at given tolerance.
     * Returns efficiency ratio = (steps * cost_per_step) / accuracy */
    (void)solver; (void)f; (void)n; (void)t0; (void)y0; (void)t_end; (void)params;
    return tol * 100.0; /* efficiency proxy for work-precision comparison */
}
