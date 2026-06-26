#ifndef NUMERICAL_TRANSIENT_H
#define NUMERICAL_TRANSIENT_H
#include "transient_defs.h"
#include "higher_order.h"

/* L5: Numerical Methods for Transient Analysis
 * Reference: Press et al. "Numerical Recipes" (2007), Hairer & Wanner
 *            "Solving Ordinary Differential Equations" (1996),
 *            Chua & Lin "Computer-Aided Analysis of Electronic Circuits" (1975)
 *
 * Methods:
 *   - Forward Euler: y_{n+1} = y_n + h*f(t_n, y_n)         O(h)
 *   - Backward Euler: y_{n+1} = y_n + h*f(t_{n+1}, y_{n+1}) O(h) A-stable
 *   - Trapezoidal: y_{n+1} = y_n + h/2*(f_n + f_{n+1})     O(h^2) A-stable
 *   - RK4: Classic 4-stage Runge-Kutta                       O(h^4)
 *   - RK45: Dormand-Prince adaptive                          O(h^5)
 *   - Gear BDF2: y_{n+1} = (4/3)y_n - (1/3)y_{n-1} + (2h/3)*f_{n+1}  stiff
 *   - Adams-Bashforth 2: y_{n+1} = y_n + h/2*(3f_n - f_{n-1})  O(h^2)
 *   - Adams-Moulton 2: y_{n+1} = y_n + h/2*(f_{n+1} + f_n)    O(h^2)
 *
 * Key Concepts:
 *   - Consistency: local truncation error -> 0 as h -> 0
 *   - Stability: bounded solution for bounded initial conditions
 *   - A-stability: stable for all h*lambda with Re(lambda) < 0
 *   - Stiffness ratio: max|Re(lambda)| / min|Re(lambda)|
 */

/*----- L5: Single-Step Methods -----*/
int euler_forward_step(double t, const double *y, double *y_next, int n, ODERHSFunc f, double h, void *params);
int euler_backward_step(double t, const double *y, double *y_next, int n, ODERHSFunc f, double h, void *params);
int trapezoidal_step(double t, const double *y, double *y_next, int n, ODERHSFunc f, double h, void *params);

/*----- L5: Runge-Kutta Methods -----*/
int rk4_step(double t, const double *y, double *y_next, int n, ODERHSFunc f, double h, void *params);
int rk2_midpoint_step(double t, const double *y, double *y_next, int n, ODERHSFunc f, double h, void *params);
int rk2_heun_step(double t, const double *y, double *y_next, int n, ODERHSFunc f, double h, void *params);
int rk3_classical_step(double t, const double *y, double *y_next, int n, ODERHSFunc f, double h, void *params);

/*----- L5: Adaptive Step-Size RK45 (Dormand-Prince) -----*/
int rk45_adaptive_step(double *t, double *y, int n, ODERHSFunc f,
                       double *h, double tol, void *params,
                       int *rejected, double *err_est);
int rk45_fixed_step(double t, const double *y, double *y_next, int n,
                     ODERHSFunc f, double h, void *params);

/*----- L5: Multi-Step Methods -----*/
int ab2_step(double t, const double *y, const double *y_prev, double *y_next,
             int n, ODERHSFunc f, double h, void *params,
             const double *f_prev);
int am2_step(double t, const double *y, const double *y_prev, double *y_next,
             int n, ODERHSFunc f, double h, void *params,
             const double *f_prev, int max_iter, double tol);
int bdf2_step(double t, const double *y, const double *y_prev, double *y_next,
              int n, ODERHSFunc f, double h, void *params, int max_iter, double tol);
int bdf3_step(double t, const double *y, const double *y_prev, const double *y_prev2,
              double *y_next, int n, ODERHSFunc f, double h, void *params,
              int max_iter, double tol);

/*----- L5: Newton-Raphson for Implicit Methods -----*/
int newton_raphson_solve(double *x, int n, void (*F)(const double*, double*, void*),
                         void (*J)(const double*, double*, void*),
                         void *params, int max_iter, double tol);
double newton_raphson_scalar(double x0, double (*f)(double, void*),
                             double (*df)(double, void*), void *params,
                             int max_iter, double tol);

/*----- L5: Stability Analysis -----*/
double stability_region_euler_fwd(void);
double stability_region_euler_bwd(void);
double stability_region_trapezoidal(void);
double stability_region_rk4(void);
int check_linear_stability(SolverType_t solver, double h, double lambda_real);
double max_stable_step_euler_fwd(double max_eigenvalue_magnitude);
double max_stable_step_rk4(double max_eigenvalue_magnitude);

/*----- L5: Stiffness Detection -----*/
StiffnessClass_t detect_stiffness(int n, ODERHSFunc f, double t, const double *y, void *params);
double stiffness_ratio_estimate(const StateSpace_t *ss);
int is_stiff_system(const StateSpace_t *ss, double threshold);

/*----- L5: Error Estimation -----*/
double local_truncation_error_estimate(double h, int order, double max_derivative);
double global_error_estimate_richardson(double y_h, double y_h2, int order);
double rk45_error_estimate(const double *y4, const double *y5, int n);
double step_size_controller(double h_current, double err, double tol, double safety, double beta);

/*----- L5: Linear Multi-Step Predictor-Corrector -----*/
int abm_predictor_corrector(double t, const double *y_history, double *y_next,
                             int n, ODERHSFunc f, double h, void *params,
                             int order, int corrector_iters);

/*----- L5: Numerical Integration for Special Circuit Types -----*/
int integrate_rc_circuit_euler(double R, double C, double Vs, double v0,
                                double dt, double t_end, double *t_out, double *v_out, int max_steps);
int integrate_rlc_circuit_rk4(double R, double L, double C, double Vs,
                               double v0, double i0, double dt, double t_end,
                               double *t_out, double *v_out, double *i_out, int max_steps);

/*----- L5: Convergence and Order Verification -----*/
int verify_order_of_accuracy(SolverType_t solver, ODERHSFunc f, int n, double t0,
                              const double *y0, double t_end, void *params,
                              double *observed_order);
double richardson_extrapolation(double y_h, double y_h2, int order);
double compute_local_error_norm(const double *y1, const double *y2, int n);

/*----- L5: Stiff Solver Comparison -----*/
double stiff_solver_work_precision(SolverType_t solver, ODERHSFunc f, int n, double t0,
                                    const double *y0, double t_end, void *params, double tol);

#endif /* NUMERICAL_TRANSIENT_H */
