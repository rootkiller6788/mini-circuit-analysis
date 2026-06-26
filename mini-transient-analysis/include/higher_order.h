#ifndef HIGHER_ORDER_H
#define HIGHER_ORDER_H
#include "transient_defs.h"

/* L3-L5: Higher-Order Systems and State-Space Methods
 * Reference: Ogata "Modern Control Engineering" (2010), Kailath "Linear Systems" (1980)
 *
 * Higher-order systems have 3+ energy storage elements.
 * State-space representation: dx/dt = A*x + B*u, y = C*x + D*u
 */

typedef struct {
    int n;
    int m;
    int p;
    double *A;
    double *B;
    double *C;
    double *D;
    double *x;
} StateSpace_t;

StateSpace_t* ss_alloc(int n, int m, int p);
void ss_free(StateSpace_t *ss);
void ss_set_identity(StateSpace_t *ss);
void ss_set_controller_canonical(StateSpace_t *ss, const double *den, const double *num, int order);
void ss_set_observer_canonical(StateSpace_t *ss, const double *den, const double *num, int order);
int ss_step_euler(StateSpace_t *ss, const double *u, double dt);
int ss_step_rk4(StateSpace_t *ss, const double *u, double dt);
void ss_compute_output(const StateSpace_t *ss, const double *u, double *y);
void ss_compute_eigenvalues(const StateSpace_t *ss, double *real, double *imag);
int ss_is_stable(const StateSpace_t *ss);
int ss_is_controllable(const StateSpace_t *ss);
int ss_is_observable(const StateSpace_t *ss);
double ss_dominant_time_constant(const StateSpace_t *ss);
void ss_to_transfer_function(const StateSpace_t *ss, double *num, double *den);
void transfer_function_to_ss(const double *num, int num_order, const double *den, int den_order, StateSpace_t *ss);

/* Modal Analysis */
int modal_decomposition(const StateSpace_t *ss, double *eigenvalues_real, double *eigenvalues_imag, double *eigenvectors);
double modal_participation_factor(const StateSpace_t *ss, int mode);
double modal_dominance_index(const StateSpace_t *ss, int mode);

/* Pole-Zero Analysis */
PoleZeroMap_t pole_zero_from_transfer_function(const double *num, int n_num, const double *den, int n_den);
int pole_zero_is_minimum_phase(const PoleZeroMap_t *pz);
double pole_zero_dominant_pole(const PoleZeroMap_t *pz);
int pole_zero_reduce_order(const PoleZeroMap_t *pz, PoleZeroMap_t *reduced, double tol);

/* Eigenvalue Sensitivity */
double eigenvalue_sensitivity_to_parameter(const StateSpace_t *ss, int row, int col, int mode);
double eigenvalue_sensitivity_to_damping(const StateSpace_t *ss, double delta_zeta, int mode);

/* Stability Margin */
double gain_margin_from_poles(const StateSpace_t *ss, double K);
double phase_margin_from_poles(const PoleZeroMap_t *pz, double omega_c);
double delay_margin(const StateSpace_t *ss);
int routh_hurwitz_stability(const double *coeffs, int order);

/* Higher-Order System Reduction */
int model_order_reduction_pade(const double *num, int n_num, const double *den, int n_den, int reduced_order, double *r_num, double *r_den);
int model_order_reduction_modal(const StateSpace_t *ss, int keep_modes, StateSpace_t *reduced);
double hankel_singular_value(const StateSpace_t *ss, int k);

/* Third-Order Transient Analysis */
double third_order_step_response(double K, double p1, double p2, double p3, double t);
double third_order_impulse_response(double K, double p1, double p2, double p3, double t);
double third_order_dominant_pole_approx(double K, double p_dominant, double t);
int third_order_valid_dominant_approx(double p_dominant, double p2, double p3);

/* Distributed Parameter Circuits */
double rc_transmission_line_step(double V0, double R_total, double C_total, double x_norm, double t);
double rc_tline_elmore_delay(double R_total, double C_total);
double rc_tline_rise_time(double R_total, double C_total);

/* State Transition Matrix */
int compute_state_transition_matrix(const StateSpace_t *ss, double dt, double *Phi);
int compute_discretized_ss(const StateSpace_t *ss_cont, StateSpace_t *ss_disc, double Ts);
int discretize_zoh(const double *A, int n, const double *B, double Ts, double *Ad, double *Bd);

/* Coupled Oscillators */
void coupled_oscillator_params(double L1, double C1, double L2, double C2, double Lm, double *omega1, double *omega2, double *omega_beat);
double coupled_oscillator_transient(double V0, double L1, double C1, double L2, double C2, double Lm, double t);

#endif /* HIGHER_ORDER_H */
