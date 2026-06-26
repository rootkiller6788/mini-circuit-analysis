#ifndef FIRST_ORDER_H
#define FIRST_ORDER_H
#include "transient_defs.h"

/*==================================================================
 * L2-L4: First-Order Transient Analysis
 *
 * L2 Core Concepts:
 *   - Time constant tau (tau = RC or tau = L/R)
 *   - Natural response: y_n(t) = y(0) * exp(-t/tau)
 *   - Forced response: steady-state response to input
 *   - Complete response = natural + forced
 *
 * L4 Fundamental Laws:
 *   - RC circuit ODE: RC * dv_c/dt + v_c = v_s(t)  (Ohm + Maxwell)
 *   - RL circuit ODE: (L/R) * di_l/dt + i_l = v_s(t)/R  (Ohm + Faraday)
 *   - Energy continuity: v_c(0+) = v_c(0-), i_l(0+) = i_l(0-)
 *   - Superposition: complete = zero-input + zero-state
 *
 * Reference: Hayt & Kemmerly Ch.8-9, Nilsson & Riedel Ch.7-8
 *==================================================================*/

/*----- L2: Time Constants -----*/
double rc_time_constant(double R, double C);
double rl_time_constant(double R, double L);
double rc_time_constant_parallel(double R1, double R2, double C);
double rl_time_constant_parallel(double R1, double R2, double L);
double thevenin_time_constant(double R_th, double C_eq);
double norton_time_constant(double R_n, double L_eq);

/*----- L4: RC Step Response (Charging)-----*/
double rc_step_charge_vc(double Vs, double v0, double R, double C, double t);
double rc_step_charge_i(double Vs, double v0, double R, double C, double t);
double rc_step_charge_vr(double Vs, double v0, double R, double C, double t);
double rc_step_charge_power_R(double Vs, double v0, double R, double C, double t);
double rc_step_charge_energy_R(double Vs, double v0, double R, double C, double t);
double rc_step_charge_energy_C(double Vs, double v0, double R, double C, double t);

/*----- L4: RC Natural Response (Discharging) -----*/
double rc_natural_vc(double v0, double R, double C, double t);
double rc_natural_i(double v0, double R, double C, double t);
double rc_natural_vr(double v0, double R, double C, double t);
double rc_natural_power_R(double v0, double R, double C, double t);
double rc_natural_energy_dissipated(double v0, double R, double C, double t);

/*----- L4: RC Complete Response (Step + Natural) -----*/
double rc_complete_vc(double Vs, double v0, double R, double C, double t);
double rc_complete_i(double Vs, double v0, double R, double C, double t);

/*----- L4: RL Step Response -----*/
double rl_step_charge_il(double Vs, double i0, double R, double L, double t);
double rl_step_charge_vl(double Vs, double i0, double R, double L, double t);
double rl_step_charge_vr(double Vs, double i0, double R, double L, double t);
double rl_step_charge_energy_L(double Vs, double i0, double R, double L, double t);
double rl_natural_il(double i0, double R, double L, double t);
double rl_natural_vl(double i0, double R, double L, double t);
double rl_complete_il(double Vs, double i0, double R, double L, double t);
double rl_complete_vl(double Vs, double i0, double R, double L, double t);

/*----- L4: RC Impulse Response -----*/
double rc_impulse_vc(double R, double C, double t);
double rc_impulse_i(double R, double C, double t);
double rl_impulse_il(double R, double L, double t);
double rl_impulse_vl(double R, double L, double t);

/*----- L4: RC Ramp Response -----*/
double rc_ramp_vc(double slope, double R, double C, double t, double v0);
double rc_ramp_i(double slope, double R, double C, double t, double v0);
double rl_ramp_il(double slope, double R, double L, double t, double i0);
double rl_ramp_vl(double slope, double R, double L, double t, double i0);

/*----- L6: RC Pulse Response -----*/
double rc_pulse_vc(double V_amplitude, double R, double C, double t, double t_pulse);
double rc_pulse_train_vc(double V_amplitude, double R, double C, double t,
                          double t_period, double duty_cycle, int n_periods);
double rl_pulse_il(double V_amplitude, double R, double L, double t, double t_pulse);

/*----- L4: First-Order Sinusoidal Response -----*/
double rc_sinusoidal_vc(double Vm, double omega, double phi, double R, double C, double t, double v0);
double rc_sinusoidal_i(double Vm, double omega, double phi, double R, double C, double t, double v0);
double rl_sinusoidal_il(double Vm, double omega, double phi, double R, double L, double t, double i0);

/*----- L4: Switching at Non-Zero Time -----*/
double rc_switched_vc(double Vs, double v0, double R, double C, double t, double t_switch);
double rl_switched_il(double Vs, double i0, double R, double L, double t, double t_switch);

/*----- L4: First-Order with Initial Conditions General -----*/
double first_order_general_response(double a, double b, double y0, double t);
double first_order_step_ss_error(double tau, double K);
double first_order_bandwidth_from_tau(double tau);
double tau_from_rise_time(double tr_10_90);

/*----- L4: Sequential Switching -----*/
double rc_sequential_vc(double V1, double V2, double R1, double R2,
                         double C, double t, double t_switch);
double rc_sequential_charge_then_discharge_vc(double Vs, double R, double C,
                                               double t, double t_hold);

/*----- L4: Multiple Time Constants (Cascaded RC) -----*/
double rc_cascaded_vc(double Vs, double R1, double C1, double R2, double C2,
                      double t, double v10, double v20);
double dominant_time_constant(double tau1, double tau2);
double effective_time_constant_cascaded(double tau1, double tau2);

/*----- L7: RC Timing Applications -----*/
double rc_delay_time(double R, double C, double V_target, double V_supply);
double rc_threshold_crossing_time(double R, double C, double v_start,
                                   double v_target, double v_final);
double rc_startup_time(double R, double C, double V_threshold, double V_supply);
double por_timeout_compute(double R, double C, double Vcc, double V_threshold);
double watchdog_timeout_compute(double R, double C, double V_threshold, double Vcc);
double debounce_settling_time_compute(double R, double C, double V_cc, double V_IH_min);

/*----- L7: RC Oscillator -----*/
double rc_oscillator_period(double R, double C, double V_high, double V_low, double V_supply);
double rc_oscillator_frequency(double R, double C, double V_high, double V_low, double V_supply);
double rc_relaxation_oscillator_freq(double R, double C, double V_th_ratio);

/*----- L4: Transient Response Bounds -----*/
double rc_max_current_instantaneous(double Vs, double v0, double R);
double rc_max_power_instantaneous(double Vs, double v0, double R);
double rl_max_voltage_instantaneous(double Vs, double i0, double R);

/*----- L7: Exponential Curve Fitting (Measure tau from data) -----*/
double estimate_tau_from_decay(double v1, double v2, double dt);
double predict_voltage_at_time(double v_measured, double t_measured,
                                double tau, double t_target);

#endif /* FIRST_ORDER_H */
