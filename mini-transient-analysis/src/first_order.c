#include "first_order.h"
#include "transient_defs.h"
#include <math.h>
#include <stdlib.h>

/*==================================================================
 * L2/L4: First-Order Transient Analysis Implementation
 *
 * The first-order ODE: tau * dy/dt + y = f(t)
 * General solution: y(t) = y_n(t) + y_f(t)
 *   y_n(t) = y(0) * exp(-t/tau)  [natural/zero-input response]
 *   y_f(t) = particular solution to forced equation [forced/zero-state]
 *
 * For step input f(t) = Vs * u(t):
 *   y(t) = y(0)*exp(-t/tau) + Vs*(1 - exp(-t/tau))
 *        = Vs + (y(0) - Vs)*exp(-t/tau)
 *==================================================================*/

/*----- L2: Time Constants -----*/

double rc_time_constant(double R, double C)
{
    if (R < 0.0 || C < 0.0) return -1.0;
    return R * C;
}

double rl_time_constant(double R, double L)
{
    if (R <= 0.0 || L < 0.0) return -1.0;
    return L / R;
}

/* Time constant for parallel RC: tau = (R1||R2)*C */
double rc_time_constant_parallel(double R1, double R2, double C)
{
    if (R1 <= 0.0 || R2 <= 0.0 || C < 0.0) return -1.0;
    double R_eq = (R1 * R2) / (R1 + R2);
    return R_eq * C;
}

/* Time constant for parallel RL: tau = L/(R1||R2) */
double rl_time_constant_parallel(double R1, double R2, double L)
{
    if (R1 <= 0.0 || R2 <= 0.0 || L < 0.0) return -1.0;
    double R_eq = (R1 * R2) / (R1 + R2);
    if (R_eq <= 0.0) return -1.0;
    return L / R_eq;
}

/* Thevenin time constant: tau = R_th * C_eq */
double thevenin_time_constant(double R_th, double C_eq)
{
    if (R_th < 0.0 || C_eq < 0.0) return -1.0;
    return R_th * C_eq;
}

/* Norton time constant: tau = L_eq / R_n */
double norton_time_constant(double R_n, double L_eq)
{
    if (R_n <= 0.0 || L_eq < 0.0) return -1.0;
    return L_eq / R_n;
}

/*----- L4: RC Step Response (Charging) -----*/

/* RC charging: v_c(t) = Vs + (v0 - Vs)*exp(-t/tau)
 * Reference: Hayt & Kemmerly, Eq. 9.12 */
double rc_step_charge_vc(double Vs, double v0, double R, double C, double t)
{
    if (t < 0.0) return v0;
    double tau = R * C;
    if (tau <= 0.0) return Vs;
    return Vs + (v0 - Vs) * exp(-t / tau);
}

/* RC charging current: i(t) = (Vs - v0)/R * exp(-t/tau)
 * Current decays exponentially from initial maximum. */
double rc_step_charge_i(double Vs, double v0, double R, double C, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0) return (Vs > v0) ? INFINITY : -INFINITY;
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    return ((Vs - v0) / R) * exp(-t / tau);
}

/* RC charging: voltage across resistor v_R(t) = (Vs - v0)*exp(-t/tau) */
double rc_step_charge_vr(double Vs, double v0, double R, double C, double t)
{
    if (t < 0.0) return Vs - v0;
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    return (Vs - v0) * exp(-t / tau);
}

/* RC charging: power dissipated in resistor P_R(t) = v_R^2 / R
 * Integrates to total energy: E_total = 0.5*C*(Vs-v0)^2 */
double rc_step_charge_power_R(double Vs, double v0, double R, double C, double t)
{
    double vr = rc_step_charge_vr(Vs, v0, R, C, t);
    if (R <= 0.0) return 0.0;
    return vr * vr / R;
}

/* RC charging: cumulated energy dissipated in resistor up to time t
 * E_R(t) = 0.5*C*(Vs-v0)^2 * (1 - exp(-2t/tau)) */
double rc_step_charge_energy_R(double Vs, double v0, double R, double C, double t)
{
    if (t < 0.0) return 0.0;
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    double E_total = 0.5 * C * (Vs - v0) * (Vs - v0);
    return E_total * (1.0 - exp(-2.0 * t / tau));
}

/* RC charging: energy stored in capacitor E_c(t) = 0.5*C*v_c(t)^2 */
double rc_step_charge_energy_C(double Vs, double v0, double R, double C, double t)
{
    double vc = rc_step_charge_vc(Vs, v0, R, C, t);
    return 0.5 * C * vc * vc;
}

/*----- L4: RC Natural Response (Discharging) -----*/

/* RC natural/discharging: v_c(t) = v0 * exp(-t/tau)
 * No external source; energy stored in capacitor dissipates through R. */
double rc_natural_vc(double v0, double R, double C, double t)
{
    if (t < 0.0) return v0;
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    return v0 * exp(-t / tau);
}

/* RC discharging current: i(t) = -(v0/R) * exp(-t/tau)
 * Current direction reverses from charging case. */
double rc_natural_i(double v0, double R, double C, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0) return 0.0;
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    return -(v0 / R) * exp(-t / tau);
}

/* RC discharging: v_R(t) = -v0 * exp(-t/tau) */
double rc_natural_vr(double v0, double R, double C, double t)
{
    return -rc_natural_vc(v0, R, C, t);
}

/* RC discharging: instantaneous power in resistor */
double rc_natural_power_R(double v0, double R, double C, double t)
{
    double vr = rc_natural_vr(v0, R, C, t);
    if (R <= 0.0) return 0.0;
    return vr * vr / R;
}

/* RC discharging: total energy dissipated up to time t
 * E_R(t) = 0.5*C*v0^2 * (1 - exp(-2t/tau)) */
double rc_natural_energy_dissipated(double v0, double R, double C, double t)
{
    if (t < 0.0) return 0.0;
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    double E_total = 0.5 * C * v0 * v0;
    return E_total * (1.0 - exp(-2.0 * t / tau));
}

/*----- L4: RC Complete Response (Step + Natural) -----*/

/* RC complete response to step: v_c(t) as above
 * Equivalent to rc_step_charge_vc */
double rc_complete_vc(double Vs, double v0, double R, double C, double t)
{
    return rc_step_charge_vc(Vs, v0, R, C, t);
}

/* RC complete response current */
double rc_complete_i(double Vs, double v0, double R, double C, double t)
{
    return rc_step_charge_i(Vs, v0, R, C, t);
}

/*----- L4: RL Step Response -----*/

/* RL step response: i_L(t) = Vs/R + (i0 - Vs/R)*exp(-R*t/L)
 * Inductor current cannot change instantaneously.
 * Final value: i_L(inf) = Vs/R (inductor acts as short at DC) */
double rl_step_charge_il(double Vs, double i0, double R, double L, double t)
{
    if (t < 0.0) return i0;
    if (R <= 0.0 || L <= 0.0) return Vs / (R + 1e-12);
    double tau = L / R;
    double i_final = Vs / R;
    return i_final + (i0 - i_final) * exp(-t / tau);
}

/* RL step: inductor voltage v_L(t) = Vs * exp(-R*t/L)
 * At t=0+, v_L jumps to Vs (inductor opposes current change). */
double rl_step_charge_vl(double Vs, double i0, double R, double L, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0 || L <= 0.0) return (t > 0.0) ? 0.0 : Vs;
    double tau = L / R;
    return (Vs - i0 * R) * exp(-t / tau);
}

/* RL step: voltage across resistor v_R(t) = Vs * (1 - exp(-t/tau)) */
double rl_step_charge_vr(double Vs, double i0, double R, double L, double t)
{
    double il = rl_step_charge_il(Vs, i0, R, L, t);
    return il * R;
}

/* RL step: energy stored in inductor E_L(t) = 0.5*L*i_L(t)^2 */
double rl_step_charge_energy_L(double Vs, double i0, double R, double L, double t)
{
    double il = rl_step_charge_il(Vs, i0, R, L, t);
    return 0.5 * L * il * il;
}

/* RL natural response (source removed): i_L(t) = i0 * exp(-R*t/L) */
double rl_natural_il(double i0, double R, double L, double t)
{
    if (t < 0.0) return i0;
    if (R <= 0.0 || L <= 0.0) return i0;
    return i0 * exp(-R * t / L);
}

/* RL natural: v_L(t) = -i0*R * exp(-R*t/L)
 * Negative voltage as magnetic field collapses. */
double rl_natural_vl(double i0, double R, double L, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0 || L <= 0.0) return 0.0;
    return -i0 * R * exp(-R * t / L);
}

/* RL complete response */
double rl_complete_il(double Vs, double i0, double R, double L, double t)
{
    return rl_step_charge_il(Vs, i0, R, L, t);
}

double rl_complete_vl(double Vs, double i0, double R, double L, double t)
{
    return rl_step_charge_vl(Vs, i0, R, L, t);
}

/*----- L4: Impulse Response -----*/

/* RC impulse response: v_c(t) = (1/(RC)) * exp(-t/RC) * u(t)
 * The impulse deposits charge instantaneously. Response is the
 * derivative of step response. */
double rc_impulse_vc(double R, double C, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0 || C <= 0.0) return 0.0;
    double tau = R * C;
    return exp(-t / tau) / tau;
}

/* RC impulse current response */
double rc_impulse_i(double R, double C, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0 || C <= 0.0) return 0.0;
    double tau = R * C;
    return -exp(-t / tau) / (R * tau);
}

/* RL impulse response: i_L(t) = (R/L) * exp(-R*t/L) * u(t) */
double rl_impulse_il(double R, double L, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0 || L <= 0.0) return 0.0;
    return (R / L) * exp(-R * t / L);
}

/* RL impulse: v_L(t) */
double rl_impulse_vl(double R, double L, double t)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0 || L <= 0.0) return 0.0;
    return R * rl_impulse_il(R, L, t);
}

/*----- L4: Ramp Response -----*/

/* RC ramp response to v_s(t) = slope * t * u(t)
 * v_c(t) = slope*(t - tau) + (v0 + slope*tau)*exp(-t/tau)
 * After transients decay, output follows ramp with delay tau. */
double rc_ramp_vc(double slope, double R, double C, double t, double v0)
{
    if (t < 0.0) return v0;
    double tau = R * C;
    if (tau <= 0.0) return slope * t;
    double vc_ss = slope * (t - tau);
    return vc_ss + (v0 - vc_ss + slope * tau) * exp(-t / tau);
}

/* RC ramp current */
double rc_ramp_i(double slope, double R, double C, double t, double v0)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0) return 0.0;
    double tau = R * C;
    if (tau <= 0.0) return slope * C;
    return slope * C * (1.0 - exp(-t / tau)) + (v0 / R) * exp(-t / tau) - (slope * tau / R) * exp(-t / tau);
}

/* RL ramp response: i_L for ramp voltage input */
double rl_ramp_il(double slope, double R, double L, double t, double i0)
{
    if (t < 0.0) return i0;
    if (R <= 0.0) return i0 + slope * t * t / (2.0 * L);
    double tau = L / R;
    double i_ss = (slope / R) * (t - tau);
    return i_ss + (i0 - i_ss + slope * tau / R) * exp(-t / tau);
}

/* RL ramp: v_L */
double rl_ramp_vl(double slope, double R, double L, double t, double i0)
{
    double il = rl_ramp_il(slope, R, L, t, i0);
    return slope * t - il * R;
}

/*----- L6: RC Pulse Response -----*/

/* RC response to a rectangular pulse of amplitude V_amplitude and width t_pulse
 * After pulse ends at t = t_pulse, circuit discharges naturally. */
double rc_pulse_vc(double V_amplitude, double R, double C, double t, double t_pulse)
{
    if (t < 0.0) return 0.0;
    double tau = R * C;
    if (tau <= 0.0) return (t <= t_pulse) ? V_amplitude : 0.0;
    if (t <= t_pulse) {
        return V_amplitude * (1.0 - exp(-t / tau));
    } else {
        double v_at_pulse = V_amplitude * (1.0 - exp(-t_pulse / tau));
        return v_at_pulse * exp(-(t - t_pulse) / tau);
    }
}

/* RC response to a periodic pulse train
 * Steady-state ripple develops after several periods. */
double rc_pulse_train_vc(double V_amplitude, double R, double C, double t,
                          double t_period, double duty_cycle, int n_periods)
{
    if (duty_cycle < 0.0 || duty_cycle > 1.0) return 0.0;
    double t_pulse = t_period * duty_cycle;
    double tau = R * C;
    if (tau <= 0.0) {
        double frac = fmod(t, t_period);
        return (frac <= t_pulse) ? V_amplitude : 0.0;
    }

    /* Compute steady-state ripple bounds */
    double a = exp(-t_pulse / tau);
    double b = exp(-(t_period - t_pulse) / tau);
    double v_high = V_amplitude * (1.0 - a) / (1.0 - a * b);
    double v_low = v_high * b;

    double frac = fmod(t, t_period);
    double n_cycle = floor(t / t_period);
    double v_initial = 0.0;

    /* For initial cycles before steady-state, compute recursively */
    if (n_cycle < (double)n_periods && n_cycle < 10) {
        for (int i = 0; i < (int)n_cycle && i < 10; i++) {
            v_initial = V_amplitude + (v_initial - V_amplitude) * a;
            v_initial = v_initial * b;
        }
    } else {
        /* Close to steady-state */
        v_initial = (frac <= t_pulse) ? v_low : v_high;
    }

    if (frac <= t_pulse) {
        return V_amplitude + (v_initial - V_amplitude) * exp(-frac / tau);
    } else {
        double v_at_end_pulse = V_amplitude + (v_initial - V_amplitude) * a;
        return v_at_end_pulse * exp(-(frac - t_pulse) / tau);
    }
}

/* RL pulse response: i_L for voltage pulse */
double rl_pulse_il(double V_amplitude, double R, double L, double t, double t_pulse)
{
    if (t < 0.0) return 0.0;
    if (R <= 0.0) return V_amplitude * t / L;
    double tau = L / R;
    double i_final = V_amplitude / R;
    if (t <= t_pulse) {
        return i_final * (1.0 - exp(-t / tau));
    } else {
        double i_at_pulse = i_final * (1.0 - exp(-t_pulse / tau));
        return i_at_pulse * exp(-(t - t_pulse) / tau);
    }
}

/*----- L4: Sinusoidal Response (Complete) -----*/

/* RC complete sinusoidal response.
 * v_s(t) = Vm * cos(omega*t + phi)
 * Complete = natural + forced sinusoidal steady-state.
 * v_c(t) = A*cos(omega*t + phi + theta) + B*exp(-t/tau) */
double rc_sinusoidal_vc(double Vm, double omega, double phi, double R, double C, double t, double v0)
{
    if (t < 0.0) return v0;
    double tau = R * C;
    if (tau <= 0.0 || omega <= 0.0) return v0;

    /* Forced sinusoidal steady-state using phasor analysis */
    double complex Zc = 1.0 / (I * omega * C);
    double complex Z = R + Zc;
    double complex Vc_phasor = Vm * cexp(I * phi * M_PI / 180.0) * Zc / Z;

    double vc_forced = cabs(Vc_phasor) * cos(omega * t + carg(Vc_phasor));
    double vc_forced_0 = cabs(Vc_phasor) * cos(carg(Vc_phasor));

    /* Natural response satisfies initial condition */
    double vc_natural = (v0 - vc_forced_0) * exp(-t / tau);

    return vc_forced + vc_natural;
}

/* RC sinusoidal current response */
double rc_sinusoidal_i(double Vm, double omega, double phi, double R, double C, double t, double v0)
{
    if (t < 0.0) return 0.0;
    double vc = rc_sinusoidal_vc(Vm, omega, phi, R, C, t, v0);
    double vs = Vm * cos(omega * t + phi * M_PI / 180.0);
    if (R <= 0.0) return 0.0;
    return (vs - vc) / R;
}

/* RL complete sinusoidal response */
double rl_sinusoidal_il(double Vm, double omega, double phi, double R, double L, double t, double i0)
{
    if (t < 0.0) return i0;
    if (R <= 0.0 || omega <= 0.0) return i0;

    double tau = L / R;
    double complex Zl = I * omega * L;
    double complex Z = R + Zl;
    double complex Il_phasor = Vm * cexp(I * phi * M_PI / 180.0) / Z;

    double il_forced = cabs(Il_phasor) * cos(omega * t + carg(Il_phasor));
    double il_forced_0 = cabs(Il_phasor) * cos(carg(Il_phasor));

    double il_natural = (i0 - il_forced_0) * exp(-t / tau);

    return il_forced + il_natural;
}

/*----- L4: Switching at Non-Zero Time -----*/

/* RC response with switch at t_switch */
double rc_switched_vc(double Vs, double v0, double R, double C, double t, double t_switch)
{
    if (t < t_switch) {
        return rc_step_charge_vc(Vs, v0, R, C, t);
    } else {
        double v_at_switch = rc_step_charge_vc(Vs, v0, R, C, t_switch);
        return rc_natural_vc(v_at_switch, R, C, t - t_switch);
    }
}

/* RL response with switch at t_switch */
double rl_switched_il(double Vs, double i0, double R, double L, double t, double t_switch)
{
    if (t < t_switch) {
        return rl_step_charge_il(Vs, i0, R, L, t);
    } else {
        double i_at_switch = rl_step_charge_il(Vs, i0, R, L, t_switch);
        return rl_natural_il(i_at_switch, R, L, t - t_switch);
    }
}

/*----- L4: General First-Order ODE Solution -----*/

/* General solution to dy/dt + a*y = b, y(0) = y0
 * y(t) = b/a + (y0 - b/a)*exp(-a*t) for a > 0 */
double first_order_general_response(double a, double b, double y0, double t)
{
    if (t < 0.0) return y0;
    if (a <= 0.0) return y0 + b * t;
    double y_ss = b / a;
    return y_ss + (y0 - y_ss) * exp(-a * t);
}

/* Steady-state error for unit step input: e_ss = 1/(1+K) for unity feedback */
double first_order_step_ss_error(double tau, double K)
{
    (void)tau;
    return 1.0 / (1.0 + K);
}

/* Bandwidth from time constant: BW (rad/s) = 1/tau
 * Reference: First-order lowpass -3dB frequency */
double first_order_bandwidth_from_tau(double tau)
{
    if (tau <= 0.0) return INFINITY;
    return 1.0 / tau;
}

/* Time constant from 10-90% rise time: tau = t_r / 2.2 */
double tau_from_rise_time(double tr_10_90)
{
    if (tr_10_90 <= 0.0) return 0.0;
    return tr_10_90 / 2.2;
}

/*----- L4: Sequential Switching -----*/

/* Two-step sequential switching: first R1 to V1, then switch to R2, V2 */
double rc_sequential_vc(double V1, double V2, double R1, double R2,
                         double C, double t, double t_switch)
{
    if (t < 0.0) return 0.0;

    if (t <= t_switch) {
        double tau1 = R1 * C;
        if (tau1 <= 0.0) return V1;
        return V1 * (1.0 - exp(-t / tau1));
    } else {
        double tau1 = R1 * C;
        double v_at_switch = (tau1 > 0.0) ? V1 * (1.0 - exp(-t_switch / tau1)) : V1;
        double tau2 = R2 * C;
        if (tau2 <= 0.0) return V2;
        return V2 + (v_at_switch - V2) * exp(-(t - t_switch) / tau2);
    }
}

/* Charge then discharge: apply Vs through R for t_hold, then disconnect */
double rc_sequential_charge_then_discharge_vc(double Vs, double R, double C,
                                               double t, double t_hold)
{
    return rc_switched_vc(Vs, 0.0, R, C, t, t_hold);
}

/*----- L4: Cascaded RC (Multiple Time Constants) -----*/

/* Two cascaded RC stages (non-interacting, buffered).
 * v2 is voltage on second capacitor.
 * Exact solution has two exponential terms. */
double rc_cascaded_vc(double Vs, double R1, double C1, double R2, double C2,
                      double t, double v10, double v20)
{
    if (t < 0.0) return v20;
    double tau1 = R1 * C1;
    double tau2 = R2 * C2;

    if (tau1 <= 0.0 && tau2 <= 0.0) return Vs;
    if (fabs(tau1 - tau2) < 1e-12) {
        /* Equal time constants: v2 = Vs * (1 - (1 + t/tau)*exp(-t/tau)) */
        return Vs * (1.0 - (1.0 + t / tau1) * exp(-t / tau1));
    }

    /* Distinct time constants */
    double A = tau1 / (tau1 - tau2);
    double B = tau2 / (tau2 - tau1);
    return Vs * (1.0 + A * exp(-t / tau1) + B * exp(-t / tau2));
}

/* Dominant time constant: the larger one dominates the response */
double dominant_time_constant(double tau1, double tau2)
{
    return (tau1 > tau2) ? tau1 : tau2;
}

/* Effective time constant for cascaded first-order systems.
 * Approximate: tau_eff ≈ tau1 + tau2 (Elmore delay) */
double effective_time_constant_cascaded(double tau1, double tau2)
{
    return tau1 + tau2;
}

/*----- L7: RC Timing Applications -----*/

/* Time to reach V_target from v_start charging toward v_final.
 * t = -tau * ln((v_final - V_target)/(v_final - v_start)) */
double rc_delay_time(double R, double C, double V_target, double V_supply)
{
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    if (V_target >= V_supply) return INFINITY;
    if (V_target <= 0.0) return 0.0;
    return -tau * log(1.0 - V_target / V_supply);
}

/* General threshold crossing time */
double rc_threshold_crossing_time(double R, double C, double v_start,
                                   double v_target, double v_final)
{
    double tau = R * C;
    if (tau <= 0.0) return 0.0;
    double numerator = v_final - v_target;
    double denominator = v_final - v_start;
    if (numerator <= 0.0 || denominator <= 0.0) return INFINITY;
    return -tau * log(numerator / denominator);
}

/* Start-up time for power supply rail */
double rc_startup_time(double R, double C, double V_threshold, double V_supply)
{
    return rc_delay_time(R, C, V_threshold, V_supply);
}

/* Power-On Reset timeout: t = R*C * ln(Vcc/(Vcc - V_threshold)) */
double por_timeout_compute(double R, double C, double Vcc, double V_threshold)
{
    if (R <= 0.0 || C <= 0.0) return 0.0;
    if (Vcc <= V_threshold || V_threshold <= 0.0) return 0.0;
    return R * C * log(Vcc / (Vcc - V_threshold));
}

/* Watchdog timeout: t = -R*C * ln(1 - V_threshold/Vcc) */
double watchdog_timeout_compute(double R, double C, double V_threshold, double Vcc)
{
    return por_timeout_compute(R, C, Vcc, V_threshold);
}

/* Debounce settling time for RC filter */
double debounce_settling_time_compute(double R, double C, double V_cc, double V_IH_min)
{
    if (R <= 0.0 || C <= 0.0) return 0.0;
    if (V_IH_min >= V_cc) return INFINITY;
    return -R * C * log(1.0 - V_IH_min / V_cc);
}

/*----- L7: RC Oscillator -----*/

/* RC relaxation oscillator period
 * Period = R*C * ln((V_supply - V_low)/(V_supply - V_high)) * 2 */
double rc_oscillator_period(double R, double C, double V_high, double V_low, double V_supply)
{
    if (R <= 0.0 || C <= 0.0) return INFINITY;
    if (V_high >= V_supply || V_low <= 0.0 || V_high <= V_low) return INFINITY;
    double charge_time = -R * C * log((V_supply - V_high) / (V_supply - V_low));
    double discharge_time = -R * C * log(V_low / V_high);
    return charge_time + discharge_time;
}

/* RC relaxation oscillator frequency */
double rc_oscillator_frequency(double R, double C, double V_high, double V_low, double V_supply)
{
    double period = rc_oscillator_period(R, C, V_high, V_low, V_supply);
    if (period <= 0.0 || isinf(period)) return 0.0;
    return 1.0 / period;
}

/* Symmetric RC oscillator using Schmitt trigger with ratio V_th/V_supply */
double rc_relaxation_oscillator_freq(double R, double C, double V_th_ratio)
{
    if (R <= 0.0 || C <= 0.0) return 0.0;
    if (V_th_ratio <= 0.0 || V_th_ratio >= 0.5) return 0.0;
    /* For symmetric oscillator: f = 1 / (2*R*C*ln((1+beta)/(1-beta))) */
    double beta = 1.0 - 2.0 * V_th_ratio;
    if (beta <= 0.0) return 0.0;
    return 1.0 / (2.0 * R * C * log((1.0 + beta) / (1.0 - beta)));
}

/*----- L4: Transient Response Bounds -----*/

/* Maximum instantaneous current in RC circuit
 * Occurs at t=0+: i_max = (Vs - v0)/R */
double rc_max_current_instantaneous(double Vs, double v0, double R)
{
    if (R <= 0.0) return INFINITY;
    return fabs(Vs - v0) / R;
}

/* Maximum instantaneous power in R */
double rc_max_power_instantaneous(double Vs, double v0, double R)
{
    double i_max = rc_max_current_instantaneous(Vs, v0, R);
    return i_max * i_max * R;
}

/* Maximum instantaneous voltage across inductor
 * At switching: v_L(0+) = Vs - i0*R */
double rl_max_voltage_instantaneous(double Vs, double i0, double R)
{
    return fabs(Vs - i0 * R);
}

/*----- L7: Exponential Curve Fitting -----*/

/* Estimate time constant from two measured points during decay.
 * tau = -dt / ln(v2/v1) */
double estimate_tau_from_decay(double v1, double v2, double dt)
{
    if (v1 <= 0.0 || v2 <= 0.0 || dt <= 0.0) return -1.0;
    if (v2 >= v1) return -1.0;
    return -dt / log(v2 / v1);
}

/* Predict voltage at time t_target given a measurement at t_measured */
double predict_voltage_at_time(double v_measured, double t_measured,
                                double tau, double t_target)
{
    if (tau <= 0.0) return v_measured;
    return v_measured * exp(-(t_target - t_measured) / tau);
}
