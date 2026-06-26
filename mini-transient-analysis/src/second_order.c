#include "second_order.h"
#include "transient_defs.h"
#include "first_order.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*==================================================================
 * L2-L4: Second-Order Transient Analysis Implementation
 *
 * The second-order ODE: d^2y/dt^2 + 2*zeta*omega_n*dy/dt + omega_n^2*y = f(t)
 *
 * Characteristic equation: s^2 + 2*zeta*omega_n*s + omega_n^2 = 0
 * Roots: s = -zeta*omega_n +/- omega_n*sqrt(zeta^2 - 1)
 *
 * Three regimes based on damping ratio zeta:
 *   Overdamped (zeta > 1): two distinct real roots
 *   Critically damped (zeta = 1): repeated real root
 *   Underdamped (0 < zeta < 1): complex conjugate roots
 *   Undamped (zeta = 0): purely imaginary roots
 *
 * Reference: Hayt & Kemmerly Ch.9, Ogata "System Dynamics" Ch.5
 *==================================================================*/

/*----- L2: RLC Series Parameters -----*/

/* Neper frequency (attenuation): alpha = R/(2L)
 * Controls the exponential decay envelope of the response. */
double rlc_series_alpha(double R, double L)
{
    if (L <= 0.0) return INFINITY;
    return R / (2.0 * L);
}

/* Undamped natural frequency: omega_n = 1/sqrt(LC)
 * The frequency at which the circuit would oscillate without damping. */
double rlc_series_omega_n(double L, double C)
{
    if (L <= 0.0 || C <= 0.0) return 0.0;
    return 1.0 / sqrt(L * C);
}

/* Damping ratio: zeta = alpha/omega_n = (R/2)*sqrt(C/L)
 * zeta < 1: underdamped (oscillatory)
 * zeta = 1: critically damped (fastest non-oscillatory)
 * zeta > 1: overdamped (sluggish, no overshoot) */
double rlc_series_zeta(double R, double L, double C)
{
    double alpha = rlc_series_alpha(R, L);
    double omega_n = rlc_series_omega_n(L, C);
    if (omega_n <= 0.0) return INFINITY;
    return alpha / omega_n;
}

/* Damped natural frequency: omega_d = omega_n * sqrt(1 - zeta^2)
 * Only real for zeta < 1. This is the observed oscillation frequency. */
double rlc_series_omega_d(double R, double L, double C)
{
    double zeta = rlc_series_zeta(R, L, C);
    double omega_n = rlc_series_omega_n(L, C);
    if (omega_n <= 0.0 || zeta >= 1.0) return 0.0;
    return omega_n * sqrt(1.0 - zeta * zeta);
}

DampingClass_t rlc_series_damping_class(double R, double L, double C)
{
    double zeta = rlc_series_zeta(R, L, C);
    return classify_damping(zeta);
}

/*----- L2: RLC Parallel Parameters -----*/

double rlc_parallel_alpha(double R, double C)
{
    if (R <= 0.0 || C <= 0.0) return INFINITY;
    return 1.0 / (2.0 * R * C);
}

double rlc_parallel_omega_n(double L, double C)
{
    if (L <= 0.0 || C <= 0.0) return 0.0;
    return 1.0 / sqrt(L * C);
}

double rlc_parallel_zeta(double R, double L, double C)
{
    double alpha = rlc_parallel_alpha(R, C);
    double omega_n = rlc_parallel_omega_n(L, C);
    if (omega_n <= 0.0) return INFINITY;
    return alpha / omega_n;
}

double rlc_parallel_omega_d(double R, double L, double C)
{
    double zeta = rlc_parallel_zeta(R, L, C);
    double omega_n = rlc_parallel_omega_n(L, C);
    if (omega_n <= 0.0 || zeta >= 1.0) return 0.0;
    return omega_n * sqrt(1.0 - zeta * zeta);
}

DampingClass_t rlc_parallel_damping_class(double R, double L, double C)
{
    double zeta = rlc_parallel_zeta(R, L, C);
    return classify_damping(zeta);
}

/*----- L4: RLC Series Step Response -----*/

/* RLC series step response for capacitor voltage v_c(t).
 * Input: DC voltage step Vs*u(t).
 * Initial conditions: v_c(0) = v0, i_L(0) = i0.
 *
 * General form: v_c(t) = v_c(inf) + transient terms
 * where v_c(inf) = Vs (capacitor is open at DC).
 *
 * The transient form depends on zeta (see sub-cases below). */
double rlc_series_step_vc(double Vs, double R, double L, double C,
                          double v0, double i0, double t)
{
    if (t < 0.0) return v0;
    if (L <= 0.0 || C <= 0.0) return Vs;

    double alpha = rlc_series_alpha(R, L);
    double omega_n = rlc_series_omega_n(L, C);
    double zeta = rlc_series_zeta(R, L, C);

    if (omega_n <= 0.0) return Vs;

    /* Steady-state: capacitor is open, v_c = Vs */
    double vc_ss = Vs;

    if (zeta > 1.0 + 1e-10) {
        /* Overdamped: v_c(t) = Vss + A1*exp(s1*t) + A2*exp(s2*t)
         * s1,2 = -alpha +/- sqrt(alpha^2 - omega_n^2) */
        double disc = sqrt(alpha * alpha - omega_n * omega_n);
        double s1 = -alpha + disc;
        double s2 = -alpha - disc;

        /* Initial conditions: v_c(0) = v0, dv_c/dt(0) = i0/C */
        double A2 = (s1 * (v0 - vc_ss) - i0 / C) / (s1 - s2);
        double A1 = v0 - vc_ss - A2;

        return vc_ss + A1 * exp(s1 * t) + A2 * exp(s2 * t);

    } else if (fabs(zeta - 1.0) < 1e-10) {
        /* Critically damped: v_c = Vss + (A1 + A2*t)*exp(-alpha*t) */
        double A1 = v0 - vc_ss;
        double A2 = i0 / C + alpha * (v0 - vc_ss);

        return vc_ss + (A1 + A2 * t) * exp(-alpha * t);

    } else if (zeta > 0.0) {
        /* Underdamped: v_c = Vss + exp(-alpha*t)*(A1*cos(wd*t) + A2*sin(wd*t)) */
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double A1 = v0 - vc_ss;
        double A2 = (i0 / C + alpha * (v0 - vc_ss)) / wd;

        return vc_ss + exp(-alpha * t) * (A1 * cos(wd * t) + A2 * sin(wd * t));

    } else {
        /* Undamped: v_c = Vss + A1*cos(omega_n*t) + A2*sin(omega_n*t) */
        double A1 = v0 - vc_ss;
        double A2 = i0 / (C * omega_n);

        return vc_ss + A1 * cos(omega_n * t) + A2 * sin(omega_n * t);
    }
}

/* RLC series step response for inductor current i_L(t).
 * At steady-state, i_L = 0 (capacitor blocks DC). */
double rlc_series_step_il(double Vs, double R, double L, double C,
                          double v0, double i0, double t)
{
    if (t < 0.0) return i0;
    if (L <= 0.0 || C <= 0.0) return 0.0;

    double alpha = rlc_series_alpha(R, L);
    double omega_n = rlc_series_omega_n(L, C);
    double zeta = rlc_series_zeta(R, L, C);

    /* i_L = C * dv_c/dt */
    double vc_ss = Vs;

    if (zeta > 1.0 + 1e-10) {
        double disc = sqrt(alpha * alpha - omega_n * omega_n);
        double s1 = -alpha + disc;
        double s2 = -alpha - disc;

        double A2 = ((v0 - vc_ss) / L + s1 * i0) / (s2 - s1);
        double A1 = i0 - A2;

        return A1 * exp(s1 * t) + A2 * exp(s2 * t);

    } else if (fabs(zeta - 1.0) < 1e-10) {
        double A1 = i0;
        double A2 = (v0 - vc_ss) / L + alpha * i0;
        return (A1 + A2 * t) * exp(-alpha * t);

    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double A1 = i0;
        double A2 = ((v0 - vc_ss) / L + alpha * i0) / wd;
        return exp(-alpha * t) * (A1 * cos(wd * t) + A2 * sin(wd * t));

    } else {
        double A1 = i0;
        double A2 = (v0 - vc_ss) / (L * omega_n);
        return A1 * cos(omega_n * t) + A2 * sin(omega_n * t);
    }
}

/* RLC series step: inductor voltage v_L = L * di_L/dt */
double rlc_series_step_vl(double Vs, double R, double L, double C,
                          double v0, double i0, double t)
{
    if (t < 0.0) return 0.0;
    if (L <= 0.0) return 0.0;

    double alpha = rlc_series_alpha(R, L);
    double omega_n = rlc_series_omega_n(L, C);
    double zeta = rlc_series_zeta(R, L, C);
    double vc_ss = Vs;

    if (zeta > 1.0 + 1e-10) {
        double disc = sqrt(alpha * alpha - omega_n * omega_n);
        double s1 = -alpha + disc;
        double s2 = -alpha - disc;
        double A2 = ((v0 - vc_ss) / L + s1 * i0) / (s2 - s1);
        double A1 = i0 - A2;
        return L * (A1 * s1 * exp(s1 * t) + A2 * s2 * exp(s2 * t));
    } else if (fabs(zeta - 1.0) < 1e-10) {
        double A1 = i0;
        double A2 = (v0 - vc_ss) / L + alpha * i0;
        return L * ((A2 - alpha * A1 - alpha * A2 * t) * exp(-alpha * t));
    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double A1 = i0;
        double A2 = ((v0 - vc_ss) / L + alpha * i0) / wd;
        double dil_dt = exp(-alpha * t) * ((-alpha * A1 + A2 * wd) * cos(wd * t) + (-alpha * A2 - A1 * wd) * sin(wd * t));
        return L * dil_dt;
    } else {
        double A1 = i0;
        double A2 = (v0 - vc_ss) / (L * omega_n);
        return L * omega_n * (-A1 * sin(omega_n * t) + A2 * cos(omega_n * t));
    }
}

/* RLC series step: resistor voltage v_R = i * R */
double rlc_series_step_vr(double Vs, double R, double L, double C,
                          double v0, double i0, double t)
{
    double il = rlc_series_step_il(Vs, R, L, C, v0, i0, t);
    return il * R;
}

/*----- L4: RLC Series Sub-type Functions (explicit forms) -----*/

double rlc_series_overdamped_vc(double Vs, double v0, double i0,
                                double s1, double s2, double C, double t)
{
    if (t < 0.0) return v0;
    double A2 = (s1 * (v0 - Vs) - i0 / C) / (s1 - s2);
    double A1 = v0 - Vs - A2;
    return Vs + A1 * exp(s1 * t) + A2 * exp(s2 * t);
}

double rlc_series_overdamped_il(double Vs, double v0, double i0,
                                double s1, double s2, double L, double t)
{
    if (t < 0.0) return i0;
    double A2 = ((v0 - Vs) / L + s1 * i0) / (s2 - s1);
    double A1 = i0 - A2;
    return A1 * exp(s1 * t) + A2 * exp(s2 * t);
}

double rlc_series_critical_vc(double Vs, double v0, double i0,
                              double alpha, double C, double t)
{
    if (t < 0.0) return v0;
    double A1 = v0 - Vs;
    double A2 = i0 / C + alpha * (v0 - Vs);
    return Vs + (A1 + A2 * t) * exp(-alpha * t);
}

double rlc_series_critical_il(double Vs, double v0, double i0,
                              double alpha, double L, double t)
{
    if (t < 0.0) return i0;
    double A1 = i0;
    double A2 = (v0 - Vs) / L + alpha * i0;
    return (A1 + A2 * t) * exp(-alpha * t);
}

double rlc_series_underdamped_vc(double Vs, double v0, double i0,
                                 double alpha, double omega_d, double C, double t)
{
    if (t < 0.0) return v0;
    double A1 = v0 - Vs;
    double A2 = (i0 / C + alpha * (v0 - Vs)) / omega_d;
    return Vs + exp(-alpha * t) * (A1 * cos(omega_d * t) + A2 * sin(omega_d * t));
}

double rlc_series_underdamped_il(double Vs, double v0, double i0,
                                 double alpha, double omega_d, double L, double t)
{
    if (t < 0.0) return i0;
    double A1 = i0;
    double A2 = ((v0 - Vs) / L + alpha * i0) / omega_d;
    return exp(-alpha * t) * (A1 * cos(omega_d * t) + A2 * sin(omega_d * t));
}

double rlc_series_undamped_vc(double Vs, double v0, double i0,
                              double omega_n, double C, double t)
{
    if (t < 0.0) return v0;
    double A1 = v0 - Vs;
    double A2 = i0 / (C * omega_n);
    return Vs + A1 * cos(omega_n * t) + A2 * sin(omega_n * t);
}

double rlc_series_undamped_il(double Vs, double v0, double i0,
                              double omega_n, double L, double t)
{
    if (t < 0.0) return i0;
    (void)Vs; /* undamped: no steady-state current through C at DC */
    double A1 = i0;
    double A2 = (v0 - Vs) / (L * omega_n);
    return A1 * cos(omega_n * t) + A2 * sin(omega_n * t);
}

/*----- L4: RLC Series Natural Response -----*/
double rlc_series_natural_vc(double R, double L, double C,
                             double v0, double i0, double t)
{
    return rlc_series_step_vc(0.0, R, L, C, v0, i0, t);
}

double rlc_series_natural_il(double R, double L, double C,
                             double v0, double i0, double t)
{
    return rlc_series_step_il(0.0, R, L, C, v0, i0, t);
}

/*----- L4: RLC Parallel Step Response -----*/

double rlc_parallel_step_v(double Is, double R, double L, double C,
                           double v0, double i0, double t)
{
    if (t < 0.0) return v0;
    if (L <= 0.0 || C <= 0.0) return 0.0;

    double alpha = rlc_parallel_alpha(R, C);
    double omega_n = rlc_parallel_omega_n(L, C);
    double zeta = rlc_parallel_zeta(R, L, C);

    if (omega_n <= 0.0) return 0.0;

    double v_ss = Is * R; /* DC steady-state: inductor short, capacitor open */

    if (zeta > 1.0 + 1e-10) {
        double disc = sqrt(alpha * alpha - omega_n * omega_n);
        double s1 = -alpha + disc;
        double s2 = -alpha - disc;
        double A2 = (s1 * (v0 - v_ss) - (Is - i0) / C) / (s1 - s2);
        double A1 = v0 - v_ss - A2;
        return v_ss + A1 * exp(s1 * t) + A2 * exp(s2 * t);
    } else if (fabs(zeta - 1.0) < 1e-10) {
        double A1 = v0 - v_ss;
        double A2 = (Is - i0) / C + alpha * (v0 - v_ss);
        return v_ss + (A1 + A2 * t) * exp(-alpha * t);
    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double A1 = v0 - v_ss;
        double A2 = ((Is - i0) / C + alpha * (v0 - v_ss)) / wd;
        return v_ss + exp(-alpha * t) * (A1 * cos(wd * t) + A2 * sin(wd * t));
    } else {
        double A1 = v0 - v_ss;
        double A2 = (Is - i0) / (C * omega_n);
        return v_ss + A1 * cos(omega_n * t) + A2 * sin(omega_n * t);
    }
}

double rlc_parallel_step_il(double Is, double R, double L, double C,
                            double v0, double i0, double t)
{
    if (t < 0.0) return i0;
    if (L <= 0.0 || C <= 0.0) return 0.0;

    double alpha = rlc_parallel_alpha(R, C);
    double omega_n = rlc_parallel_omega_n(L, C);
    double zeta = rlc_parallel_zeta(R, L, C);

    if (omega_n <= 0.0) return 0.0;

    double il_ss = Is;

    if (zeta > 1.0 + 1e-10) {
        double disc = sqrt(alpha * alpha - omega_n * omega_n);
        double s1 = -alpha + disc;
        double s2 = -alpha - disc;
        double A2 = (s1 * (i0 - il_ss) + v0 / L) / (s1 - s2);
        double A1 = i0 - il_ss - A2;
        return il_ss + A1 * exp(s1 * t) + A2 * exp(s2 * t);
    } else if (fabs(zeta - 1.0) < 1e-10) {
        double A1 = i0 - il_ss;
        double A2 = v0 / L + alpha * (i0 - il_ss);
        return il_ss + (A1 + A2 * t) * exp(-alpha * t);
    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double A1 = i0 - il_ss;
        double A2 = (v0 / L + alpha * (i0 - il_ss)) / wd;
        return il_ss + exp(-alpha * t) * (A1 * cos(wd * t) + A2 * sin(wd * t));
    } else {
        double A1 = i0 - il_ss;
        double A2 = v0 / (L * omega_n);
        return il_ss + A1 * cos(omega_n * t) + A2 * sin(omega_n * t);
    }
}

double rlc_parallel_step_ic(double Is, double R, double L, double C,
                            double v0, double i0, double t)
{
    double v = rlc_parallel_step_v(Is, R, L, C, v0, i0, t);
    double il = rlc_parallel_step_il(Is, R, L, C, v0, i0, t);
    return Is - v / R - il;
}

/*----- L4: RLC Parallel Natural Response -----*/

double rlc_parallel_natural_v(double R, double L, double C,
                              double v0, double i0, double t)
{
    return rlc_parallel_step_v(0.0, R, L, C, v0, i0, t);
}

double rlc_parallel_natural_il(double R, double L, double C,
                               double v0, double i0, double t)
{
    return rlc_parallel_step_il(0.0, R, L, C, v0, i0, t);
}

/*----- L4: RLC Impulse Response -----*/

double rlc_series_impulse_vc(double R, double L, double C, double t)
{
    /* Impulse response = derivative of step response (zero initial conditions)
     * Interpret as response to initial current impulse: i0 = 1/L, v0 = 0 */
    return rlc_series_step_il(0.0, R, L, C, 0.0, 1.0 / L, t);
}

double rlc_series_impulse_il(double R, double L, double C, double t)
{
    if (t < 0.0) return 0.0;
    double zeta = rlc_series_zeta(R, L, C);
    double omega_n = rlc_series_omega_n(L, C);
    double alpha = rlc_series_alpha(R, L);

    if (omega_n <= 0.0) return 0.0;

    if (zeta > 1.0) {
        double disc = sqrt(alpha * alpha - omega_n * omega_n);
        return (1.0 / (2.0 * disc)) * (exp((-alpha + disc) * t) - exp((-alpha - disc) * t));
    } else if (fabs(zeta - 1.0) < 1e-10) {
        return t * exp(-alpha * t);
    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        return exp(-alpha * t) * sin(wd * t) / wd;
    } else {
        return sin(omega_n * t) / omega_n;
    }
}

double rlc_parallel_impulse_v(double R, double L, double C, double t)
{
    if (t < 0.0) return 0.0;
    double zeta = rlc_parallel_zeta(R, L, C);
    double omega_n = rlc_parallel_omega_n(L, C);
    double alpha = rlc_parallel_alpha(R, C);

    if (omega_n <= 0.0) return 0.0;

    if (zeta > 1.0) {
        double disc = sqrt(alpha * alpha - omega_n * omega_n);
        return (1.0 / (2.0 * C * disc)) * (exp((-alpha + disc) * t) - exp((-alpha - disc) * t));
    } else if (fabs(zeta - 1.0) < 1e-10) {
        return (t / C) * exp(-alpha * t);
    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        return exp(-alpha * t) * sin(wd * t) / (C * wd);
    } else {
        return sin(omega_n * t) / (C * omega_n);
    }
}

/*----- L4: RLC Sinusoidal Response (Complete) -----*/

double rlc_series_sinusoidal_vc(double Vm, double omega, double phi,
                                double R, double L, double C,
                                double v0, double i0, double t)
{
    if (t < 0.0) return v0;
    if (L <= 0.0 || C <= 0.0) return v0;

    /* Forced response using phasor analysis */
    double complex Zl = I * omega * L;
    double complex Zc = 1.0 / (I * omega * C);
    double complex Z = R + Zl + Zc;
    double complex Vc_phasor = Vm * cexp(I * phi * M_PI / 180.0) * Zc / Z;

    double vc_forced = cabs(Vc_phasor) * cos(omega * t + carg(Vc_phasor));
    double vc_forced_0 = cabs(Vc_phasor) * cos(carg(Vc_phasor));

    /* Natural response: same form as natural response from initial conditions */
    double vc_natural = rlc_series_natural_vc(R, L, C, v0 - vc_forced_0, i0, t);

    return vc_forced + vc_natural;
}

double rlc_parallel_sinusoidal_v(double Im, double omega, double phi,
                                  double R, double L, double C,
                                  double v0, double i0, double t)
{
    if (t < 0.0) return v0;
    if (L <= 0.0 || C <= 0.0) return v0;

    double complex Yr = 1.0 / R;
    double complex Yl = 1.0 / (I * omega * L);
    double complex Yc = I * omega * C;
    double complex Y = Yr + Yl + Yc;
    double complex V_phasor = Im * cexp(I * phi * M_PI / 180.0) / Y;

    double v_forced = cabs(V_phasor) * cos(omega * t + carg(V_phasor));
    double v_forced_0 = cabs(V_phasor) * cos(carg(V_phasor));

    double v_natural = rlc_parallel_natural_v(R, L, C, v0 - v_forced_0, i0, t);

    return v_forced + v_natural;
}

/*----- L6: Second-Order Performance Metrics -----*/

double second_order_overshoot_pct(double zeta)
{
    return overshoot_from_zeta(zeta);
}

double second_order_peak_time(double omega_n, double zeta)
{
    return peak_time_from_params(omega_n, zeta);
}

double second_order_settling_time_2pct(double zeta, double omega_n)
{
    return settling_time_2pct(zeta, omega_n);
}

double second_order_settling_time_5pct(double zeta, double omega_n)
{
    return settling_time_5pct(zeta, omega_n);
}

double second_order_rise_time_10_90(double omega_n, double zeta)
{
    return rise_time_10_90_from_params(omega_n, zeta);
}

/* L6: Delay time: t_d ≈ (1 + 0.7*zeta) / omega_n */
double second_order_delay_time(double omega_n, double zeta)
{
    if (omega_n <= 0.0) return INFINITY;
    return (1.0 + 0.7 * zeta) / omega_n;
}

/* L6: Period of oscillation: T_d = 2*pi/omega_d */
double second_order_period_of_oscillation(double omega_d)
{
    if (omega_d <= 0.0) return INFINITY;
    return 2.0 * M_PI / omega_d;
}

/* L6: Decay ratio: ratio of successive peaks = exp(-2*pi*zeta/sqrt(1-zeta^2)) */
double second_order_decay_ratio(double zeta)
{
    if (zeta >= 1.0 || zeta <= 0.0) return 1.0;
    return exp(-2.0 * M_PI * zeta / sqrt(1.0 - zeta * zeta));
}

/* L6: Number of oscillations before settling within 2%
 * N ≈ ln(0.02) / (-2*pi*zeta/sqrt(1-zeta^2)) */
double second_order_num_oscillations_to_settle(double zeta)
{
    if (zeta >= 1.0 || zeta <= 0.0) return 0.0;
    double log_dec = -2.0 * M_PI * zeta / sqrt(1.0 - zeta * zeta);
    if (log_dec >= 0.0) return INFINITY;
    return log(0.02) / log_dec;
}

/*----- L4: General Second-Order Characteristic Roots -----*/

/* Compute roots of a*s^2 + b*s + c = 0
 * Standard form: s^2 + (b/a)*s + (c/a) = 0
 * Roots: s = (-b +/- sqrt(b^2 - 4ac)) / (2a) */
void second_order_characteristic_roots(double a, double b, double c,
                                        double *r1_real, double *r1_imag,
                                        double *r2_real, double *r2_imag)
{
    if (fabs(a) < 1e-15) {
        /* Degenerate to first order */
        if (r1_real) *r1_real = (fabs(b) > 1e-15) ? -c / b : 0.0;
        if (r1_imag) *r1_imag = 0.0;
        if (r2_real) *r2_real = 0.0;
        if (r2_imag) *r2_imag = 0.0;
        return;
    }

    double disc = b * b - 4.0 * a * c;

    if (disc > 0.0) {
        double sqrt_disc = sqrt(disc);
        if (r1_real) *r1_real = (-b + sqrt_disc) / (2.0 * a);
        if (r1_imag) *r1_imag = 0.0;
        if (r2_real) *r2_real = (-b - sqrt_disc) / (2.0 * a);
        if (r2_imag) *r2_imag = 0.0;
    } else if (fabs(disc) < 1e-15) {
        double r = -b / (2.0 * a);
        if (r1_real) *r1_real = r;
        if (r1_imag) *r1_imag = 0.0;
        if (r2_real) *r2_real = r;
        if (r2_imag) *r2_imag = 0.0;
    } else {
        double real_part = -b / (2.0 * a);
        double imag_part = sqrt(-disc) / (2.0 * a);
        if (r1_real) *r1_real = real_part;
        if (r1_imag) *r1_imag = imag_part;
        if (r2_real) *r2_real = real_part;
        if (r2_imag) *r2_imag = -imag_part;
    }
}

double second_order_discriminant(double a, double b, double c)
{
    return b * b - 4.0 * a * c;
}

DampingClass_t second_order_classify_by_roots(double discriminant)
{
    if (discriminant > 1e-10) return DAMPING_OVERDAMPED;
    if (fabs(discriminant) < 1e-10) return DAMPING_CRITICALLY_DAMPED;
    return DAMPING_UNDERDAMPED;
}

/*----- L4: Energy in RLC Circuits -----*/

double rlc_series_energy_c(double C, double vc)
{
    return 0.5 * C * vc * vc;
}

double rlc_series_energy_l(double L, double il)
{
    return 0.5 * L * il * il;
}

double rlc_series_total_energy(double C, double L, double vc, double il)
{
    return 0.5 * C * vc * vc + 0.5 * L * il * il;
}

double rlc_series_energy_dissipated_R(double R, double il, double dt)
{
    return il * il * R * dt;
}

/*----- L6: Pole-Zero Map for Second-Order Systems -----*/

PoleZeroMap_t second_order_pole_zero_map(double zeta, double omega_n, double K)
{
    PoleZeroMap_t pz;
    memset(&pz, 0, sizeof(pz));

    pz.dc_gain = K;
    pz.n_poles = 2;
    pz.n_zeros = 0;

    if (zeta >= 1.0) {
        double alpha = zeta * omega_n;
        double disc = omega_n * sqrt(zeta * zeta - 1.0);
        pz.poles[0].sigma = -alpha + disc;
        pz.poles[0].omega = 0.0;
        pz.poles[1].sigma = -alpha - disc;
        pz.poles[1].omega = 0.0;
    } else if (zeta > 0.0) {
        double alpha = zeta * omega_n;
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        pz.poles[0].sigma = -alpha;
        pz.poles[0].omega = wd;
        pz.poles[1].sigma = -alpha;
        pz.poles[1].omega = -wd;
    } else {
        pz.poles[0].sigma = 0.0;
        pz.poles[0].omega = omega_n;
        pz.poles[1].sigma = 0.0;
        pz.poles[1].omega = -omega_n;
    }

    return pz;
}

/*----- L4: Normalized Responses -----*/

double normalized_step_response(double zeta, double omega_n, double t)
{
    if (t < 0.0) return 0.0;
    if (omega_n <= 0.0) return 1.0;

    if (zeta > 1.0 + 1e-10) {
        double disc = omega_n * sqrt(zeta * zeta - 1.0);
        double s1 = -zeta * omega_n + disc;
        double s2 = -zeta * omega_n - disc;
        double A2 = s1 / (s1 - s2);
        double A1 = -s2 / (s1 - s2);
        return 1.0 + A1 * exp(s1 * t) + A2 * exp(s2 * t);
    } else if (fabs(zeta - 1.0) < 1e-10) {
        return 1.0 - (1.0 + omega_n * t) * exp(-omega_n * t);
    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double alpha = zeta * omega_n;
        return 1.0 - exp(-alpha * t) * (cos(wd * t) + (alpha / wd) * sin(wd * t));
    } else {
        return 1.0 - cos(omega_n * t);
    }
}

double normalized_impulse_response(double zeta, double omega_n, double t)
{
    if (t < 0.0) return 0.0;
    if (omega_n <= 0.0) return 0.0;

    if (zeta > 1.0 + 1e-10) {
        double disc = omega_n * sqrt(zeta * zeta - 1.0);
        double s1 = -zeta * omega_n + disc;
        double s2 = -zeta * omega_n - disc;
        return (omega_n * omega_n / (2.0 * disc)) * (exp(s1 * t) - exp(s2 * t));
    } else if (fabs(zeta - 1.0) < 1e-10) {
        return omega_n * omega_n * t * exp(-omega_n * t);
    } else if (zeta > 0.0) {
        double wd = omega_n * sqrt(1.0 - zeta * zeta);
        double alpha = zeta * omega_n;
        return (omega_n * omega_n / wd) * exp(-alpha * t) * sin(wd * t);
    } else {
        return omega_n * sin(omega_n * t);
    }
}

/*----- L6: Design Formulas -----*/

double design_R_for_zeta_series(double zeta, double L, double C)
{
    if (L <= 0.0 || C <= 0.0) return -1.0;
    if (zeta <= 0.0) return 0.0;
    return 2.0 * zeta * sqrt(L / C);
}

double design_R_for_zeta_parallel(double zeta, double L, double C)
{
    if (L <= 0.0 || C <= 0.0) return -1.0;
    if (zeta <= 0.0) return INFINITY;
    return sqrt(L / C) / (2.0 * zeta);
}

double design_L_for_omega_n(double omega_n, double C)
{
    if (omega_n <= 0.0 || C <= 0.0) return -1.0;
    return 1.0 / (omega_n * omega_n * C);
}

double design_C_for_omega_n(double omega_n, double L)
{
    if (omega_n <= 0.0 || L <= 0.0) return -1.0;
    return 1.0 / (omega_n * omega_n * L);
}

/*----- L7: Second-Order Filter Transient -----*/

double lowpass_rc2_step_response(double Vs, double R1, double C1,
                                  double R2, double C2, double t)
{
    return rc_cascaded_vc(Vs, R1, C1, R2, C2, t, 0.0, 0.0);
}

/* Sallen-Key lowpass step response approximation
 * For Q <= 0.5 (two real poles): overdamped
 * For Q > 0.5: underdamped second-order response */
double sallen_key_lp_step_response(double Vs, double Q, double omega_n, double t)
{
    if (t < 0.0) return 0.0;
    double zeta = 1.0 / (2.0 * Q);
    return Vs * normalized_step_response(zeta, omega_n, t);
}

/*----- L8: Time-Varying RLC Parameters -----*/

double rlc_tv_step_vc(double Vs, double R0, double dRdt, double L,
                      double C, double v0, double i0, double t)
{
    if (t < 0.0) return v0;
    if (L <= 0.0 || C <= 0.0) return Vs;

    double R_t = R0 + dRdt * t;
    if (R_t < 0.0) R_t = 0.0;

    /* Approximate using instantaneous parameters
     * (exact solution requires numerical integration for TV systems) */
    double alpha_t = (L > 0.0) ? R_t / (2.0 * L) : INFINITY;
    double omega_n = 1.0 / sqrt(L * C);
    double zeta_t = (omega_n > 0.0) ? alpha_t / omega_n : INFINITY;

    if (zeta_t > 1.0) {
        double disc = sqrt(alpha_t * alpha_t - omega_n * omega_n);
        double s1 = -alpha_t + disc;
        double s2 = -alpha_t - disc;
        double A2 = (s1 * (v0 - Vs) - i0 / C) / (s1 - s2);
        double A1 = v0 - Vs - A2;
        return Vs + A1 * exp(s1 * t) + A2 * exp(s2 * t);
    } else {
        double wd = (zeta_t < 1.0) ? omega_n * sqrt(1.0 - zeta_t * zeta_t) : 0.0;
        if (wd > 0.0) {
            double A1 = v0 - Vs;
            double A2 = (i0 / C + alpha_t * (v0 - Vs)) / wd;
            return Vs + exp(-alpha_t * t) * (A1 * cos(wd * t) + A2 * sin(wd * t));
        } else {
            double A1 = v0 - Vs;
            double A2 = i0 / C + alpha_t * (v0 - Vs);
            return Vs + (A1 + A2 * t) * exp(-alpha_t * t);
        }
    }
}
