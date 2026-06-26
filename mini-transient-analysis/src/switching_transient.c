#include "switching_transient.h"
#include "transient_defs.h"
#include <math.h>
#include <stdlib.h>

/*==================================================================
 * L6-L7: Switching Transient Applications
 *
 * This module implements transient analysis for semiconductor
 * switching circuits, power converters, relay/solenoid drivers,
 * motor drivers, ESD protection, transmission lines, and PLLs.
 *
 * Reference: Mohan, Undeland & Robbins "Power Electronics" (2003),
 *           Erickson & Maksimovic "Fundamentals of Power Electronics" (2001),
 *           Baliga "Power Semiconductor Devices" (2008),
 *           Sedra & Smith "Microelectronic Circuits" (2020)
 *==================================================================*/

/*----- L6: BJT Switching Transients -----*/

double bjt_turn_on_delay(double Vcc, double Rc, double Rb, double Vbe_on, double Cbe, double Cbc, double beta)
{
    if (Rb <= 0.0 || Cbe <= 0.0) return 0.0;
    double tau = Rb * (Cbe + Cbc * (1.0 + beta * Rc / Rb));
    if (tau <= 0.0) return 0.0;
    double Vb_final = Vcc * Rb / (Rb + Rc * beta);
    return tau * log(Vb_final / (Vb_final - Vbe_on + 1e-12));
}

double bjt_rise_time_collector(double Vcc, double Rc, double Ic_max, double Cbc, double beta, double f_T)
{
    if (f_T <= 0.0 || Ic_max <= 0.0) return INFINITY;
    double tau_eff = beta / (2.0 * M_PI * f_T);
    double Vce_change = Vcc - 0.2;
    return tau_eff * log((Ic_max * Rc) / (Ic_max * Rc - Vce_change + 1e-12));
}

double bjt_storage_time(double Ic_sat, double Ib_fwd, double Ib_rev, double tau_s)
{
    if (tau_s <= 0.0) return 0.0;
    if (Ib_fwd + Ib_rev <= 0.0) return INFINITY;
    return tau_s * log((Ib_fwd + Ib_rev) / (Ic_sat / 10.0 + Ib_rev + 1e-12));
}

double bjt_fall_time_collector(double Vcc, double Rc, double Ic_sat, double Cbc, double f_T)
{
    if (f_T <= 0.0) return INFINITY;
    double tau = 1.0 / (2.0 * M_PI * f_T);
    return tau * log((Vcc / Rc) / (0.1 * Ic_sat + 1e-12));
}

double bjt_switching_loss_per_cycle(double Vcc, double Ic, double t_on, double t_off, double f_sw)
{
    return 0.5 * Vcc * Ic * (t_on + t_off) * f_sw;
}

double bjt_base_drive_current(double Ic, double beta_forced)
{
    if (beta_forced <= 0.0) return INFINITY;
    return Ic / beta_forced;
}

double bjt_anti_saturation_diode_current(double Ic, double Vce_sat, double Vd, double Rb)
{
    if (Rb <= 0.0) return 0.0;
    return (Vd - Vce_sat) / Rb;
}

double bjt_baker_clamp_design(double Ic, double beta, double Vbe, double Vd, double Rb)
{
    if (beta <= 0.0 || Rb <= 0.0) return 0.0;
    double Ib = Ic / beta;
    return Ib + (Vd - Vbe) / Rb;
}

/*----- L6: MOSFET Switching Transients -----*/

double mosfet_turn_on_delay(double Vth, double Vdrive, double Rg, double Ciss)
{
    if (Rg <= 0.0 || Ciss <= 0.0 || Vdrive <= Vth) return INFINITY;
    double tau = Rg * Ciss;
    return tau * log(Vdrive / (Vdrive - Vth + 1e-12));
}

double mosfet_rise_time_drain(double Vdd, double Rdson, double Id, double Rg, double Cgd, double Vdrive, double Vth)
{
    if (Rg <= 0.0 || Cgd <= 0.0) return INFINITY;
    double Vgp = Vth + Id / (2.0 * 1.0); /* transconductance approximated */
    double ig = (Vdrive - Vgp) / Rg;
    if (ig <= 0.0) return INFINITY;
    return Cgd * (Vdd - Id * Rdson) / ig;
}

double mosfet_fall_time_drain(double Vdd, double Id, double Rg, double Cgd, double Vth)
{
    if (Rg <= 0.0 || Cgd <= 0.0) return INFINITY;
    double ig = Vth / Rg;
    if (ig <= 0.0) return INFINITY;
    return Cgd * Vdd / ig;
}

double mosfet_gate_charge_total(double Qgs, double Qgd, double Qg_remainder)
{
    return Qgs + Qgd + Qg_remainder;
}

double mosfet_gate_drive_current(double Qg_total, double t_sw_desired)
{
    if (t_sw_desired <= 0.0) return INFINITY;
    return Qg_total / t_sw_desired;
}

double mosfet_gate_resistor(double Vdrive, double Ig_peak)
{
    if (Ig_peak <= 0.0) return INFINITY;
    return Vdrive / Ig_peak;
}

double mosfet_miller_plateau_duration(double Qgd, double Vdrive, double Vth, double Rg)
{
    if (Rg <= 0.0) return INFINITY;
    double Ig_miller = (Vdrive - Vth) / Rg;
    if (Ig_miller <= 0.0) return INFINITY;
    return Qgd / Ig_miller;
}

double mosfet_switching_energy(double Vds, double Id, double t_on, double t_off)
{
    return 0.5 * Vds * Id * (t_on + t_off);
}

double mosfet_output_capacitance_loss(double Coss, double Vds, double f_sw)
{
    return 0.5 * Coss * Vds * Vds * f_sw;
}

/*----- L6: Diode Reverse Recovery -----*/

double diode_reverse_recovery_charge(double Irr, double t_rr)
{
    return 0.5 * Irr * t_rr;
}

double diode_reverse_recovery_peak(double di_dt, double Q_rr)
{
    if (di_dt <= 0.0) return 0.0;
    return sqrt(2.0 * Q_rr * di_dt);
}

double diode_reverse_recovery_time(double Irr, double di_dt)
{
    if (di_dt <= 0.0) return INFINITY;
    return Irr / di_dt;
}

double diode_recovery_softness_factor(double ta, double tb)
{
    if (tb <= 0.0) return 0.0;
    return ta / tb;
}

double diode_forward_recovery_voltage(double Vf_max, double t_fr, double t)
{
    if (t_fr <= 0.0) return Vf_max;
    if (t >= t_fr) return 0.7;
    return 0.7 + (Vf_max - 0.7) * (1.0 - t / t_fr);
}

double diode_recovery_loss_per_cycle(double Vr, double Irr, double t_rr, double f_sw)
{
    return 0.25 * Vr * Irr * t_rr * f_sw;
}

/*----- L7: Buck Converter Transients -----*/

double buck_startup_inrush(double Vin, double Cout, double L, double t)
{
    if (L <= 0.0 || Cout <= 0.0) return 0.0;
    double omega = 1.0 / sqrt(L * Cout);
    return Vin * sqrt(Cout / L) * sin(omega * t);
}

double buck_load_step_undershoot(double delta_I, double Cout, double L, double f_sw, double duty)
{
    if (Cout <= 0.0) return INFINITY;
    double t_response = 1.0 / f_sw;
    return delta_I * t_response / Cout + 0.1 * delta_I * duty / Cout;
}

double buck_load_step_settling_time(double L, double Cout, double zeta_closed_loop)
{
    if (zeta_closed_loop <= 0.0) return INFINITY;
    double omega_n = 1.0 / sqrt(L * Cout);
    return 4.0 / (zeta_closed_loop * omega_n);
}

double buck_inductor_current_ripple(double Vin, double Vout, double L, double f_sw)
{
    if (L <= 0.0 || f_sw <= 0.0) return INFINITY;
    double duty = Vout / Vin;
    return Vin * duty * (1.0 - duty) / (L * f_sw);
}

double buck_output_voltage_ripple(double delta_IL, double Cout, double ESR, double f_sw)
{
    if (Cout <= 0.0 || f_sw <= 0.0) return INFINITY;
    double delta_v_cap = delta_IL / (8.0 * Cout * f_sw);
    double delta_v_esr = delta_IL * ESR;
    return delta_v_cap + delta_v_esr;
}

double buck_critical_inductance(double Vin, double Vout, double Iout_min, double f_sw)
{
    if (f_sw <= 0.0 || Iout_min <= 0.0) return INFINITY;
    double duty = Vout / Vin;
    return Vin * duty * (1.0 - duty) / (2.0 * Iout_min * f_sw);
}

double buck_duty_cycle_ideal(double Vout, double Vin)
{
    if (Vin <= 0.0) return 0.0;
    return Vout / Vin;
}

double buck_startup_soft_start_time(double Cout, double I_limit, double Vout)
{
    if (I_limit <= 0.0) return INFINITY;
    return Cout * Vout / I_limit;
}

/*----- L7: Boost Converter Transients -----*/

double boost_startup_current(double Vin, double L, double R_load, double t)
{
    if (L <= 0.0) return 0.0;
    double tau = L / (R_load + 1e-6);
    double I_final = Vin / (R_load + 1e-6);
    return I_final * (1.0 - exp(-t / tau));
}

double boost_duty_cycle_ideal(double Vout, double Vin)
{
    if (Vout <= Vin) return 0.0;
    return 1.0 - Vin / Vout;
}

double boost_inductor_current_ripple(double Vin, double L, double duty, double f_sw)
{
    if (L <= 0.0 || f_sw <= 0.0) return INFINITY;
    return Vin * duty / (L * f_sw);
}

double boost_output_voltage_ripple(double Vout, double duty, double Cout, double R_load, double f_sw)
{
    if (Cout <= 0.0 || f_sw <= 0.0 || R_load <= 0.0) return INFINITY;
    return Vout * duty / (R_load * Cout * f_sw);
}

double boost_right_half_plane_zero(double Vout, double L, double duty, double R_load)
{
    if (L <= 0.0 || R_load <= 0.0) return 0.0;
    return R_load * (1.0 - duty) * (1.0 - duty) / L;
}

/*----- L7: Gate Driver Transients -----*/

double gate_driver_peak_current(double Vdrive, double Rg_total, double Rg_internal)
{
    double R_total = Rg_total + Rg_internal;
    if (R_total <= 0.0) return INFINITY;
    return Vdrive / R_total;
}

double gate_driver_power_dissipation(double Qg_total, double Vdrive, double f_sw)
{
    return Qg_total * Vdrive * f_sw;
}

double gate_driver_bootstrap_capacitor(double Qg_total, double delta_V_allowed, double Q_leakage, double f_sw)
{
    if (delta_V_allowed <= 0.0) return INFINITY;
    return (Qg_total + Q_leakage / f_sw) / delta_V_allowed;
}

double gate_driver_dead_time(double t_off_max, double t_on_min, double safety_margin)
{
    return (t_off_max - t_on_min) * safety_margin;
}

/*----- L7: Snubber Circuit Design -----*/

double snubber_rc_R(double V_peak, double I_peak, double safety_factor)
{
    if (I_peak <= 0.0) return INFINITY;
    return V_peak / (I_peak * safety_factor);
}

double snubber_rc_C(double I_peak, double t_rise, double V_peak)
{
    if (V_peak <= 0.0 || t_rise <= 0.0) return 0.0;
    return I_peak * t_rise / V_peak;
}

double snubber_power_loss(double C, double V_peak, double f_sw)
{
    return 0.5 * C * V_peak * V_peak * f_sw;
}

double snubber_rc_turn_off_capacitor(double I_load, double t_fall, double V_max_allowed)
{
    if (V_max_allowed <= 0.0) return INFINITY;
    return I_load * t_fall / V_max_allowed;
}

double snubber_diode_recovery(double L_stray, double Irr, double V_clamp)
{
    if (V_clamp <= 0.0) return 0.0;
    return 0.5 * L_stray * Irr * Irr / V_clamp;
}

/*----- L7: Relay and Solenoid Transients -----*/

double relay_coil_time_constant(double L_coil, double R_coil)
{
    if (R_coil <= 0.0) return INFINITY;
    return L_coil / R_coil;
}

double relay_pull_in_time(double V_drive, double V_pull_in, double tau)
{
    if (V_drive <= V_pull_in || tau <= 0.0) return INFINITY;
    return tau * log(V_drive / (V_drive - V_pull_in + 1e-12));
}

double relay_drop_out_time(double V_hold, double V_drop_out, double tau)
{
    if (V_drop_out <= 0.0 || tau <= 0.0) return INFINITY;
    return tau * log(V_hold / V_drop_out);
}

double relay_flyback_voltage(double I_coil, double L_coil, double C_snub)
{
    if (C_snub <= 0.0) return INFINITY;
    return I_coil * sqrt(L_coil / C_snub);
}

double relay_release_voltage_spike(double Vs, double L_coil, double R_coil, double t_off)
{
    double tau = relay_coil_time_constant(L_coil, R_coil);
    if (tau <= 0.0 || R_coil <= 0.0) return INFINITY;
    double I_at_off = Vs / R_coil;
    double di_dt = -I_at_off / tau;
    return -L_coil * di_dt;
}

double solenoid_peak_current(double V_drive, double R_coil)
{
    if (R_coil <= 0.0) return INFINITY;
    return V_drive / R_coil;
}

double solenoid_pwm_hold_duty(double V_drive, double I_hold, double I_peak)
{
    if (I_peak <= 0.0) return 0.0;
    return I_hold / I_peak;
}

double solenoid_response_time(double L_coil, double R_coil, double I_target, double V_drive)
{
    double tau = relay_coil_time_constant(L_coil, R_coil);
    if (tau <= 0.0 || R_coil <= 0.0) return INFINITY;
    double I_final = V_drive / R_coil;
    if (I_target >= I_final) return INFINITY;
    return tau * log(I_final / (I_final - I_target + 1e-12));
}

/*----- L7: DC Motor Driver Transients -----*/

double dc_motor_startup_current(double V_supply, double R_winding)
{
    if (R_winding <= 0.0) return INFINITY;
    return V_supply / R_winding;
}

double dc_motor_electrical_tau(double L_winding, double R_winding)
{
    return relay_coil_time_constant(L_winding, R_winding);
}

double dc_motor_mechanical_tau(double J_rotor, double R_winding, double Kt, double Ke)
{
    double D = Kt * Ke / R_winding;
    if (D <= 0.0) return INFINITY;
    return J_rotor / D;
}

double dc_motor_current_ripple_pwm(double V_supply, double L_winding, double f_pwm, double duty)
{
    if (L_winding <= 0.0 || f_pwm <= 0.0) return INFINITY;
    return V_supply * duty * (1.0 - duty) / (L_winding * f_pwm);
}

double dc_motor_back_emf(double Ke, double rpm)
{
    return Ke * rpm;
}

double h_bridge_shoot_through_dead_time(double t_rise, double t_fall, double safety_margin)
{
    return (t_rise + t_fall) * safety_margin;
}

/*----- L7: ESD and Protection Transients -----*/

double hbm_esd_current_peak(double V_esd, double R_hbm)
{
    if (R_hbm <= 0.0) return 0.0;
    return V_esd / R_hbm;
}

double hbm_esd_time_constant(double C_hbm, double R_hbm)
{
    return C_hbm * R_hbm;
}

double hbm_esd_energy(double C_hbm, double V_esd)
{
    return 0.5 * C_hbm * V_esd * V_esd;
}

double cdm_esd_current_peak(double C_device, double V_esd, double R_discharge)
{
    if (R_discharge <= 0.0) return 0.0;
    return V_esd / R_discharge;
}

double tvs_clamping_voltage_transient(double V_br, double I_pp, double R_dyn)
{
    return V_br + I_pp * R_dyn;
}

double tvs_response_time(double C_junction, double R_dyn)
{
    return C_junction * R_dyn;
}

/*----- L7: Power Sequencing and Supervision -----*/

double power_sequencing_delay(double R, double C, double V_threshold, double V_supply)
{
    if (R <= 0.0 || C <= 0.0 || V_threshold >= V_supply) return INFINITY;
    return -R * C * log(1.0 - V_threshold / V_supply);
}

double voltage_supervisor_timeout(double R, double C, double Vcc, double V_threshold)
{
    return power_sequencing_delay(R, C, V_threshold, Vcc);
}

double inrush_current_limiter_energy(double C_bulk, double Vin)
{
    return 0.5 * C_bulk * Vin * Vin;
}

double soft_start_inrush_reduction_pct(double I_inrush, double I_rated)
{
    if (I_inrush <= 0.0) return 0.0;
    return (1.0 - I_rated / I_inrush) * 100.0;
}

/*----- L7: Transmission Line Transients -----*/

double tl_reflection_coefficient(double ZL, double Z0)
{
    if (ZL + Z0 == 0.0) return 1.0;
    return (ZL - Z0) / (ZL + Z0);
}

double tl_propagation_delay(double length_m, double velocity_factor)
{
    if (velocity_factor <= 0.0) return INFINITY;
    return length_m / (3.0e8 * velocity_factor);
}

double tl_voltage_at_load_initial(double V_incident, double ZL, double Z0)
{
    return V_incident * (1.0 + tl_reflection_coefficient(ZL, Z0));
}

double tl_source_reflection_coefficient(double ZS, double Z0)
{
    return tl_reflection_coefficient(ZS, Z0);
}

double tl_load_reflection_coefficient(double ZL, double Z0)
{
    return tl_reflection_coefficient(ZL, Z0);
}

double tl_voltage_bounce_diagram(double V0, double ZS, double ZL, double Z0, double t, double Td)
{
    if (Td <= 0.0) return V0 * ZL / (ZL + ZS);
    double Gamma_S = tl_source_reflection_coefficient(ZS, Z0);
    double Gamma_L = tl_load_reflection_coefficient(ZL, Z0);
    double T = 1.0 + Gamma_L;

    int n_lattice = (int)(t / Td);
    double V = 0.0;
    double V_forward = V0 * Z0 / (Z0 + ZS);

    for (int i = 0; i <= n_lattice && i < 100; i++) {
        double V_inc = V_forward * pow(Gamma_S * Gamma_L, (double)i);
        if (i > 0) V_inc *= Gamma_L;
        double t_arrival = (2.0 * i + 1.0) * Td;
        if (t >= t_arrival) V += V_inc * T;
    }
    return V;
}

double tl_ringing_frequency(double length_m, double velocity_factor)
{
    double Td = tl_propagation_delay(length_m, velocity_factor);
    if (Td <= 0.0) return 0.0;
    return 1.0 / (4.0 * Td);
}

double tl_critical_length(double t_rise, double velocity_factor)
{
    if (velocity_factor <= 0.0) return INFINITY;
    return 3.0e8 * velocity_factor * t_rise / 2.0;
}

int tl_is_lumped_or_distributed(double length_m, double t_rise, double velocity_factor)
{
    double L_crit = tl_critical_length(t_rise, velocity_factor);
    return (length_m > L_crit) ? 1 : 0; /* 0=lumped, 1=distributed */
}

/*----- L7: Transformer Inrush -----*/

double transformer_inrush_current_peak(double V_peak, double R_winding, double L_winding, double Br, double Bs)
{
    if (R_winding <= 0.0) return INFINITY;
    double phi_offset = (Bs - Br) / Bs;
    return V_peak * phi_offset / R_winding;
}

double transformer_inrush_decay_tau(double L_winding, double R_winding)
{
    return relay_coil_time_constant(L_winding, R_winding);
}

double transformer_inrush_current_at_t(double I_peak, double tau, double t)
{
    if (tau <= 0.0) return I_peak;
    return I_peak * exp(-t / tau);
}

/*----- L7/L8: PLL Transient -----*/

double pll_natural_frequency(double Kpd, double Kvco, double N)
{
    if (N <= 0.0) return 0.0;
    return sqrt(Kpd * Kvco / N);
}

double pll_damping_factor(double Kpd, double Kvco, double N, double R1, double C1)
{
    double omega_n = pll_natural_frequency(Kpd, Kvco, N);
    if (omega_n <= 0.0) return 0.0;
    return 0.5 * omega_n * R1 * C1;
}

double pll_lock_time_estimate(double wn, double zeta, double freq_step, double tol)
{
    if (wn <= 0.0 || zeta <= 0.0 || tol <= 0.0) return INFINITY;
    return -log(tol / freq_step) / (zeta * wn);
}

double pll_phase_error_step_response(double delta_phi, double wn, double zeta, double t)
{
    if (t < 0.0) return delta_phi;
    if (zeta >= 1.0) {
        double alpha = zeta * wn;
        double disc = wn * sqrt(zeta * zeta - 1.0);
        double s1 = -alpha + disc, s2 = -alpha - disc;
        double A1 = delta_phi * s2 / (s2 - s1);
        double A2 = delta_phi * s1 / (s1 - s2);
        return A1 * exp(s1 * t) + A2 * exp(s2 * t);
    } else {
        double wd = wn * sqrt(1.0 - zeta * zeta);
        double alpha = zeta * wn;
        return delta_phi * exp(-alpha * t) * (cos(wd * t) + (alpha / wd) * sin(wd * t));
    }
}

double pll_frequency_step_response(double delta_f, double wn, double zeta, double t)
{
    if (t < 0.0) return 0.0;
    if (zeta < 1.0) {
        double wd = wn * sqrt(1.0 - zeta * zeta);
        return delta_f * (1.0 - exp(-zeta * wn * t) * (cos(wd * t) + (zeta * wn / wd) * sin(wd * t)));
    }
    return delta_f * (1.0 - exp(-zeta * wn * t) * (1.0 + zeta * wn * t));
}

double pll_loop_filter_corner(double R1, double C1)
{
    if (R1 <= 0.0 || C1 <= 0.0) return 0.0;
    return 1.0 / (2.0 * M_PI * R1 * C1);
}

/*----- L7: Capacitor Inrush and Lifetime -----*/

double capacitor_inrush_current(double Vin, double ESR, double ESL_cap)
{
    double R_eff = ESR;
    if (R_eff <= 0.0) R_eff = 0.01;
    return Vin / R_eff;
}

double capacitor_charge_time_to_voltage(double C, double R, double Vin, double V_target)
{
    if (R <= 0.0 || C <= 0.0 || V_target >= Vin) return INFINITY;
    double tau = R * C;
    return -tau * log(1.0 - V_target / Vin);
}

double capacitor_self_heating_temp_rise(double I_rms, double ESR, double Rth)
{
    return I_rms * I_rms * ESR * Rth;
}

double capacitor_lifetime_derating(double L0, double T_op, double T_rated, double V_op, double V_rated, double n)
{
    if (V_rated <= 0.0 || T_rated <= T_op) return 0.0;
    double V_factor = pow(V_op / V_rated, n);
    double T_factor = pow(2.0, (T_rated - T_op) / 10.0);
    return L0 * V_factor * T_factor;
}
