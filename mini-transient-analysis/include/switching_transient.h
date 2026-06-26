#ifndef SWITCHING_TRANSIENT_H
#define SWITCHING_TRANSIENT_H
#include "transient_defs.h"

/* L6-L7: Switching Transient Analysis
 * Reference: Mohan, Undeland & Robbins "Power Electronics" (2003),
 *           Erickson & Maksimovic "Fundamentals of Power Electronics" (2001),
 *           Baliga "Fundamentals of Power Semiconductor Devices" (2008)
 */

/*----- L6: BJT Switching Transients -----*/
double bjt_turn_on_delay(double Vcc, double Rc, double Rb, double Vbe_on, double Cbe, double Cbc, double beta);
double bjt_rise_time_collector(double Vcc, double Rc, double Ic_max, double Cbc, double beta, double f_T);
double bjt_storage_time(double Ic_sat, double Ib_fwd, double Ib_rev, double tau_s);
double bjt_fall_time_collector(double Vcc, double Rc, double Ic_sat, double Cbc, double f_T);
double bjt_switching_loss_per_cycle(double Vcc, double Ic, double t_on, double t_off, double f_sw);
double bjt_base_drive_current(double Ic, double beta_forced);
double bjt_anti_saturation_diode_current(double Ic, double Vce_sat, double Vd, double Rb);
double bjt_baker_clamp_design(double Ic, double beta, double Vbe, double Vd, double Rb);

/*----- L6: MOSFET Switching Transients -----*/
double mosfet_turn_on_delay(double Vth, double Vdrive, double Rg, double Ciss);
double mosfet_rise_time_drain(double Vdd, double Rdson, double Id, double Rg, double Cgd, double Vdrive, double Vth);
double mosfet_fall_time_drain(double Vdd, double Id, double Rg, double Cgd, double Vth);
double mosfet_gate_charge_total(double Qgs, double Qgd, double Qg_remainder);
double mosfet_gate_drive_current(double Qg_total, double t_sw_desired);
double mosfet_gate_resistor(double Vdrive, double Ig_peak);
double mosfet_miller_plateau_duration(double Qgd, double Vdrive, double Vth, double Rg);
double mosfet_switching_energy(double Vds, double Id, double t_on, double t_off);
double mosfet_output_capacitance_loss(double Coss, double Vds, double f_sw);

/*----- L6: Diode Reverse Recovery -----*/
double diode_reverse_recovery_charge(double Irr, double t_rr);
double diode_reverse_recovery_peak(double di_dt, double Q_rr);
double diode_reverse_recovery_time(double Irr, double di_dt);
double diode_recovery_softness_factor(double ta, double tb);
double diode_forward_recovery_voltage(double Vf_max, double t_fr, double t);
double diode_recovery_loss_per_cycle(double Vr, double Irr, double t_rr, double f_sw);

/*----- L7: Buck Converter Transients -----*/
double buck_startup_inrush(double Vin, double Cout, double L, double t);
double buck_load_step_undershoot(double delta_I, double Cout, double L, double f_sw, double duty);
double buck_load_step_settling_time(double L, double Cout, double zeta_closed_loop);
double buck_inductor_current_ripple(double Vin, double Vout, double L, double f_sw);
double buck_output_voltage_ripple(double delta_IL, double Cout, double ESR, double f_sw);
double buck_critical_inductance(double Vin, double Vout, double Iout_min, double f_sw);
double buck_duty_cycle_ideal(double Vout, double Vin);
double buck_startup_soft_start_time(double Cout, double I_limit, double Vout);

/*----- L7: Boost Converter Transients -----*/
double boost_startup_current(double Vin, double L, double R_load, double t);
double boost_duty_cycle_ideal(double Vout, double Vin);
double boost_inductor_current_ripple(double Vin, double L, double duty, double f_sw);
double boost_output_voltage_ripple(double Vout, double duty, double Cout, double R_load, double f_sw);
double boost_right_half_plane_zero(double Vout, double L, double duty, double R_load);

/*----- L7: Gate Driver Transients -----*/
double gate_driver_peak_current(double Vdrive, double Rg_total, double Rg_internal);
double gate_driver_power_dissipation(double Qg_total, double Vdrive, double f_sw);
double gate_driver_bootstrap_capacitor(double Qg_total, double delta_V_allowed, double Q_leakage, double f_sw);
double gate_driver_dead_time(double t_off_max, double t_on_min, double safety_margin);

/*----- L7: Snubber Circuit Design -----*/
double snubber_rc_R(double V_peak, double I_peak, double safety_factor);
double snubber_rc_C(double I_peak, double t_rise, double V_peak);
double snubber_power_loss(double C, double V_peak, double f_sw);
double snubber_rc_turn_off_capacitor(double I_load, double t_fall, double V_max_allowed);
double snubber_diode_recovery(double L_stray, double Irr, double V_clamp);

/*----- L7: Relay and Solenoid Transients -----*/
double relay_coil_time_constant(double L_coil, double R_coil);
double relay_pull_in_time(double V_drive, double V_pull_in, double tau);
double relay_drop_out_time(double V_hold, double V_drop_out, double tau);
double relay_flyback_voltage(double I_coil, double L_coil, double C_snub);
double relay_release_voltage_spike(double Vs, double L_coil, double R_coil, double t_off);
double solenoid_peak_current(double V_drive, double R_coil);
double solenoid_pwm_hold_duty(double V_drive, double I_hold, double I_peak);
double solenoid_response_time(double L_coil, double R_coil, double I_target, double V_drive);

/*----- L7: DC Motor Driver Transients -----*/
double dc_motor_startup_current(double V_supply, double R_winding);
double dc_motor_electrical_tau(double L_winding, double R_winding);
double dc_motor_mechanical_tau(double J_rotor, double R_winding, double Kt, double Ke);
double dc_motor_current_ripple_pwm(double V_supply, double L_winding, double f_pwm, double duty);
double dc_motor_back_emf(double Ke, double rpm);
double h_bridge_shoot_through_dead_time(double t_rise, double t_fall, double safety_margin);

/*----- L7: ESD and Protection Transients -----*/
double hbm_esd_current_peak(double V_esd, double R_hbm);
double hbm_esd_time_constant(double C_hbm, double R_hbm);
double hbm_esd_energy(double C_hbm, double V_esd);
double cdm_esd_current_peak(double C_device, double V_esd, double R_discharge);
double tvs_clamping_voltage_transient(double V_br, double I_pp, double R_dyn);
double tvs_response_time(double C_junction, double R_dyn);

/*----- L7: Power Sequencing and Supervision -----*/
double power_sequencing_delay(double R, double C, double V_threshold, double V_supply);
double voltage_supervisor_timeout(double R, double C, double Vcc, double V_threshold);
double inrush_current_limiter_energy(double C_bulk, double Vin);
double soft_start_inrush_reduction_pct(double I_inrush, double I_rated);

/*----- L7: Transmission Line Transients -----*/
double tl_reflection_coefficient(double ZL, double Z0);
double tl_propagation_delay(double length_m, double velocity_factor);
double tl_voltage_at_load_initial(double V_incident, double ZL, double Z0);
double tl_source_reflection_coefficient(double ZS, double Z0);
double tl_load_reflection_coefficient(double ZL, double Z0);
double tl_voltage_bounce_diagram(double V0, double ZS, double ZL, double Z0, double t, double Td);
double tl_ringing_frequency(double length_m, double velocity_factor);
double tl_critical_length(double t_rise, double velocity_factor);
int tl_is_lumped_or_distributed(double length_m, double t_rise, double velocity_factor);

/*----- L7: Transformer Inrush -----*/
double transformer_inrush_current_peak(double V_peak, double R_winding, double L_winding, double Br, double Bs);
double transformer_inrush_decay_tau(double L_winding, double R_winding);
double transformer_inrush_current_at_t(double I_peak, double tau, double t);

/*----- L7/L8: PLL Transient -----*/
double pll_natural_frequency(double Kpd, double Kvco, double N);
double pll_damping_factor(double Kpd, double Kvco, double N, double R1, double C1);
double pll_lock_time_estimate(double wn, double zeta, double freq_step, double tol);
double pll_phase_error_step_response(double delta_phi, double wn, double zeta, double t);
double pll_frequency_step_response(double delta_f, double wn, double zeta, double t);
double pll_loop_filter_corner(double R1, double C1);

/*----- L7: Capacitor Inrush and Lifetime -----*/
double capacitor_inrush_current(double Vin, double ESR, double ESL_cap);
double capacitor_charge_time_to_voltage(double C, double R, double Vin, double V_target);
double capacitor_self_heating_temp_rise(double I_rms, double ESR, double Rth);
double capacitor_lifetime_derating(double L0, double T_op, double T_rated, double V_op, double V_rated, double n);

#endif /* SWITCHING_TRANSIENT_H */
