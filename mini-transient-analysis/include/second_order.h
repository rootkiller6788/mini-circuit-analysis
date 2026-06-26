#ifndef SECOND_ORDER_H
#define SECOND_ORDER_H
#include "transient_defs.h"

double rlc_series_alpha(double R, double L);
double rlc_series_omega_n(double L, double C);
double rlc_series_zeta(double R, double L, double C);
double rlc_series_omega_d(double R, double L, double C);
DampingClass_t rlc_series_damping_class(double R, double L, double C);

double rlc_parallel_alpha(double R, double C);
double rlc_parallel_omega_n(double L, double C);
double rlc_parallel_zeta(double R, double L, double C);
double rlc_parallel_omega_d(double R, double L, double C);
DampingClass_t rlc_parallel_damping_class(double R, double L, double C);

double rlc_series_step_vc(double Vs, double R, double L, double C,
                          double v0, double i0, double t);
double rlc_series_step_il(double Vs, double R, double L, double C,
                          double v0, double i0, double t);
double rlc_series_step_vl(double Vs, double R, double L, double C,
                          double v0, double i0, double t);
double rlc_series_step_vr(double Vs, double R, double L, double C,
                          double v0, double i0, double t);

double rlc_series_overdamped_vc(double Vs, double v0, double i0,
                                double s1, double s2, double C, double t);
double rlc_series_overdamped_il(double Vs, double v0, double i0,
                                double s1, double s2, double L, double t);

double rlc_series_critical_vc(double Vs, double v0, double i0,
                              double alpha, double C, double t);
double rlc_series_critical_il(double Vs, double v0, double i0,
                              double alpha, double L, double t);

double rlc_series_underdamped_vc(double Vs, double v0, double i0,
                                 double alpha, double omega_d, double C, double t);
double rlc_series_underdamped_il(double Vs, double v0, double i0,
                                 double alpha, double omega_d, double L, double t);

double rlc_series_undamped_vc(double Vs, double v0, double i0,
                              double omega_n, double C, double t);
double rlc_series_undamped_il(double Vs, double v0, double i0,
                              double omega_n, double L, double t);

double rlc_series_natural_vc(double R, double L, double C,
                             double v0, double i0, double t);
double rlc_series_natural_il(double R, double L, double C,
                             double v0, double i0, double t);

double rlc_parallel_step_v(double Is, double R, double L, double C,
                           double v0, double i0, double t);
double rlc_parallel_step_il(double Is, double R, double L, double C,
                            double v0, double i0, double t);
double rlc_parallel_step_ic(double Is, double R, double L, double C,
                            double v0, double i0, double t);

double rlc_parallel_natural_v(double R, double L, double C,
                              double v0, double i0, double t);
double rlc_parallel_natural_il(double R, double L, double C,
                               double v0, double i0, double t);

double rlc_series_impulse_vc(double R, double L, double C, double t);
double rlc_series_impulse_il(double R, double L, double C, double t);
double rlc_parallel_impulse_v(double R, double L, double C, double t);

double rlc_series_sinusoidal_vc(double Vm, double omega, double phi,
                                double R, double L, double C,
                                double v0, double i0, double t);
double rlc_parallel_sinusoidal_v(double Im, double omega, double phi,
                                  double R, double L, double C,
                                  double v0, double i0, double t);

double second_order_overshoot_pct(double zeta);
double second_order_peak_time(double omega_n, double zeta);
double second_order_settling_time_2pct(double zeta, double omega_n);
double second_order_settling_time_5pct(double zeta, double omega_n);
double second_order_rise_time_10_90(double omega_n, double zeta);
double second_order_delay_time(double omega_n, double zeta);
double second_order_period_of_oscillation(double omega_d);
double second_order_decay_ratio(double zeta);
double second_order_num_oscillations_to_settle(double zeta);

void second_order_characteristic_roots(double a, double b, double c,
                                        double *r1_real, double *r1_imag,
                                        double *r2_real, double *r2_imag);
double second_order_discriminant(double a, double b, double c);
DampingClass_t second_order_classify_by_roots(double discriminant);

double rlc_series_energy_c(double C, double vc);
double rlc_series_energy_l(double L, double il);
double rlc_series_total_energy(double C, double L, double vc, double il);
double rlc_series_energy_dissipated_R(double R, double il, double dt);

PoleZeroMap_t second_order_pole_zero_map(double zeta, double omega_n, double K);

double normalized_step_response(double zeta, double omega_n, double t);
double normalized_impulse_response(double zeta, double omega_n, double t);

double design_R_for_zeta_series(double zeta, double L, double C);
double design_R_for_zeta_parallel(double zeta, double L, double C);
double design_L_for_omega_n(double omega_n, double C);
double design_C_for_omega_n(double omega_n, double L);

double lowpass_rc2_step_response(double Vs, double R1, double C1,
                                  double R2, double C2, double t);
double sallen_key_lp_step_response(double Vs, double Q, double omega_n, double t);

double rlc_tv_step_vc(double Vs, double R0, double dRdt, double L,
                      double C, double v0, double i0, double t);

#endif /* SECOND_ORDER_H */
