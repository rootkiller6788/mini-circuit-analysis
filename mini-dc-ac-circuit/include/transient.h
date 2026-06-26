#ifndef TRANSIENT_H
#define TRANSIENT_H
#include "circuit_elements.h"
typedef struct { double vC_init; double iL_init; } InitialConditions_t;
typedef enum { SOLVER_ANALYTIC=0, SOLVER_EULER=1, SOLVER_RK4=2, SOLVER_TRAP=3 } SolverType_t;
double rc_time_constant(double R, double C);
double rl_time_constant(double R, double L);
double rc_charging_vc(double Vf, double v0, double tau, double t);
double rc_charging_i(double Vf, double v0, double R, double tau, double t);
double rc_discharging_vc(double v0, double tau, double t);
double rc_discharging_i(double v0, double R, double tau, double t);
double rl_charging_il(double If, double i0, double tau, double t);
double rl_charging_vl(double Vs, double i0, double R, double tau, double t);
double rl_discharging_il(double i0, double tau, double t);
double capacitor_energy(double C, double v);
double inductor_energy(double L, double i);
typedef enum { RLC_OVERDAMPED=0,RLC_CRITICALLY_DAMPED=1,RLC_UNDERDAMPED=2,RLC_UNDAMPED=3 } RLC_Regime_t;
double rlc_series_alpha(double R, double L);
double rlc_series_omega0(double L, double C);
double rlc_series_zeta(double R, double L, double C);
double rlc_series_omega_d(double R, double L, double C);
RLC_Regime_t rlc_series_regime(double R, double L, double C);
double rlc_series_step_vc(double Vs,double R,double L,double C,double v0,double i0,double t);
double rlc_series_step_il(double Vs,double R,double L,double C,double v0,double i0,double t);
double rlc_series_natural_vc(double R,double L,double C,double v0,double i0,double t);
double rlc_series_natural_il(double R,double L,double C,double v0,double i0,double t);
double rlc_parallel_step_il(double Is,double R,double L,double C,double i0,double v0,double t);
double rlc_parallel_step_vc(double Is,double R,double L,double C,double i0,double v0,double t);
double rlc_parallel_alpha(double R, double C);
double rlc_parallel_zeta(double R, double L, double C);
typedef struct { int n_states; double *v_caps; double *i_inductors; double t; double dt; SolverType_t solver; } TransientState_t;
int transient_euler_step(TransientState_t *st,void (*f)(const TransientState_t*,double*,double),double t);
int transient_rk4_step(TransientState_t *st,void (*f)(const TransientState_t*,double*,double),double t);
int transient_simulate(TransientState_t *st,void (*f)(const TransientState_t*,double*,double),double ts,double te,double dt,void (*cb)(double,const TransientState_t*,void*),void *ud);
double rc_ac_complete_vc(double Vm,double w,double ph,double R,double C,double v0,double t);
double rl_ac_complete_il(double Vm,double w,double ph,double R,double L,double i0,double t);
double timer_555_astable_freq(double R1,double R2,double C);
double timer_555_astable_duty(double R1,double R2);
double timer_555_monostable_pulse(double R,double C);
double rc_oscillator_freq_symmetric(double R,double C);
double rc_slew_rate(double V_step,double R,double C);
double rc_rise_time_10_90(double R,double C);
double rc_fall_time_90_10(double R,double C);
double rc_delay_time_50(double R,double C);
double varactor_capacitance(double Cj0,double Vr,double Vbi,double m);
double saturable_inductance(double L0,double I_cur,double Isat);
double temp_dependent_resistance(double R0,double T0,double T,double alpha);
#endif
