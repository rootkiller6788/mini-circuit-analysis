#include "transient.h"
#include "ac_analysis.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef I
#define I _Complex_I
#endif

/* Time Constants L2 L4 */
double rc_time_constant(double R,double C){return R*C;}
double rl_time_constant(double R,double L){if(R<=0.0)return INFINITY;return L/R;}

/* RC Charging (Step Response) L6 */
double rc_charging_vc(double Vf,double v0,double tau,double t){
  if(tau<=0.0)return Vf;return Vf+(v0-Vf)*exp(-t/tau);}
double rc_charging_i(double Vf,double v0,double R,double tau,double t){
  if(R==0.0||tau<=0.0)return 0.0;return((Vf-v0)/R)*exp(-t/tau);}
double rc_discharging_vc(double v0,double tau,double t){
  if(tau<=0.0)return 0.0;return v0*exp(-t/tau);}
double rc_discharging_i(double v0,double R,double tau,double t){
  if(R==0.0||tau<=0.0)return 0.0;return-(v0/R)*exp(-t/tau);}

/* RL Charging L6 */
double rl_charging_il(double If,double i0,double tau,double t){
  if(tau<=0.0)return If;return If+(i0-If)*exp(-t/tau);}
double rl_charging_vl(double Vs,double i0,double R,double tau,double t){
  if(tau<=0.0)return 0.0;return(Vs-i0*R)*exp(-t/tau);}
double rl_discharging_il(double i0,double tau,double t){
  if(tau<=0.0)return 0.0;return i0*exp(-t/tau);}

/* Energy Storage L1 */
double capacitor_energy(double C,double v){return 0.5*C*v*v;}
double inductor_energy(double L,double i){return 0.5*L*i*i;}

/* RLC Series Parameters L6 */
double rlc_series_alpha(double R,double L){if(L<=0.0)return INFINITY;return R/(2.0*L);}
double rlc_series_omega0(double L,double C){if(L<=0.0||C<=0.0)return 0.0;return 1.0/sqrt(L*C);}
double rlc_series_zeta(double R,double L,double C){double a=rlc_series_alpha(R,L);double w0=rlc_series_omega0(L,C);if(w0<=0.0)return INFINITY;return a/w0;}
double rlc_series_omega_d(double R,double L,double C){double z=rlc_series_zeta(R,L,C);double w0=rlc_series_omega0(L,C);if(z>=1.0)return 0.0;return w0*sqrt(1.0-z*z);}
RLC_Regime_t rlc_series_regime(double R,double L,double C){
  double z=rlc_series_zeta(R,L,C);if(z>1.0)return RLC_OVERDAMPED;
  else if(fabs(z-1.0)<1e-9)return RLC_CRITICALLY_DAMPED;
  else if(z>0.0)return RLC_UNDERDAMPED;else return RLC_UNDAMPED;}

/* RLC Series Step Response L6 */
double rlc_series_step_vc(double Vs,double R,double L,double C,double v0,double i0,double t){
  if(t<0.0)return v0;
  double z=rlc_series_zeta(R,L,C);double w0=rlc_series_omega0(L,C);double a=rlc_series_alpha(R,L);
  double vc_ss=Vs;
  if(z>1.0){
    double s1=-a+sqrt(a*a-w0*w0);double s2=-a-sqrt(a*a-w0*w0);
    double A2=((s1*(v0-Vs))-i0/C)/(s1-s2);double A1=v0-Vs-A2;
    return vc_ss+A1*exp(s1*t)+A2*exp(s2*t);}
  else if(fabs(z-1.0)<1e-9){
    double A1=v0-Vs;double A2=i0/C+a*(v0-Vs);
    return vc_ss+(A1+A2*t)*exp(-a*t);}
  else{
    double wd=w0*sqrt(1.0-z*z);double A1=v0-Vs;double A2=(i0/C+a*(v0-Vs))/wd;
    return vc_ss+exp(-a*t)*(A1*cos(wd*t)+A2*sin(wd*t));}}
double rlc_series_step_il(double Vs,double R,double L,double C,double v0,double i0,double t){
  if(t<0.0)return i0;
  double z=rlc_series_zeta(R,L,C);double w0=rlc_series_omega0(L,C);double a=rlc_series_alpha(R,L);
  if(z>1.0){
    double s1=-a+sqrt(a*a-w0*w0);double s2=-a-sqrt(a*a-w0*w0);
    double A2=((v0-Vs)/L+s1*i0)/(s2-s1);double A1=i0-A2;
    return A1*exp(s1*t)+A2*exp(s2*t);}
  else if(fabs(z-1.0)<1e-9){
    double A1=i0;double A2=(v0-Vs)/L+a*i0;
    return(A1+A2*t)*exp(-a*t);}
  else{
    double wd=w0*sqrt(1.0-z*z);double A1=i0;double A2=((v0-Vs)/L+a*i0)/wd;
    return exp(-a*t)*(A1*cos(wd*t)+A2*sin(wd*t));}}

/* RLC Natural Response L6 */
double rlc_series_natural_vc(double R,double L,double C,double v0,double i0,double t){
  return rlc_series_step_vc(0.0,R,L,C,v0,i0,t);}
double rlc_series_natural_il(double R,double L,double C,double v0,double i0,double t){
  return rlc_series_step_il(0.0,R,L,C,v0,i0,t);}

/* Parallel RLC L6 */
double rlc_parallel_alpha(double R,double C){if(R<=0.0||C<=0.0)return INFINITY;return 1.0/(2.0*R*C);}
double rlc_parallel_zeta(double R,double L,double C){double a=rlc_parallel_alpha(R,C);double w0=1.0/sqrt(L*C);if(w0<=0.0)return INFINITY;return a/w0;}
double rlc_parallel_step_il(double Is,double R,double L,double C,double i0,double v0,double t){
  if(t<0.0)return i0;double a=rlc_parallel_alpha(R,C);double w0=1.0/sqrt(L*C);
  double z=rlc_parallel_zeta(R,L,C);double il_ss=Is;
  if(z>1.0){double s1=-a+sqrt(a*a-w0*w0);double s2=-a-sqrt(a*a-w0*w0);
    double A2=((s1*(i0-Is))+v0/L)/(s1-s2);double A1=i0-Is-A2;
    return il_ss+A1*exp(s1*t)+A2*exp(s2*t);}
  else if(fabs(z-1.0)<1e-9){double A1=i0-Is;double A2=v0/L+a*(i0-Is);
    return il_ss+(A1+A2*t)*exp(-a*t);}
  else{double wd=w0*sqrt(1.0-z*z);double A1=i0-Is;double A2=(v0/L+a*(i0-Is))/wd;
    return il_ss+exp(-a*t)*(A1*cos(wd*t)+A2*sin(wd*t));}}
double rlc_parallel_step_vc(double Is,double R,double L,double C,double i0,double v0,double t){
  if(t<0.0)return v0;double il=rlc_parallel_step_il(Is,R,L,C,i0,v0,t);
  double a=rlc_parallel_alpha(R,C);double w0=1.0/sqrt(L*C);
  double il_t=rlc_parallel_step_il(Is,R,L,C,i0,v0,t);
  double dil_dt=(t>0)?(il_t-rlc_parallel_step_il(Is,R,L,C,i0,v0,t*0.999))/ (t*0.001):0;
  return L*dil_dt;}

/* Numerical Integration L5 */
int transient_euler_step(TransientState_t *st,void (*f)(const TransientState_t*,double*,double),double t){
  if(!st||!f||st->n_states<=0)return -1;double *dx=(double*)malloc((size_t)st->n_states*sizeof(double));
  if(!dx)return -1;f(st,dx,t);
  for(int i=0;i<st->n_states;i++){if(i<st->n_states)st->v_caps[i]+=st->dt*dx[i];}
  st->t=t+st->dt;free(dx);return 0;}
int transient_rk4_step(TransientState_t *st,void (*f)(const TransientState_t*,double*,double),double t){
  if(!st||!f||st->n_states<=0)return -1;int ns=st->n_states;double dt=st->dt;
  double *k1=(double*)malloc((size_t)(4*ns)*sizeof(double));if(!k1)return -1;
  double *k2=k1+ns,*k3=k2+ns,*k4=k3+ns;
  f(st,k1,t);
  TransientState_t st2=*st;for(int i=0;i<ns;i++)st2.v_caps[i]=st->v_caps[i]+0.5*dt*k1[i];f(&st2,k2,t+0.5*dt);
  TransientState_t st3=*st;for(int i=0;i<ns;i++)st3.v_caps[i]=st->v_caps[i]+0.5*dt*k2[i];f(&st3,k3,t+0.5*dt);
  TransientState_t st4=*st;for(int i=0;i<ns;i++)st4.v_caps[i]=st->v_caps[i]+dt*k3[i];f(&st4,k4,t+dt);
  for(int i=0;i<ns;i++)st->v_caps[i]+=(dt/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  st->t=t+dt;free(k1);return 0;}
int transient_simulate(TransientState_t *st,void (*f)(const TransientState_t*,double*,double),double ts,double te,double dt,void (*cb)(double,const TransientState_t*,void*),void *ud){
  if(!st||!f||ts>=te||dt<=0)return -1;st->t=ts;st->dt=dt;int steps=0;
  while(st->t<te){if(cb)cb(st->t,st,ud);
    if(st->solver==SOLVER_RK4)transient_rk4_step(st,f,st->t);else transient_euler_step(st,f,st->t);steps++;
    if(steps>1000000)break;}return steps;}

/* Complete Response L6 Natural+Forced */
double rc_ac_complete_vc(double Vm,double w,double ph,double R,double C,double v0,double t){
  double tau=R*C;double ph_rad=ph*M_PI/180.0;
  double complex Z=R+1.0/(I*w*C);double complex Vc_ss=Vm*cpx_from_polar(1.0,0)/ (I*w*C*Z);
  double vc_forced=cabs(Vc_ss)*cos(w*t+carg(Vc_ss));
  double vc_natural=(v0-cabs(Vc_ss)*cos(carg(Vc_ss)))*exp(-t/tau);
  return vc_forced+vc_natural;}
double rl_ac_complete_il(double Vm,double w,double ph,double R,double L,double i0,double t){
  double tau=L/R;double ph_rad=ph*M_PI/180.0;
  double complex Z=R+I*w*L;double complex Il_ss=Vm*cpx_from_polar(1.0,ph_rad)/Z;
  double il_forced=cabs(Il_ss)*cos(w*t+carg(Il_ss));
  double il_natural=(i0-cabs(Il_ss)*cos(carg(Il_ss)))*exp(-t/tau);
  return il_forced+il_natural;}

/* Timer Applications L7 */
double timer_555_astable_freq(double R1,double R2,double C){if(R1<=0||R2<=0||C<=0)return 0;return 1.44/((R1+2*R2)*C);}
double timer_555_astable_duty(double R1,double R2){double denom=R1+2*R2;if(denom<=0)return 0;return(R1+R2)/denom;}
double timer_555_monostable_pulse(double R,double C){if(R<=0||C<=0)return 0;return 1.1*R*C;}
double rc_oscillator_freq_symmetric(double R,double C){if(R<=0||C<=0)return 0;return 1.0/(2.0*R*C*log(3.0));}

/* Nonlinear Models L8 */
double varactor_capacitance(double Cj0,double Vr,double Vbi,double m){
  if(Vbi<=0)return Cj0;double ratio=1.0+Vr/Vbi;if(ratio<=0)return Cj0;return Cj0/pow(ratio,m);}
double saturable_inductance(double L0,double I_cur,double Isat){if(Isat<=0)return L0;return L0/(1.0+(I_cur/Isat)*(I_cur/Isat));}
double temp_dependent_resistance(double R0,double T0,double T,double alpha){return R0*(1.0+alpha*(T-T0));}

/* Step Response for arbitrary RLC L6 */
double step_response_universal(double t,double tau,double order){
  if(t<0)return 0;if(order==1)return 1.0-exp(-t/tau);
  return 1.0-exp(-t/tau)*(1.0+t/tau);}
double impulse_response_rc(double t,double R,double C){if(t<0)return 0;if(R<=0||C<=0)return 0;return exp(-t/(R*C))/(R*C);}
double impulse_response_rl(double t,double R,double L){if(t<0)return 0;if(R<=0||L<=0)return 0;return(R/L)*exp(-(R/L)*t);}

/* Ramp Response L6 */
double rc_ramp_response_vc(double slope,double R,double C,double t,double v0){
  double tau=R*C;double vc_ss=slope*(t-tau);return vc_ss+(v0-vc_ss+slope*tau)*exp(-t/tau);}

/* Pulse Response L6 */
double rc_pulse_response_vc(double V_amplitude,double R,double C,double t,double t_pulse){
  double tau=R*C;if(t<0)return 0;
  if(t<=t_pulse)return V_amplitude*(1.0-exp(-t/tau));
  double v_at_pulse=V_amplitude*(1.0-exp(-t_pulse/tau));return v_at_pulse*exp(-(t-t_pulse)/tau);}

/* Sinusoidal Steady State in Time Domain L6 */
double rc_sinusoidal_response_vc(double Vm,double omega,double phi,double R,double C,double t){
  double tau=R*C;
  double complex Z=R+1.0/(I*omega*C);
  double complex Vc_phasor=Vm*(1.0/(I*omega*C))/Z;
  double mag=cabs(Vc_phasor);double ph=carg(Vc_phasor);
  return mag*cos(omega*t+phi+ph);}
double rl_sinusoidal_response_il(double Vm,double omega,double phi,double R,double L,double t){
  double tau=L/R;double complex Z=R+I*omega*L;double complex Il=Vm/Z;
  return cabs(Il)*cos(omega*t+phi+carg(Il));}

/* Laplace Transform Based Analysis L3 L8 */
double inverse_laplace_first_order(double pole,double t){if(t<0)return 0;return exp(-pole*t);}
double inverse_laplace_second_order_real(double p1,double p2,double t,double *out){
  if(t<0)return 0;double d=p2-p1;if(fabs(d)<1e-10)return t*exp(-p1*t);return(exp(-p1*t)-exp(-p2*t))/d;}
double inverse_laplace_second_order_complex(double alpha,double wd,double t){
  if(t<0)return 0;return exp(-alpha*t)*sin(wd*t)/wd;}

/* Numerical Transient (Euler) for RC circuit L5 */
void rc_transient_euler_simulate(double R,double C,double V0,double dt,double t_end,double *v_out){
  if(!v_out||dt<=0)return;int steps=(int)(t_end/dt);double v=V0;
  for(int i=0;i<=steps&&i<1000;i++){v_out[i]=v;double i_c=(v)/R;v-=i_c*dt/C;}}

/* Slew Rate and Rise Time L7 */
double rc_slew_rate(double V_step,double R,double C){if(R<=0||C<=0)return INFINITY;return V_step/(R*C);}
double rc_rise_time_10_90(double R,double C){return 2.2*R*C;}
double rc_fall_time_90_10(double R,double C){return 2.2*R*C;}
double rc_delay_time_50(double R,double C){return 0.693*R*C;}

/* Power-on Reset Timing L7 */
double por_timeout(double R,double C,double Vcc,double V_threshold){
  if(Vcc<=V_threshold||V_threshold<=0)return 0;return R*C*log(Vcc/(Vcc-V_threshold));}

/* State-Space Representation L8 */
typedef struct{int n;double *A;double *B;double *C;double *D;double *x;}StateSpace_t;
StateSpace_t* ss_alloc(int n,int m,int p){StateSpace_t*ss=malloc(sizeof(StateSpace_t));
  ss->n=n;ss->A=calloc(n*n,sizeof(double));ss->B=calloc(n*m,sizeof(double));
  ss->C=calloc(p*n,sizeof(double));ss->D=calloc(p*m,sizeof(double));ss->x=calloc(n,sizeof(double));return ss;}
void ss_free(StateSpace_t*ss){if(!ss)return;free(ss->A);free(ss->B);free(ss->C);free(ss->D);free(ss->x);free(ss);}
void ss_step_euler(StateSpace_t*ss,double*u,double dt){
  if(!ss)return;double*x=ss->x;int n=ss->n;
  double *dx=calloc(n,sizeof(double));for(int i=0;i<n;i++){dx[i]=0;for(int j=0;j<n;j++)dx[i]+=ss->A[i*n+j]*x[j];for(int j=0;j<1;j++)dx[i]+=ss->B[i]*u[j];}
  for(int i=0;i<n;i++)x[i]+=dx[i]*dt;free(dx);}

/* Transmission Line Transient L8 */
double tl_reflection_coefficient(double ZL,double Z0){if(ZL+Z0==0)return 1.0;return(ZL-Z0)/(ZL+Z0);}
double tl_transmission_coefficient(double ZL,double Z0){return 1.0+tl_reflection_coefficient(ZL,Z0);}
double tl_voltage_at_load(double V_incident,double ZL,double Z0){return V_incident*tl_transmission_coefficient(ZL,Z0);}
double tl_propagation_delay(double length_m,double velocity_factor){return length_m/(3e8*velocity_factor);}
double tl_ringing_period(double length_m,double velocity_factor){return 4.0*tl_propagation_delay(length_m,velocity_factor);}
double tl_termination_resistor_parallel(double Z0,double C_load,double trise){
  if(trise<=0)return Z0;return Z0;}

/* Snubber Circuit Design (L7, L8) */
double snubber_rc_resistor(double V_peak,double I_peak,double safety_factor){return V_peak/(I_peak*safety_factor);}
double snubber_rc_capacitor(double I_peak,double trise,double V_peak){if(V_peak<=0||trise<=0)return 0;return I_peak*trise/V_peak;}
double snubber_power_dissipation(double C,double V_peak,double f_sw){return 0.5*C*V_peak*V_peak*f_sw;}

/* Gate Driver Transient L7 */
double gate_drive_current(double Qg_total,double t_sw){if(t_sw<=0)return INFINITY;return Qg_total/t_sw;}
double gate_resistor_min(double V_drive,double I_peak_gate){if(I_peak_gate<=0)return INFINITY;return V_drive/I_peak_gate;}
double miller_plateau_charge(double Qgd,double V_plateau,double V_drive,double R_gate){
  if(R_gate<=0)return 0;return Qgd/((V_drive-V_plateau)/R_gate);}

/* Zener Diode Transient L7 */
double zener_regulation_voltage_variation(double I_min,double I_max,double Rz,double Vz_nom){
  if(Vz_nom<=0)return 0;return Rz*(I_max-I_min)/Vz_nom*100.0;}
double zener_power_dissipation(double Vz,double Iz){return Vz*Iz;}
double zener_series_resistor(double Vin_min,double Vin_max,double Vz,double Iz_min,double Iz_max){
  double Rs_max=(Vin_min-Vz)/(Iz_min+1e-6);double Rs_min=(Vin_max-Vz)/(Iz_max);
  if(Rs_min>Rs_max)return -1;return(Rs_min+Rs_max)/2.0;}

/* Inrush Current L7 */
double inrush_current_peak(double V_in,double R_series,double C_bulk){if(R_series<=0)return INFINITY;return V_in/R_series;}
double inrush_energy_joules(double C,double V){return 0.5*C*V*V;}
double ntc_thermistor_cold_resistance(double R25,double B,double T_cold){return R25*exp(B*(1.0/(T_cold+273.15)-1.0/298.15));}

/* Bootstrap Capacitor L7 */
double bootstrap_capacitor_min(double Qg_total,double dV_allowed,double Q_leakage,double f_sw){
  if(dV_allowed<=0)return 0;return(Qg_total+Q_leakage/f_sw)/dV_allowed;}

/* Phase-Locked Loop Transient L8 */
double pll_natural_frequency(double Kpd,double Kvco,double N_div,double R1,double C1){
  if(N_div<=0||C1<=0)return 0;return sqrt(Kpd*Kvco/(N_div*C1));}
double pll_damping_factor(double wn,double R1,double C1){return 0.5*wn*R1*C1;}
double pll_lock_time_estimate(double wn,double zeta,double freq_step,double tol_hz){
  if(wn<=0||tol_hz<=0)return INFINITY;return -log(tol_hz/freq_step)/(zeta*wn);}
double pll_loop_bandwidth(double wn,double zeta){return wn*sqrt(1+2*zeta*zeta+sqrt(4*zeta*zeta*zeta*zeta+4*zeta*zeta+2));}

/* Electrostatic Discharge (ESD) Models L7 */
double hbm_esd_current_peak(double V_esd,double R_hbm){if(R_hbm<=0)return 0;return V_esd/R_hbm;}
double hbm_esd_energy(double C_hbm,double V_esd){return 0.5*C_hbm*V_esd*V_esd;}
double cdm_esd_current_peak(double C_device,double V_esd,double R_discharge){if(R_discharge<=0)return 0;return V_esd/R_discharge;}

/* Soft Start Circuit L7 */
double soft_start_ramp_time(double C_ss,double I_charge,double V_ref){if(I_charge<=0)return INFINITY;return C_ss*V_ref/I_charge;}
double soft_start_inrush_reduction(double I_inrush,double I_rated){return(1.0-I_rated/I_inrush)*100.0;}

/* Load Step Response L7 */
double load_step_voltage_droop(double delta_I,double C_out,double t_response,double L_parasitic,double trise){
  double dv_cap=delta_I*t_response/C_out;double dv_esl=L_parasitic*delta_I/trise;return dv_cap+dv_esl;}
double load_step_settling_time(double L,double C,double zeta){double wn=1/sqrt(L*C);return 4.0/(zeta*wn);}
double load_transient_recovery_cycles(double t_settle,double f_sw){return t_settle*f_sw;}

/* Ringing and Overshoot L6 */
double overshoot_percent(double V_peak,double V_steady){if(V_steady<=0)return 0;return(V_peak/V_steady-1.0)*100.0;}
double overshoot_from_zeta(double zeta){if(zeta>=1)return 0;return 100.0*exp(-M_PI*zeta/sqrt(1-zeta*zeta));}
double ringing_frequency(double fn,double zeta){if(zeta>=1)return 0;return fn*sqrt(1-zeta*zeta);}
double ringing_decay_envelope(double t,double zeta,double wn){return exp(-zeta*wn*t);}
double Q_from_overshoot(double overshoot_pct){if(overshoot_pct<=0)return INFINITY;double d=log(100.0/overshoot_pct)/M_PI;return sqrt(d*d+1)/(2*d);}

/* Sequenced Power Supply L7 */
double power_sequencing_delay(double C_delay,double R_delay,double V_threshold,double V_supply){
  return -R_delay*C_delay*log(1.0-V_threshold/V_supply);}
double power_supply_tracking_ratio(double V_tracked,double V_master){if(V_master<=0)return 0;return V_tracked/V_master;}

/* Watchdog Timer L7 */
double watchdog_timeout_rc(double R,double C,double V_threshold,double Vcc){return -R*C*log(1.0-V_threshold/Vcc);}
int watchdog_ok_delay_ms(double R,double C){return(int)(1.1*R*C*1000);}

/* Debounce Circuit L7 */
double debounce_rc_time_constant(double R_pullup,double C_filter){return R_pullup*C_filter;}
double debounce_settling_time(double R,double C,double V_cc,double V_IH_min){return -R*C*log(1.0-V_IH_min/V_cc);}
int debounce_samples_needed(double bounce_time_ms,double sample_period_ms){if(sample_period_ms<=0)return 1;return(int)(bounce_time_ms/sample_period_ms)+1;}

/* Motor Driver Transient L7 */
double motor_startup_current(double V_supply,double R_winding,double R_ds_on){return V_supply/(R_winding+R_ds_on);}
double motor_back_emf(double Ke_V_per_rpm,double rpm){return Ke_V_per_rpm*rpm;}
double motor_current_ripple(double V_supply,double L_winding,double f_pwm,double duty){if(L_winding<=0||f_pwm<=0)return 0;return V_supply*duty*(1-duty)/(L_winding*f_pwm);}
double motor_mechanical_time_constant(double J,double R,double Kt,double Ke){double D=Kt*Ke/R;if(D<=0)return INFINITY;return J/D;}
double motor_electrical_time_constant(double L,double R){if(R<=0)return INFINITY;return L/R;}
double h_bridge_shoot_through_dead_time_ns(double tr,double tf,double safety_margin){return(tr+tf)*safety_margin;}

/* Relay Coil Transient L7 */
double relay_coil_time_constant(double L_coil,double R_coil){return L_coil/R_coil;}
double relay_pull_in_time_ms(double V_drive,double V_pull_in,double tau){return tau*log(V_drive/(V_drive-V_pull_in))*1000;}
double relay_drop_out_time_ms(double V_hold,double V_drop_out,double tau){return tau*log(V_hold/V_drop_out)*1000;}
double relay_flyback_diode_peak_current(double I_coil){return I_coil;}
double relay_snubber_voltage_peak(double I_coil,double L_coil,double C_snubber){if(C_snubber<=0)return INFINITY;return I_coil*sqrt(L_coil/C_snubber);}

/* Solenoid Driver L7 */
double solenoid_peak_current(double V_drive,double R_coil){if(R_coil<=0)return 0;return V_drive/R_coil;}
double solenoid_hold_current(double V_hold,double R_coil){if(R_coil<=0)return 0;return V_hold/R_coil;}
double solenoid_pwm_hold_voltage(double V_drive,double I_hold,double I_peak){return V_drive*I_hold/I_peak;}
double solenoid_energy_joules(double L,double I_cur){return 0.5*L*I_cur*I_cur;}

/* Piezo Driver L7 */
double piezo_capacitance(double epsilon,double area_m2,double thickness_m){if(thickness_m<=0)return 0;return epsilon*area_m2/thickness_m;}
double piezo_resonant_frequency(double L_mech,double C_mech){return 1.0/(2*M_PI*sqrt(L_mech*C_mech));}
double piezo_driver_voltage_swing(double V_pp,double C_piezo,double f){return V_pp*C_piezo*f;}
double piezo_driver_power(double V_pp,double C_piezo,double f){return C_piezo*V_pp*V_pp*f;}
