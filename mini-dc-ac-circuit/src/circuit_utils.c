#include "circuit_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Unit Conversions L1 */
double deg_to_rad(double deg){return deg*M_PI/180.0;}
double rad_to_deg(double rad){return rad*180.0/M_PI;}
double hz_to_radps(double hz){return 2.0*M_PI*hz;}
double radps_to_hz(double radps){return radps/(2.0*M_PI);}
double celsius_to_kelvin(double c){return c+273.15;}
double db_to_linear(double db){return pow(10.0,db/10.0);}
double linear_to_db(double lin){if(lin<=0)return -INFINITY;return 10*log10(lin);}
double dbm_to_watts(double dbm){return pow(10.0,(dbm-30.0)/10.0);}
double watts_to_dbm(double watts){if(watts<=0)return -INFINITY;return 10*log10(watts)+30.0;}
double dbm_to_volts_rms(double dbm,double Z0){if(Z0<=0)return 0;return sqrt(pow(10.0,(dbm-30.0)/10.0)*Z0);}
double volts_rms_to_dbm(double Vrms,double Z0){if(Vrms<=0||Z0<=0)return -INFINITY;return 10*log10(Vrms*Vrms/Z0)+30.0;}

/* SI Prefix L1 */
double pico_to_base(double p){return p*1e-12;}
double nano_to_base(double n){return n*1e-9;}
double micro_to_base(double u){return u*1e-6;}
double milli_to_base(double m){return m*1e-3;}
double kilo_to_base(double k){return k*1e3;}
double mega_to_base(double M){return M*1e6;}
double giga_to_base(double G){return G*1e9;}
double base_to_milli(double v){return v*1e3;}
double base_to_micro(double v){return v*1e6;}
double base_to_nano(double v){return v*1e9;}
double base_to_pico(double v){return v*1e12;}

/* Standard Resistor Values E12 E24 L1 */
static const double e12_vals[12]={1.0,1.2,1.5,1.8,2.2,2.7,3.3,3.9,4.7,5.6,6.8,8.2};
static const double e24_vals[24]={1.0,1.1,1.2,1.3,1.5,1.6,1.8,2.0,2.2,2.4,2.7,3.0,3.3,3.6,3.9,4.3,4.7,5.1,5.6,6.2,6.8,7.5,8.2,9.1};
static double nearest_series(double R,const double *vals,int n){
  if(R<=0||n<=0)return 0;double decade=pow(10.0,floor(log10(R)));double mantissa=R/decade;
  if(mantissa<1.0){mantissa*=10;decade/=10;}if(mantissa>10.0){mantissa/=10;decade*=10;}
  double best=vals[0];double best_diff=fabs(mantissa-best);
  for(int i=1;i<n;i++){double diff=fabs(mantissa-vals[i]);if(diff<best_diff){best_diff=diff;best=vals[i];}}
  if(best>=10.0)return best*decade/10.0;return best*decade;}
double nearest_e12_value(double R){return nearest_series(R,e12_vals,12);}
double nearest_e24_value(double R){return nearest_series(R,e24_vals,24);}
static int is_in_series(double R,const double *vals,int n,double tol){
  double near=nearest_series(R,vals,n);if(near<=0)return 0;return fabs(R-near)/near<=tol;}
int is_standard_e12(double R,double tol){return is_in_series(R,e12_vals,12,tol);}
int is_standard_e24(double R,double tol){return is_in_series(R,e24_vals,24,tol);}

/* Rounding L1 */
double round_to_sigdig(double val,int n){if(n<=0)return val;if(val==0)return 0;
  double scale=pow(10.0,ceil(log10(fabs(val)))-n);return round(val/scale)*scale;}
double round_to_decimals(double val,int n){double scale=pow(10.0,n);return round(val*scale)/scale;}

/* String Utilities */
const char* unit_to_string(Unit_t unit){switch(unit){case UNIT_VOLT:return"V";case UNIT_AMPERE:return"A";case UNIT_OHM:return"Ohm";case UNIT_FARAD:return"F";case UNIT_HENRY:return"H";case UNIT_WATT:return"W";case UNIT_VAR:return"var";case UNIT_VA:return"VA";case UNIT_RAD_PER_S:return"rad/s";case UNIT_HERTZ:return"Hz";case UNIT_SECOND:return"s";case UNIT_DEGREE:return"deg";case UNIT_DB:return"dB";case UNIT_NEPER:return"Np";default:return"";}}
const char* element_type_to_string(ElementType_t t){switch(t){case ELEM_RESISTOR:return"Resistor";case ELEM_CAPACITOR:return"Capacitor";case ELEM_INDUCTOR:return"Inductor";case ELEM_DC_VOLTAGE_SRC:return"DC_Vsrc";case ELEM_DC_CURRENT_SRC:return"DC_Isrc";case ELEM_AC_VOLTAGE_SRC:return"AC_Vsrc";case ELEM_AC_CURRENT_SRC:return"AC_Isrc";case ELEM_VCVS:return"VCVS";case ELEM_VCCS:return"VCCS";case ELEM_CCVS:return"CCVS";case ELEM_CCCS:return"CCCS";case ELEM_OPAMP_IDEAL:return"OpAmp";case ELEM_TRANSFORMER_IDEAL:return"Transformer";case ELEM_GYRATOR:return"Gyrator";default:return"Unknown";}}

/* Print Utilities */
void print_impedance(Impedance_t Z,FILE *fp){if(!fp)fp=stdout;fprintf(fp,"Z = %.4f + j%.4f Ohm",Z.re,Z.im);}
void print_complex_power(ComplexPower_t cp,FILE *fp){
  if(!fp)fp=stdout;fprintf(fp,"S = %.4f VA P = %.4f W Q = %.4f var pf = %.4f %s",cp.S,cp.P,cp.Q,cp.pf,cp.leading?"leading":"lagging");}
void print_resonance(Resonance_t r,FILE *fp){
  if(!fp)fp=stdout;fprintf(fp,"f0=%.4f Hz Q=%.4f BW=%.4f Hz zeta=%.4f",r.f0_hz,r.Q,r.BW_hz,r.zeta);}

/* Error Codes */
const char* error_string(ErrorCode_t e){switch(e){
  case ERR_NONE:return"No error";case ERR_NULL_POINTER:return"Null pointer";
  case ERR_DIVISION_BY_ZERO:return"Division by zero";case ERR_SINGULAR_MATRIX:return"Singular matrix";
  case ERR_NO_CONVERGENCE:return"No convergence";case ERR_INVALID_ARGUMENT:return"Invalid argument";
  case ERR_OUT_OF_MEMORY:return"Out of memory";case ERR_CIRCUIT_NOT_CONNECTED:return"Circuit not connected";
  case ERR_FLOATING_NODE:return"Floating node";case ERR_VOLTAGE_LOOP:return"Voltage loop";
  case ERR_CURRENT_CUTSET:return"Current cutset";case ERR_MAX_ELEMENTS:return"Max elements exceeded";
  case ERR_MAX_NODES:return"Max nodes exceeded";case ERR_UNKNOWN_ELEMENT:return"Unknown element";
  default:return"Unknown error";}}

/* Random Numbers for Monte Carlo L8 */
static unsigned long _rand_seed=123456789;
double random_uniform(double a,double b){_rand_seed=(_rand_seed*1103515245+12345)&0x7fffffff;return a+(b-a)*(_rand_seed/2147483647.0);}
double random_normal(double mean,double stddev){
  double u1=random_uniform(0,1);double u2=random_uniform(0,1);
  double z=sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);return mean+stddev*z;}
double resistor_with_tolerance(double nominal,double tol_pct){
  double stddev=nominal*tol_pct/(3.0*100.0);return random_normal(nominal,stddev);}

/* Double Comparison Utilities */
int double_near(double a,double b,double rel_tol){if(a==b)return 1;double diff=fabs(a-b);double scale=fmax(fabs(a),fabs(b));if(scale==0)return 1;return diff/scale<=rel_tol;}
int double_near_abs(double a,double b,double abs_tol){return fabs(a-b)<=abs_tol;}
int is_valid_double(double x){return !isnan(x)&&!isinf(x);}

/* CRC-16 for Netlist Verification */
unsigned short crc16_netlist(const Circuit_t *ckt){
  if(!ckt)return 0;unsigned short crc=0xFFFF;unsigned char *p=(unsigned char*)ckt;
  for(size_t i=0;i<sizeof(*ckt);i++){crc^=p[i];for(int j=0;j<8;j++){if(crc&1)crc=(crc>>1)^0xA001;else crc>>=1;}}return crc;}
size_t safe_strcpy(char *dst,const char *src,size_t dst_size){
  if(!dst||!src||dst_size==0)return 0;size_t i=0;while(i<dst_size-1&&src[i]){dst[i]=src[i];i++;}dst[i]=0;return i;}
double parallel_equiv(double a,double b){if(a+b==0)return 0;return(a*b)/(a+b);}
double parallel_equiv_3(double a,double b,double c){double ab=parallel_equiv(a,b);return parallel_equiv(ab,c);}

/* Advanced Unit Conversions L1 */
double celsius_to_fahrenheit(double c){return c*9.0/5.0+32.0;}
double fahrenheit_to_celsius(double f){return(f-32.0)*5.0/9.0;}
double kelvin_to_celsius(double k){return k-273.15;}
double inches_to_meters(double inches){return inches*0.0254;}
double mils_to_meters(double mils){return mils*2.54e-5;}
double awg_to_diameter_mm(int awg){if(awg<0||awg>40)return 0;return 0.127*pow(92.0,(36.0-awg)/39.0);}
double awg_to_resistance_per_meter(int awg){double d=awg_to_diameter_mm(awg);if(d<=0)return 0;double area=M_PI*d*d/4.0*1e-6;return 1.68e-8/area;}
int nearest_awg_from_diameter_mm(double dia_mm){
  if(dia_mm<=0)return -1;int best=0;double best_diff=1e9;
  for(int awg=0;awg<=40;awg++){double d=awg_to_diameter_mm(awg);double diff=fabs(d-dia_mm);if(diff<best_diff){best_diff=diff;best=awg;}}return best;}

/* Thermal Calculations L7 */
double resistor_temp_rise(double power_w,double thermal_resistance_cw){return power_w*thermal_resistance_cw;}
double resistor_derating(double rated_power,double amb_temp,double max_temp,double derate_factor){
  if(amb_temp>=max_temp)return 0;return rated_power*(1.0-derate_factor*(amb_temp-25.0)/(max_temp-25.0));}
double junction_temperature(double p_diss,double rth_ja,double t_amb){return t_amb+p_diss*rth_ja;}
double heatsink_thermal_resistance(double tj_max,double t_amb,double p_max,double rth_jc,double rth_cs){
  double rth_total=(tj_max-t_amb)/p_max;double rth_sa=rth_total-rth_jc-rth_cs;if(rth_sa<0)return 0;return rth_sa;}

/* PCB Trace Calculations L7 */
double trace_resistance(double resistivity,double length_m,double width_m,double thickness_m){
  if(width_m<=0||thickness_m<=0)return INFINITY;return resistivity*length_m/(width_m*thickness_m);}
double trace_inductance_uh(double length_mm,double width_mm){
  if(width_mm<=0)return 0;return 0.0002*length_mm*(log(2.0*length_mm/width_mm)+0.5+0.2235*width_mm/length_mm);}
double trace_capacitance_pf(double length_mm,double width_mm,double height_mm,double er){
  if(height_mm<=0)return 0;return 0.00885*er*length_mm*width_mm/height_mm;}
double trace_characteristic_impedance(double L_per_m,double C_per_m){if(C_per_m<=0)return INFINITY;return sqrt(L_per_m/C_per_m);}
double trace_max_current_a(double width_mm,double thickness_oz,double temp_rise_c){
  double area_mil2=width_mm/0.0254*(thickness_oz*1.378);double k=0.048;return k*pow(temp_rise_c,0.44)*pow(area_mil2,0.725);}

/* Capacitor related L7 */
double capacitor_reactance(double C,double freq){if(C<=0||freq<=0)return INFINITY;return 1.0/(2*M_PI*freq*C);}
double inductor_reactance(double L,double freq){if(L<=0)return 0;return 2*M_PI*freq*L;}
double capacitor_impedance_magnitude(double C,double freq){return capacitor_reactance(C,freq);}
double capacitor_esr_power_loss(double I_rms,double esr){return I_rms*I_rms*esr;}
double capacitor_ripple_voltage(double I_load,double freq,double C,double esr){
  if(C<=0||freq<=0)return INFINITY;double dv_cap=I_load/(freq*C);double dv_esr=I_load*esr;return dv_cap+dv_esr;}
double capacitor_discharge_time(double C,double V0,double V_final,double R_load){
  if(C<=0||R_load<=0||V0<=V_final||V_final<=0)return 0;return R_load*C*log(V0/V_final);}

/* Battery Related L7 */
double battery_life_hours(double capacity_ah,double current_a){if(current_a<=0)return INFINITY;return capacity_ah/current_a;}
double battery_voltage_sag(double Voc,double I_load,double R_internal){return Voc-I_load*R_internal;}
double peukert_capacity(double rated_ah,double rated_hours,double actual_current,double peukert_const){
  double rated_current=rated_ah/rated_hours;double t=rated_hours*pow(rated_current/actual_current,peukert_const);return actual_current*t;}

/* Decibel Calculations L7 */
double voltage_gain_db(double vout,double vin){if(vin<=0)return 0;return 20*log10(vout/vin);}
double current_gain_db(double iout,double iin){if(iin<=0)return 0;return 20*log10(iout/iin);}
double power_gain_db(double pout,double pin){if(pin<=0)return -INFINITY;return 10*log10(pout/pin);}
double db_to_voltage_ratio(double db){return pow(10.0,db/20.0);}
double db_to_power_ratio(double db){return pow(10.0,db/10.0);}

/* Impedance in series/parallel L3 */
Impedance_t impedance_add_three(Impedance_t Z1,Impedance_t Z2,Impedance_t Z3){
  Impedance_t Z;Z.re=Z1.re+Z2.re+Z3.re;Z.im=Z1.im+Z2.im+Z3.im;return Z;}
Impedance_t impedance_parallel_three(Impedance_t Z1,Impedance_t Z2,Impedance_t Z3){
  return impedance_parallel(impedance_parallel(Z1,Z2),Z3);}
double impedance_magnitude(Impedance_t Z){return sqrt(Z.re*Z.re+Z.im*Z.im);}
double impedance_phase_deg(Impedance_t Z){return atan2(Z.im,Z.re)*180.0/M_PI;}
Impedance_t resistance_to_impedance(double R){Impedance_t Z;Z.re=R;Z.im=0;return Z;}

/* Bridge Circuits L6 */
double maxwell_bridge_inductance(double R1,double R2,double R3,double C1){return R1*R2*C1;}
double maxwell_bridge_resistance(double R1,double R2,double R3){return R1*R2/R3;}
double hay_bridge_inductance(double R1,double R2,double R3,double C1,double w){double w2=w*w;return R1*R2*C1/(1+w2*R1*R1*C1*C1);}
double schering_bridge_capacitance(double R3,double R4,double C2){return C2*R4/R3;}
double schering_bridge_dissipation(double R3,double R4,double C2,double Cx){return 2*M_PI*1000*R4*C2;}

/* Three-Phase Power L7 */
double three_phase_active_power(double V_ll_rms,double I_line_rms,double pf){return sqrt(3.0)*V_ll_rms*I_line_rms*pf;}
double three_phase_reactive_power(double V_ll_rms,double I_line_rms,double pf){return sqrt(3.0)*V_ll_rms*I_line_rms*sin(acos(pf));}
double three_phase_apparent_power(double V_ll_rms,double I_line_rms){return sqrt(3.0)*V_ll_rms*I_line_rms;}
double line_to_phase_voltage_wye(double V_ll){return V_ll/sqrt(3.0);}
double phase_to_line_voltage_wye(double V_ph){return V_ph*sqrt(3.0);}
double line_current_delta(double I_phase){return I_phase*sqrt(3.0);}
double phase_current_delta(double I_line){return I_line/sqrt(3.0);}

/* Noise Calculations L7 */
double thermal_noise_voltage_rms(double R,double T_kelvin,double BW_hz){return sqrt(4.0*1.380649e-23*T_kelvin*R*BW_hz);}
double thermal_noise_power_dbm(double T_kelvin,double BW_hz){return 10*log10(1.380649e-23*T_kelvin*BW_hz*1000);}
double shot_noise_current_rms(double I_dc,double BW_hz){return sqrt(2.0*1.602176634e-19*I_dc*BW_hz);}
double flicker_noise_corner_freq(double Kf,double f_L,double f_H){return exp((log(f_H)+log(f_L))/2);}
double total_noise_voltage_rms(const double *noise_sources,int n){double sum=0;for(int i=0;i<n;i++)sum+=noise_sources[i]*noise_sources[i];return sqrt(sum);}

/* Bode Plot Helpers L6 */
double bode_asymptote_mag_db(double f,double f_corner,int order){if(f<=0||f_corner<=0)return 0;if(f<f_corner)return 0;return order*20*log10(f/f_corner);}
double bode_asymptote_phase_deg(double f,double f_corner,int order){if(f<=0||f_corner<=0)return 0;double decade=f/f_corner;if(decade<0.1)return 0;if(decade>10)return order*90.0;return order*45.0*(log10(decade)+1);}

/* Component Value Codes L1 */
int resistor_color_code_value(const int *bands,int n_bands){
  static const int digit_vals[10]={0,1,2,3,4,5,6,7,8,9};
  if(n_bands==4)return(bands[0]*10+bands[1])*(int)pow(10,bands[2]);
  if(n_bands==5)return(bands[0]*100+bands[1]*10+bands[2])*(int)pow(10,bands[3]);
  return(bands[0]*10+bands[1])*(int)pow(10,bands[2]);}
int resistor_color_code_tolerance(int band){static const int tol[11]={0,1,2,0,0,0,5,0,0,0,0};return(band>=0&&band<11)?tol[band]:20;}
const char* resistor_color_name(int digit){
  static const char* names[10]={"Black","Brown","Red","Orange","Yellow","Green","Blue","Violet","Grey","White"};
  return(digit>=0&&digit<10)?names[digit]:"Unknown";}

/* SMD Resistor Code Parser L1 */
double smd_resistor_code_3digit(int code){int val=code/10;int exp=code%10;return val*pow(10.0,exp);}
double smd_resistor_code_4digit(int code){int val=code/10;int exp=code%10;return val*pow(10.0,exp);}
double smd_resistor_code_eia96(int code){static const double e96[96]={1.00};return 1.0;}

/* Capacitor Code Parser L1 */
double capacitor_code_3digit_pf(int code){int val=code/10;int exp=code%10;return val*pow(10.0,exp);}
double capacitor_code_4digit_pf(int code){int val=code/10;int exp=code%10;return val*pow(10.0,exp);}

/* Power Supply Analysis L7 */
double transformer_va_rating(double V_secondary,double I_secondary){return V_secondary*I_secondary;}
double rectifier_peak_voltage(double V_rms,double Vf_diode){return V_rms*sqrt(2.0)-2*Vf_diode;}
double smoothing_capacitor_min(double I_load,double V_ripple,double f_line){if(f_line<=0||V_ripple<=0)return 0;return I_load/(2*f_line*V_ripple);}
double linear_regulator_efficiency(double V_in,double V_out,double I_load,double I_q){return(V_out*I_load)/(V_in*(I_load+I_q))*100.0;}
double linear_regulator_power_dissipation(double V_in,double V_out,double I_load){return(V_in-V_out)*I_load;}
double switching_regulator_duty_buck(double V_out,double V_in,double efficiency){if(V_in<=0)return 0;return V_out/(efficiency*V_in);}
double switching_regulator_duty_boost(double V_out,double V_in){if(V_out<=V_in)return 0;return 1.0-V_in/V_out;}
double buck_inductor_ripple_current(double V_in,double V_out,double L,double f_sw,double duty){
  if(L<=0||f_sw<=0)return 0;return(V_in-V_out)*duty/(L*f_sw);}
double buck_output_ripple_voltage(double delta_IL,double f_sw,double C_out,double ESR){
  if(f_sw<=0||C_out<=0)return 0;return delta_IL*(ESR+1.0/(8*f_sw*C_out));}

/* PWM and Duty Cycle L7 */
double pwm_duty_cycle_percent(double t_on,double t_period){if(t_period<=0)return 0;return t_on/t_period*100.0;}
double pwm_average_voltage(double V_high,double duty_pct){return V_high*duty_pct/100.0;}
double pwm_frequency(double t_period){if(t_period<=0)return 0;return 1.0/t_period;}
double pwm_resolution_bits(double f_clk,double f_pwm){if(f_pwm<=0)return 0;return log2(f_clk/f_pwm);}
double servo_pulse_width_ms(double angle_deg,double min_ms,double max_ms,double max_deg){
  if(max_deg<=0)return min_ms;return min_ms+(max_ms-min_ms)*angle_deg/max_deg;}

/* ADC/DAC Quantization L7 */
double adc_resolution_volts(double V_ref,int bits){if(bits<=0)return 0;return V_ref/pow(2.0,bits);}
double adc_quantization_noise_rms(double V_lsb){return V_lsb/sqrt(12.0);}
double adc_snr_db(int bits){return 6.02*bits+1.76;}
double adc_enob(double snr_db){return(snr_db-1.76)/6.02;}
int adc_voltage_to_code(double V_in,double V_ref,int bits){if(V_ref<=0)return 0;if(V_in>=V_ref)return(1<<bits)-1;if(V_in<=0)return 0;return(int)(V_in/V_ref*(1<<bits));}
double adc_code_to_voltage(int code,double V_ref,int bits){if(bits<=0)return 0;return code*V_ref/(1<<bits);}

/* Pull-up/Pull-down Resistor Selection L7 */
double pullup_resistor_max(double Vcc,double V_IH_min,double I_leakage){if(I_leakage<=0)return 1e6;return(Vcc-V_IH_min)/I_leakage;}
double pullup_resistor_min(double Vcc,double V_OL_max,double I_OL_max){if(I_OL_max<=0)return 100;return(Vcc-V_OL_max)/I_OL_max;}
double pulldown_resistor_max(double V_IL_max,double I_leakage){if(I_leakage<=0)return 1e6;return V_IL_max/I_leakage;}

/* LED Current Limiting L7 */
double led_current_limit_resistor(double V_supply,double V_forward,double I_forward){if(I_forward<=0)return INFINITY;return(V_supply-V_forward)/I_forward;}
double led_power_dissipation(double Vf,double If){return Vf*If;}
double led_efficiency_lumens_per_watt(double lumens,double P_watts){if(P_watts<=0)return 0;return lumens/P_watts;}
double led_series_resistor_7segment(double V_supply,double Vf,double If,int n_segments){return led_current_limit_resistor(V_supply,Vf,If);}

/* Operational Amplifier Circuits L7 */
double opamp_inverting_gain(double Rf,double Rin){if(Rin<=0)return -INFINITY;return -Rf/Rin;}
double opamp_non_inverting_gain(double Rf,double Rg){if(Rg<=0)return 0;return 1.0+Rf/Rg;}
double opamp_differential_gain(double Rf,double R1,double R2,double Rg){return Rf/R1;}
double opamp_integrator_fc(double R,double C){return 1.0/(2*M_PI*R*C);}
double opamp_differentiator_fc(double R,double C){return 1.0/(2*M_PI*R*C);}
double opamp_sallen_key_lp_fc(double R1,double R2,double C1,double C2){return 1.0/(2*M_PI*sqrt(R1*R2*C1*C2));}
double opamp_sallen_key_lp_Q(double R1,double R2,double C1,double C2){return sqrt(R1*R2*C1*C2)/(C2*(R1+R2));}
double opamp_offset_voltage_output(double V_os,double gain){return V_os*gain;}
double opamp_input_bias_current_error(double I_bias,double R_feedback){return I_bias*R_feedback;}
double opamp_cmrr_error_db(double CMRR_db,double V_cm,double V_diff){return CMRR_db+20*log10(V_cm/V_diff);}
double opamp_slew_rate_limit_freq(double SR_V_us,double V_peak){if(V_peak<=0)return INFINITY;return SR_V_us/(2e6*M_PI*V_peak);}
double opamp_gain_bandwidth_limit_freq(double GBW_MHz,double gain){if(gain<=0)return 0;return GBW_MHz*1e6/gain;}
double opamp_total_output_noise_uv_rms(const double *noise_sources,int n){double sum=0;for(int i=0;i<n;i++)sum+=noise_sources[i]*noise_sources[i];return sqrt(sum);}

/* Comparator Circuits L7 */
double comparator_hysteresis_voltage(double Vcc,double R1,double R2){return Vcc*R2/(R1+R2);}
double comparator_threshold_high(double Vref,double Vhyst){return Vref+Vhyst/2.0;}
double comparator_threshold_low(double Vref,double Vhyst){return Vref-Vhyst/2.0;}
double window_comparator_R_values(double V_ref_high,double V_ref_low,double Vcc,double I_bias){
  if(I_bias<=0)return 0;return(V_ref_high-V_ref_low)/I_bias;}

/* Current Sensing L7 */
double shunt_resistor_value(double V_full_scale,double I_max){if(I_max<=0)return 0;return V_full_scale/I_max;}
double shunt_power_dissipation(double I_max,double R_shunt){return I_max*I_max*R_shunt;}
double current_sense_amplifier_gain(double V_out_fs,double V_shunt_fs){if(V_shunt_fs<=0)return 0;return V_out_fs/V_shunt_fs;}
double high_side_current_sense_CMV_error(double Vcm,double CMRR_db){return Vcm/pow(10.0,CMRR_db/20.0);}
double kelvin_connection_resistance(double R_trace_per_m,double length_m){return R_trace_per_m*length_m;}

/* Protection Circuits L7 */
double fuse_i2t_rating(double I_fault,double t_clear){return I_fault*I_fault*t_clear;}
double tvs_clamping_voltage(double V_breakdown,double I_peak,double R_dynamic){return V_breakdown+I_peak*R_dynamic;}
double esd_protection_resistor(double V_esd,double I_max){if(I_max<=0)return INFINITY;return V_esd/I_max;}
double polyfuse_hold_current_derating(double I_hold_25C,double T_ambient){return I_hold_25C*(1.0-0.005*(T_ambient-25.0));}
double crowbar_scr_gate_trigger_voltage(double V_overvoltage,double R1,double R2){return V_overvoltage*R2/(R1+R2);}

/* Voltage Reference Circuits L7 */
double zener_reference_voltage(double Vz,double Iz,double Rz){return Vz+Iz*Rz;}
double bandgap_reference_voltage(double Vbe,double delta_Vbe,double R2_R1_ratio){return Vbe+delta_Vbe*R2_R1_ratio;}
double reference_tempco_ppm_per_c(double V_ref_T1,double V_ref_T2,double T1,double T2){if(T1==T2)return 0;return(V_ref_T2-V_ref_T1)/V_ref_T1*1e6/(T2-T1);}
double reference_line_regulation(double dV_out,double dV_in){if(dV_in==0)return 0;return dV_out/dV_in;}
double reference_load_regulation(double dV_out,double dI_load){if(dI_load==0)return 0;return dV_out/dI_load;}

/* Analog Switch / Multiplexer L7 */
double analog_switch_on_resistance(double V_supply,double I_leakage){return V_supply/I_leakage;}
double mux_settling_time(double R_on,double C_load,int bits){return R_on*C_load*bits*log(2.0);}
double mux_charge_injection_error(double Q_injection,double C_hold){if(C_hold<=0)return INFINITY;return Q_injection/C_hold;}
double mux_crosstalk_db(double C_ch_ch,double C_ch_gnd){if(C_ch_gnd<=0)return -INFINITY;return 20*log10(C_ch_ch/C_ch_gnd);}

/* Oscillator Start-up L7 */
double crystal_startup_time_ms(double Q,double f0,double gm_crit_ratio){if(f0<=0)return 0;return Q/(M_PI*f0)*log(10/gm_crit_ratio)*1000;}
double crystal_load_capacitance(double C1,double C2,double C_stray){return(C1*C2)/(C1+C2)+C_stray;}
double crystal_pullability_ppm_per_pF(double C_motional,double C0,double C_L){if(C_L<=0)return 0;return C_motional*1e6/(2*(C0+C_L)*(C0+C_L));}
double rc_oscillator_startup_time(double R,double C,double gain){return R*C*3*gain;}

/* Battery Charging L7 */
double li_ion_charge_current_c_rate(double capacity_mAh,double C_rate){return capacity_mAh*C_rate/1000.0;}
double li_ion_cc_cv_transition_voltage(double V_full,int n_cells){return V_full*n_cells;}
double li_ion_state_of_charge_pct(double V_cell){if(V_cell<=3.0)return 0;if(V_cell>=4.2)return 100;return(V_cell-3.0)/1.2*100.0;}
double nimh_delta_v_termination_mV_per_cell(void){return -5.0;}
double battery_capacity_temperature_derating(double cap_25C,double T_celsius){return cap_25C*(1.0-0.01*fmax(T_celsius-25,0));}

/* Automotive Electrical L7 */
double load_dump_voltage_clamping(double V_alternator,double Z_source){return V_alternator*1.5;}
double iso7637_pulse1_voltage(double V_battery,double R_source){return -100;}
double reverse_battery_protection_diode_drop(double I_load,double Vf_diode){return Vf_diode;}
double cranking_voltage_droop_min(double V_battery_nominal,double R_internal,double I_starter){return fmax(V_battery_nominal-I_starter*R_internal,6.0);}
double can_bus_termination_resistance(void){return 120.0;}
double lin_bus_pullup_resistance(double V_bat,double V_diode,double I_pullup){if(I_pullup<=0)return 1000;return(V_bat-V_diode)/I_pullup;}
