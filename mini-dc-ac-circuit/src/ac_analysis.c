#include "ac_analysis.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>

#ifndef I
#define I _Complex_I
#endif

/* Complex Number Utilities L3 */
double complex cpx_from_polar(double mag,double ph_deg){double ph=ph_deg*M_PI/180.0;return mag*cos(ph)+I*mag*sin(ph);}
double cpx_magnitude(double complex z){return cabs(z);}
double cpx_phase_deg(double complex z){return carg(z)*180.0/M_PI;}
double cpx_phase_rad(double complex z){return carg(z);}
double complex cpx_conjugate(double complex z){return conj(z);}

/* AC Impedance Functions L4 V=Z*I */
double complex ac_impedance_resistor(double R){return R+0*I;}
double complex ac_impedance_capacitor(double C,double omega){
  if(C<=0.0||omega<=0.0)return INFINITY+0*I;return 0-I/(omega*C);}
double complex ac_impedance_inductor(double L,double omega){
  if(L<=0.0||omega<=0.0)return 0+0*I;return 0+I*omega*L;}
double complex ac_impedance_capacitor_real(double C,double omega,double esr,double esl){
  if(C<=0.0||omega<=0.0)return esr-I/(omega*C);return esr+I*(omega*esl-1.0/(omega*C));}
double complex ac_impedance_inductor_real(double L,double omega,double dcr,double rc){
  double complex Z_ser=dcr+I*omega*L;if(rc<=0.0)return Z_ser;return Z_ser*rc/(Z_ser+rc);}
double complex ac_impedance_rlc_series(double R,double L,double C,double omega){
  double Xl=omega*L;double Xc=(omega*C>0)?(-1.0/(omega*C)):(-INFINITY);return R+I*(Xl+Xc);}
double complex ac_impedance_rlc_parallel(double R,double L,double C,double omega){
  double G=(R>0.0)?1.0/R:INFINITY;double Bc=omega*C;double Bl=(omega*L>0.0)?(-1.0/(omega*L)):(-INFINITY);
  double complex Y=G+I*(Bc+Bl);return 1.0/Y;}
double complex ac_equivalent_impedance(const Circuit_t *ckt,int n1,int n2,double omega){
  if(!ckt)return INFINITY+0*I;
  double complex Z_eq=0+0*I;
  for(int k=0;k<ckt->n_branches;k++){const Branch_t *b=&ckt->branches[k];
    if((b->node_from==n1&&b->node_to==n2)||(b->node_from==n2&&b->node_to==n1)){
      if(b->elem_type==ELEM_RESISTOR)Z_eq=b->value+0*I;
      else if(b->elem_type==ELEM_CAPACITOR)Z_eq=ac_impedance_capacitor(b->value,omega);
      else if(b->elem_type==ELEM_INDUCTOR)Z_eq=ac_impedance_inductor(b->value,omega);}}
  return Z_eq;}

/* Complex Gaussian Elimination L5 */
static int complex_gaussian_elim(double complex *A,const double complex *b,double complex *x,int n){
  if(!A||!b||!x||n<=0)return -1;
  for(int i=0;i<n;i++)x[i]=b[i];
  for(int col=0;col<n;col++){
    int max_row=col;double max_val=cabs(A[col*n+col]);
    for(int row=col+1;row<n;row++){double val=cabs(A[row*n+col]);if(val>max_val){max_val=val;max_row=row;}}
    if(max_val<1e-15)return -1;
    if(max_row!=col){for(int j=col;j<n;j++){double complex t=A[col*n+j];A[col*n+j]=A[max_row*n+j];A[max_row*n+j]=t;}double complex t=x[col];x[col]=x[max_row];x[max_row]=t;}
    double complex pivot=A[col*n+col];
    for(int row=col+1;row<n;row++){double complex factor=A[row*n+col]/pivot;A[row*n+col]=0;for(int j=col+1;j<n;j++)A[row*n+j]-=factor*A[col*n+j];x[row]-=factor*x[col];}}
  for(int i=n-1;i>=0;i--){double complex sum=x[i];for(int j=i+1;j<n;j++)sum-=A[i*n+j]*x[j];x[i]=sum/A[i*n+i];}
  return 0;}

/* AC Nodal Analysis L5 */
int ac_nodal_solve(int n,const double complex *Y,const double complex *I_src,double complex *V){
  if(!Y||!I_src||!V||n<=0)return -1;int sz=n*n;
  double complex *A=(double complex*)malloc((size_t)sz*sizeof(double complex));if(!A)return -1;
  for(int i=0;i<sz;i++)A[i]=Y[i];int r=complex_gaussian_elim(A,I_src,V,n);free(A);return r;}

/* Build AC Admittance Matrix L5 */
int ac_nodal_build_Y(const Circuit_t *ckt,double omega,double complex *Y,double complex *I_src){
  if(!ckt||!Y||!I_src)return -1;if(ckt->n_nodes<2)return -1;
  int n=ckt->n_nodes-1;int ns=n*n;for(int i=0;i<ns;i++)Y[i]=0;for(int i=0;i<n;i++)I_src[i]=0;
  int *nm=(int*)malloc((size_t)ckt->n_nodes*sizeof(int));if(!nm)return -1;
  int mi=0;for(int i=0;i<ckt->n_nodes;i++){if(i==ckt->ground_node)nm[i]=-1;else nm[i]=mi++;}
  for(int k=0;k<ckt->n_resistors;k++){const Resistor_t *r=&ckt->resistors[k];int ni=r->t[0].node_id,nj=r->t[1].node_id;
    if(r->R<=0.0)continue;double complex y=1.0/r->R;
    if(nm[ni]>=0)Y[nm[ni]*n+nm[ni]]+=y;if(nm[nj]>=0)Y[nm[nj]*n+nm[nj]]+=y;
    if(nm[ni]>=0&&nm[nj]>=0){Y[nm[ni]*n+nm[nj]]-=y;Y[nm[nj]*n+nm[ni]]-=y;}}
  for(int k=0;k<ckt->n_capacitors;k++){const Capacitor_t *c=&ckt->capacitors[k];int ni=c->t[0].node_id,nj=c->t[1].node_id;
    if(c->C<=0.0)continue;double complex y=I*omega*c->C;
    if(nm[ni]>=0)Y[nm[ni]*n+nm[ni]]+=y;if(nm[nj]>=0)Y[nm[nj]*n+nm[nj]]+=y;
    if(nm[ni]>=0&&nm[nj]>=0){Y[nm[ni]*n+nm[nj]]-=y;Y[nm[nj]*n+nm[ni]]-=y;}}
  for(int k=0;k<ckt->n_inductors;k++){const Inductor_t *l=&ckt->inductors[k];int ni=l->t[0].node_id,nj=l->t[1].node_id;
    if(l->L<=0.0)continue;double complex y=1.0/(I*omega*l->L);
    if(nm[ni]>=0)Y[nm[ni]*n+nm[ni]]+=y;if(nm[nj]>=0)Y[nm[nj]*n+nm[nj]]+=y;
    if(nm[ni]>=0&&nm[nj]>=0){Y[nm[ni]*n+nm[nj]]-=y;Y[nm[nj]*n+nm[ni]]-=y;}}
  for(int k=0;k<ckt->n_ac_isrc;k++){const ACCurrentSource_t *s=&ckt->ac_isrc[k];
    int ni=s->t[0].node_id,nj=s->t[1].node_id;
    double complex Is=cpx_from_polar(s->I_rms,s->phase_deg);
    if(nm[ni]>=0)I_src[nm[ni]]+=Is;if(nm[nj]>=0)I_src[nm[nj]]-=Is;}
  for(int k=0;k<ckt->n_ac_vsrc;k++){const ACVoltageSource_t *s=&ckt->ac_vsrc[k];
    int np=s->t[0].node_id,nm_s=s->t[1].node_id;double complex Vs=cpx_from_polar(s->V_rms,s->phase_deg);
    if(nm_s==ckt->ground_node&&nm[np]>=0){int idx=nm[np];for(int j=0;j<n;j++)Y[idx*n+j]=0;Y[idx*n+idx]=1;I_src[idx]=Vs;}
    else if(np==ckt->ground_node&&nm[nm_s]>=0){int idx=nm[nm_s];for(int j=0;j<n;j++)Y[idx*n+j]=0;Y[idx*n+idx]=1;I_src[idx]=-Vs;}}
  free(nm);return n;}

/* AC Mesh Analysis L5 */
int ac_mesh_solve(int m,const double complex *Z,const double complex *V,double complex *I_mesh){
  if(!Z||!V||!I_mesh||m<=0)return -1;int sz=m*m;
  double complex *A=(double complex*)malloc((size_t)sz*sizeof(double complex));if(!A)return -1;
  for(int i=0;i<sz;i++)A[i]=Z[i];int r=complex_gaussian_elim(A,V,I_mesh,m);free(A);return r;}

/* Frequency Response L5 L6 */
int freq_response_rc_lowpass(double R,double C,double fs,double fe,int np,FreqPoint_t *resp){
  if(!resp||np<=0||fs<=0||fe<=fs)return -1;double fc=1.0/(2*M_PI*R*C);
  for(int i=0;i<np;i++){double f=fs+(fe-fs)*i/(np-1);resp[i].freq_hz=f;
    resp[i].magnitude=1.0/sqrt(1.0+(f/fc)*(f/fc));resp[i].phase_deg=-atan2(f/fc,1.0)*180.0/M_PI;}return 0;}
int freq_response_rc_highpass(double R,double C,double fs,double fe,int np,FreqPoint_t *resp){
  if(!resp||np<=0||fs<=0||fe<=fs)return -1;double fc=1.0/(2*M_PI*R*C);
  for(int i=0;i<np;i++){double f=fs+(fe-fs)*i/(np-1);resp[i].freq_hz=f;
    double ratio=f/fc;resp[i].magnitude=ratio/sqrt(1.0+ratio*ratio);resp[i].phase_deg=90.0-atan2(ratio,1.0)*180.0/M_PI;}return 0;}

/* Series RLC Frequency Response */
int freq_response_rlc_series(double R,double L,double C,RLC_Output_t out,double fs,double fe,int np,FreqPoint_t *resp){
  if(!resp||np<=0||fs<=0||fe<=fs)return -1;
  double w0=1.0/sqrt(L*C);double f0=w0/(2*M_PI);
  for(int i=0;i<np;i++){double f=fs+(fe-fs)*i/(np-1);resp[i].freq_hz=f;double w=2*M_PI*f;
    double complex Z=R+I*(w*L-1.0/(w*C));
    double complex H;if(out==FR_RLC_ACROSS_R)H=R/Z;
    else if(out==FR_RLC_ACROSS_L)H=I*w*L/Z;else H=(1.0/(I*w*C))/Z;
    resp[i].magnitude=cabs(H);resp[i].phase_deg=carg(H)*180.0/M_PI;}return 0;}
int freq_response_rlc_parallel(double R,double L,double C,double fs,double fe,int np,FreqPoint_t *resp){
  if(!resp||np<=0||fs<=0||fe<=fs)return -1;
  for(int i=0;i<np;i++){double f=fs+(fe-fs)*i/(np-1);resp[i].freq_hz=f;double w=2*M_PI*f;
    double complex Y=1.0/R+I*(w*C-1.0/(w*L));double complex H=1.0/(1.0+R*Y);
    resp[i].magnitude=cabs(H);resp[i].phase_deg=carg(H)*180.0/M_PI;}return 0;}

/* Complex Power L2 L4 S=P+jQ */
ComplexPower_t complex_power_from_vi(double Vr,double Vp,double Ir,double Ip){
  ComplexPower_t cp;double ph=(Vp-Ip)*M_PI/180.0;cp.phi_deg=Vp-Ip;cp.pf=cos(ph);
  cp.P=Vr*Ir*cos(ph);cp.Q=Vr*Ir*sin(ph);cp.S=Vr*Ir;
  cp.leading=(cp.Q<0)?1:0;return cp;}
ComplexPower_t complex_power_from_v_z(double Vr,double complex Z){
  ComplexPower_t cp;double I2=(Vr*Vr)/cabs(Z*conj(Z));
  cp.S=Vr*Vr/cabs(Z);cp.P=I2*creal(Z);cp.Q=I2*cimag(Z);
  cp.pf=(cp.S>0)?cp.P/cp.S:1.0;cp.phi_deg=atan2(cp.Q,cp.P)*180.0/M_PI;cp.leading=(cp.Q<0)?1:0;return cp;}
ComplexPower_t complex_power_from_i_z(double Ir,double complex Z){
  ComplexPower_t cp;cp.P=Ir*Ir*creal(Z);cp.Q=Ir*Ir*cimag(Z);
  cp.S=Ir*Ir*cabs(Z);cp.pf=(cp.S>0)?cp.P/cp.S:1.0;cp.phi_deg=atan2(cp.Q,cp.P)*180.0/M_PI;cp.leading=(cp.Q<0)?1:0;return cp;}

/* Power Factor Correction L6 L7 */
double pf_correction_capacitance(double P,double V,double f,double pf_old,double pf_target){
  double phi_old=acos(pf_old);double phi_new=acos(pf_target);
  double omega=2*M_PI*f;if(omega<=0||V<=0)return 0.0;
  return P*(tan(phi_old)-tan(phi_new))/(omega*V*V);}
double pf_correction_unity(double P,double V,double f,double pf_old,int lag){
  if(!lag)return 0.0;double phi=acos(pf_old);double omega=2*M_PI*f;
  if(omega<=0||V<=0)return 0.0;return P*tan(phi)/(omega*V*V);}

/* Resonance Analysis L6 */
Resonance_t resonance_series_analyze(double R,double L,double C){
  Resonance_t res;res.omega0=1.0/sqrt(L*C);res.f0_hz=res.omega0/(2*M_PI);
  res.Q=res.omega0*L/R;res.zeta=1.0/(2*res.Q);res.BW_radps=R/L;res.BW_hz=res.BW_radps/(2*M_PI);
  res.f_lower_hz=res.f0_hz-res.BW_hz/2;res.f_upper_hz=res.f0_hz+res.BW_hz/2;res.Z_at_res=R;
  if(res.f_lower_hz<0)res.f_lower_hz=0;return res;}
Resonance_t resonance_parallel_analyze(double R,double L,double C){
  Resonance_t res;res.omega0=1.0/sqrt(L*C);res.f0_hz=res.omega0/(2*M_PI);
  res.Q=R/(res.omega0*L);res.zeta=1.0/(2*res.Q);res.BW_radps=1.0/(R*C);res.BW_hz=res.BW_radps/(2*M_PI);
  res.f_lower_hz=res.f0_hz-res.BW_hz/2;res.f_upper_hz=res.f0_hz+res.BW_hz/2;res.Z_at_res=R;
  if(res.f_lower_hz<0)res.f_lower_hz=0;if(res.f_upper_hz>res.f0_hz*10)res.f_upper_hz=res.f0_hz*10;return res;}

/* Filter Transfer Functions L6 */
double rc_lowpass_cutoff(double R,double C){return 1.0/(2*M_PI*R*C);}
double complex rc_lowpass_h(double R,double C,double complex s){return 1.0/(1.0+s*R*C);}
double rc_lowpass_magnitude(double R,double C,double f){double rc=R*C;double w=2*M_PI*f;return 1.0/sqrt(1.0+w*w*rc*rc);}
double rc_lowpass_phase(double R,double C,double f){double w=2*M_PI*f;return -atan2(w*R*C,1.0)*180.0/M_PI;}
double rc_highpass_cutoff(double R,double C){return 1.0/(2*M_PI*R*C);}
double complex rc_highpass_h(double R,double C,double complex s){double complex num=s*R*C;return num/(1.0+num);}
double rc_highpass_magnitude(double R,double C,double f){double w=2*M_PI*f;double x=w*R*C;return x/sqrt(1.0+x*x);}
double sallen_key_lp_cutoff(double R1,double R2,double C1,double C2){return 1.0/(2*M_PI*sqrt(R1*R2*C1*C2));}
double sallen_key_lp_q(double R1,double R2,double C1,double C2){double n=sqrt(R1*R2*C1*C2);return n/(C2*(R1+R2));}
double complex sallen_key_lp_h(double R1,double R2,double C1,double C2,double complex s){
  double complex denom=1.0+s*C2*(R1+R2)+s*s*R1*R2*C1*C2;return 1.0/denom;}
double sallen_key_hp_cutoff(double R1,double R2,double C1,double C2){return 1.0/(2*M_PI*sqrt(R1*R2*C1*C2));}
double bp_center_freq(double L,double C){return 1.0/(2*M_PI*sqrt(L*C));}
double bp_bandwidth(double R,double L){return R/(2*M_PI*L);}
double bp_quality(double f0,double bw){if(bw<=0)return INFINITY;return f0/bw;}
double notch_attenuation_db(double R,double L,double C,double f){
  double w=2*M_PI*f;double complex Y=1.0/R+I*(w*C-1.0/(w*L));double complex Z=1.0/Y;
  double complex H=Z/(Z+R);return 20*log10(cabs(H));}

/* Bode Plot Generator L5 L6 */
static double complex eval_poly(const double *c,int deg,double complex s){
  double complex val=0;for(int i=deg;i>=0;i--)val=val*s+c[i];return val;}
int bode_compute(const double *b,int m,const double *a,int n,double fs,double fe,int np,double *mag,double *ph){
  if(!b||!a||!mag||!ph||np<=0||fs<=0||fe<=fs)return -1;
  for(int i=0;i<np;i++){double f=fs+(fe-fs)*i/(np-1);double w=2*M_PI*f;double complex s=I*w;
    double complex H=eval_poly(b,m,s)/eval_poly(a,n,s);mag[i]=20*log10(cabs(H));ph[i]=carg(H)*180.0/M_PI;}return 0;}
int bode_single_pole(double fp,double fs,double fe,int np,double *mag,double *ph){
  double a[2]={1.0,1.0/(2*M_PI*fp)};double b[1]={1.0};return bode_compute(b,0,a,1,fs,fe,np,mag,ph);}
int bode_single_zero(double fz,double fs,double fe,int np,double *mag,double *ph){
  double b[2]={0.0,1.0/(2*M_PI*fz)};double a[1]={1.0};return bode_compute(b,1,a,0,fs,fe,np,mag,ph);}
int bode_complex_poles(double wn,double zeta,double fs,double fe,int np,double *mag,double *ph){
  double a[3]={wn*wn,2*zeta*wn,1.0};double b[1]={wn*wn};return bode_compute(b,0,a,2,fs,fe,np,mag,ph);}

/* Impedance Matching L6 L8 */
ImpedanceMatch_t impedance_match_l_network(double Rs,double Rl,double f){
  ImpedanceMatch_t m;m.valid=0;
  if(Rs<=0||Rl<=0||f<=0){m.L_henries=m.C_farads=0;return m;}
  double w=2*M_PI*f;
  if(Rl>Rs){double Q=sqrt(Rl/Rs-1.0);m.C_farads=Q/(w*Rl);m.L_henries=Q*Rs/w;m.type=MATCH_LP;m.valid=1;}
  else{double Q=sqrt(Rs/Rl-1.0);m.L_henries=Q*Rl/w;m.C_farads=Q/(w*Rs);m.type=MATCH_HP;m.valid=1;}
  return m;}
PiMatch_t impedance_match_pi_network(double Rs,double Rl,double f,double Q){
  PiMatch_t m;m.valid=0;if(Rs<=0||Rl<=0||f<=0||Q<=0){m.C1=m.C2=m.L=0;return m;}
  double w=2*M_PI*f;double Rv=(Rs<Rl)?Rs:Rl;
  m.L=Rv*Q/w;m.C1=Q/(w*Rs);m.C2=Q/(w*Rl);m.valid=1;return m;}
TMatch_t impedance_match_t_network(double Rs,double Rl,double f,double Q){
  TMatch_t m;m.valid=0;if(Rs<=0||Rl<=0||f<=0||Q<=0){m.L1=m.L2=m.C=0;return m;}
  double w=2*M_PI*f;double Rv=(Rs>Rl)?Rs:Rl;
  m.C=Q/(w*Rv);m.L1=Q*Rs/w;m.L2=Q*Rl/w;m.valid=1;return m;}

/* RF/Microwave Functions L7 L8 */
double complex reflection_coefficient(double complex ZL,double Z0){
  if(Z0<=0)return 1+0*I;return(ZL-Z0)/(ZL+Z0);}
double return_loss_db(double complex ZL,double Z0){
  double complex G=reflection_coefficient(ZL,Z0);double mag=cabs(G);
  if(mag<=0)return INFINITY;return -20*log10(mag);}
double vswr(double complex ZL,double Z0){
  double complex G=reflection_coefficient(ZL,Z0);double mag=cabs(G);
  if(mag>=1.0)return INFINITY;return(1.0+mag)/(1.0-mag);}
double mismatch_loss_db(double complex ZL,double Z0){
  double complex G=reflection_coefficient(ZL,Z0);double mag2=cabs(G)*cabs(G);
  if(mag2>=1.0)return -INFINITY;return -10*log10(1.0-mag2);}
double unilateral_transducer_gain(double s21,double s11,double s22,double complex gs,double complex gl){
  double gs_mag2=cabs(gs)*cabs(gs);double gl_mag2=cabs(gl)*cabs(gl);
  double Gs=(1.0-gs_mag2)/cabs(1.0-s11*gs)/cabs(1.0-s11*gs);
  double G0=s21*s21;
  double Gl=(1.0-gl_mag2)/cabs(1.0-s22*gl)/cabs(1.0-s22*gl);
  return Gs*G0*Gl;}
double rollett_stability_k(double s11m,double s12m,double s21m,double s22m){
  double delta=s11m*s22m-s12m*s21m;double num=1.0-s11m*s11m-s22m*s22m+delta*delta;
  double denom=2.0*s12m*s21m;if(denom<=0)return INFINITY;return num/denom;}
double stability_delta_mag(double complex s11,double complex s12,double complex s21,double complex s22){
  return cabs(s11*s22-s12*s21);}
double cascaded_noise_figure_linear(const double *nf,const double *g,int n){
  if(!nf||!g||n<=0)return 0;double NF=nf[0];double Gcum=g[0];
  for(int i=1;i<n;i++){NF+=(nf[i]-1.0)/Gcum;Gcum*=g[i];}return NF;}
double cascaded_noise_figure_db(const double *nf_db,const double *g_db,int n){
  if(!nf_db||!g_db||n<=0)return 0;
  double *nfl=(double*)malloc((size_t)n*sizeof(double));
  double *gl=(double*)malloc((size_t)n*sizeof(double));
  if(!nfl||!gl){free(nfl);free(gl);return 0;}
  for(int i=0;i<n;i++){nfl[i]=pow(10.0,nf_db[i]/10.0);gl[i]=pow(10.0,g_db[i]/10.0);}
  double NF=cascaded_noise_figure_linear(nfl,gl,n);free(nfl);free(gl);return 10*log10(NF);}

/* dB Conversion L7 */


double p1db_from_psat(double p_sat_dbm){return p_sat_dbm-10.0;}
double iip3_from_two_tone(double pin,double dp){return pin+dp/2.0;}
double oip3_from_iip3(double iip3,double g){return iip3+g;}

/* Additional Filter Types L6 */
double rl_lowpass_cutoff(double R,double L){return R/(2*M_PI*L);}
double rl_highpass_cutoff(double R,double L){return R/(2*M_PI*L);}
double lc_lowpass_cutoff(double L,double C){return 1.0/(2*M_PI*sqrt(L*C));}

/* Filter Order Estimation L6 */
int butterworth_order_estimate(double fs,double fp,double As,double Ap){
  if(fs<=fp||As<=Ap)return -1;double N=log10((pow(10,As/10)-1)/(pow(10,Ap/10)-1))/(2*log10(fs/fp));return(int)ceil(N);}

/* Attenuator Networks L6 L7 */
void pi_attenuator(double Z0,double atten_db,double *R1,double *R2,double *R3){
  double a=pow(10.0,atten_db/20.0);double k=(a-1)/(a+1);if(R1)*R1=Z0*(a*a-1)/(2*a);
  if(R2)*R2=Z0*k;if(R3)*R3=*R1;}
void t_attenuator(double Z0,double atten_db,double *R1,double *R2,double *R3){
  double a=pow(10.0,atten_db/20.0);if(R1)*R1=Z0*(a-1)/(a+1);
  if(R2)*R2=2*Z0*a/(a*a-1);if(R3)*R3=*R1;}

/* LC Resonant Tank L6 */
double lc_tank_impedance_at_resonance(double L,double C,double R_parasitic){
  double w0=1.0/sqrt(L*C);double Q=w0*L/R_parasitic;return Q*Q*R_parasitic;}
double lc_tank_3db_bandwidth(double f0,double Q){if(Q<=0)return INFINITY;return f0/Q;}

/* Coupling Capacitor Selection L7 */
double coupling_capacitor_min(double f_low,double R_load){if(f_low<=0||R_load<=0)return 0;return 1.0/(2*M_PI*f_low*R_load);}
double bypass_capacitor_min(double f_low,double R_equiv){if(f_low<=0||R_equiv<=0)return 0;return 1.0/(2*M_PI*f_low*R_equiv);}

/* Transformer Impedance Reflection L4 */
double transformer_reflected_impedance(double Z_load,double n_turns_ratio){return Z_load/(n_turns_ratio*n_turns_ratio);}
double transformer_primary_impedance(double Z_load,double n,double Lm,double omega){
  double complex Z_reflected=transformer_reflected_impedance(Z_load,n);
  double complex Z_mag=I*omega*Lm;return cabs(Z_reflected*Z_mag/(Z_reflected+Z_mag));}

/* dBm to Voltage Conversions L7 */
double dbm_to_vpeak(double dbm,double Z0){double Vrms=sqrt(pow(10.0,(dbm-30.0)/10.0)*Z0);return Vrms*sqrt(2.0);}
double vpeak_to_dbm(double Vpeak,double Z0){double Vrms=Vpeak/sqrt(2.0);return 10*log10(Vrms*Vrms/Z0)+30.0;}
double dbuv_to_dbm(double dbuv){return dbuv-107.0;}
double dbm_to_dbuv(double dbm){return dbm+107.0;}

/* Antenna Impedance Matching Helper L7 */
double quarter_wave_transformer_Z0(double Z_load,double Z_source){return sqrt(Z_load*Z_source);}
double quarter_wave_length_m(double freq_hz,double velocity_factor){return 3e8*velocity_factor/(4*freq_hz);}
double stub_length_open_c(double freq_hz,double Z0,double target_Z){
  double w=2*M_PI*freq_hz;double B=1.0/target_Z;return atan(B*Z0)/w*3e8;}

/* Multi-stage Amplifier Analysis L7 */
double multistage_gain_linear(const double *gains_linear,int n){double g=1;for(int i=0;i<n;i++)g*=gains_linear[i];return g;}
double multistage_gain_db(const double *gains_db,int n){double g=0;for(int i=0;i<n;i++)g+=gains_db[i];return g;}
double multistage_bandwidth_cascaded(double *bw_hz,int n,double exponent){
  if(n<=0)return INFINITY;double sum=0;for(int i=0;i<n;i++){if(bw_hz[i]<=0)continue;sum+=pow(1.0/bw_hz[i],exponent);}
  if(sum<=0)return INFINITY;return pow(1.0/sum,1.0/exponent);}
double cascaded_risetime_ns(const double *tr_ns,int n){double sum=0;for(int i=0;i<n;i++)sum+=tr_ns[i]*tr_ns[i];return sqrt(sum);}

/* Intermodulation Distortion L8 */
double imd3_level(double fund_dbm,double im3_dbm){return fund_dbm-im3_dbm;}
static double fund_dbm(double p);
double imd3_intercept_point(double p1_dbm,double p2_dbm,double im3_1_dbm,double im3_2_dbm){
  double slope=(im3_2_dbm-im3_1_dbm)/(p2_dbm-p1_dbm);return p1_dbm+(fund_dbm(p1_dbm)-im3_1_dbm)/(3.0-slope);}
static double fund_dbm(double p){return p;}
double sfdr_db(double iip3_dbm,double noise_floor_dbm){return(2.0*iip3_dbm+noise_floor_dbm)/3.0-noise_floor_dbm;}

/* Error Vector Magnitude L7 */
double evm_percent(const double complex *ideal,const double complex *measured,int n){
  if(n<=0)return 0;double num=0,den=0;
  for(int i=0;i<n;i++){double err=cabs(measured[i]-ideal[i]);num+=err*err;den+=cabs(ideal[i])*cabs(ideal[i]);}
  if(den<=0)return 0;return sqrt(num/den)*100.0;}
double evm_db(double evm_pct){if(evm_pct<=0)return -INFINITY;return 20*log10(evm_pct/100.0);}

/* Phase Noise to Jitter L7 */
double phase_noise_to_rms_jitter_sec(double phase_noise_dbc_hz,double freq_hz,double offset_hz){
  if(freq_hz<=0||offset_hz<=0)return 0;return sqrt(2*pow(10,phase_noise_dbc_hz/10)*offset_hz)/(2*M_PI*freq_hz);}
double rms_jitter_to_phase_noise(double rms_jitter_sec,double freq_hz,double offset_hz){
  if(freq_hz<=0||offset_hz<=0||rms_jitter_sec<=0)return -INFINITY;
  double A=2*M_PI*freq_hz*rms_jitter_sec;return 10*log10(A*A/(2*offset_hz));}

/* Transmission Line Basics L7 */
double txline_propagation_delay_ps_per_m(double er_effective){return 1000.0*sqrt(er_effective)/0.3;}
double txline_wavelength_m(double freq_hz,double er_effective){return 3e8/(freq_hz*sqrt(er_effective));}
double txline_impedance_microstrip(double h_mm,double w_mm,double er){
  if(w_mm<=0||h_mm<=0)return 0;double wh=w_mm/h_mm;
  if(wh<=1){double ee=(er+1)/2+(er-1)/2*(1/sqrt(1+12/wh)+0.04*(1-wh)*(1-wh));return 60*log(8/wh+wh/4)/sqrt(ee);}
  else{double ee=(er+1)/2+(er-1)/2*(1/sqrt(1+12/wh));return 120*M_PI/(sqrt(ee)*(wh+1.393+0.667*log(wh+1.444)));}}
double txline_loss_db_per_m(double freq_hz,double Z0,double loss_tangent,double er_effective){
  return 27.3*sqrt(er_effective)*loss_tangent*freq_hz/(3e8);}

/* Skin Effect L8 */
double skin_depth_m(double freq_hz,double conductivity_sm,double mu_r){
  if(freq_hz<=0||conductivity_sm<=0)return INFINITY;return 1.0/sqrt(M_PI*freq_hz*4e-7*M_PI*mu_r*conductivity_sm);}
double skin_effect_ac_resistance(double R_dc,double freq_hz,double radius_m,double conductivity_sm,double mu_r){
  double delta=skin_depth_m(freq_hz,conductivity_sm,mu_r);if(delta<=0||radius_m<=0)return R_dc;
  return R_dc*radius_m/(2*delta);}

/* Mixer and Frequency Conversion L8 */
double mixer_conversion_loss_db(double P_RF_dbm,double P_IF_dbm){return P_RF_dbm-P_IF_dbm;}
double mixer_noise_figure_db(double conversion_loss_db,double NF_IF_amp_db){
  double cl=pow(10,conversion_loss_db/10);double nf_if=pow(10,NF_IF_amp_db/10);return 10*log10(cl*nf_if);}
double image_frequency(double f_RF,double f_LO,int low_side){return low_side?f_LO*2-f_RF:f_RF-2*f_LO;}
double image_rejection_ratio_db(double f_RF,double f_IF,double Q){if(Q<=0)return 0;return 20*log10(Q*f_IF/f_RF);}

/* Oscillator Phase Noise L8 */
double leeson_phase_noise_dbc_hz(double f_m,double f0,double QL,double F,double kT,double P_sig){
  if(f_m<=0||f0<=0||QL<=0||P_sig<=0)return 0;return 10*log10(2*F*kT/P_sig*(1+(f0/(2*QL*f_m))*(f0/(2*QL*f_m))));}
double phase_noise_to_jitter_rad(double L_dbc_hz,double f_start,double f_stop){
  if(f_start>=f_stop)return 0;return sqrt(2*pow(10,L_dbc_hz/10)*(f_stop-f_start));}

/* IQ Modulation L8 */
void iq_modulate(double I_data,double Q_data,double carrier_freq,double t,double *rf_out){
  *rf_out=I_data*cos(2*M_PI*carrier_freq*t)-Q_data*sin(2*M_PI*carrier_freq*t);}
void iq_demodulate(double rf_signal,double carrier_freq,double t,double *I_out,double *Q_out){
  *I_out=rf_signal*cos(2*M_PI*carrier_freq*t);*Q_out=-rf_signal*sin(2*M_PI*carrier_freq*t);}
double iq_imbalance_image_rejection(double gain_imbalance_db,double phase_imbalance_deg){
  double g=pow(10,gain_imbalance_db/20);double p=phase_imbalance_deg*M_PI/180;
  return 10*log10((g*g+1+2*g*cos(p))/(g*g+1-2*g*cos(p)));}

/* Smith Chart helper L8 */
double smith_chart_gamma_to_impedance(double gamma_mag,double gamma_ang_deg,double Z0){
  double complex G=cpx_from_polar(gamma_mag,gamma_ang_deg);return cabs(Z0*(1.0+G)/(1.0-G));}
double smith_chart_impedance_to_gamma_mag(double Z,double Z0){return fabs((Z-Z0)/(Z+Z0));}
double smith_chart_vswr_circle_radius(double VSWR){return(VSWR-1)/(VSWR+1);}

/* S-Parameter Conversions L8 */
void s_to_z_params(double s11_mag,double s11_ang,double s12_mag,double s12_ang,
                    double s21_mag,double s21_ang,double s22_mag,double s22_ang,double Z0,
                    double *z11,double *z12,double *z21,double *z22){
  double complex S11=cpx_from_polar(s11_mag,s11_ang);double complex S12=cpx_from_polar(s12_mag,s12_ang);
  double complex S21=cpx_from_polar(s21_mag,s21_ang);double complex S22=cpx_from_polar(s22_mag,s22_ang);
  double complex det_S=S11*S22-S12*S21;
  double complex denom=(1.0-S11)*(1.0-S22)-S12*S21;
  *z11=Z0*cabs(((1+S11)*(1-S22)+S12*S21)/denom);*z12=Z0*cabs((2*S12)/denom);
  *z21=Z0*cabs((2*S21)/denom);*z22=Z0*cabs(((1-S11)*(1+S22)+S12*S21)/denom);}

/* Wilkinson Power Divider L7 */
double wilkinson_Z_quarterwave(double Z0){return Z0*sqrt(2.0);}
double wilkinson_isolation_resistor(double Z0){return 2.0*Z0;}

/* Directional Coupler L8 */
double coupler_coupling_db(double P_input,double P_coupled){return 10*log10(P_input/P_coupled);}
double coupler_directivity_db(double P_coupled,double P_isolated){if(P_isolated<=0)return INFINITY;return 10*log10(P_coupled/P_isolated);}
double coupler_insertion_loss_db(double P_input,double P_output){return 10*log10(P_input/P_output);}

/* Link Budget L7 */
double link_budget_rx_power_dbm(double tx_power_dbm,double tx_gain_dbi,double path_loss_db,double rx_gain_dbi){return tx_power_dbm+tx_gain_dbi-path_loss_db+rx_gain_dbi;}
double free_space_path_loss_db(double distance_m,double freq_hz){return 20*log10(distance_m)+20*log10(freq_hz)-147.55;}
double friis_rx_power_watts(double tx_power_w,double tx_gain,double rx_gain,double wavelength_m,double distance_m){
  if(distance_m<=0)return 0;return tx_power_w*tx_gain*rx_gain*wavelength_m*wavelength_m/(16*M_PI*M_PI*distance_m*distance_m);}
double link_margin_db(double rx_power_dbm,double rx_sensitivity_dbm){return rx_power_dbm-rx_sensitivity_dbm;}
double rx_sensitivity_dbm(double noise_figure_db,double bandwidth_hz,double required_snr_db){return -174+10*log10(bandwidth_hz)+noise_figure_db+required_snr_db;}
double maximum_range_m(double tx_power_dbm,double tx_gain_dbi,double rx_gain_dbi,double rx_sensitivity_dbm,double freq_hz){
  double margin_db=tx_power_dbm+tx_gain_dbi+rx_gain_dbi-rx_sensitivity_dbm;return pow(10,(margin_db+147.55-20*log10(freq_hz))/20);}
