#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "ac_analysis.h"

static int tr=0,tp=0;
#define T(n) do{tr++;printf("  %s: ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

static void t1(void){T("Complex");double complex z=cpx_from_polar(1,90);assert(fabs(cpx_magnitude(z)-1)<1e-9);P();}
static void t2(void){T("AC Z");double complex Zc=ac_impedance_capacitor(1e-6,1000);assert(cimag(Zc)<0);P();}
static void t3(void){T("RLC Z");double complex Z=ac_impedance_rlc_series(50,1e-3,1e-6,1000);assert(fabs(creal(Z)-50)<1e-9);P();}
static void t4(void){T("Power");ComplexPower_t cp=complex_power_from_vi(120,0,10,-30);assert(cp.P>0);P();}
static void t5(void){T("Resonance");Resonance_t r=resonance_series_analyze(10,1e-3,1e-6);assert(r.f0_hz>0);P();}
static void t6(void){T("Filter");double fc=rc_lowpass_cutoff(1000,1e-6);double m=rc_lowpass_magnitude(1000,1e-6,fc);assert(fabs(m-0.7071)<0.002);P();}
static void t7(void){T("VSWR");double v=vswr(75,50);assert(v>=1);P();}
static void t8(void){T("NF");double nf[]={2,3};double g[]={10,10};double NF=cascaded_noise_figure_linear(nf,g,2);assert(NF>2);P();}

int main(void){printf("=== test_ac ===\n");t1();t2();t3();t4();t5();t6();t7();t8();printf("Results: %d/%d\n",tp,tr);return(tp==tr)?0:1;}
