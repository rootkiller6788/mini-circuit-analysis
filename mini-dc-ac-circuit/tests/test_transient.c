#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "transient.h"

static int tr=0,tp=0;
#define T(n) do{tr++;printf("  %s: ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

static void t1(void){T("RC tau");double tau=rc_time_constant(1000,1e-6);assert(fabs(tau-0.001)<1e-9);P();}
static void t2(void){T("RL tau");double tau=rl_time_constant(100,1e-3);assert(fabs(tau-1e-5)<1e-12);P();}
static void t3(void){T("RC charge");double vc=rc_charging_vc(10,0,0.001,0.001);assert(fabs(vc-10*(1-exp(-1)))<0.01);P();}
static void t4(void){T("RC discharge");double vc=rc_discharging_vc(10,0.001,0.001);assert(fabs(vc-10*exp(-1))<0.01);P();}
static void t5(void){T("Energy");assert(fabs(capacitor_energy(1e-6,10)-5e-5)<1e-12);P();}
static void t6(void){T("RLC regime");RLC_Regime_t reg=rlc_series_regime(200,1e-3,1e-6);assert(reg==RLC_OVERDAMPED);P();}
static void t7(void){T("555 Timer");double f=timer_555_astable_freq(1000,10000,1e-6);assert(f>0);P();}
static void t8(void){T("Temp coeff");double R=temp_dependent_resistance(1000,25,50,0.00393);assert(fabs(R-1098.25)<0.1);P();}

int main(void){printf("=== test_transient ===\n");t1();t2();t3();t4();t5();t6();t7();t8();printf("Results: %d/%d\n",tp,tr);return(tp==tr)?0:1;}
