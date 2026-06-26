#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "dc_analysis.h"

static int tests_run = 0, tests_passed = 0;
#define TEST(n) do { tests_run++; printf("  %s: ", n); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void t_ohm(void){TEST("Ohm");assert(fabs(ohm_law_voltage(0.01,1000)-10)<1e-9);PASS();}
static void t_joule(void){TEST("Joule");assert(fabs(joule_power_from_vi(10,0.01)-0.1)<1e-9);PASS();}
static void t_kcl(void){TEST("KCL");double c[]={1,2,-3};assert(fabs(kcl_residual(c,3))<1e-12);PASS();}
static void t_kvl(void){TEST("KVL");double v[]={5,-3,-2};assert(fabs(kvl_residual(v,3))<1e-12);PASS();}
static void t_vdiv(void){TEST("Vdiv");assert(fabs(voltage_divider(10,1000,1000)-5)<1e-9);PASS();}
static void t_cdiv(void){TEST("Cdiv");assert(fabs(current_divider_2r(1,100,100)-0.5)<1e-9);PASS();}
static void t_eqr(void){TEST("EqR");double R[]={100,200,300};assert(fabs(resistance_series(R,3)-600)<1e-9);assert(fabs(resistance_parallel_2(100,100)-50)<1e-9);PASS();}
static void t_thev(void){TEST("Thev");TheveninEquiv_t th=thevenin_from_oc_sc(10,1);assert(fabs(th.V_th-10)<1e-9);assert(fabs(th.R_th-10)<1e-9);PASS();}
static void t_maxp(void){TEST("MaxP");assert(fabs(max_power_transfer_pmax(10,10)-2.5)<1e-9);PASS();}
static void t_dw(void){TEST("D-W");double r1,r2,r3;delta_to_wye(300,300,300,&r1,&r2,&r3);assert(fabs(r1-100)<1e-9);PASS();}
static void t_wb(void){TEST("WB");assert(wheatstone_is_balanced(1000,1000,1000,1000,1e-6));PASS();}

int main(void){printf("=== test_dc ===\n");t_ohm();t_joule();t_kcl();t_kvl();t_vdiv();t_cdiv();t_eqr();t_thev();t_maxp();t_dw();t_wb();printf("Results: %d/%d\n",tests_passed,tests_run);return(tests_passed==tests_run)?0:1;}
