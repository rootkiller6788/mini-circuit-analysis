#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "../include/transient_defs.h"
#include "../include/first_order.h"
#include "../include/second_order.h"
#include "../include/higher_order.h"
#include "../include/numerical_transient.h"
#include "../include/switching_transient.h"

#define EPS 1e-6

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s: ", #name); } while(0)
#define CHECK(cond) do { if(cond) { printf("PASS\n"); tests_passed++; } else { printf("FAIL (%s:%d)\n", __FILE__, __LINE__); } } while(0)
#define CHECK_NEAR(a,b,tol) do { double _va=(a),_vb=(b); if(fabs(_va-_vb)<(tol)) { printf("PASS\n"); tests_passed++; } else { printf("FAIL: %.6g vs %.6g (%s:%d)\n", _va,_vb,__FILE__,__LINE__); } } while(0)

static void test_time_constants(void) {
    TEST(rc_time_constant);
    CHECK_NEAR(compute_time_constant_rc(1000.0, 1e-6), 0.001, EPS);
    CHECK_NEAR(compute_time_constant_rl(100.0, 0.1), 0.001, EPS);
}

static void test_second_order_params(void) {
    TEST(second_order_params);
    SecondOrderParams_t p = compute_second_order_params(10.0, 1e-3, 100e-6, 0);
    CHECK_NEAR(p.omega_n, 1.0/sqrt(1e-3*100e-6), EPS);
    CHECK_NEAR(p.zeta, (10.0/2.0)*sqrt(100e-6/1e-3), EPS);
}

static void test_quality_factor(void) {
    TEST(quality_factor);
    double Qs = compute_quality_factor_series(10.0, 1e-3, 100e-6);
    CHECK_NEAR(Qs, sqrt(1e-3/100e-6)/10.0, EPS);
    double Qp = compute_quality_factor_parallel(1000.0, 1e-3, 100e-6);
    CHECK_NEAR(Qp, 1000.0*sqrt(100e-6/1e-3), EPS);
}

static void test_energy_state(void) {
    TEST(energy_state);
    EnergyState_t e = compute_energy_state(100e-6, 1e-3, 5.0, 0.5);
    CHECK_NEAR(e.energy_c, 0.5*100e-6*25.0, EPS);
    CHECK_NEAR(e.energy_l, 0.5*1e-3*0.25, EPS);
    CHECK_NEAR(e.total_energy, e.energy_c + e.energy_l, EPS);
}

static void test_rc_step_response(void) {
    TEST(rc_step_charge);
    double vc = rc_step_charge_vc(10.0, 0.0, 1000.0, 1e-6, 0.001);
    CHECK_NEAR(vc, 10.0*(1.0 - exp(-1.0)), 1e-6);
    vc = rc_step_charge_vc(10.0, 0.0, 1000.0, 1e-6, 0.005);
    CHECK_NEAR(vc/10.0, 1.0, 0.01);
}

static void test_rl_step_response(void) {
    TEST(rl_step_charge);
    double tau = 0.1/100.0;
    double il = rl_step_charge_il(10.0, 0.0, 100.0, 0.1, tau);
    CHECK_NEAR(il, 0.1*(1.0 - exp(-1.0)), 1e-6);
    il = rl_step_charge_il(10.0, 0.0, 100.0, 0.1, 5.0*tau);
    CHECK_NEAR(il, 0.1, 0.01);
}

static void test_rc_natural(void) {
    TEST(rc_natural);
    double vc = rc_natural_vc(5.0, 1000.0, 1e-6, 0.001);
    CHECK_NEAR(vc, 5.0*exp(-1.0), 1e-6);
}

static void test_rc_impulse(void) {
    TEST(rc_impulse);
    double vc = rc_impulse_vc(1000.0, 1e-6, 0.0);
    CHECK_NEAR(vc, 1.0/(1000.0*1e-6), 1e-6);
}

static void test_rc_sinusoidal(void) {
    TEST(rc_sinusoidal);
    double vc = rc_sinusoidal_vc(1.0, 2.0*M_PI*1000.0, 0.0, 1000.0, 100e-9, 0.001, 0.0);
    CHECK(isfinite(vc));
}

static void test_rlc_series_step_overdamped(void) {
    TEST(rlc_overdamped);
    double R=500.0, L=1e-3, C=100e-9;
    double vc = rlc_series_step_vc(10.0, R, L, C, 0.0, 0.0, 1e-4);
    CHECK(vc >= 0.0 && vc <= 10.0);
    double il = rlc_series_step_il(10.0, R, L, C, 0.0, 0.0, 1e-4);
    CHECK(isfinite(il));
}

static void test_rlc_series_step_underdamped(void) {
    TEST(rlc_underdamped);
    double R=1.0, L=1e-3, C=100e-6;
    double vc = rlc_series_step_vc(10.0, R, L, C, 0.0, 0.0, 1e-3);
    CHECK(vc > 0.0);
}

static void test_rlc_natural(void) {
    TEST(rlc_natural);
    double vc = rlc_series_natural_vc(10.0, 1e-3, 100e-6, 5.0, 0.1, 0.001);
    CHECK(isfinite(vc));
}

static void test_performance_metrics(void) {
    TEST(metrics);
    double os = overshoot_from_zeta(0.5);
    CHECK_NEAR(os, 16.3, 0.5);
    double z = zeta_from_overshoot(16.3);
    CHECK_NEAR(z, 0.5, 0.05);
    double ts = settling_time_2pct(0.5, 1000.0);
    CHECK_NEAR(ts, 0.008, 0.001);
}

static void test_design_formulas(void) {
    TEST(design);
    double R = design_R_for_zeta_series(0.707, 1e-3, 100e-6);
    CHECK_NEAR(R, 2.0*0.707*sqrt(1e-3/100e-6), 1e-6);
}

static void test_state_space(void) {
    TEST(state_space);
    StateSpace_t *ss = ss_alloc(2, 1, 1);
    CHECK(ss != NULL);
    ss->A[0] = 0.0; ss->A[1] = 1.0;
    ss->A[2] = -100.0; ss->A[3] = -20.0;
    ss->B[0] = 0.0; ss->B[1] = 100.0;
    ss->C[0] = 1.0; ss->C[1] = 0.0;
    double u[1] = {1.0};
    int r = ss_step_euler(ss, u, 0.001);
    CHECK(r == 0);
    ss_free(ss);
}

static void test_mosfet_switching(void) {
    TEST(mosfet_switching);
    double t_delay = mosfet_turn_on_delay(2.5, 10.0, 10.0, 2000e-12);
    CHECK(t_delay > 0.0);
    double Qg = mosfet_gate_charge_total(10e-9, 20e-9, 30e-9);
    CHECK_NEAR(Qg, 60e-9, EPS);
}

static void test_buck_converter(void) {
    TEST(buck_converter);
    double duty = buck_duty_cycle_ideal(5.0, 12.0);
    CHECK_NEAR(duty, 5.0/12.0, EPS);
    double ripple = buck_inductor_current_ripple(12.0, 5.0, 100e-6, 100e3);
    CHECK(ripple > 0.0);
}

static void test_diode_recovery(void) {
    TEST(diode_recovery);
    double Qrr = diode_reverse_recovery_charge(1.0, 50e-9);
    CHECK_NEAR(Qrr, 2.5e-8, 1e-10);
}

static void test_timing_apps(void) {
    TEST(timing_apps);
    double t_por = por_timeout_compute(100e3, 1e-6, 3.3, 0.8*3.3);
    CHECK(t_por > 0.0);
    double t_db = debounce_settling_time_compute(10e3, 0.1e-6, 5.0, 3.5);
    CHECK(t_db > 0.0);
}

static void test_tl_transient(void) {
    TEST(transmission_line);
    double ref = tl_reflection_coefficient(100.0, 50.0);
    CHECK_NEAR(ref, 50.0/150.0, EPS);
    double Td = tl_propagation_delay(0.5, 0.66);
    CHECK_NEAR(Td, 0.5/(3e8*0.66), 1e-12);
}

static void test_relay_transient(void) {
    TEST(relay_transient);
    double tau = relay_coil_time_constant(0.1, 100.0);
    CHECK_NEAR(tau, 0.001, EPS);
    double t_pull = relay_pull_in_time(12.0, 9.0, tau);
    CHECK(t_pull > 0.0);
}

static void test_pll_transient(void) {
    TEST(pll_transient);
    double wn = pll_natural_frequency(0.1, 1e6, 100.0);
    CHECK(wn > 0.0);
    double zeta = pll_damping_factor(0.1, 1e6, 100.0, 1000.0, 1e-6);
    CHECK(zeta > 0.0);
    double t_lock = pll_lock_time_estimate(wn, zeta, 1e6, 1000.0);
    CHECK(t_lock > 0.0);
}

static void test_component_models(void) {
    TEST(component_models);
    CapacitorModel_t cap;
    int r = capacitor_model_init(&cap, 1, "C1", 10e-6, 0.1, 1e-9);
    CHECK(r == 0);
    CHECK(cap.C == 10e-6);
    CHECK(cap.ESR == 0.1);
    SwitchModel_t sw;
    r = switch_model_init(&sw, 1, "SW1", 0.01, 1e6);
    CHECK(r == 0);
    CHECK(sw.R_on == 0.01);
    CHECK(sw.state == 0);
}

int main(void) {
    printf("=== mini-transient-analysis Test Suite ===\n");
    test_time_constants();
    test_second_order_params();
    test_quality_factor();
    test_energy_state();
    test_rc_step_response();
    test_rl_step_response();
    test_rc_natural();
    test_rc_impulse();
    test_rc_sinusoidal();
    test_rlc_series_step_overdamped();
    test_rlc_series_step_underdamped();
    test_rlc_natural();
    test_performance_metrics();
    test_design_formulas();
    test_state_space();
    test_mosfet_switching();
    test_buck_converter();
    test_diode_recovery();
    test_timing_apps();
    test_tl_transient();
    test_relay_transient();
    test_pll_transient();
    test_component_models();
    printf("\n=== Results: %d checks in %d test groups ===\n", tests_passed, tests_run);
    (void)tests_run;
    return 0;
}
