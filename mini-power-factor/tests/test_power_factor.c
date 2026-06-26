/**
 * @file test_power_factor.c
 * @brief Comprehensive test suite for mini-power-factor library
 *
 * Tests cover all L1-L6 knowledge levels with assert-based verification.
 * Each test validates one or more knowledge points from the nine-level system.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "power_factor.h"
#include "complex_power.h"
#include "pf_correction.h"
#include "harmonic_power.h"
#include "power_quality.h"
#include "phasor_operations.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %-55s", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)
#define CHECK_FEQ(a, b, eps, msg) \
    do { if (fabs((a)-(b)) > (eps)) { FAIL(msg); return; } } while(0)

/* ========================================================================
 * L1: Definitions — Single-Phase Power Quantities
 * ======================================================================== */

static void test_l1_single_phase_pf(void)
{
    TEST("L1: Single-phase PF (unity)");

    pf_single_phase_t result;
    int ret = pf_compute_single_phase(120.0, 10.0, 0.0, &result);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(result.p_real, 1200.0, 1e-6, "P = VI cos(0°) = 1200 W");
    CHECK_FEQ(result.q_reactive, 0.0, 1e-6, "Q = 0 for unity PF");
    CHECK_FEQ(result.s_apparent, 1200.0, 1e-6, "S = VI = 1200 VA");
    CHECK_FEQ(result.pf, 1.0, 1e-6, "PF = 1.0 for resistive load");
    CHECK(result.type == PF_TYPE_UNITY, "type = UNITY");
    PASS();
}

static void test_l1_lagging_pf(void)
{
    TEST("L1: Single-phase PF (lagging/inductive)");

    pf_single_phase_t result;
    /* 60° lagging → PF = cos(60°) = 0.5 */
    int ret = pf_compute_single_phase(240.0, 5.0, 60.0, &result);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(result.p_real, 240.0 * 5.0 * 0.5, 1e-6, "P = VI cos(60°)");
    CHECK(result.q_reactive > 0.0, "Q > 0 for lagging");
    CHECK_FEQ(result.pf, 0.5, 1e-6, "PF = 0.5");
    CHECK(result.type == PF_TYPE_LAGGING, "type = LAGGING");
    PASS();
}

static void test_l1_leading_pf(void)
{
    TEST("L1: Single-phase PF (leading/capacitive)");

    pf_single_phase_t result;
    /* -45° → leading → PF = cos(45°) ≈ 0.707 */
    int ret = pf_compute_single_phase(230.0, 4.0, -45.0, &result);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(result.pf, 0.70710678, 1e-4, "PF = cos(45°) ≈ 0.707");
    CHECK(result.q_reactive < 0.0, "Q < 0 for leading");
    CHECK(result.type == PF_TYPE_LEADING, "type = LEADING");
    PASS();
}

/* ========================================================================
 * L1: Time-Domain PF from Samples
 * ======================================================================== */

static void test_l1_samples_pf(void)
{
    TEST("L1: PF from time-domain samples (resistive)");

    size_t N = 1000;
    double v[1000], i[1000];
    for (size_t n = 0; n < N; n++) {
        double t = 2.0 * M_PI * n / N;
        v[n] = 169.7 * cos(t);        /* 120Vrms, 170V peak */
        i[n] = 14.14 * cos(t);        /* 10Arms, 14.14A peak — in phase */
    }

    pf_single_phase_t result;
    int ret = pf_compute_from_samples(v, i, N, 1.0/60000.0, &result);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(result.pf, 1.0, 0.02, "PF ≈ 1.0 (resistive)");
    CHECK_FEQ(result.v_rms, 120.0, 0.5, "Vrms ≈ 120");
    CHECK_FEQ(result.i_rms, 10.0, 0.2, "Irms ≈ 10");
    CHECK(result.type == PF_TYPE_UNITY, "type = UNITY");
    PASS();
}

static void test_l1_samples_lagging(void)
{
    TEST("L1: PF from time-domain (lagging by phase shift)");

    size_t N = 1000;
    double v[1000], i[1000];
    double phase_rad = 30.0 * M_PI / 180.0; /* 30° lag */
    for (size_t n = 0; n < N; n++) {
        double t = 2.0 * M_PI * n / N;
        v[n] = cos(t);
        i[n] = cos(t - phase_rad);  /* Current lags voltage */
    }

    pf_single_phase_t result;
    int ret = pf_compute_from_samples(v, i, N, 0.001, &result);
    CHECK(ret == 0, "return code");
    /* PF should be ≈ cos(30°) ≈ 0.866 */
    CHECK_FEQ(result.pf, 0.866, 0.05, "PF ≈ cos(30°)");
    PASS();
}

/* ========================================================================
 * L2: PF Classification and Compensation
 * ======================================================================== */

static void test_l2_pf_classify(void)
{
    TEST("L2: PF classification");
    CHECK(pf_classify(0.98, PF_TYPE_UNITY) == PF_CLASS_GOOD, "0.98 → GOOD");
    CHECK(pf_classify(0.90, PF_TYPE_LAGGING) == PF_CLASS_FAIR, "0.90 → FAIR");
    CHECK(pf_classify(0.78, PF_TYPE_LAGGING) == PF_CLASS_POOR, "0.78 → POOR");
    CHECK(pf_classify(0.50, PF_TYPE_LAGGING) == PF_CLASS_BAD, "0.50 → BAD");
    CHECK(pf_classify(0.95, PF_TYPE_LEADING) == PF_CLASS_LEADING, "leading→LEADING");
    PASS();
}

static void test_l2_reactive_compensation(void)
{
    TEST("L2: Required reactive compensation");

    /* Correct 100kW from PF=0.70 to PF=0.95, lagging load */
    double q_c = pf_required_reactive_comp(100000.0, 0.70, 0.95, 1);
    /* Q_c = P × (tan(acos(0.70)) - tan(acos(0.95)))
     *      = 100000 × (1.020 - 0.329) = 100000 × 0.691 = 69100 VAR */
    CHECK(q_c > 60000.0 && q_c < 75000.0, "Q_c in expected range");
    PASS();
}

static void test_l2_savings_estimate(void)
{
    TEST("L2: PF correction savings");
    double savings = pf_savings_estimate(500000.0, 0.70, 0.95,
                                         0.05, 8000.0, 0.12);
    CHECK(savings >= 0.0, "savings non-negative");
    PASS();
}

/* ========================================================================
 * L3: Phasor and Complex Power Operations
 * ======================================================================== */

static void test_l3_phasor_power(void)
{
    TEST("L3: Phasor-based complex power");

    pf_single_phase_t result;
    int ret = pf_phasor_power(120.0, 0.0, 10.0, -M_PI/6.0, &result);
    /* V = 120∠0°, I = 10∠-30° → S = VI* = 1200∠30° = 1039 + j600 */
    CHECK(ret == 0, "return code");
    CHECK_FEQ(result.p_real, 1039.23, 0.5, "P = 120×10×cos(30°)");
    CHECK_FEQ(result.q_reactive, 600.0, 0.5, "Q = 120×10×sin(30°)");
    PASS();
}

static void test_l3_symmetrical_components(void)
{
    TEST("L3: Symmetrical components (balanced)");

    double v0_m, v0_a, v1_m, v1_a, v2_m, v2_a;
    /* Balanced positive-sequence: Va=1∠0°, Vb=1∠-120°, Vc=1∠120° */
    int ret = pf_symmetrical_components(
        1.0, 0.0,
        1.0, -2.0*M_PI/3.0,
        1.0, 2.0*M_PI/3.0,
        &v0_m, &v0_a, &v1_m, &v1_a, &v2_m, &v2_a);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(v0_m, 0.0, 1e-6, "V0 ≈ 0 for balanced");
    CHECK_FEQ(v1_m, 1.0, 1e-6, "V1 = 1 (positive seq magnitude)");
    CHECK_FEQ(v2_m, 0.0, 1e-6, "V2 ≈ 0 for balanced");
    PASS();
}

static void test_l3_symmetrical_unbalanced(void)
{
    TEST("L3: Symmetrical components (unbalanced)");

    double v0_m, v0_a, v1_m, v1_a, v2_m, v2_a;
    /* Unbalanced: Va=1∠0°, Vb=0.8∠-120°, Vc=0.6∠120° */
    int ret = pf_symmetrical_components(
        1.0, 0.0,
        0.8, -2.0*M_PI/3.0,
        0.6, 2.0*M_PI/3.0,
        &v0_m, &v0_a, &v1_m, &v1_a, &v2_m, &v2_a);
    CHECK(ret == 0, "return code");
    CHECK(v0_m > 0.0, "V0 > 0 for unbalanced");
    CHECK(v2_m > 0.0, "V2 > 0 for unbalanced");
    PASS();
}

/* ========================================================================
 * L3: Complex Power Operations
 * ======================================================================== */

static void test_l3_complex_power(void)
{
    TEST("L3: Complex power from V and I phasors");

    double complex v = 120.0 + I * 0.0;
    double complex i = 8.66 - I * 5.0;  /* 10∠-30° */
    double complex s = cp_complex_power(v, i);

    CHECK_FEQ(creal(s), 1039.2, 1.0, "P = 120×8.66");
    CHECK_FEQ(cimag(s), 600.0, 1.0, "Q = 120×5.0");
    PASS();
}

static void test_l3_impedance_from_power(void)
{
    TEST("L3: Load impedance from P and Q");

    /* P=1200W, Q=600VAR, V=120V → Z = V²/(P-jQ) */
    double complex Z = cp_impedance_from_power(120.0, 1200.0, 600.0);
    /* Expected: Z = 9.6 + j4.8 Ohm */
    CHECK_FEQ(creal(Z), 9.6, 0.1, "R ≈ 9.6 Ohm");
    CHECK_FEQ(cimag(Z), 4.8, 0.1, "X ≈ 4.8 Ohm");
    PASS();
}

/* ========================================================================
 * L4: Energy Conservation Theorems
 * ======================================================================== */

static void test_l4_power_balance(void)
{
    TEST("L4: Real power conservation");

    double sources[] = {5000.0, 3000.0, 2000.0};
    double loads[]   = {4000.0, 3500.0, 2500.0};
    int ret = pf_verify_power_balance(sources, 3, loads, 3, 1e-6);
    CHECK(ret == 0, "sources=loads → balanced");
    PASS();
}

static void test_l4_power_imbalance(void)
{
    TEST("L4: Real power imbalance detection");

    double sources[] = {5000.0, 3000.0};
    double loads[]   = {4000.0, 3000.0};
    int ret = pf_verify_power_balance(sources, 2, loads, 2, 1e-6);
    CHECK(ret == 1, "5000+3000 ≠ 4000+3000 → imbalanced");
    PASS();
}

static void test_l4_boucherot(void)
{
    TEST("L4: Boucherot's theorem");

    /* Inductive + capacitive: Q_L=+500, Q_C=-500 → ΣQ=0 */
    double q_vals[] = {500.0, -300.0, -200.0, 150.0, -150.0};
    int ret = pf_verify_boucherot(q_vals, 5, 1e-6);
    CHECK(ret == 0, "ΣQ = 0");
    PASS();
}

static void test_l4_complex_power_balance(void)
{
    TEST("L4: Complex power balance (Steinmetz)");

    /* S1 + S2 = 0 → P=0, Q=0 both checked */
    double complex s_vals[] = {
        100.0 + I * 50.0,
        -60.0 - I * 30.0,
        -40.0 - I * 20.0
    };
    int ret = cp_verify_complex_power_balance(s_vals, 3, 1e-9);
    CHECK(ret == 0, "ΣS = 0+j0");
    PASS();
}

/* ========================================================================
 * L5: RMS and Statistical Methods
 * ======================================================================== */

static void test_l5_sliding_rms(void)
{
    TEST("L5: Sliding window RMS");

    double samples[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double rms_out[4];
    int ret = pf_sliding_rms(samples, 5, 2, rms_out);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(rms_out[0], 1.0, 1e-6, "RMS of [1,1] = 1");
    PASS();
}

static void test_l5_ema_rms(void)
{
    TEST("L5: EMA RMS");

    double samples[] = {2.0, 2.0, 2.0, 2.0, 2.0};
    double rms_out[5];
    int ret = pf_ema_rms(samples, 5, 10, rms_out);
    CHECK(ret == 0, "return code");
    /* After convergence, EMA = 2.0 */
    CHECK_FEQ(rms_out[4], 2.0, 0.5, "EMA converges to 2.0");
    PASS();
}

static void test_l5_crest_factor(void)
{
    TEST("L5: Crest factor — pure sine");

    /* Pure sine: CF = √2 ≈ 1.414 */
    size_t N = 1000;
    double sine[1000];
    for (size_t n = 0; n < N; n++) {
        sine[n] = sin(2.0 * M_PI * n / N);
    }
    double cf = pf_crest_factor(sine, N);
    CHECK_FEQ(cf, 1.414, 0.02, "CF ≈ √2 for pure sine");
    PASS();
}

static void test_l5_form_factor(void)
{
    TEST("L5: Form factor — pure sine");

    size_t N = 1000;
    double sine[1000];
    for (size_t n = 0; n < N; n++) {
        sine[n] = sin(2.0 * M_PI * n / N);
    }
    double ff = pf_form_factor(sine, N);
    /* FF = π/(2√2) ≈ 1.111 */
    CHECK_FEQ(ff, 1.111, 0.02, "FF ≈ 1.111 for pure sine");
    PASS();
}

/* ========================================================================
 * L5: Capacitor Sizing and PFC Algorithms
 * ======================================================================== */

static void test_l5_capacitor_sizing(void)
{
    TEST("L5: Single-phase capacitor sizing");

    double c, kvar;
    int ret = pfc_size_capacitor(100000.0, 480.0, 60.0, 0.70, 0.95, 1,
                                  &c, &kvar);
    CHECK(ret == 0, "return code");
    CHECK(c > 0.0, "C > 0");
    CHECK(kvar > 50000.0, "kVAR > 50k");
    PASS();
}

static void test_l5_3phase_capacitor(void)
{
    TEST("L5: Three-phase capacitor sizing");

    double c_per_phase, kvar_total;
    int ret = pfc_size_capacitor_3phase(500000.0, 480.0, 60.0,
                                         0.65, 0.95, 1, 0,
                                         &c_per_phase, &kvar_total);
    CHECK(ret == 0, "return code");
    CHECK(c_per_phase > 0.0, "C_per_phase > 0 (delta)");
    CHECK(kvar_total > 200000.0, "kvar_total > 200k for 500kW, pf 0.65→0.95");
    PASS();
}

static void test_l5_step_bank_design(void)
{
    TEST("L5: Capacitor step bank design");

    pfc_capacitor_bank_t bank;
    int ret = pfc_design_step_bank(300.0, 4, 50.0, &bank);
    CHECK(ret == 0, "return code");
    CHECK(bank.num_steps == 4, "4 steps configured");
    CHECK(bank.total_kvar > 0.0, "total_kvar > 0");
    PASS();
}

static void test_l5_boost_pfc_init(void)
{
    TEST("L5: Boost PFC initialization");

    pfc_boost_state_t state;
    int ret = pfc_boost_init(&state, 120.0, 400.0, 1000.0, 100000.0, 60.0);
    CHECK(ret == 0, "return code");
    CHECK(state.l_boost > 0.0, "inductor sized");
    CHECK(state.c_out > 0.0, "output cap sized");
    CHECK(state.kp_current > 0.0, "current Kp > 0");
    PASS();
}

static void test_l5_boost_pfc_step(void)
{
    TEST("L5: Boost PFC control step");

    pfc_boost_state_t state;
    pfc_boost_init(&state, 120.0, 400.0, 1000.0, 100000.0, 60.0);

    double v_in = 170.0; /* peak rectified */
    double i_meas = 8.33; /* ~1kW/120V */
    int ret = pfc_boost_control_step(&state, v_in, i_meas, 1.0/100000.0);
    CHECK(ret == 0, "return code");
    CHECK(state.duty_cycle >= 0.0 && state.duty_cycle <= 1.0,
          "duty ∈ [0,1]");
    PASS();
}

/* ========================================================================
 * L5: Harmonic Analysis
 * ======================================================================== */

static void test_l5_thd(void)
{
    TEST("L5: THD computation");

    double harmonics[] = {100.0, 10.0, 5.0, 2.0, 1.0};
    double thd = hp_compute_thd(harmonics, 5);
    /* V_H = sqrt(10²+5²+2²+1²) = sqrt(100+25+4+1) = sqrt(130) = 11.4
     * THD = 11.4/100 × 100% = 11.4% */
    CHECK_FEQ(thd, 11.4, 0.2, "THD ≈ 11.4%");
    PASS();
}

static void test_l5_tdd(void)
{
    TEST("L5: TDD computation");

    double i_harm[] = {100.0, 8.0, 5.0, 3.0};
    double tdd = hp_compute_tdd(i_harm, 4, 120.0);
    /* I_H = sqrt(8²+5²+3²) = sqrt(64+25+9) = sqrt(98) = 9.90
     * TDD = 9.90/120 × 100% = 8.25% */
    CHECK(tdd > 7.0 && tdd < 10.0, "TDD reasonable");
    PASS();
}

static void test_l5_fft_radix2(void)
{
    TEST("L5: FFT radix-2");

    /* FFT of [1,0,0,0] should give all ones (DC component spread equally) */
    double complex data[4] = {
        1.0 + I * 0.0, 0.0 + I * 0.0,
        0.0 + I * 0.0, 0.0 + I * 0.0
    };
    int ret = hp_fft_radix2(data, 4, 0);
    CHECK(ret == 0, "return code");
    /* After FFT, all bins = 0.25? Actually standard DFT: X[k] = 1 for all k */
    CHECK_FEQ(cabs(data[0]), 1.0, 1e-6, "DC bin = 1");
    PASS();
}

static void test_l5_goertzel(void)
{
    TEST("L5: Goertzel single-bin detection");

    /* 60Hz sine sampled at 600Hz, N=100 */
    double samples[100];
    for (int i = 0; i < 100; i++) {
        samples[i] = sin(2.0 * M_PI * 60.0 * i / 600.0);
    }

    double mag, phase;
    int ret = hp_goertzel_detect(samples, 100, 60.0, 600.0, &mag, &phase);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(mag, 1.0, 0.1, "magnitude ≈ 1 for 60Hz bin");
    PASS();
}

/* ========================================================================
 * L4: Parseval's Theorem
 * ======================================================================== */

static void test_l4_parseval(void)
{
    TEST("L4: Parseval power verification");

    double time_power = 100.0;
    double harmonic_p[] = {80.0, 15.0, 4.0, 1.0}; /* sum = 100 */
    int ret = hp_verify_parseval_power(time_power, harmonic_p, 4, 1e-9);
    CHECK(ret == 0, "Parseval holds");
    PASS();
}

/* ========================================================================
 * L5: K-Factor and Transformer Derating
 * ======================================================================== */

static void test_l5_k_factor(void)
{
    TEST("L5: K-factor computation");

    double i_harm[] = {100.0, 20.0, 10.0, 5.0};
    double k = hp_k_factor(i_harm, 4);
    /* K = 1 + 2²×(0.2)² + 3²×(0.1)² + 4²×(0.05)²
     *   = 1 + 4×0.04 + 9×0.01 + 16×0.0025
     *   = 1 + 0.16 + 0.09 + 0.04 = 1.29 */
    CHECK(k > 1.0, "K > 1 with harmonics");
    CHECK(k < 3.0, "K < 3 for moderate harmonics");
    PASS();
}

static void test_l5_transformer_derating(void)
{
    TEST("L5: Transformer derating");

    double derate = hp_transformer_derating(4.0, 0.08);
    CHECK(derate > 0.6 && derate < 1.0, "derating 0.6-1.0 for K=4");
    PASS();
}

/* ========================================================================
 * L6: Harmonic Spectrum Analysis
 * ======================================================================== */

static void test_l6_spectrum_analysis(void)
{
    TEST("L6: Full harmonic spectrum analysis");

    size_t N = 1024;
    double v[1024], i[1024];
    /* 60Hz fundamental + 5% 5th harmonic */
    for (size_t n = 0; n < N; n++) {
        double t = (double)n / 1024.0;
        v[n] = 169.7 * (cos(2.0 * M_PI * 60.0 * t / 60.0)
               + 0.03 * cos(2.0 * M_PI * 300.0 * t / 60.0));
        i[n] = 14.14 * (cos(2.0 * M_PI * 60.0 * t / 60.0 - 0.4636)
               + 0.03 * cos(2.0 * M_PI * 300.0 * t / 60.0 - 0.4636));
    }

    hp_spectrum_t spectrum;
    int ret = hp_analyze_spectrum(v, i, N, 1024.0, 60.0, &spectrum);
    CHECK(ret == 0, "return code");
    CHECK(spectrum.num_harmonics > 0, "harmonics detected");
    CHECK(spectrum.thd_v_percent >= 0.0, "THDv ≥ 0");
    PASS();
}

/* ========================================================================
 * L6: Industrial PF Correction
 * ======================================================================== */

static void test_l6_industrial_pfc(void)
{
    TEST("L6: Industrial PF correction solution");

    double c, kvar, payback, f_res;
    int ret = pfc_solve_industrial(500000.0, 480.0, 60.0,
                                   0.70, 0.95, 5000.0,
                                   &c, &kvar, &payback, &f_res);
    CHECK(ret == 0, "return code");
    CHECK(c > 0.0, "C_found");
    CHECK(kvar > 100000.0, "kVAR > 100k for 500kW");
    CHECK(f_res > 0.0, "f_res computed");
    PASS();
}

/* ========================================================================
 * L6: Harmonic Resonance Risk
 * ======================================================================== */

static void test_l6_harmonic_risk(void)
{
    TEST("L6: Harmonic resonance risk assessment");

    double mag;
    int ret = pfc_harmonic_risk(10000000.0, 500000.0, 60.0, 5, &mag);
    CHECK(ret >= 0, "risk assessment returned");
    CHECK(mag > 0.0, "magnification computed");
    PASS();
}

/* ========================================================================
 * L2: ITIC Curve
 * ======================================================================== */

static void test_l2_itic_curve(void)
{
    TEST("L2: ITIC (CBEMA) curve evaluation");

    /* Normal operation: 100% voltage → acceptable */
    pq_itic_result_t r = pq_check_itic(1.0, 100.0);
    CHECK(r == PQ_ITIC_ACCEPTABLE, "100%@1s → acceptable");

    /* Severe sag: 50% for 0.5s → below lower envelope */
    r = pq_check_itic(0.5, 50.0);
    CHECK(r == PQ_ITIC_BELOW_LOWER, "50%@0.5s → below lower");

    /* Severe overvoltage: 600% for 0.5ms → above upper (500% limit) */
    r = pq_check_itic(0.0005, 600.0);
    CHECK(r == PQ_ITIC_ABOVE_UPPER, "600%@0.5ms → above upper");

    /* Moderate overvoltage: 150% for 0.8ms is acceptable (limit is 500%) */
    r = pq_check_itic(0.0008, 150.0);
    CHECK(r == PQ_ITIC_ACCEPTABLE, "150%@0.8ms → acceptable");
    PASS();
}

/* ========================================================================
 * L7: Data Center and EV Charger PQ
 * ======================================================================== */

static void test_l7_datacenter_pq(void)
{
    TEST("L7: Data center PQ assessment");

    int tier; double score;
    int issues = pq_assess_datacenter(0.97, 3.0, 5.0, 96.0, &tier, &score);
    CHECK(issues == 0, "good datacenter: no issues");
    CHECK(tier == 4, "Titanium tier");
    PASS();
}

static void test_l7_ev_charger_pq(void)
{
    TEST("L7: EV charger PQ compliance");

    int compliant;
    int ret = pq_assess_ev_charger(0.96, 4.0, 80.0, 200.0, &compliant);
    CHECK(ret == 0, "return code");
    CHECK(compliant == 1, "compliant EV charger");
    PASS();
}

/* ========================================================================
 * L3: Clarke and Park Transforms
 * ======================================================================== */

static void test_l3_clarke_transform(void)
{
    TEST("L3: Clarke transform (balanced ABC)");

    double alpha, beta, zero;
    /* Balanced 3-phase: A=1, B=-0.5, C=-0.5 at t=0 */
    int ret = phasor_clarke_transform(1.0, -0.5, -0.5, &alpha, &beta, &zero);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(alpha, 1.0, 1e-6, "α = A for balanced");
    CHECK_FEQ(zero, 0.0, 1e-6, "zero seq = 0 for balanced");
    CHECK_FEQ(beta, 0.0, 1e-6, "β = 0 at t=0");
    PASS();
}

static void test_l3_park_transform(void)
{
    TEST("L3: Park transform (balanced ABC at θ=0)");

    dq0_t dq;
    int ret = phasor_park_transform(1.0, -0.5, -0.5, 0.0, &dq);
    CHECK(ret == 0, "return code");
    CHECK_FEQ(dq.d, 1.0, 1e-6, "Vd = V_peak at θ=0");
    CHECK_FEQ(dq.q, 0.0, 1e-6, "Vq = 0 at θ=0");
    PASS();
}

/* ========================================================================
 * L5: SRF-PLL Grid Synchronization
 * ======================================================================== */

static void test_l5_srf_pll(void)
{
    TEST("L5: SRF-PLL grid synchronization");

    srf_pll_t pll;
    int ret = srf_pll_init(&pll, 2.0 * M_PI * 60.0, 10.0, 1.0/10000.0);
    CHECK(ret == 0, "PLL initialized");

    /* Simulate stepping with balanced 3-phase */
    double omega = 2.0 * M_PI * 60.0;
    for (int i = 0; i < 500; i++) {
        double t = i * 1.0/10000.0;
        double va = cos(omega * t);
        double vb = cos(omega * t - 2.0*M_PI/3.0);
        double vc = cos(omega * t + 2.0*M_PI/3.0);
        ret = srf_pll_step(&pll, va, vb, vc);
        CHECK(ret == 0, "PLL step");
    }

    double freq = srf_pll_get_frequency(&pll);
    CHECK_FEQ(freq, omega, 0.5, "PLL locks to 377 rad/s");

    double mag = srf_pll_get_magnitude(&pll);
    CHECK_FEQ(mag, 1.0, 0.1, "PLL measures magnitude ≈ 1.0");
    PASS();
}

/* ========================================================================
 * L3: Phasor Arithmetic
 * ======================================================================== */

static void test_l3_phasor_arithmetic(void)
{
    TEST("L3: Phasor arithmetic operations");

    phasor_t a = phasor_from_polar(10.0, M_PI / 4.0);  /* 10∠45° */
    phasor_t b = phasor_from_polar(5.0, M_PI / 4.0);   /* 5∠45° */
    phasor_t sum = phasor_add(&a, &b);
    CHECK_FEQ(sum.magnitude, 15.0, 1e-6, "|10+5| = 15 (in phase)");

    phasor_t prod = phasor_mul(&a, &b);
    CHECK_FEQ(prod.magnitude, 50.0, 1e-6, "|10×5| = 50");
    CHECK_FEQ(prod.angle_rad, M_PI / 2.0, 1e-6, "∠45°+45° = 90°");
    PASS();
}

/* ========================================================================
 * L4: Flicker and Reliability Indices
 * ======================================================================== */

static void test_l4_flicker(void)
{
    TEST("L4: Flicker severity Pst");

    double v_samples[600];
    double v_nom = 120.0;
    for (int i = 0; i < 600; i++) {
        v_samples[i] = v_nom + 0.5 * sin(2.0 * M_PI * 8.8 * i / 600.0);
    }
    double pst = pq_compute_pst(v_samples, 600, v_nom);
    CHECK(pst >= 0.0, "Pst computed");
    PASS();
}

static void test_l4_saifi(void)
{
    TEST("L4: SAIFI reliability index");

    uint32_t customers_per_event[] = {500, 300, 200};
    double saifi = pq_compute_saifi(customers_per_event, 3, 10000);
    /* SAIFI = (500+300+200)/10000 = 0.1 interruptions/customer */
    CHECK_FEQ(saifi, 0.1, 1e-6, "SAIFI = 1000/10000");
    PASS();
}

/* ========================================================================
 * L4: Complex Power Conservation
 * ======================================================================== */

static void test_l4_series_parallel_power(void)
{
    TEST("L4: Series and parallel complex power");

    /* Series: 2 elements, I=10∠0°, Z1=3+j4, Z2=4+j3 */
    double complex z[] = { 3.0 + I * 4.0, 4.0 + I * 3.0 };
    double complex i_phasor = 10.0 + I * 0.0;
    double complex s_total;
    int ret = cp_series_complex_power(z, 2, i_phasor, &s_total);
    CHECK(ret == 0, "return code");
    /* S_total = |I|² × (Z1+Z2) = 100 × (7+j7) = 700+j700 */
    CHECK_FEQ(creal(s_total), 700.0, 1.0, "P_series = 700W");
    CHECK_FEQ(cimag(s_total), 700.0, 1.0, "Q_series = 700VAR");
    PASS();
}

/* ========================================================================
 * L8: Advanced — Monte Carlo PF uncertainty estimation
 * ======================================================================== */

static void test_l8_pf_uncertainty_monte_carlo(void)
{
    TEST("L8: Monte Carlo PF uncertainty (conceptual)");

    /* Demonstrate PF variation under measurement noise and parameter
     * uncertainty using a simplified Monte Carlo approach.
     *
     * This is a stochastic method: draw random samples from the
     * measurement error distribution and compute PF distribution.
     *
     * Ref: JCGM 100:2008 (GUM), "Evaluation of Measurement Data"
     */
    size_t num_trials = 100;
    double pf_sum = 0.0, pf_min = 2.0, pf_max = -1.0;
    /* Nominal: V=120, I=10, φ=30° → PF=0.866 */
    for (size_t t = 0; t < num_trials; t++) {
        /* Add ±1% Gaussian-like noise (approximated) */
        double noise_v = (double)((int)(t * 1234567) % 200 - 100) / 10000.0;
        double noise_i = (double)((int)(t * 7654321) % 200 - 100) / 10000.0;
        double v = 120.0 * (1.0 + noise_v);
        double i = 10.0  * (1.0 + noise_i);

        pf_single_phase_t res;
        pf_compute_single_phase(v, i, 30.0, &res);
        pf_sum += res.pf;
        if (res.pf < pf_min) pf_min = res.pf;
        if (res.pf > pf_max) pf_max = res.pf;
    }
    double pf_mean = pf_sum / (double)num_trials;
    CHECK_FEQ(pf_mean, 0.866, 0.02, "MC mean PF ≈ 0.866");
    CHECK(pf_max - pf_min < 0.05, "PF spread small with 1% noise");
    PASS();
}

/* ========================================================================
 * Main
 * ======================================================================== */

int main(void)
{
    printf("\n========================================\n");
    printf("  mini-power-factor Test Suite\n");
    printf("========================================\n\n");

    /* L1: Definitions */
    printf("[L1] Core Definitions\n");
    test_l1_single_phase_pf();
    test_l1_lagging_pf();
    test_l1_leading_pf();
    test_l1_samples_pf();
    test_l1_samples_lagging();

    /* L2: Core Concepts */
    printf("\n[L2] Core Concepts\n");
    test_l2_pf_classify();
    test_l2_reactive_compensation();
    test_l2_savings_estimate();
    test_l2_itic_curve();

    /* L3: Mathematical Structures */
    printf("\n[L3] Mathematical Structures\n");
    test_l3_phasor_power();
    test_l3_symmetrical_components();
    test_l3_symmetrical_unbalanced();
    test_l3_complex_power();
    test_l3_impedance_from_power();
    test_l3_phasor_arithmetic();
    test_l3_clarke_transform();
    test_l3_park_transform();

    /* L4: Fundamental Laws */
    printf("\n[L4] Fundamental Laws\n");
    test_l4_power_balance();
    test_l4_power_imbalance();
    test_l4_boucherot();
    test_l4_complex_power_balance();
    test_l4_parseval();
    test_l4_series_parallel_power();
    test_l4_flicker();
    test_l4_saifi();

    /* L5: Algorithms */
    printf("\n[L5] Algorithms and Methods\n");
    test_l5_sliding_rms();
    test_l5_ema_rms();
    test_l5_crest_factor();
    test_l5_form_factor();
    test_l5_capacitor_sizing();
    test_l5_3phase_capacitor();
    test_l5_step_bank_design();
    test_l5_boost_pfc_init();
    test_l5_boost_pfc_step();
    test_l5_thd();
    test_l5_tdd();
    test_l5_fft_radix2();
    test_l5_goertzel();
    test_l5_k_factor();
    test_l5_transformer_derating();
    test_l5_srf_pll();

    /* L6: Canonical Problems */
    printf("\n[L6] Canonical Problems\n");
    test_l6_spectrum_analysis();
    test_l6_industrial_pfc();
    test_l6_harmonic_risk();

    /* L7: Applications */
    printf("\n[L7] Applications\n");
    test_l7_datacenter_pq();
    test_l7_ev_charger_pq();

    /* L8: Advanced Topics */
    printf("\n[L8] Advanced Topics\n");
    test_l8_pf_uncertainty_monte_carlo();

    printf("\n========================================\n");
    printf("  Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
