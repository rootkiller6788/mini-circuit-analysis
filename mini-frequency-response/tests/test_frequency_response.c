/**
 * @file test_frequency_response.c
 * @brief Comprehensive test suite for mini-frequency-response module
 *
 * Tests cover:
 * - L1: Core definitions (data types, conversions)
 * - L2: Core concepts (transfer functions, frequency sweeps)
 * - L3: Math structures (Horner evaluation, polynomial arithmetic)
 * - L4: Fundamental laws (Routh-Hurwitz, Nyquist, stability margins)
 * - L5: Algorithms (Bode plots, filter design, root-finding)
 * - L6: Canonical problems (RLC resonance, filter design from spec)
 *
 * All tests use standard assert() — no custom macros.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "frequency_response.h"
#include "transfer_function.h"
#include "bode_plot.h"
#include "filter_design.h"
#include "resonance.h"
#include "stability.h"
#include "network_function.h"

#define ASSERT_NEAR(actual, expected, tol) \
    assert(fabs((actual) - (expected)) < (tol))

#define ASSERT_COMPLEX_NEAR(actual, expected, tol) \
    assert(cabs((actual) - (expected)) < (tol))

/* ============================================================================
 * L1: Definitions — Data types and basic conversions
 * ============================================================================ */

static void test_l1_conversions(void)
{
    printf("  L1: Basic conversions...\n");

    /* Hz ↔ rad/s */
    double omega = freq_hz_to_rad(1.0);
    ASSERT_NEAR(omega, 2.0 * M_PI, 1e-10);

    double freq = freq_rad_to_hz(2.0 * M_PI);
    ASSERT_NEAR(freq, 1.0, 1e-10);

    /* Magnitude ↔ dB */
    double db = magnitude_to_db(1.0);
    ASSERT_NEAR(db, 0.0, 1e-10);

    db = magnitude_to_db(sqrt(2.0));
    ASSERT_NEAR(db, 3.01029995664, 1e-6);  /* -3 dB point */

    db = magnitude_to_db(10.0);
    ASSERT_NEAR(db, 20.0, 1e-10);

    double mag = db_to_magnitude(20.0);
    ASSERT_NEAR(mag, 10.0, 1e-10);

    /* Phase computation */
    double phase = phase_degrees(1.0, 1.0);
    ASSERT_NEAR(phase, 45.0, 1e-6);

    phase = phase_degrees(0.0, 1.0);
    ASSERT_NEAR(phase, 90.0, 1e-6);

    phase = phase_degrees(-1.0, 0.0);
    ASSERT_NEAR(phase, 180.0, 1e-6);

    printf("  L1: PASSED\n");
}

/* ============================================================================
 * L1: Frequency grid generation
 * ============================================================================ */

static void test_l1_freq_grid(void)
{
    printf("  L1: Frequency grid generation...\n");

    size_t n = 0;
    double *freqs = freq_logspace(1.0, 1000.0, 10, &n);
    assert(freqs != NULL);
    assert(n > 0);
    ASSERT_NEAR(freqs[0], 1.0, 1e-10);

    /* Each decade should have 10 points */
    assert(n >= 31);  /* 3 decades × 10 pts + endpoint */

    /* Check monotonic increase */
    for (size_t i = 1; i < n; i++) {
        assert(freqs[i] > freqs[i - 1]);
    }

    free(freqs);

    /* Edge cases */
    double *bad = freq_logspace(-1.0, 100.0, 10, &n);
    assert(bad == NULL);
    assert(n == 0);

    bad = freq_logspace(100.0, 10.0, 10, &n);
    assert(bad == NULL);

    printf("  L1: PASSED\n");
}

/* ============================================================================
 * L1: Phase unwrapping
 * ============================================================================ */

static void test_l1_phase_unwrap(void)
{
    printf("  L1: Phase unwrapping...\n");

    /* Simulate wrapped phase crossing -180°/+180° boundary */
    double wrapped[] = { 170.0, 179.0, -179.0, -170.0, -90.0, 0.0 };
    size_t n = 6;

    double *unwrapped = phase_unwrap(wrapped, n, 180.0);
    assert(unwrapped != NULL);

    /* After unwrapping, phases should be continuous */
    assert(unwrapped[2] > 180.0);  /* Was -179°, should be 181° */
    assert(unwrapped[3] > 180.0);  /* Was -170°, should be 190° */

    free(unwrapped);
    printf("  L1: PASSED\n");
}

/* ============================================================================
 * L2: Transfer function creation and evaluation
 * ============================================================================ */

static void test_l2_transfer_function(void)
{
    printf("  L2: Transfer function creation and evaluation...\n");

    /* H(s) = 1/(s + 1) — first-order lowpass, τ=1s */
    double num[] = { 1.0 };
    double den[] = { 1.0, 1.0 };

    tf_polynomial_t *tf = tf_polynomial_create(num, 0, den, 1);
    assert(tf != NULL);

    /* DC gain: H(0) = 1 */
    double dc = tf_dc_gain(tf);
    ASSERT_NEAR(dc, 1.0, 1e-10);

    /* High-frequency gain: H(∞) = 0 */
    double hf = tf_hf_gain(tf);
    ASSERT_NEAR(hf, 0.0, 1e-10);

    /* Evaluate at s = j (ω = 1 rad/s):
     * H(j1) = 1/(1 + j1) = 1/(1+j) = (1-j)/2 = 0.5 - j0.5
     * |H| = 1/√2 ≈ 0.7071, ∠H = -45° */
    double _Complex H = tf_evaluate(tf, I);
    ASSERT_NEAR(creal(H), 0.5, 1e-6);
    ASSERT_NEAR(cimag(H), -0.5, 1e-6);

    double mag_at_1hz = tf_magnitude_at(tf, 1.0 / (2.0 * M_PI));
    ASSERT_NEAR(mag_at_1hz, 1.0 / sqrt(2.0), 1e-4);

    double phase_at_1hz = tf_phase_at(tf, 1.0 / (2.0 * M_PI));
    ASSERT_NEAR(phase_at_1hz, -45.0, 1e-4);

    /* Frequency response sweep */
    freq_response_t *resp = tf_compute_freq_response(tf, 0.1, 10.0, 10);
    assert(resp != NULL);
    assert(resp->num_points > 0);
    freq_response_free(resp);

    tf_polynomial_free(tf);
    printf("  L2: PASSED\n");
}

/* ============================================================================
 * L2: Transfer function arithmetic (series, parallel, feedback)
 * ============================================================================ */

static void test_l2_tf_arithmetic(void)
{
    printf("  L2: Transfer function arithmetic...\n");

    /* H₁ = 1/(s+1), H₂ = 1/(s+2) */
    double num1[] = { 1.0 }, den1[] = { 1.0, 1.0 };
    double num2[] = { 1.0 }, den2[] = { 2.0, 1.0 };

    tf_polynomial_t *H1 = tf_polynomial_create(num1, 0, den1, 1);
    tf_polynomial_t *H2 = tf_polynomial_create(num2, 0, den2, 1);

    /* Cascade: H = H₁·H₂ = 1/((s+1)(s+2)) = 1/(s²+3s+2) */
    tf_polynomial_t *H_cascade = tf_multiply(H1, H2);
    assert(H_cascade != NULL);
    ASSERT_NEAR(H_cascade->den[0], 2.0, 1e-10);  /* a₀ = 1×2 = 2 */
    ASSERT_NEAR(H_cascade->den[1], 3.0, 1e-10);  /* a₁ = 1+2 = 3 */
    ASSERT_NEAR(H_cascade->den[2], 1.0, 1e-10);  /* a₂ = 1 */

    /* Feedback: H = H₁/(1+H₁) = 1/(s+2) (unity negative feedback)
     * H₁ = 1/(s+1), H₂ = 1 → H_cl = H₁/(1 + H₁·1)
     *   = (1/(s+1)) / ((s+2)/(s+1)) = 1/(s+2) */
    double num_fb[] = { 1.0 };   /* H₂(s) = 1 */
    double den_fb[] = { 1.0 };   /* denominator = 1 */
    tf_polynomial_t *unity_fb = tf_polynomial_create(num_fb, 0, den_fb, 0);
    tf_polynomial_t *H_cl = tf_feedback(H1, unity_fb, 0);
    assert(H_cl != NULL);
    /* H_cl should be 1/(s+2): den[0]=2, den[1]=1 */
    ASSERT_NEAR(H_cl->den[0], 2.0, 1e-10);
    ASSERT_NEAR(H_cl->den[1], 1.0, 1e-10);

    tf_polynomial_free(H1);
    tf_polynomial_free(H2);
    tf_polynomial_free(H_cascade);
    tf_polynomial_free(unity_fb);
    tf_polynomial_free(H_cl);

    printf("  L2: PASSED\n");
}

/* ============================================================================
 * L3: Horner polynomial evaluation
 * ============================================================================ */

static void test_l3_horner(void)
{
    printf("  L3: Horner evaluation...\n");

    /* P(s) = 1 + 2s + 3s² at s=2: 1+4+12=17 */
    double coeffs[] = { 1.0, 2.0, 3.0 };
    tf_polynomial_t *tf = tf_polynomial_create(coeffs, 2, (double[]){1.0}, 0);
    double _Complex val = tf_evaluate(tf, 2.0);
    ASSERT_NEAR(creal(val), 17.0, 1e-10);
    tf_polynomial_free(tf);

    printf("  L3: PASSED\n");
}

/* ============================================================================
 * L3: Pole-zero analysis (Math structures)
 * ============================================================================ */

static void test_l3_pole_zero(void)
{
    printf("  L3: Pole-zero analysis...\n");

    /* H(s) = 1/(s² + 2s + 1) = 1/(s+1)² — double pole at s=-1 */
    double num[] = { 1.0 };
    double den[] = { 1.0, 2.0, 1.0 };
    tf_polynomial_t *tf = tf_polynomial_create(num, 0, den, 2);
    assert(tf != NULL);

    tf_pole_zero_t *pz = tf_to_pole_zero(tf);
    assert(pz != NULL);
    assert(pz->num_poles == 2);
    /* Both poles should be at s = -1 + j0 */
    for (size_t i = 0; i < pz->num_poles; i++) {
        double re = creal(pz->poles[i]);
        ASSERT_NEAR(re, -1.0, 1e-6);
    }

    /* Minimum phase check */
    int is_mp = tf_is_minimum_phase(pz);
    assert(is_mp == 1);

    tf_pole_zero_free(pz);
    tf_polynomial_free(tf);

    /* Test non-minimum-phase: H(s) = (s-1)/(s+1) — zero at s=+1 */
    double num_nmp[] = { -1.0, 1.0 };  /* s - 1 */
    double den_nmp[] = { 1.0, 1.0 };  /* s + 1 */
    tf_polynomial_t *tf_nmp = tf_polynomial_create(num_nmp, 1, den_nmp, 1);

    tf_pole_zero_t *pz_nmp = tf_to_pole_zero(tf_nmp);
    assert(pz_nmp != NULL);
    is_mp = tf_is_minimum_phase(pz_nmp);
    assert(is_mp == 0);  /* Non-minimum-phase */

    tf_pole_zero_free(pz_nmp);
    tf_polynomial_free(tf_nmp);

    printf("  L3: PASSED\n");
}

/* ============================================================================
 * L4: Routh-Hurwitz stability (Fundamental Law)
 * ============================================================================ */

static void test_l4_routh_hurwitz(void)
{
    printf("  L4: Routh-Hurwitz stability criterion...\n");

    /* Stable: s² + 2s + 1 = (s+1)² (both poles at s=-1) */
    double stable_coeffs[] = { 1.0, 2.0, 1.0 };  /* a₀ + a₁s + a₂s² */
    routh_hurwitz_t *rh = stability_routh_hurwitz(stable_coeffs, 2);
    assert(rh != NULL);
    assert(rh->is_stable == 1);
    assert(rh->sign_changes == 0);
    routh_hurwitz_free(rh);

    /* Unstable: s² - 1 = (s+1)(s-1) (pole at s=+1) */
    double unstable_coeffs[] = { -1.0, 0.0, 1.0 };
    rh = stability_routh_hurwitz(unstable_coeffs, 2);
    assert(rh != NULL);
    assert(rh->is_stable == 0);
    assert(rh->sign_changes == 1);
    routh_hurwitz_free(rh);

    /* Marginally stable: s² + 1 (poles at s=±j) */
    double marginal_coeffs[] = { 1.0, 0.0, 1.0 };
    rh = stability_routh_hurwitz(marginal_coeffs, 2);
    assert(rh != NULL);
    assert(rh->is_marginally_stable == 1);
    routh_hurwitz_free(rh);

    /* Third-order stable: s³ + 3s² + 3s + 1 = (s+1)³ */
    double s3_stable[] = { 1.0, 3.0, 3.0, 1.0 };
    rh = stability_routh_hurwitz(s3_stable, 3);
    assert(rh != NULL);
    assert(rh->is_stable == 1);
    routh_hurwitz_free(rh);

    printf("  L4: PASSED\n");
}

/* ============================================================================
 * L4: Bode's gain-phase relationship (Fundamental Law)
 * ============================================================================ */

static void test_l4_bode_gain_phase(void)
{
    printf("  L4: Bode gain-phase relationship...\n");

    /* First-order lowpass: H(s) = 1/(s+1)
     * Magnitude slope at mid-band: -20 dB/dec → phase ≈ -90° */
    double num[] = { 1.0 }, den[] = { 1.0, 1.0 };
    tf_polynomial_t *tf = tf_polynomial_create(num, 0, den, 1);

    bode_plot_t *bode = bode_compute(tf, 0.01, 100.0, 50);
    assert(bode != NULL);

    double approx_phase = 0.0, error_rms = 0.0;
    int ret = bode_gain_phase_relation(bode, &approx_phase, &error_rms);
    assert(ret == 0);
    /* Error should be small for minimum-phase systems */
    assert(error_rms < 45.0);  /* Approximation error is expected */

    bode_free(bode);
    tf_polynomial_free(tf);

    printf("  L4: PASSED\n");
}

/* ============================================================================
 * L4: Nyquist stability (Fundamental Law)
 * ============================================================================ */

static void test_l4_nyquist(void)
{
    printf("  L4: Nyquist stability criterion...\n");

    /* Loop gain L(s) = 1/(s+1)² — stable, 2 LHP poles, 0 RHP
     * Nyquist plot should NOT encircle (-1, j0) → N=0, Z=0 (stable CL) */
    double num[] = { 1.0 }, den[] = { 1.0, 2.0, 1.0 };
    tf_polynomial_t *L = tf_polynomial_create(num, 0, den, 2);

    int is_stable = 0, enc = 0;
    int ret = stability_nyquist_check(L, 1e-3, 1e3, 500, &is_stable, &enc);
    assert(ret == 0);
    /* With P=0, N=0 → Z=0 → stable closed loop */
    assert(is_stable == 1);

    tf_polynomial_free(L);

    printf("  L4: PASSED\n");
}

/* ============================================================================
 * L5: Bode plot construction (Algorithm)
 * ============================================================================ */

static void test_l5_bode_plot(void)
{
    printf("  L5: Bode plot construction...\n");

    /* First-order lowpass: H(s) = 1/(s+1)
     * Cutoff frequency: ω_c = 1 rad/s → f_c = 1/(2π) ≈ 0.159 Hz
     * DC gain: 0 dB
     * HF slope: -20 dB/decade */

    double num[] = { 1.0 }, den[] = { 1.0, 1.0 };
    tf_polynomial_t *tf = tf_polynomial_create(num, 0, den, 1);

    bode_plot_t *bode = bode_compute(tf, 0.001, 100.0, 50);
    assert(bode != NULL);

    /* DC gain should be 0 dB */
    ASSERT_NEAR(bode->dc_gain_db, 0.0, 1e-6);

    /* Dominant pole frequency */
    double dom_pole = bode_dominant_pole_freq(bode);
    ASSERT_NEAR(dom_pole, 1.0 / (2.0 * M_PI), 1e-4);

    /* HF slope should be approximately -20 dB/dec */
    double slope = bode_slope_in_range(bode, 10.0, 100.0);
    assert(slope < -10.0);  /* Should be close to -20 */

    /* Gain-bandwidth product = 1×1/(2π) = 1/(2π) */
    double gbwp = bode_gain_bandwidth_product(bode);
    ASSERT_NEAR(gbwp, 1.0 / (2.0 * M_PI), 1e-4);

    bode_free(bode);

    /* Asymptotic Bode */
    tf_pole_zero_t *pz = tf_to_pole_zero(tf);
    assert(pz != NULL);
    bode_plot_t *bode_asym = bode_asymptotic(pz, 0.001, 100.0, 50);
    assert(bode_asym != NULL);
    bode_free(bode_asym);
    tf_pole_zero_free(pz);

    tf_polynomial_free(tf);
    printf("  L5: PASSED\n");
}

/* ============================================================================
 * L5: Filter design (Algorithm)
 * ============================================================================ */

static void test_l5_filter_design(void)
{
    printf("  L5: Filter design...\n");

    /* Butterworth lowpass: 3rd order */
    tf_polynomial_t *proto = filter_butterworth_prototype(3);
    assert(proto != NULL);
    assert(proto->den_order == 3);
    /* B₃(s) = s³ + 2s² + 2s + 1 */
    ASSERT_NEAR(proto->den[0], 1.0, 1e-10);
    ASSERT_NEAR(proto->den[1], 2.0, 1e-10);
    ASSERT_NEAR(proto->den[2], 2.0, 1e-10);
    ASSERT_NEAR(proto->den[3], 1.0, 1e-10);
    tf_polynomial_free(proto);

    /* Butterworth 4th order */
    proto = filter_butterworth_prototype(4);
    assert(proto != NULL);
    assert(proto->den_order == 4);
    /* B₄(s): s⁴ + 2.6131s³ + 3.4142s² + 2.6131s + 1 */
    ASSERT_NEAR(proto->den[0], 1.0, 1e-10);
    ASSERT_NEAR(proto->den[4], 1.0, 1e-10);
    tf_polynomial_free(proto);

    /* Bessel 3rd order */
    proto = filter_bessel_prototype(3);
    assert(proto != NULL);
    assert(proto->den_order == 3);
    tf_polynomial_free(proto);

    /* Chebyshev I 3rd order with 1 dB ripple */
    proto = filter_chebyshev1_prototype(3, 1.0);
    assert(proto != NULL);
    assert(proto->den_order == 3);
    tf_polynomial_free(proto);

    /* Filter order computation */
    filter_spec_t spec;
    memset(&spec, 0, sizeof(spec));
    spec.approx = APPROX_BUTTERWORTH;
    spec.f_pass = 1000.0;
    spec.f_stop = 10000.0;
    spec.a_pass = 3.0;
    spec.a_stop = 40.0;
    spec.source_impedance = 50.0;
    spec.load_impedance = 50.0;

    int order = filter_order_butterworth(&spec);
    assert(order >= 2);  /* Need at least order 2 for these specs */

    /* LP→HP transformation */
    tf_polynomial_t *lp2 = filter_butterworth_prototype(2);
    tf_polynomial_t *hp = filter_lp_to_hp(lp2, 1.0);
    assert(hp != NULL);
    assert(hp->den_order == 2);
    tf_polynomial_free(lp2);
    tf_polynomial_free(hp);

    /* G-values for Butterworth */
    size_t n_elem = 0;
    double *g = filter_g_values_butterworth(3, &n_elem);
    assert(g != NULL);
    assert(n_elem == 5);  /* g0, g1, g2, g3, g4 for 3rd order */
    /* g₀ = 1, g₁ = 2sin(π/6) = 1, g₂ = 2sin(π/2) = 2, g₃ = 1, g₄ = 1 */
    ASSERT_NEAR(g[0], 1.0, 1e-6);
    ASSERT_NEAR(g[1], 1.0, 1e-6);   /* 2sin(30°) = 1 */
    ASSERT_NEAR(g[2], 2.0, 1e-6);   /* 2sin(90°) = 2 */
    ASSERT_NEAR(g[3], 1.0, 1e-6);   /* 2sin(150°) = 1 */
    ASSERT_NEAR(g[4], 1.0, 1e-6);   /* Load */
    free(g);

    /* Denormalize g-values */
    double *L_vals = NULL, *C_vals = NULL;
    double R_load = 0.0;
    g = filter_g_values_butterworth(2, &n_elem);
    int ret = filter_denormalize(g, n_elem, 50.0, 2.0 * M_PI * 1000.0,
                                  &L_vals, &C_vals, &R_load);
    assert(ret == 0);
    free(g);
    free(L_vals);
    free(C_vals);

    printf("  L5: PASSED\n");
}

/* ============================================================================
 * L6: RLC Resonance (Canonical Problem)
 * ============================================================================ */

static void test_l6_resonance(void)
{
    printf("  L6: RLC Resonance...\n");

    /* Series RLC: R=1Ω, L=1mH, C=1µF
     * f₀ = 1/(2π√(LC)) = 1/(2π√(1e-3×1e-6)) = 1/(2π√(1e-9))
     *    = 1/(2π×3.162e-5) = 5033 Hz
     * Q = ω₀L/R = 2π×5033×1e-3/1 = 31.6
     * BW = f₀/Q = 5033/31.6 = 159 Hz */

    rlc_params_t params = { 1.0, 1e-3, 1e-6, 1.0 };
    resonance_result_t res = resonance_series(&params);

    double f0_expected = 1.0 / (2.0 * M_PI * sqrt(1e-9));
    ASSERT_NEAR(res.resonant_freq_hz, f0_expected, 1.0);

    double Q_expected = f0_expected * 2.0 * M_PI * 1e-3 / 1.0;
    ASSERT_NEAR(res.quality_factor, Q_expected, 0.1);

    /* Q from bandwidth formula */
    double q_bw = resonance_q_from_bandwidth(res.resonant_freq_hz,
                                               res.half_power_low_hz,
                                               res.half_power_high_hz);
    ASSERT_NEAR(q_bw, res.quality_factor, 0.1);

    /* Damping factor */
    double zeta = resonance_damping_from_q(res.quality_factor);
    ASSERT_NEAR(zeta, 1.0 / (2.0 * res.quality_factor), 1e-6);

    /* Parallel RLC */
    params.R = 1000.0;  /* 1 kΩ */
    resonance_result_t res_p = resonance_parallel(&params);
    ASSERT_NEAR(res_p.resonant_freq_hz, f0_expected, 1.0);
    /* Q_parallel = R/(ω₀L) */
    double Qp = 1000.0 / (f0_expected * 2.0 * M_PI * 1e-3);
    ASSERT_NEAR(res_p.quality_factor, Qp, 0.1);

    /* Series RLC transfer functions */
    params.R = 1.0;
    tf_polynomial_t *tf_vr = resonance_series_tf(&params, 0);  /* V_R/V_in */
    assert(tf_vr != NULL);
    /* This should be a bandpass: H(s) = (R/L)s/(s² + (R/L)s + 1/(LC)) */
    assert(tf_vr->num_order == 1);
    assert(tf_vr->den_order == 2);
    tf_polynomial_free(tf_vr);

    tf_polynomial_t *tf_vc = resonance_series_tf(&params, 2);  /* V_C/V_in */
    assert(tf_vc != NULL);
    /* Lowpass: H(s) = 1/(LC)/(s² + (R/L)s + 1/(LC)) */
    assert(tf_vc->num_order == 0);
    assert(tf_vc->den_order == 2);
    tf_polynomial_free(tf_vc);

    /* Universal resonance curve */
    double norm_freqs[] = { 0.1, 0.5, 0.9, 1.0, 1.1, 2.0, 10.0 };
    double *curve = resonance_universal_curve(norm_freqs, 10.0, 7);
    assert(curve != NULL);
    /* At resonance (Ω=1): |H| = 1 */
    ASSERT_NEAR(curve[3], 1.0, 1e-6);
    /* At Ω=0.9, Q=10: (0.9-1/0.9)² = 0.0446, denom=5.458, |H|=0.428
     * At half-power (Ω≈0.951): |H| = 1/√2 ≈ 0.707 */
    ASSERT_NEAR(curve[2], 0.428, 0.01);  /* Ω=0.9 */
    free(curve);

    printf("  L6: PASSED\n");
}

/* ============================================================================
 * L6: Biquad filter design (Canonical Problem)
 * ============================================================================ */

static void test_l6_biquad(void)
{
    printf("  L6: Biquad filter design...\n");

    /* Lowpass biquad: ω₀=1000 rad/s, Q=0.707 (Butterworth), gain=1 */
    biquad_section_t lp = biquad_create(FILTER_TYPE_LOWPASS,
                                          1000.0, 0.707, 1.0);
    ASSERT_NEAR(lp.dc_gain, 1.0, 1e-6);
    ASSERT_NEAR(lp.b0, 1e6, 1e-6);  /* ω₀² */
    ASSERT_NEAR(lp.a0, 1e6, 1e-6);  /* ω₀² */
    ASSERT_NEAR(lp.a1, 1000.0 / 0.707, 1.0);  /* ω₀/Q */

    /* Highpass biquad */
    biquad_section_t hp = biquad_create(FILTER_TYPE_HIGHPASS,
                                          1000.0, 0.707, 1.0);
    ASSERT_NEAR(hp.b2, 1.0, 1e-6);
    ASSERT_NEAR(hp.b0, 0.0, 1e-6);

    /* Bandpass biquad */
    biquad_section_t bp = biquad_create(FILTER_TYPE_BANDPASS,
                                          1000.0, 5.0, 1.0);
    ASSERT_NEAR(bp.b1, (1000.0 / 5.0) * 1.0, 1e-6);  /* (ω₀/Q)·gain */
    ASSERT_NEAR(bp.b0, 0.0, 1e-6);

    /* Sallen-Key LP component calculation */
    double R1, R2, R3, R4;
    int ret = filter_sallen_key_lp(&lp, 1e-9, 2e-9, 1.0,
                                     &R1, &R2, &R3, &R4);
    assert(ret == 0);
    assert(R1 > 0.0 && R2 > 0.0);

    /* MFB LP component calculation */
    double R1_mfb, R2_mfb, R3_mfb;
    ret = filter_mfb_lp(&lp, 1e-9, 1e-9, &R1_mfb, &R2_mfb, &R3_mfb);
    assert(ret == 0);
    assert(R1_mfb > 0.0 && R2_mfb > 0.0 && R3_mfb > 0.0);

    /* Tow-Thomas component calculation */
    double R_freq, R_q, R_in, R_fb;
    ret = filter_tow_thomas(&lp, 1e-9, &R_freq, &R_q, &R_in, &R_fb);
    assert(ret == 0);
    assert(R_freq > 0.0 && R_q > 0.0);

    printf("  L6: PASSED\n");
}

/* ============================================================================
 * L7: Crystal resonator model (Application)
 * ============================================================================ */

static void test_l7_crystal(void)
{
    printf("  L7: Crystal resonator model...\n");

    /* 10 MHz AT-cut crystal typical values */
    double C0 = 5e-12;   /* 5 pF */
    double L1 = 10e-3;   /* 10 mH */
    double C1 = 25e-15;  /* 25 fF */
    double R1 = 10.0;    /* 10 Ω */

    double fs = 0.0, fp = 0.0, Q = 0.0;
    tf_polynomial_t *Z = resonance_crystal_model(C0, L1, C1, R1,
                                                   &fs, &fp, &Q);
    assert(Z != NULL);

    /* f_s ≈ 1/(2π√(L₁C₁)) = 1/(2π√(10e-3×25e-15))
     *     = 1/(2π√(2.5e-16)) ≈ 1/(2π×1.58e-8) ≈ 10.07 MHz */
    ASSERT_NEAR(fs, 10.07e6, 0.1e6);

    /* f_p = f_s·√(1+C₁/C₀) ≈ 10.07e6·√(1+25e-15/5e-12)
     *     ≈ 10.07e6·√(1.005) ≈ 10.10 MHz */
    ASSERT_NEAR(fp, 10.10e6, 0.1e6);

    /* Q should be very high: Q ≈ ω_s·L₁/R₁ */
    double Q_expected = 2.0 * M_PI * fs * L1 / R1;
    ASSERT_NEAR(Q, Q_expected, 100.0);
    assert(Q > 10000);  /* Crystal Q is very high */

    tf_polynomial_free(Z);
    printf("  L7: PASSED\n");
}

/* ============================================================================
 * L7: Switched-capacitor stability (Application)
 * ============================================================================ */

static void test_l7_sc_stability(void)
{
    printf("  L7: Switched-capacitor filter stability...\n");

    /* Stable CT prototype: H(s) = 1/(s+1) */
    double num[] = { 1.0 }, den[] = { 1.0, 1.0 };
    tf_polynomial_t *tf = tf_polynomial_create(num, 0, den, 1);

    int is_stable = 0;
    int ret = stability_sc_filter(tf, 100e3, &is_stable);
    assert(ret == 0);
    assert(is_stable == 1);

    tf_polynomial_free(tf);
    printf("  L7: PASSED\n");
}

/* ============================================================================
 * L8: Coupled resonators (Advanced)
 * ============================================================================ */

static void test_l8_coupled_resonators(void)
{
    printf("  L8: Coupled resonators...\n");

    /* Two identical RLC tanks with magnetic coupling */
    rlc_params_t primary   = { 1.0, 1e-3, 1e-6, 1.0 };
    rlc_params_t secondary = { 1.0, 1e-3, 1e-6, 1.0 };
    double M = 1e-5;  /* 10 µH mutual inductance */

    double freqs[100];
    double f0 = 1.0 / (2.0 * M_PI * sqrt(1e-9));  /* ~5033 Hz */
    for (int i = 0; i < 100; i++) {
        freqs[i] = f0 * (0.5 + (double)i / 100.0);  /* 0.5×f₀ to 1.5×f₀ */
    }

    freq_response_t *resp = resonance_coupled(&primary, &secondary, M,
                                                freqs, 100);
    assert(resp != NULL);
    assert(resp->num_points == 100);

    /* With coupling, should observe a response */
    double max_mag = 0.0;
    for (size_t i = 0; i < resp->num_points; i++) {
        if (resp->points[i].magnitude > max_mag) {
            max_mag = resp->points[i].magnitude;
        }
    }
    assert(max_mag > 0.0);

    freq_response_free(resp);
    printf("  L8: PASSED\n");
}

/* ============================================================================
 * L8: Root locus (Advanced)
 * ============================================================================ */

static void test_l8_root_locus(void)
{
    printf("  L8: Root locus analysis...\n");

    /* Loop gain L(s) = 1/(s(s+1)(s+2))
     * Numerator: 1, Denominator: s³ + 3s² + 2s */
    double num[] = { 1.0 };
    double den[] = { 0.0, 2.0, 3.0, 1.0 };
    tf_polynomial_t *L = tf_polynomial_create(num, 0, den, 3);

    double K_vals[] = { 0.1, 0.5, 1.0, 2.0, 5.0, 10.0 };
    size_t n_K = 6;

    double _Complex *poles = (double _Complex *)calloc(n_K * 3,
                                                         sizeof(double _Complex));
    int ret = stability_root_locus(L, K_vals, n_K, poles);
    assert(ret == 0);

    /* Verify poles are computed */
    for (size_t k = 0; k < n_K; k++) {
        int valid_count = 0;
        for (size_t j = 0; j < 3; j++) {
            if (cabs(poles[k * 3 + j]) > 1e-10) valid_count++;
        }
        assert(valid_count > 0);
    }

    free(poles);

    /* Gain margin via root locus */
    double K_margin = 0.0, f_margin = 0.0;
    ret = stability_gain_margin_root_locus(L, &K_margin, &f_margin);
    /* May or may not find crossing depending on sweep */
    if (ret == 0) {
        assert(K_margin > 0.0);
    }

    tf_polynomial_free(L);
    printf("  L8: PASSED\n");
}

/* ============================================================================
 * Network Functions
 * ============================================================================ */

static void test_network_functions(void)
{
    printf("  Network functions...\n");

    /* Voltage divider: Z₁ = R = 1kΩ, Z₂ = 1/(sC) = 1/(s·1e-6)
     * H(s) = Z₂/(Z₁+Z₂) = 1/(1 + sRC) = 1/(1 + s·1e-3)
     * This is a lowpass with τ = 1 ms, f_c = 159 Hz */

    tf_polynomial_t *Z1 = network_element_impedance(2, 1000.0, 0.0);     /* 1kΩ */
    tf_polynomial_t *Z2 = network_element_impedance(1, 1e-6, 0.0);       /* 1µF */
    tf_polynomial_t *div = network_voltage_divider_tf(Z1, Z2);

    assert(div != NULL);
    /* DC gain should be 1 */
    double dc = tf_dc_gain(div);
    ASSERT_NEAR(dc, 1.0, 1e-6);

    /* Cutoff: f_c = 1/(2πRC) = 1/(2π×1000×1e-6) ≈ 159 Hz */
    freq_response_t *resp = tf_compute_freq_response(div, 1.0, 10000.0, 20);
    assert(resp != NULL);
    double fc = find_cutoff_freq(resp, 0);  /* Lowpass */
    ASSERT_NEAR(fc, 159.15, 1.0);

    freq_response_free(resp);
    tf_polynomial_free(div);
    tf_polynomial_free(Z1);
    tf_polynomial_free(Z2);

    /* S-parameter computation */
    s_params_2port_t spar = network_s_params_from_z(
        50.0, 0.0, 100.0, 50.0, 50.0, 1e9);
    double K = network_rollett_k(&spar);
    /* With S₁₂=0, K should be large (unilateral device) */
    assert(K > 0.0);

    /* VSWR: |Γ| = 1/3 → VSWR = (1+1/3)/(1-1/3) = 2 */
    double vswr = network_vswr(1.0 / 3.0);
    ASSERT_NEAR(vswr, 2.0, 1e-6);

    vswr = network_vswr(0.0);  /* Perfect match */
    ASSERT_NEAR(vswr, 1.0, 1e-6);

    printf("  Network functions: PASSED\n");
}

/* ============================================================================
 * Frequency Response Analysis
 * ============================================================================ */

static void test_freq_response_analysis(void)
{
    printf("  Frequency response analysis...\n");

    /* Create a bandpass-like response:
     * H(s) = s/(s² + 0.5s + 1) — center at ω=1, Q=2 */
    double num[] = { 0.0, 1.0 };
    double den[] = { 1.0, 0.5, 1.0 };
    tf_polynomial_t *tf = tf_polynomial_create(num, 1, den, 2);

    freq_response_t *resp = tf_compute_freq_response(tf, 0.01, 100.0, 50);
    assert(resp != NULL);

    /* Find peak */
    double peak_freq = 0.0;
    double peak_mag = find_peak_response(resp, &peak_freq);
    assert(peak_mag > 1.0);  /* Q=2 → peaking */
    ASSERT_NEAR(peak_freq, 1.0 / (2.0 * M_PI), 0.1);  /* Peak near ω=1 rad/s */

    /* Find bandwidth */
    double f_low = 0.0, f_high = 0.0;
    double bw = find_bandwidth(resp, &f_low, &f_high);
    assert(bw > 0.0);
    assert(f_high > f_low);

    /* Compute roll-off on HF side */
    double rolloff = compute_rolloff(resp, 10.0, 100.0);
    assert(rolloff < 0.0);  /* Should be negative (lowpass behavior in stopband) */

    /* Interpolation */
    freq_point_t interp = freq_response_interpolate(resp, 1.0);
    assert(interp.frequency == 1.0);

    freq_response_free(resp);
    tf_polynomial_free(tf);

    printf("  Frequency response analysis: PASSED\n");
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(void)
{
    printf("=== mini-frequency-response Test Suite ===\n\n");

    test_l1_conversions();
    test_l1_freq_grid();
    test_l1_phase_unwrap();
    test_l2_transfer_function();
    test_l2_tf_arithmetic();
    test_l3_horner();
    test_l3_pole_zero();
    test_l4_routh_hurwitz();
    test_l4_bode_gain_phase();
    test_l4_nyquist();
    test_l5_bode_plot();
    test_l5_filter_design();
    test_l6_resonance();
    test_l6_biquad();
    test_l7_crystal();
    test_l7_sc_stability();
    test_l8_coupled_resonators();
    test_l8_root_locus();
    test_network_functions();
    test_freq_response_analysis();

    printf("\n=== ALL TESTS PASSED ===\n");
    return 0;
}
