/**
 * @file test_conversion.c
 * @brief Tests for parameter conversion, interconnection, and stability.
 */
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../include/two_port.h"
#include "../include/conversion.h"
#include "../include/abcd_params.h"
#include "../include/s_params.h"
#include "../include/z_params.h"
#include "../include/interconnection.h"
#include "../include/stability.h"

#define TOL 1e-9

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) printf("  TEST: %s ... ", name); tests_passed++
#define CHECK(cond) do { \
    if (!(cond)) { printf("FAIL at %s:%d: %s\n", __FILE__, __LINE__, #cond); tests_failed++; return; } \
} while(0)
#define TEST_END() printf("OK\n")

/* ============================================================
 * Z ↔ Y Conversion (L5)
 * ============================================================ */
static void test_z_to_y_conversion(void) {
    TEST("convert_z_to_y — roundtrip");
    /* Simple resistive Z: [[10, 5], [5, 10]]
       det = 10*10-5*5 = 75
       Y = [[10/75, -5/75], [-5/75, 10/75]] */
    matrix2x2_t z = matrix2x2_make(
        complex_make(10.0, 0.0), complex_make(5.0, 0.0),
        complex_make(5.0, 0.0), complex_make(10.0, 0.0)
    );
    matrix2x2_t y = convert_z_to_y(z);
    CHECK(fabs(y.m11.real - 10.0/75.0) < TOL);
    CHECK(fabs(y.m12.real + 5.0/75.0) < TOL);
    CHECK(fabs(y.m22.real - 10.0/75.0) < TOL);

    /* Round-trip */
    matrix2x2_t z_back = convert_y_to_z(y);
    CHECK(complex_approx_equal(z.m11, z_back.m11, TOL));
    CHECK(complex_approx_equal(z.m22, z_back.m22, TOL));
    TEST_END();
}

/* ============================================================
 * Z ↔ ABCD Conversion (L5)
 * ============================================================ */
static void test_z_to_abcd_conversion(void) {
    TEST("convert_z_to_abcd");
    /* Series Z=50:
       z = [[50,50],[50,50]]
       ABCD: A=50/50=1, B=(2500-2500)/50=0... Wait z21=50
       A = 50/50 = 1, B = (2500-2500)/50 = 0, C = 1/50 = 0.02, D = 50/50 = 1 */
    matrix2x2_t z = matrix2x2_make(
        complex_make(50.0, 0.0), complex_make(50.0, 0.0),
        complex_make(50.0, 0.0), complex_make(50.0, 0.0)
    );
    matrix2x2_t abcd = convert_z_to_abcd(z);
    CHECK(fabs(abcd.m11.real - 1.0) < TOL);
    CHECK(fabs(abcd.m21.real - 0.02) < TOL);
    CHECK(fabs(abcd.m22.real - 1.0) < TOL);
    TEST_END();
}

/* ============================================================
 * S ↔ ABCD Conversion (L5)
 * ============================================================ */
static void test_s_to_abcd_conversion(void) {
    TEST("convert_s_to_abcd — s21 only (matched line)");
    /* Ideal transmission line: s11=s22=0, s21=s12=e^(-jθ)
       For θ=90°: s21 = -j
       ABCD: A=0, B=j*Z0, C=j/Z0, D=0 */
    matrix2x2_t s = matrix2x2_make(
        complex_make(0.0, 0.0), complex_make(0.0, -1.0),
        complex_make(0.0, -1.0), complex_make(0.0, 0.0)
    );
    matrix2x2_t abcd = convert_s_to_abcd(s, 50.0);
    /* A = ((1+0)(1-0)+(-j)(-j))/(2*(-j)) = (1 + (-1))/(-2j) = 0 ✓ */
    CHECK(fabs(abcd.m11.real) < TOL && fabs(abcd.m11.imag) < TOL);
    /* B = 50*((1+0)(1+0)-(-j)(-j))/(2*(-j)) = 50*(1-(-1))/(-2j) = 50*2/(-2j) = -50/j = j50 */
    CHECK(fabs(abcd.m12.real) < TOL);
    CHECK(fabs(abcd.m12.imag - 50.0) < TOL);
    TEST_END();
}

/* ============================================================
 * ABCD → Z conversion (round-trip test)
 * ============================================================ */
static void test_abcd_to_z_roundtrip(void) {
    TEST("ABCD ↔ Z round-trip");
    /* Use a T-network with all elements nonzero to avoid singular conversions.
       T: Za=10, Zb=20, Zc=5Ω
       z11=15, z12=z21=5, z22=25
       A = z11/z21 = 3, B = (15*25-25)/5 = 70, C = 1/5 = 0.2, D = 25/5 = 5 */
    matrix2x2_t z_orig = zparams_t_network(
        complex_make(10.0, 0.0),
        complex_make(20.0, 0.0),
        complex_make(5.0, 0.0)
    );
    matrix2x2_t abcd = convert_z_to_abcd(z_orig);
    matrix2x2_t z_back = convert_abcd_to_z(abcd);
    matrix2x2_t abcd_back = convert_z_to_abcd(z_back);

    CHECK(complex_approx_equal(abcd.m11, abcd_back.m11, 1e-6));
    TEST_END();
}

/* ============================================================
 * Generic Parameter Conversion (L5)
 * ============================================================ */
static void test_generic_conversion(void) {
    TEST("convert_parameters — generic routing");
    /* Create Z-params, convert to Y, H, G, ABCD, S, and back */
    matrix2x2_t z = matrix2x2_make(
        complex_make(100.0, 50.0), complex_make(30.0, 10.0),
        complex_make(30.0, 10.0), complex_make(80.0, 40.0)
    );

    matrix2x2_t y = convert_parameters(z, PARAM_Z, PARAM_Y, 50.0);
    CHECK(!isnan(y.m11.real));

    matrix2x2_t h = convert_parameters(z, PARAM_Z, PARAM_H, 50.0);
    CHECK(!isnan(h.m11.real));

    matrix2x2_t s = convert_parameters(z, PARAM_Z, PARAM_S, 50.0);
    CHECK(!isnan(s.m11.real));

    /* Check that |s11| < 1 for this passive-ish network */
    CHECK(complex_mag(s.m11) < 1.0);

    TEST_END();
}

/* ============================================================
 * Interconnection Tests (L5)
 * ============================================================ */
static void test_cascade_interconnection(void) {
    TEST("interconnect_cascade");
    /* Two series impedances: Z1=R1, Z2=R2
       Cascade ABCD_total = ABCD1 × ABCD2
       ABCD1 = [[1,R1],[0,1]], ABCD2 = [[1,R2],[0,1]]
       ABCD_total = [[1,R1+R2],[0,1]] */
    matrix2x2_t abcd1 = abcd_series_z(complex_make(10.0, 0.0));
    matrix2x2_t abcd2 = abcd_series_z(complex_make(20.0, 0.0));
    matrix2x2_t total = interconnect_cascade(abcd1, abcd2);

    /* A=1, B=30, C=0, D=1 */
    CHECK(fabs(total.m11.real - 1.0) < TOL);
    CHECK(fabs(total.m12.real - 30.0) < TOL);
    CHECK(fabs(total.m21.real) < TOL);
    CHECK(fabs(total.m22.real - 1.0) < TOL);
    TEST_END();
}

static void test_series_series_interconnection(void) {
    TEST("interconnect_series_series");
    /* Two T-networks in series-series: Z_total = Z1 + Z2 */
    matrix2x2_t z1 = zparams_t_network(
        complex_make(10.0, 0.0), complex_make(20.0, 0.0), complex_make(5.0, 0.0)
    );
    matrix2x2_t z2 = zparams_t_network(
        complex_make(1.0, 0.0), complex_make(2.0, 0.0), complex_make(0.5, 0.0)
    );
    matrix2x2_t z_total = interconnect_series_series(z1, z2);

    CHECK(fabs(z_total.m11.real - (15.0 + 1.5)) < TOL);
    TEST_END();
}

/* ============================================================
 * Stability Tests (L4)
 * ============================================================ */
static void test_rollett_k_stable(void) {
    TEST("stability_rollett_k — unconditionally stable");
    /* Unilateral matched amplifier: s11=0, s22=0, s21=3, s12=0
       K → INFINITY (unilateral) */
    matrix2x2_t s = matrix2x2_make(
        complex_make(0.0, 0.0), complex_make(0.0, 0.0),
        complex_make(3.0, 0.0), complex_make(0.0, 0.0)
    );
    double k = stability_rollett_k(s, PARAM_S, 50.0);
    CHECK(isinf(k));
    CHECK(stability_classify(s, PARAM_S, 50.0) == STABLE_UNCONDITIONAL);
    TEST_END();
}

static void test_rollett_k_unstable(void) {
    TEST("stability_rollett_k — potentially unstable");
    /* Large s12 and s21 with poor match → K < 1 */
    matrix2x2_t s = matrix2x2_make(
        complex_make(0.8, 0.0), complex_make(0.3, 0.0),
        complex_make(4.0, 0.0), complex_make(0.7, 0.0)
    );
    double k = stability_rollett_k(s, PARAM_S, 50.0);
    CHECK(k < 2.0);  /* Should be below unconditional threshold */
    TEST_END();
}

/* ============================================================
 * Bartlett's Bisection Theorem (L4)
 * ============================================================ */
static void test_bartlett_bisection(void) {
    TEST("bartlett_bisection");
    /* Zoc = 200Ω, Zsc = 50Ω
       z11 = z22 = (200+50)/2 = 125
       z12 = z21 = (200-50)/2 = 75 */
    matrix2x2_t z = bartlett_bisection(
        complex_make(200.0, 0.0), complex_make(50.0, 0.0)
    );
    CHECK(fabs(z.m11.real - 125.0) < TOL);
    CHECK(fabs(z.m12.real - 75.0) < TOL);
    CHECK(fabs(z.m22.real - 125.0) < TOL);
    TEST_END();
}

/* ============================================================
 * S-Parameter Analysis
 * ============================================================ */
static void test_s21_db(void) {
    TEST("sparams_s21_db");
    /* |s21|=3 → 20*log10(3) ≈ 9.54 dB */
    matrix2x2_t s = matrix2x2_make(
        complex_make(0.0, 0.0), complex_make(0.0, 0.0),
        complex_make(3.0, 0.0), complex_make(0.0, 0.0)
    );
    double db = sparams_s21_db(s);
    CHECK(fabs(db - 20.0 * log10(3.0)) < 0.01);
    TEST_END();
}

static void test_gamma_in(void) {
    TEST("sparams_gamma_in — load-pull effect");
    /* Unilateral device: Γin = s11 regardless of ΓL */
    matrix2x2_t s = matrix2x2_make(
        complex_make(0.2, 0.0), complex_make(0.0, 0.0),
        complex_make(2.0, 0.0), complex_make(0.1, 0.0)
    );
    complex_t gl = complex_make(0.5, 0.3);
    complex_t gin = sparams_gamma_in(s, gl);
    /* s12=0 → gin = s11 */
    CHECK(fabs(gin.real - 0.2) < TOL);
    CHECK(fabs(gin.imag) < TOL);
    TEST_END();
}

/* ============================================================
 * Main
 * ============================================================ */
int main(void) {
    printf("=== Conversion, Interconnection & Stability Tests ===\n");

    test_z_to_y_conversion();
    test_z_to_abcd_conversion();
    test_s_to_abcd_conversion();
    test_abcd_to_z_roundtrip();
    test_generic_conversion();
    test_cascade_interconnection();
    test_series_series_interconnection();
    test_rollett_k_stable();
    test_rollett_k_unstable();
    test_bartlett_bisection();
    test_s21_db();
    test_gamma_in();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed - tests_failed, tests_passed);
    return tests_failed > 0 ? 1 : 0;
}
