/**
 * @file test_core.c
 * @brief Tests for core two-port operations: complex math, matrices, Z/Y/H/G params.
 *
 * All tests use standard assert(). Each test validates a specific
 * mathematical fact or circuit theory identity.
 */
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../include/two_port.h"
#include "../include/z_params.h"
#include "../include/y_params.h"
#include "../include/h_params.h"
#include "../include/g_params.h"

#define TOL 1e-9

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    printf("  TEST: %s ... ", name); \
    tests_passed++; \
} while(0)

#define CHECK(cond) do { \
    if (!(cond)) { \
        printf("FAIL at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        tests_failed++; \
        return; \
    } \
} while(0)

#define TEST_END() printf("OK\n")

/* ================================================================
 * Complex Arithmetic (L3)
 * ================================================================ */
static void test_complex_add(void) {
    TEST("complex_add");
    complex_t a = complex_make(3.0, 4.0);
    complex_t b = complex_make(1.0, -2.0);
    complex_t c = complex_add(a, b);
    CHECK(fabs(c.real - 4.0) < TOL);
    CHECK(fabs(c.imag - 2.0) < TOL);
    TEST_END();
}

static void test_complex_mul(void) {
    TEST("complex_mul");
    /* (3+j4)*(1+j2) = 3+j6+j4-8 = -5+j10 */
    complex_t a = complex_make(3.0, 4.0);
    complex_t b = complex_make(1.0, 2.0);
    complex_t c = complex_mul(a, b);
    CHECK(fabs(c.real + 5.0) < TOL);
    CHECK(fabs(c.imag - 10.0) < TOL);
    TEST_END();
}

static void test_complex_div(void) {
    TEST("complex_div");
    /* (4+j6)/(2+j2) = (4+j6)(2-j2)/(4+4) = (8+12+j12-j8)/8 = (20+j4)/8 = 2.5+j0.5 */
    complex_t a = complex_make(4.0, 6.0);
    complex_t b = complex_make(2.0, 2.0);
    complex_t c = complex_div(a, b);
    CHECK(fabs(c.real - 2.5) < TOL);
    CHECK(fabs(c.imag - 0.5) < TOL);

    /* Division by zero */
    complex_t zero = complex_make(0.0, 0.0);
    complex_t d = complex_div(a, zero);
    CHECK(isnan(d.real));
    TEST_END();
}

static void test_complex_mag(void) {
    TEST("complex_mag");
    /* |3+j4| = 5 */
    complex_t z = complex_make(3.0, 4.0);
    double mag = complex_mag(z);
    CHECK(fabs(mag - 5.0) < TOL);
    TEST_END();
}

static void test_complex_arg(void) {
    TEST("complex_arg");
    /* arg(1+j0) = 0, arg(j) = π/2, arg(-1) = π */
    CHECK(fabs(complex_arg(complex_make(1.0, 0.0))) < TOL);
    CHECK(fabs(complex_arg(complex_make(0.0, 1.0)) - M_PI/2.0) < TOL);
    CHECK(fabs(complex_arg(complex_make(-1.0, 0.0)) - M_PI) < TOL);
    TEST_END();
}

static void test_complex_conj(void) {
    TEST("complex_conj");
    complex_t z = complex_make(5.0, -3.0);
    complex_t zc = complex_conj(z);
    CHECK(fabs(zc.real - 5.0) < TOL);
    CHECK(fabs(zc.imag - 3.0) < TOL);
    TEST_END();
}

static void test_complex_exp(void) {
    TEST("complex_exp (Euler's formula)");
    /* e^(jπ) = cos(π)+j*sin(π) = -1+j0 */
    complex_t z = complex_make(0.0, M_PI);
    complex_t ez = complex_exp(z);
    CHECK(fabs(ez.real + 1.0) < 1e-9);
    CHECK(fabs(ez.imag) < 1e-9);

    /* e^(1+j0) = e */
    complex_t z2 = complex_make(1.0, 0.0);
    complex_t ez2 = complex_exp(z2);
    CHECK(fabs(ez2.real - 2.718281828) < 1e-9);
    CHECK(fabs(ez2.imag) < 1e-9);
    TEST_END();
}

static void test_complex_sqrt(void) {
    TEST("complex_sqrt");
    /* sqrt(-1) = j */
    complex_t z = complex_make(-1.0, 0.0);
    complex_t sz = complex_sqrt(z);
    CHECK(fabs(sz.real) < 1e-9);
    CHECK(fabs(sz.imag - 1.0) < 1e-9);

    /* sqrt(4) = 2 */
    complex_t z2 = complex_make(4.0, 0.0);
    complex_t sz2 = complex_sqrt(z2);
    CHECK(fabs(sz2.real - 2.0) < 1e-9);
    CHECK(fabs(sz2.imag) < 1e-9);
    TEST_END();
}

/* ================================================================
 * Matrix Operations (L3)
 * ================================================================ */
static void test_matrix_det(void) {
    TEST("matrix2x2_det");
    complex_t one = complex_make(1.0, 0.0);
    complex_t two = complex_make(2.0, 0.0);
    complex_t three = complex_make(3.0, 0.0);
    complex_t four = complex_make(4.0, 0.0);

    /* det([[1,2],[3,4]]) = 1*4 - 2*3 = -2 */
    matrix2x2_t m = matrix2x2_make(one, two, three, four);
    complex_t d = matrix2x2_det(m);
    CHECK(fabs(d.real + 2.0) < TOL);
    CHECK(fabs(d.imag) < TOL);
    TEST_END();
}

static void test_matrix_inv(void) {
    TEST("matrix2x2_inv");
    complex_t one = complex_make(1.0, 0.0);
    complex_t two = complex_make(2.0, 0.0);
    complex_t three = complex_make(3.0, 0.0);
    complex_t four = complex_make(4.0, 0.0);

    /* Inverse: [[4,-2],[-3,1]] / (1*4-2*3) = [[-2,1],[1.5,-0.5]] */
    matrix2x2_t m = matrix2x2_make(one, two, three, four);
    matrix2x2_t inv = matrix2x2_inv(m);
    CHECK(fabs(inv.m11.real + 2.0) < TOL);
    CHECK(fabs(inv.m12.real - 1.0) < TOL);
    CHECK(fabs(inv.m21.real - 1.5) < TOL);
    CHECK(fabs(inv.m22.real + 0.5) < TOL);

    /* Verify: m * inv = I */
    matrix2x2_t prod = matrix2x2_mul(m, inv);
    CHECK(fabs(prod.m11.real - 1.0) < TOL);
    CHECK(fabs(prod.m22.real - 1.0) < TOL);
    TEST_END();
}

static void test_matrix_mul(void) {
    TEST("matrix2x2_mul — cascade property");
    /* ABCD for series Z=2: [[1,2],[0,1]]
       ABCD for shunt Y=0.5: [[1,0],[0.5,1]]
       Cascade: [[1,2],[0,1]] × [[1,0],[0.5,1]] = [[2,2],[0.5,1]] */
    complex_t z_el = complex_make(2.0, 0.0);
    complex_t y_el = complex_make(0.5, 0.0);
    complex_t zero = complex_make(0.0, 0.0);
    complex_t one = complex_make(1.0, 0.0);

    matrix2x2_t series_z = matrix2x2_make(one, z_el, zero, one);
    matrix2x2_t shunt_y = matrix2x2_make(one, zero, y_el, one);
    matrix2x2_t cascade = matrix2x2_mul(series_z, shunt_y);

    CHECK(fabs(cascade.m11.real - 2.0) < TOL);
    CHECK(fabs(cascade.m12.real - 2.0) < TOL);
    CHECK(fabs(cascade.m21.real - 0.5) < TOL);
    CHECK(fabs(cascade.m22.real - 1.0) < TOL);
    TEST_END();
}

/* ================================================================
 * Z-Parameter Tests (L1, L2)
 * ================================================================ */
static void test_zparams_t_network(void) {
    TEST("zparams_t_network");
    /* T-network: Za=10, Zb=20, Zc=5
       z11 = 10+5=15, z12=z21=5, z22=20+5=25 */
    complex_t za = complex_make(10.0, 0.0);
    complex_t zb = complex_make(20.0, 0.0);
    complex_t zc = complex_make(5.0, 0.0);
    matrix2x2_t z = zparams_t_network(za, zb, zc);

    CHECK(fabs(z.m11.real - 15.0) < TOL);
    CHECK(fabs(z.m12.real - 5.0) < TOL);
    CHECK(fabs(z.m21.real - 5.0) < TOL);
    CHECK(fabs(z.m22.real - 25.0) < TOL);
    TEST_END();
}

static void test_zparams_rlc(void) {
    TEST("zparams_series_rlc — resonance");
    /* Series RLC: R=10, L=1e-3, C=1e-6
       ω0 = 1/sqrt(LC) = 1/sqrt(1e-9) = 31623 rad/s
       At ω0: Z = R = 10Ω */
    double omega0 = 1.0 / sqrt(1e-3 * 1e-6);
    matrix2x2_t z = zparams_series_rlc(10.0, 1e-3, 1e-6, omega0);
    CHECK(fabs(z.m11.real - 10.0) < 1e-6);
    CHECK(fabs(z.m11.imag) < 1e-6);
    TEST_END();
}

static void test_zparams_input_impedance(void) {
    TEST("zparams_input_impedance");
    /* T-network Za=Zb=50, Zc=100 → z11=150, z12=z21=100, z22=150
       Terminate with ZL=50:
       Zin = 150 - 10000/(150+50) = 150 - 50 = 100 Ω */
    matrix2x2_t z = zparams_t_network(
        complex_make(50.0, 0.0),
        complex_make(50.0, 0.0),
        complex_make(100.0, 0.0)
    );
    complex_t zl = complex_make(50.0, 0.0);
    complex_t zin = zparams_input_impedance(z, zl);
    CHECK(fabs(zin.real - 100.0) < 1e-6);
    TEST_END();
}

static void test_zparams_to_t_network(void) {
    TEST("zparams_to_t_network — round-trip");
    /* Create T-network, extract, verify */
    complex_t za = complex_make(10.0, 5.0);
    complex_t zb = complex_make(20.0, -3.0);
    complex_t zc = complex_make(7.0, 2.0);
    matrix2x2_t z = zparams_t_network(za, zb, zc);

    complex_t ea, eb, ec;
    int ret = zparams_to_t_network(z, &ea, &eb, &ec);
    CHECK(ret == 0);
    CHECK(complex_approx_equal(ea, za, TOL));
    CHECK(complex_approx_equal(eb, zb, TOL));
    CHECK(complex_approx_equal(ec, zc, TOL));
    TEST_END();
}

/* ================================================================
 * Y-Parameter Tests
 * ================================================================ */
static void test_yparams_pi_network(void) {
    TEST("yparams_pi_network");
    /* π-network: Ya=0.01, Yb=0.02, Yc=0.03 S
       y11 = 0.01+0.02=0.03, y12=-0.02, y21=-0.02, y22=0.02+0.03=0.05 */
    complex_t ya = complex_make(0.01, 0.0);
    complex_t yb = complex_make(0.02, 0.0);
    complex_t yc = complex_make(0.03, 0.0);
    matrix2x2_t y = yparams_pi_network(ya, yb, yc);

    CHECK(fabs(y.m11.real - 0.03) < TOL);
    CHECK(fabs(y.m12.real + 0.02) < TOL);
    CHECK(fabs(y.m21.real + 0.02) < TOL);
    CHECK(fabs(y.m22.real - 0.05) < TOL);
    TEST_END();
}

static void test_yparams_mosfet_hf(void) {
    TEST("yparams_mosfet_hf — fT calculation");
    double gm = 0.1;    /* 100 mS */
    double cgs = 1e-12; /* 1 pF */
    double cgd = 0.1e-12; /* 0.1 pF */
    double ft = yparams_f_t_mosfet(gm, cgs, cgd);
    /* fT ≈ 0.1 / (2π * 1.1e-12) ≈ 14.5 GHz */
    double ft_expected = gm / (2.0 * M_PI * (cgs + cgd));
    CHECK(fabs(ft - ft_expected) < 1e3);
    CHECK(ft > 1e9);  /* Should be in GHz range */
    TEST_END();
}

/* ================================================================
 * H-Parameter Tests
 * ================================================================ */
static void test_hparams_bjt_ce(void) {
    TEST("hparams_bjt_ce");
    /* Typical BJT: rbb=50, rpi=2.5k, beta=100, ro=50k */
    matrix2x2_t h = hparams_bjt_ce(50.0, 2500.0, 100.0, 50000.0);

    /* h_ie ≈ rbb + rpi = 2550 Ω */
    CHECK(fabs(h.m11.real - 2550.0) < 1.0);
    /* h_fe = 100 */
    CHECK(fabs(h.m21.real - 100.0) < TOL);
    /* h_oe = 1/50k = 20 µS */
    CHECK(fabs(h.m22.real - 2e-5) < 1e-9);
    TEST_END();
}

static void test_hparams_ce_to_cb(void) {
    TEST("hparams_ce_to_cb — alpha ≈ -0.99");
    matrix2x2_t h_ce = hparams_bjt_ce(50.0, 2500.0, 100.0, 50000.0);
    matrix2x2_t h_cb = hparams_ce_to_cb(h_ce);

    /* h_fb ≈ -α = -β/(1+β) = -100/101 ≈ -0.9901 */
    double alpha = 100.0 / 101.0;
    CHECK(fabs(h_cb.m21.real + alpha) < 0.001);
    TEST_END();
}

static void test_hparams_beta_cutoff(void) {
    TEST("hparams_beta_cutoff_freq");
    /* rπ=2.5k, Cπ=10pF → fβ = 1/(2π*2500*1e-11) ≈ 6.37 MHz */
    double fbeta = hparams_beta_cutoff_freq(2500.0, 10e-12);
    double expected = 1.0 / (2.0 * M_PI * 2500.0 * 10e-12);
    CHECK(fabs(fbeta - expected) < 1000.0);
    TEST_END();
}

/* ================================================================
 * G-Parameter Tests
 * ================================================================ */
static void test_gparams_fet_cs_gain(void) {
    TEST("gparams_fet_cs intrinsic gain");
    /* gm=1mS, ro=50k → |gain| = 1e-3*50000 = 50 */
    matrix2x2_t g = gparams_fet_cs(1e-3, 50000.0, 1e6);
    double a_int = gparams_intrinsic_gain(g);
    CHECK(fabs(a_int - 50.0) < 0.1);
    TEST_END();
}

/* ================================================================
 * Reflection Coefficient & VSWR (L2)
 * ================================================================ */
static void test_reflection_coefficient(void) {
    TEST("reflection_coefficient");
    /* Z = 50Ω, Z0=50 → Γ = 0 */
    complex_t z = complex_make(50.0, 0.0);
    complex_t g = reflection_coefficient(z, 50.0);
    CHECK(complex_is_zero(g, TOL));

    /* Z = 100Ω, Z0=50 → Γ = (100-50)/(100+50) = 1/3 */
    complex_t z2 = complex_make(100.0, 0.0);
    complex_t g2 = reflection_coefficient(z2, 50.0);
    CHECK(fabs(g2.real - 1.0/3.0) < TOL);

    /* Z = 0 (short), Z0=50 → Γ = -1 */
    complex_t z3 = complex_make(0.0, 0.0);
    complex_t g3 = reflection_coefficient(z3, 50.0);
    CHECK(fabs(g3.real + 1.0) < TOL);
    TEST_END();
}

static void test_vswr(void) {
    TEST("VSWR from reflection coefficient");
    /* |Γ| = 1/3 → VSWR = (1+1/3)/(1-1/3) = (4/3)/(2/3) = 2 */
    complex_t g = complex_make(1.0/3.0, 0.0);
    double v = vswr_from_gamma(g);
    CHECK(fabs(v - 2.0) < TOL);

    /* |Γ| = 0 → VSWR = 1 */
    complex_t g2 = complex_make(0.0, 0.0);
    double v2 = vswr_from_gamma(g2);
    CHECK(fabs(v2 - 1.0) < TOL);
    TEST_END();
}

/* ================================================================
 * Maximum Power Transfer (L4)
 * ================================================================ */
static void test_conjugate_match(void) {
    TEST("conjugate_match — maximum power transfer");
    complex_t zs = complex_make(50.0, 30.0);  /* RS=50, XS=30 (inductive) */
    complex_t zl_opt = conjugate_match(zs);
    CHECK(fabs(zl_opt.real - 50.0) < TOL);
    CHECK(fabs(zl_opt.imag + 30.0) < TOL);  /* -30 (capacitive) */
    TEST_END();
}

/* ================================================================
 * Miller's Theorem (L2)
 * ================================================================ */
static void test_miller_split(void) {
    TEST("miller_split_impedance");
    /* Zf = j100 (capacitor), Av = -10
       Z1 = j100/(1-(-10)) = j100/11 ≈ j9.09
       Z2 = j100*(-10)/(-10-1) = -j1000/-11 ≈ j90.9 */
    complex_t zf = complex_make(0.0, 100.0);
    complex_t av = complex_make(-10.0, 0.0);
    complex_t z1, z2;
    miller_split_impedance(zf, av, &z1, &z2);

    CHECK(fabs(z1.imag - 100.0/11.0) < 0.01);
    CHECK(fabs(z2.imag - 1000.0/11.0) < 0.01);
    TEST_END();
}

/* ================================================================
 * Main
 * ================================================================ */
int main(void) {
    printf("=== Core Two-Port Network Tests ===\n");

    /* Complex arithmetic */
    test_complex_add();
    test_complex_mul();
    test_complex_div();
    test_complex_mag();
    test_complex_arg();
    test_complex_conj();
    test_complex_exp();
    test_complex_sqrt();

    /* Matrix operations */
    test_matrix_det();
    test_matrix_inv();
    test_matrix_mul();

    /* Z-parameters */
    test_zparams_t_network();
    test_zparams_rlc();
    test_zparams_input_impedance();
    test_zparams_to_t_network();

    /* Y-parameters */
    test_yparams_pi_network();
    test_yparams_mosfet_hf();

    /* H-parameters */
    test_hparams_bjt_ce();
    test_hparams_ce_to_cb();
    test_hparams_beta_cutoff();

    /* G-parameters */
    test_gparams_fet_cs_gain();

    /* Reflection and VSWR */
    test_reflection_coefficient();
    test_vswr();

    /* Maximum power transfer */
    test_conjugate_match();

    /* Miller's theorem */
    test_miller_split();

    printf("\n=== Results: %d/%d tests passed ===\n",
           tests_passed - tests_failed, tests_passed);
    return tests_failed > 0 ? 1 : 0;
}
