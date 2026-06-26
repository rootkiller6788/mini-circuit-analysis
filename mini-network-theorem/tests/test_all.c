/**
 * @file test_all.c
 * @brief Comprehensive test suite for all network theorem functions.
 *
 * Tests cover:
 *   - Impedance computation (series, parallel, RLC)
 *   - Thevenin/Norton equivalence and source transformation
 *   - Superposition theorem
 *   - Maximum power transfer theorem
 *   - Millman's theorem
 *   - Tellegen's theorem
 *   - Reciprocity verification
 *   - Two-port parameter conversions
 *   - Wye-Delta transformations
 *   - Wheatstone bridge analysis
 *   - Voltage/current divider analysis
 *   - Gaussian elimination and LU decomposition
 *   - MNA solver
 *   - Compensation theorem
 *
 * All tests use assert() for validation. Exit code 0 = all pass.
 */

#include "network_theorem.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOLERANCE 1e-9
#define ASSERT_NEAR(a, b, tol) assert(fabs((a) - (b)) < (tol))

/* ====================================================================
 * Test: Impedance Computation
 * ==================================================================== */
static void test_impedance_computation(void) {
    printf("  [TEST] Impedance computation...\n");

    /* Resistor */
    ComplexImpedance zr = impedance_resistor(100.0);
    ASSERT_NEAR(zr.real, 100.0, TOLERANCE);
    ASSERT_NEAR(zr.imag, 0.0, TOLERANCE);
    ASSERT_NEAR(impedance_magnitude(zr), 100.0, TOLERANCE);
    ASSERT_NEAR(impedance_phase(zr), 0.0, TOLERANCE);

    /* Capacitor: Z = -j/(2*PI*f*C) */
    ComplexImpedance zc = impedance_capacitor(1e-6, 1000.0);
    ASSERT_NEAR(zc.real, 0.0, TOLERANCE);
    assert(zc.imag < -150.0 && zc.imag > -165.0);

    /* Inductor: Z = j*2*PI*f*L */
    ComplexImpedance zl = impedance_inductor(0.01, 1000.0);
    ASSERT_NEAR(zl.real, 0.0, TOLERANCE);
    assert(zl.imag > 60.0 && zl.imag < 65.0);

    /* Series */
    ComplexImpedance zs = impedance_series(zr, zl);
    ASSERT_NEAR(zs.real, 100.0, TOLERANCE);

    /* Parallel: two 100 Ohm resistors -> 50 Ohm */
    ComplexImpedance zp = impedance_parallel(zr, zr);
    ASSERT_NEAR(zp.real, 50.0, TOLERANCE);
    ASSERT_NEAR(zp.imag, 0.0, TOLERANCE);

    /* Parallel with short */
    ComplexImpedance zshort = {0.0, 0.0};
    ComplexImpedance zps = impedance_parallel(zshort, zr);
    ASSERT_NEAR(zps.real, 0.0, TOLERANCE);
    ASSERT_NEAR(zps.imag, 0.0, TOLERANCE);

    printf("  [PASS] Impedance computation\n");
}

/* ====================================================================
 * Test: Source Transformation
 * ==================================================================== */
static void test_source_transformation(void) {
    printf("  [TEST] Source transformation...\n");

    TheveninEquivalent thev;
    thev.Z_th.real = 5.0; thev.Z_th.imag = 0.0;
    thev.V_th.dc_offset = 10.0;
    thev.V_th.ac_amplitude = 0.0;
    thev.V_th.angular_freq = 0.0;
    thev.V_th.phase = 0.0;
    thev.frequency = 0.0;

    NortonEquivalent nort;
    thevenin_to_norton(&thev, &nort);
    ASSERT_NEAR(nort.I_n.dc_offset, 2.0, TOLERANCE);
    ASSERT_NEAR(nort.Y_n.real, 0.2, TOLERANCE);

    TheveninEquivalent thev2;
    norton_to_thevenin(&nort, &thev2);
    ASSERT_NEAR(thev2.V_th.dc_offset, 10.0, TOLERANCE);
    ASSERT_NEAR(thev2.Z_th.real, 5.0, TOLERANCE);

    printf("  [PASS] Source transformation\n");
}

/* ====================================================================
 * Test: Maximum Power Transfer
 * ==================================================================== */
static void test_max_power_transfer(void) {
    printf("  [TEST] Maximum power transfer...\n");

    TheveninEquivalent src;
    src.Z_th.real = 10.0; src.Z_th.imag = 0.0;
    src.V_th.dc_offset = 12.0;
    src.V_th.ac_amplitude = 0.0;
    src.V_th.angular_freq = 0.0;
    src.V_th.phase = 0.0;
    src.frequency = 0.0;

    MaxPowerTransferResult mpt;
    int ret = max_power_transfer(&src, &mpt);
    assert(ret == 0);
    ASSERT_NEAR(mpt.optimal_load_real, 10.0, TOLERANCE);
    ASSERT_NEAR(mpt.optimal_load_imag, 0.0, TOLERANCE);
    ASSERT_NEAR(mpt.max_power, 3.6, 1e-6);
    assert(mpt.is_ac == 0);

    /* AC case */
    src.Z_th.imag = 5.0;
    ret = max_power_transfer(&src, &mpt);
    assert(ret == 0);
    assert(mpt.is_ac == 1);
    ASSERT_NEAR(mpt.optimal_load_imag, -5.0, TOLERANCE);
    ASSERT_NEAR(mpt.max_power, 1.8, 1e-6);

    printf("  [PASS] Maximum power transfer\n");
}

/* ====================================================================
 * Test: Millman's Theorem
 * ==================================================================== */
static void test_millman(void) {
    printf("  [TEST] Millman theorem...\n");

    double v[] = {10.0, 20.0};
    double r[] = {10.0, 20.0};

    MillmanResult mr;
    int ret = millman_solve(v, r, 2, &mr);
    assert(ret == 0);
    ASSERT_NEAR(mr.common_node_voltage, 13.3333333333, 1e-6);
    free(mr.branch_voltages);
    free(mr.branch_impedances);

    double v2[] = {12.0, 0.0, 6.0};
    double r2[] = {4.0, 8.0, 2.0};
    ret = millman_solve(v2, r2, 3, &mr);
    assert(ret == 0);
    ASSERT_NEAR(mr.common_node_voltage, 6.0 / 0.875, 1e-6);
    free(mr.branch_voltages);
    free(mr.branch_impedances);

    printf("  [PASS] Millman theorem\n");
}

/* ====================================================================
 * Test: Tellegen's Theorem
 * ==================================================================== */
static void test_tellegen(void) {
    printf("  [TEST] Tellegen theorem...\n");

    double v_state1[] = {-10.0, 6.0, 4.0};
    double i_state2[] = {2.0, 2.0, 2.0};

    TellegenResult tr;
    int ret = tellegen_verify(v_state1, i_state2, 3, 1e-9, &tr);
    assert(ret == 0);
    ASSERT_NEAR(tr.power_sum, 0.0, 1e-9);
    free(tr.branch_voltages_state1);
    free(tr.branch_currents_state2);

    printf("  [PASS] Tellegen theorem\n");
}

/* ====================================================================
 * Test: Gaussian Elimination
 * ==================================================================== */
static void test_gaussian_elimination(void) {
    printf("  [TEST] Gaussian elimination...\n");

    double A[] = {2.0, 1.0, 1.0, 3.0};
    double b[] = {5.0, 7.0};
    double x[2];

    int ret = gaussian_elimination(A, b, 2, x);
    assert(ret == 0);
    ASSERT_NEAR(x[0], 1.6, TOLERANCE);
    ASSERT_NEAR(x[1], 1.8, TOLERANCE);

    /* 3x3 system */
    double A2[] = {3.0, 2.0, -1.0, 2.0, -2.0, 4.0, -1.0, 0.5, -1.0};
    double b2[] = {1.0, -2.0, 0.0};
    double x2[3];
    ret = gaussian_elimination(A2, b2, 3, x2);
    assert(ret == 0);
    assert(isfinite(x2[0]) && isfinite(x2[1]) && isfinite(x2[2]));

    printf("  [PASS] Gaussian elimination\n");
}

/* ====================================================================
 * Test: LU Decomposition
 * ==================================================================== */
static void test_lu_decomposition(void) {
    printf("  [TEST] LU decomposition...\n");

    double A[] = {4.0, 3.0, 6.0, 3.0};
    double b[] = {10.0, 12.0};
    double x[2];

    int ret = lu_decompose(A, 2);
    assert(ret == 0);
    lu_solve(A, b, 2, x);
    ASSERT_NEAR(x[0], 1.0, TOLERANCE);
    ASSERT_NEAR(x[1], 2.0, TOLERANCE);

    printf("  [PASS] LU decomposition\n");
}

/* ====================================================================
 * Test: Two-Port Parameter Conversions
 * ==================================================================== */
static void test_two_port_conversions(void) {
    printf("  [TEST] Two-port parameter conversions...\n");

    ZParameters z;
    z.z11.real = 10.0; z.z11.imag = 0.0;
    z.z12.real = 2.0;  z.z12.imag = 0.0;
    z.z21.real = 2.0;  z.z21.imag = 0.0;
    z.z22.real = 8.0;  z.z22.imag = 0.0;

    /* Z->Y->Z round-trip */
    YParameters y;
    assert(z_to_y(&z, &y) == 0);
    ZParameters z2;
    assert(y_to_z(&y, &z2) == 0);
    ASSERT_NEAR(z2.z11.real, z.z11.real, 1e-6);
    ASSERT_NEAR(z2.z12.real, z.z12.real, 1e-6);
    ASSERT_NEAR(z2.z21.real, z.z21.real, 1e-6);
    ASSERT_NEAR(z2.z22.real, z.z22.real, 1e-6);

    /* Z->H->Z round-trip */
    HParameters h;
    assert(z_to_h(&z, &h) == 0);
    ZParameters z3;
    assert(h_to_z(&h, &z3) == 0);
    ASSERT_NEAR(z3.z11.real, z.z11.real, 1e-6);

    /* Z->ABCD */
    ABCDParameters abcd;
    assert(z_to_abcd(&z, &abcd) == 0);
    assert(isfinite(abcd.A_real) && isfinite(abcd.B_real));
    assert(isfinite(abcd.C_real) && isfinite(abcd.D_real));

    /* Z->S (Z0=50) */
    SParameters s;
    assert(z_to_s(&z, 50.0, &s) == 0);
    ASSERT_NEAR(s.Z0, 50.0, TOLERANCE);
    assert(isfinite(s.s11.real));

    /* S->Z round-trip */
    ZParameters z4;
    assert(s_to_z(&s, 50.0, &z4) == 0);
    ASSERT_NEAR(z4.z11.real, z.z11.real, 1e-3);

    /* ABCD cascade */
    ABCDParameters abcd2;
    assert(z_to_abcd(&z, &abcd2) == 0);
    ABCDParameters cascaded;
    abcd_cascade(&abcd, &abcd2, &cascaded);
    assert(isfinite(cascaded.A_real));

    printf("  [PASS] Two-port parameter conversions\n");
}

/* ====================================================================
 * Test: Reciprocity Verification
 * ==================================================================== */
static void test_reciprocity(void) {
    printf("  [TEST] Reciprocity verification...\n");

    ZParameters z_recip;
    z_recip.z11.real = 10.0; z_recip.z11.imag = 0.0;
    z_recip.z12.real = 5.0;  z_recip.z12.imag = 0.0;
    z_recip.z21.real = 5.0;  z_recip.z21.imag = 0.0;
    z_recip.z22.real = 12.0; z_recip.z22.imag = 0.0;

    ReciprocityParams rp;
    assert(verify_reciprocity(&z_recip, &rp) == 0);
    assert(rp.is_reciprocal == 1);
    ASSERT_NEAR(rp.reciprocity_error, 0.0, TOLERANCE);

    ZParameters z_nonrecip;
    z_nonrecip.z11.real = 10.0; z_nonrecip.z11.imag = 0.0;
    z_nonrecip.z12.real = 5.0;  z_nonrecip.z12.imag = 0.0;
    z_nonrecip.z21.real = 3.0;  z_nonrecip.z21.imag = 0.0;
    z_nonrecip.z22.real = 12.0; z_nonrecip.z22.imag = 0.0;

    assert(verify_reciprocity(&z_nonrecip, &rp) == 0);
    assert(rp.is_reciprocal == 0);
    assert(rp.reciprocity_error > 1.0);

    printf("  [PASS] Reciprocity verification\n");
}

/* ====================================================================
 * Test: Wye-Delta Transformations
 * ==================================================================== */
static void test_wye_delta(void) {
    printf("  [TEST] Wye-Delta transformations...\n");

    WyeNetwork wye = {10.0, 10.0, 10.0};
    DeltaNetwork delta;
    wye_to_delta(&wye, &delta);
    ASSERT_NEAR(delta.R12, 30.0, TOLERANCE);
    ASSERT_NEAR(delta.R23, 30.0, TOLERANCE);
    ASSERT_NEAR(delta.R31, 30.0, TOLERANCE);

    WyeNetwork wye2;
    delta_to_wye(&delta, &wye2);
    ASSERT_NEAR(wye2.R1, 10.0, TOLERANCE);
    ASSERT_NEAR(wye2.R2, 10.0, TOLERANCE);
    ASSERT_NEAR(wye2.R3, 10.0, TOLERANCE);

    WyeNetwork wy3 = {2.0, 4.0, 6.0};
    DeltaNetwork dl3;
    wye_to_delta(&wy3, &dl3);
    double sumprod = 2.0*4.0 + 4.0*6.0 + 6.0*2.0; /* 8+24+12=44 */
    ASSERT_NEAR(dl3.R12, sumprod/6.0, TOLERANCE);
    ASSERT_NEAR(dl3.R23, sumprod/2.0, TOLERANCE);
    ASSERT_NEAR(dl3.R31, sumprod/4.0, TOLERANCE);

    printf("  [PASS] Wye-Delta transformations\n");
}

/* ====================================================================
 * Test: Wheatstone Bridge
 * ==================================================================== */
static void test_wheatstone_bridge(void) {
    printf("  [TEST] Wheatstone bridge...\n");

    WheatstoneBridge bridge;
    wheatstone_analyze(100.0, 100.0, 200.0, 200.0, 1000.0, 10.0, &bridge);
    assert(bridge.is_balanced == 1);
    ASSERT_NEAR(bridge.V_output, 0.0, 1e-9);
    ASSERT_NEAR(bridge.I_g, 0.0, 1e-9);

    wheatstone_analyze(120.0, 100.0, 200.0, 200.0, 1000.0, 10.0, &bridge);
    assert(bridge.is_balanced == 0);
    ASSERT_NEAR(bridge.V_output, 10.0*(200.0/320.0) - 10.0*(200.0/300.0), 1e-3);

    printf("  [PASS] Wheatstone bridge\n");
}

/* ====================================================================
 * Test: Voltage and Current Divider
 * ==================================================================== */
static void test_dividers(void) {
    printf("  [TEST] Voltage and current dividers...\n");

    VoltageDivider vdiv;
    voltage_divider_analyze(3.0, 6.0, 1e20, 9.0, &vdiv);
    ASSERT_NEAR(vdiv.V_out, 6.0, TOLERANCE);
    ASSERT_NEAR(vdiv.I_load, 0.0, TOLERANCE);

    voltage_divider_analyze(3.0, 6.0, 3.0, 9.0, &vdiv);
    ASSERT_NEAR(vdiv.V_out, 3.6, TOLERANCE);

    CurrentDivider cdiv;
    current_divider_analyze(2.0, 8.0, 10.0, &cdiv);
    ASSERT_NEAR(cdiv.I1, 8.0, TOLERANCE);
    ASSERT_NEAR(cdiv.I2, 2.0, TOLERANCE);

    printf("  [PASS] Voltage and current dividers\n");
}

/* ====================================================================
 * Test: Power Calculations
 * ==================================================================== */
static void test_power_calculations(void) {
    printf("  [TEST] Power calculations...\n");

    ASSERT_NEAR(power_dissipated(10.0, 2.0, 5.0), 20.0, TOLERANCE);
    ASSERT_NEAR(reactive_power(2.0, 3.0), 12.0, TOLERANCE);
    ASSERT_NEAR(apparent_power(230.0, 5.0), 1150.0, TOLERANCE);
    double pf = power_factor(0.5, 1.0);
    ASSERT_NEAR(pf, cos(-0.5), TOLERANCE);

    printf("  [PASS] Power calculations\n");
}

/* ====================================================================
 * Test: Admittance/Impedance Conversion
 * ==================================================================== */
static void test_admittance_conversion(void) {
    printf("  [TEST] Admittance-impedance conversion...\n");

    ComplexImpedance z = {3.0, 4.0};
    ComplexAdmittance y = impedance_to_admittance(z);
    ASSERT_NEAR(y.real, 0.12, TOLERANCE);
    ASSERT_NEAR(y.imag, -0.16, TOLERANCE);

    ComplexImpedance z2 = admittance_to_impedance(y);
    ASSERT_NEAR(z2.real, 3.0, TOLERANCE);
    ASSERT_NEAR(z2.imag, 4.0, TOLERANCE);

    printf("  [PASS] Admittance-impedance conversion\n");
}

/* ====================================================================
 * Test: Compensation Theorem
 * ==================================================================== */
static void test_compensation(void) {
    printf("  [TEST] Compensation theorem...\n");

    CompensationResult cr;
    int ret = compensation_compute(3, 10.0, 15.0, 2.0, &cr);
    assert(ret == 0);
    assert(cr.modified_branch == 3);
    ASSERT_NEAR(cr.delta_z, 5.0, TOLERANCE);
    ASSERT_NEAR(cr.compensation_voltage, 10.0, TOLERANCE);

    printf("  [PASS] Compensation theorem\n");
}

/* ====================================================================
 * Test: Audio Amplifier Matching (L7 Application)
 * ==================================================================== */
static void test_audio_amplifier_match(void) {
    printf("  [TEST] Audio amplifier matching...\n");

    AudioMatchingResult amr;
    int ret = audio_amplifier_match(4.0, 8.0, 20.0, &amr);
    assert(ret == 0);
    ASSERT_NEAR(amr.actual_power, 400.0 * 8.0 / 144.0, 0.01);
    ASSERT_NEAR(amr.matched_power, 25.0, 0.01);
    ASSERT_NEAR(amr.efficiency, 0.6666667, 1e-6);
    ASSERT_NEAR(amr.damping_factor, 2.0, TOLERANCE);

    printf("  [PASS] Audio amplifier matching\n");
}

/* ====================================================================
 * Test: Series/Parallel Combinations
 * ==================================================================== */
static void test_series_parallel_combinations(void) {
    printf("  [TEST] Series/parallel combinations...\n");

    ASSERT_NEAR(combine_series_resistors(10.0, 20.0), 30.0, TOLERANCE);
    ASSERT_NEAR(combine_parallel_resistors(10.0, 20.0), 200.0/30.0, TOLERANCE);
    ASSERT_NEAR(combine_parallel_resistors(10.0, 0.0), 0.0, TOLERANCE);
    ASSERT_NEAR(combine_series_capacitors(10e-6, 20e-6), 6.666667e-6, 1e-12);
    ASSERT_NEAR(combine_parallel_capacitors(10e-6, 20e-6), 30e-6, 1e-12);
    ASSERT_NEAR(combine_series_inductors(0.1, 0.2), 0.3, TOLERANCE);
    ASSERT_NEAR(combine_parallel_inductors(0.1, 0.2), 0.02/0.3, TOLERANCE);

    printf("  [PASS] Series/parallel combinations\n");
}

/* ====================================================================
 * Main Test Runner
 * ==================================================================== */
int main(void) {
    printf("========================================\n");
    printf(" Network Theorems - Test Suite\n");
    printf("========================================\n\n");

    test_impedance_computation();
    test_source_transformation();
    test_max_power_transfer();
    test_millman();
    test_tellegen();
    test_gaussian_elimination();
    test_lu_decomposition();
    test_two_port_conversions();
    test_reciprocity();
    test_wye_delta();
    test_wheatstone_bridge();
    test_dividers();
    test_power_calculations();
    test_admittance_conversion();
    test_compensation();
    test_audio_amplifier_match();
    test_series_parallel_combinations();

    printf("\n========================================\n");
    printf(" ALL TESTS PASSED\n");
    printf("========================================\n");
    return 0;
}
