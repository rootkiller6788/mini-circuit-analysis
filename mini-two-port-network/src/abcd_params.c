/**
 * @file abcd_params.c
 * @brief ABCD (Transmission) Parameter Implementation
 *
 * ABCD parameters are the cascade-multiplicative representation,
 * essential for transmission line, filter, and RF chain analysis.
 */

#include "../include/abcd_params.h"
#include <float.h>

/* ============================================================================
 * L1: ABCD Creation Functions
 * ============================================================================ */

matrix2x2_t abcd_create(complex_t a, complex_t b, complex_t c, complex_t d) {
    return matrix2x2_make(a, b, c, d);
}

/**
 * Series impedance Z: ABCD = [[1, Z], [0, 1]].
 *
 * Verification:
 *   V1 = 1*V2 + Z*(-I2) = V2 - Z*I2  (voltage drop across Z)
 *   I1 = 0*V2 + 1*(-I2) = -I2        (same current through series element)
 * ✓
 */
matrix2x2_t abcd_series_z(complex_t z) {
    return matrix2x2_make(
        complex_make(1.0, 0.0), z,
        complex_make(0.0, 0.0), complex_make(1.0, 0.0)
    );
}

/**
 * Shunt admittance Y: ABCD = [[1, 0], [Y, 1]].
 *
 * Verification:
 *   V1 = 1*V2 + 0*(-I2) = V2        (same voltage)
 *   I1 = Y*V2 + 1*(-I2) = Y*V1 - I2  (current divides)
 *   I1 + I2 = Y*V1                   (the difference goes through Y)
 * ✓
 */
matrix2x2_t abcd_shunt_y(complex_t y) {
    return matrix2x2_make(
        complex_make(1.0, 0.0), complex_make(0.0, 0.0),
        y, complex_make(1.0, 0.0)
    );
}

/**
 * Ideal transformer: ABCD = [[1/n, 0], [0, n]].
 *
 *   V1 = V2/n  (step-down if n < 1)
 *   I1 = -n*I2  (current steps up if voltage steps down)
 */
matrix2x2_t abcd_ideal_transformer(double n) {
    if (n < 1e-30) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }
    return matrix2x2_make(
        complex_make(1.0 / n, 0.0), complex_make(0.0, 0.0),
        complex_make(0.0, 0.0), complex_make(n, 0.0)
    );
}

/**
 * Transmission line ABCD parameters.
 *
 * General (lossy) line:
 *   γ = α + jβ  (propagation constant)
 *   ABCD = [[cosh(γl), Z0*sinh(γl)], [sinh(γl)/Z0, cosh(γl)]]
 *
 * cosh(a+jb) = cosh(a)cos(b) + j*sinh(a)sin(b)
 * sinh(a+jb) = sinh(a)cos(b) + j*cosh(a)sin(b)
 *
 * Lossless case (α = 0):
 *   cosh(jβl) = cos(βl)
 *   sinh(jβl) = j*sin(βl)
 *   ABCD = [[cos(βl), jZ0*sin(βl)], [j*sin(βl)/Z0, cos(βl)]]
 *
 * Special cases:
 *   l = λ/4 (βl = π/2): ABCD = [[0, jZ0], [j/Z0, 0]]
 *     → Zin = Z0²/ZL (quarter-wave transformer)
 *   l = λ/2 (βl = π): ABCD = [[-1, 0], [0, -1]]
 *     → Zin = ZL (half-wave repeats impedance)
 *
 * Course: ETH 227-0455 — Transmission line ABCD
 * Ref: Pozar §2.1
 */
matrix2x2_t abcd_transmission_line(double z0, double alpha,
                                    double beta, double length) {
    double alpha_l = alpha * length;
    double beta_l = beta * length;

    /* cosh(αl + jβl) = cosh(αl)cos(βl) + j*sinh(αl)sin(βl) */
    double cosh_al = cosh(alpha_l);
    double sinh_al = sinh(alpha_l);
    double cos_bl = cos(beta_l);
    double sin_bl = sin(beta_l);

    complex_t cosh_gl = complex_make(cosh_al * cos_bl, sinh_al * sin_bl);
    complex_t sinh_gl = complex_make(sinh_al * cos_bl, cosh_al * sin_bl);

    complex_t a = cosh_gl;
    complex_t b = complex_make(z0 * sinh_gl.real, z0 * sinh_gl.imag);
    complex_t c = complex_make(sinh_gl.real / z0, sinh_gl.imag / z0);
    complex_t d = cosh_gl;

    return matrix2x2_make(a, b, c, d);
}

/**
 * L-network: ABCD = series × shunt.
 *
 *   ABCD = [[1, Z_series], [0, 1]] × [[1, 0], [Y_shunt, 1]]
 *        = [[1 + Z_series*Y_shunt, Z_series], [Y_shunt, 1]]
 */
matrix2x2_t abcd_l_network(complex_t z_series, complex_t y_shunt) {
    matrix2x2_t series_abcd = abcd_series_z(z_series);
    matrix2x2_t shunt_abcd = abcd_shunt_y(y_shunt);
    return matrix2x2_mul(series_abcd, shunt_abcd);
}

/**
 * π-network: shunt-sha → series → shunt-shc.
 *
 * ABCD = [[1,0],[Ya,1]] × [[1,1/Yb],[0,1]] × [[1,0],[Yc,1]]
 *
 * Full expansion:
 *   A = 1 + Yc/Yb
 *   B = 1/Yb
 *   C = Ya + Yc + Ya*Yc/Yb
 *   D = 1 + Ya/Yb
 */
matrix2x2_t abcd_pi_network(complex_t ya, complex_t yb, complex_t yc) {
    if (complex_is_zero(yb, 1e-30)) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }

    complex_t yb_inv = complex_div(complex_make(1.0, 0.0), yb);

    complex_t a = complex_add(complex_make(1.0, 0.0), complex_mul(yc, yb_inv));
    complex_t b = yb_inv;
    complex_t c = complex_add(
        complex_add(ya, yc),
        complex_div(complex_mul(ya, yc), yb)
    );
    complex_t d = complex_add(complex_make(1.0, 0.0), complex_mul(ya, yb_inv));

    return matrix2x2_make(a, b, c, d);
}

/**
 * T-network: series-Za → shunt-Zc → series-Zb.
 *
 * ABCD = [[1,Za],[0,1]] × [[1,0],[1/Zc,1]] × [[1,Zb],[0,1]]
 *
 * Full expansion:
 *   A = 1 + Za/Zc
 *   B = Za + Zb + Za*Zb/Zc
 *   C = 1/Zc
 *   D = 1 + Zb/Zc
 */
matrix2x2_t abcd_t_network(complex_t za, complex_t zb, complex_t zc) {
    if (complex_is_zero(zc, 1e-30)) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }

    complex_t zc_inv = complex_div(complex_make(1.0, 0.0), zc);

    complex_t a = complex_add(complex_make(1.0, 0.0), complex_mul(za, zc_inv));
    complex_t b = complex_add(
        complex_add(za, zb),
        complex_div(complex_mul(za, zb), zc)
    );
    complex_t c = zc_inv;
    complex_t d = complex_add(complex_make(1.0, 0.0), complex_mul(zb, zc_inv));

    return matrix2x2_make(a, b, c, d);
}

/* ============================================================================
 * L3: ABCD Analysis Functions
 * ============================================================================ */

/**
 * Voltage transfer: H = V2/V1 = ZL / (A*ZL + B).
 *
 * Derivation:
 *   V1 = A*V2 + B*(-I2)
 *   I1 = C*V2 + D*(-I2)
 *   V2 = -ZL*I2 → I2 = -V2/ZL
 *
 *   V1 = A*V2 + B*V2/ZL = V2*(A + B/ZL)
 *   V2/V1 = 1/(A + B/ZL) = ZL/(A*ZL + B) ✓
 */
complex_t abcd_voltage_transfer(matrix2x2_t abcd, complex_t zl) {
    complex_t denom = complex_add(
        complex_mul(abcd.m11, zl),
        abcd.m12
    );
    return complex_div(zl, denom);
}

/**
 * Input impedance: Zin = (A*ZL + B) / (C*ZL + D).
 *
 * Derivation:
 *   Zin = V1/I1 = (A*V2 - B*I2) / (C*V2 - D*I2)
 *   Substitute V2 = -ZL*I2:
 *   Zin = (A*(-ZL*I2) - B*I2) / (C*(-ZL*I2) - D*I2)
 *       = (-A*ZL*I2 - B*I2) / (-C*ZL*I2 - D*I2)
 *       = (A*ZL + B) / (C*ZL + D)  ✓
 *
 * This is the fundamental impedance transformation formula.
 * It encompasses all transmission line and matching network behavior.
 */
complex_t abcd_input_impedance(matrix2x2_t abcd, complex_t zl) {
    complex_t num = complex_add(
        complex_mul(abcd.m11, zl),
        abcd.m12
    );
    complex_t denom = complex_add(
        complex_mul(abcd.m21, zl),
        abcd.m22
    );
    return complex_div(num, denom);
}

/**
 * Image impedance at port 1: Zi1 = sqrt(A*B / (C*D)).
 *
 * The image impedance is defined such that if port 2 is terminated in Zi2,
 * the input impedance at port 1 is Zi1.
 *
 * Zi1² = (A*B) / (C*D)
 * Zi2² = (D*B) / (C*A)
 *
 * Image parameters form the basis of classical filter design (constant-k
 * and m-derived filters). The image impedance is the iterative impedance
 * for an infinite cascade.
 */
complex_t abcd_image_impedance_1(matrix2x2_t abcd) {
    complex_t num = complex_mul(abcd.m11, abcd.m12);
    complex_t denom = complex_mul(abcd.m21, abcd.m22);
    complex_t ratio = complex_div(num, denom);
    return complex_sqrt(ratio);
}

complex_t abcd_image_impedance_2(matrix2x2_t abcd) {
    complex_t num = complex_mul(abcd.m22, abcd.m12);
    complex_t denom = complex_mul(abcd.m21, abcd.m11);
    complex_t ratio = complex_div(num, denom);
    return complex_sqrt(ratio);
}

/**
 * Image transfer constant: θ = acosh(sqrt(A*D)).
 *
 * For a symmetric network (A = D): θ = acosh(A).
 * cosh(θ) = (A+D)/2 for the cascade.
 *
 * θ = α + jβ where:
 *   α = attenuation constant (Np/section)
 *   β = phase constant (rad/section)
 *
 * Filter passband: α = 0, β ∈ (0, π).
 * Filter stopband: α > 0, β = 0 or π.
 */
complex_t abcd_image_transfer_constant(matrix2x2_t abcd) {
    complex_t ad = complex_mul(abcd.m11, abcd.m22);
    complex_t sqrt_ad = complex_sqrt(ad);
    /* cosh(θ) = sqrt(A*D) for image-matched cascade */
    /* θ = acosh(sqrt(A*D)) = ln(sqrt(A*D) + sqrt(A*D-1)) */
    complex_t ad_minus_1 = complex_sub(sqrt_ad, complex_make(1.0, 0.0));
    complex_t sqrt_ad_m1 = complex_sqrt(ad_minus_1);
    complex_t sum = complex_add(sqrt_ad, sqrt_ad_m1);

    /* ln(z) = ln|z| + j*arg(z) */
    return complex_make(log(complex_mag(sum)), complex_arg(sum));
}

complex_t abcd_propagation_constant(matrix2x2_t abcd) {
    /* cosh(γ) = (A+D)/2 */
    complex_t half_sum = complex_make(
        (abcd.m11.real + abcd.m22.real) / 2.0,
        (abcd.m11.imag + abcd.m22.imag) / 2.0
    );
    /* acosh(z) = ln(z + sqrt(z²-1)) */
    complex_t z_sq = complex_mul(half_sum, half_sum);
    complex_t z_sq_m1 = complex_sub(z_sq, complex_make(1.0, 0.0));
    complex_t sqrt_term = complex_sqrt(z_sq_m1);
    complex_t sum = complex_add(half_sum, sqrt_term);
    return complex_make(log(complex_mag(sum)), complex_arg(sum));
}

int abcd_is_reciprocal(matrix2x2_t abcd, double tolerance) {
    complex_t det = matrix2x2_det(abcd);
    return complex_approx_equal(det, complex_make(1.0, 0.0), tolerance);
}

/**
 * Insertion loss from ABCD in a Z0 system.
 *
 * S21 = 2 / (A + B/Z0 + C*Z0 + D)
 * IL(dB) = -20*log10(|S21|)
 */
double abcd_insertion_loss_db(matrix2x2_t abcd, double z0) {
    double z0_inv = 1.0 / z0;
    complex_t a_plus_bz0 = complex_add(
        abcd.m11,
        complex_make(abcd.m12.real * z0_inv, abcd.m12.imag * z0_inv)
    );
    complex_t cz0_plus_d = complex_add(
        complex_make(abcd.m21.real * z0, abcd.m21.imag * z0),
        abcd.m22
    );
    complex_t denom = complex_add(a_plus_bz0, cz0_plus_d);
    complex_t s21 = complex_div(complex_make(2.0, 0.0), denom);
    double mag_s21 = complex_mag(s21);
    if (mag_s21 < 1e-30) return INFINITY;
    return -20.0 * log10(mag_s21);
}

/* ============================================================================
 * L4: Cascade Analysis
 * ============================================================================ */

/**
 * Cascade n networks: ABCD_total = ABCD_1 × ABCD_2 × ... × ABCD_n.
 *
 * This is the unique power of ABCD parameters — cascading is just
 * matrix multiplication. For n stages, this performs n-1 matrix
 * multiplications (O(n)).
 *
 * Application: A complete RF chain with n components.
 */
matrix2x2_t abcd_cascade_n(const matrix2x2_t *stages, int n) {
    if (n <= 0) {
        /* Identity matrix = direct connection (wire) */
        return matrix2x2_make(
            complex_make(1.0, 0.0), complex_make(0.0, 0.0),
            complex_make(0.0, 0.0), complex_make(1.0, 0.0)
        );
    }

    matrix2x2_t result = stages[0];
    for (int i = 1; i < n; i++) {
        result = matrix2x2_mul(result, stages[i]);
    }
    return result;
}

/**
 * De-embedding: DUT = left⁻¹ × total × right⁻¹.
 *
 * This is how VNA measurements are processed. The total measured
 * ABCD matrix includes the DUT plus cables, connectors, and fixtures.
 * By measuring the left and right fixtures separately (using calibration
 * standards: Short, Open, Load, Thru — SOLT), we can mathematically
 * remove their effects.
 *
 * More sophisticated calibration methods (TRL, LRM) use different
 * standards but the same de-embedding principle.
 *
 * Course: ETH 227-0455 — VNA calibration and de-embedding
 */
matrix2x2_t abcd_deembed(matrix2x2_t abcd_total,
                          matrix2x2_t abcd_left,
                          matrix2x2_t abcd_right) {
    matrix2x2_t left_inv = matrix2x2_inv(abcd_left);
    matrix2x2_t right_inv = matrix2x2_inv(abcd_right);
    matrix2x2_t temp = matrix2x2_mul(left_inv, abcd_total);
    return matrix2x2_mul(temp, right_inv);
}
