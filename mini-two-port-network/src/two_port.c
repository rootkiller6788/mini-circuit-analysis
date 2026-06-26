/**
 * @file two_port.c
 * @brief Core Two-Port Network Implementation
 *
 * Implements complex arithmetic, 2x2 matrix operations, and fundamental
 * two-port analysis functions (input/output impedance, gain, stability
 * property checks, Miller's theorem, maximum power transfer).
 *
 * All functions are designed to teach individual concepts from the
 * nine-school circuit analysis curriculum.
 */

#include "../include/two_port.h"
#include <stdio.h>
#include <string.h>
#include <float.h>

/* ============================================================================
 * L3: Complex Number Arithmetic
 * ============================================================================ */

/**
 * Complex addition: z1 + z2 = (a1+a2) + j(b1+b2)
 * Reference: Kreyszig, "Advanced Engineering Mathematics", Ch. 13
 * Complexity: O(1)
 */
complex_t complex_make(double real, double imag) {
    complex_t z;
    z.real = real;
    z.imag = imag;
    return z;
}

complex_t complex_add(complex_t a, complex_t b) {
    return complex_make(a.real + b.real, a.imag + b.imag);
}

complex_t complex_sub(complex_t a, complex_t b) {
    return complex_make(a.real - b.real, a.imag - b.imag);
}

/**
 * Complex multiplication using the FOIL method:
 * (a+jb)(c+jd) = ac + jad + jbc + j²bd = (ac-bd) + j(ad+bc)
 * The sign inversion on bd comes from j² = -1.
 */
complex_t complex_mul(complex_t a, complex_t b) {
    return complex_make(
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real
    );
}

/**
 * Complex division via rationalization:
 * (a+jb)/(c+jd) = (a+jb)(c-jd) / (c²+d²)
 *               = (ac+bd)/(c²+d²) + j(bc-ad)/(c²+d²)
 *
 * Special case: denominator zero → return NAN
 */
complex_t complex_div(complex_t a, complex_t b) {
    double denom = b.real * b.real + b.imag * b.imag;
    if (denom < 1e-30) {
        /* Division by zero — return NaN to signal error */
        return complex_make(NAN, NAN);
    }
    return complex_make(
        (a.real * b.real + a.imag * b.imag) / denom,
        (a.imag * b.real - a.real * b.imag) / denom
    );
}

/**
 * Complex magnitude: |z| = sqrt(a² + b²)
 * This is the Euclidean norm in the complex plane.
 * For engineering: |Z| gives the magnitude of impedance, |S21| gives gain.
 */
double complex_mag(complex_t z) {
    return sqrt(z.real * z.real + z.imag * z.imag);
}

/**
 * Complex phase angle: φ = atan2(imag, real)
 * Returns the principal argument in [-π, π].
 * Key properties:
 *   - atan2(0, positive) = 0 (purely resistive)
 *   - atan2(positive, 0) = π/2 (purely inductive)
 *   - atan2(negative, 0) = -π/2 (purely capacitive)
 */
double complex_arg(complex_t z) {
    return atan2(z.imag, z.real);
}

complex_t complex_conj(complex_t z) {
    return complex_make(z.real, -z.imag);
}

/**
 * Complex exponential: e^(a+jb) = e^a * (cos(b) + j sin(b))
 * This is Euler's formula combined with real exponentiation.
 * Applications:
 *   - Time-harmonic fields: E = E0 * e^(jωt)
 *   - Phasor rotation: e^(jθ) advances phase by θ
 *   - Damped oscillations: e^(σt) * e^(jωt) = e^((σ+jω)t)
 *
 * Reference: Euler, "Introductio in analysin infinitorum" (1748)
 * Course: MIT 6.003 — Euler's formula in signal processing
 */
complex_t complex_exp(complex_t z) {
    double mag = exp(z.real);
    return complex_make(mag * cos(z.imag), mag * sin(z.imag));
}

/**
 * Complex from polar coordinates: r∠θ → r*cos(θ) + j*r*sin(θ)
 * This is the inverse of the magnitude/phase representation.
 */
complex_t complex_polar(double r, double theta) {
    return complex_make(r * cos(theta), r * sin(theta));
}

/**
 * Complex square root using the algebraic method.
 * For z = a + jb:
 *   sqrt(z) = ±(sqrt((|z|+a)/2) + j*sign(b)*sqrt((|z|-a)/2))
 *
 * We return the principal branch (positive real part, or positive
 * imaginary part if real part is zero).
 *
 * Verification: (x+jy)² = (x²-y²) + j(2xy) = a + jb
 * Solving: x²-y²=a, 2xy=b → x = sqrt((|z|+a)/2), y = sign(b)*sqrt((|z|-a)/2)
 *
 * Application: Reflection coefficient → impedance conversion involves sqrt.
 */
complex_t complex_sqrt(complex_t z) {
    double mag = complex_mag(z);
    if (mag < 1e-30) return complex_make(0.0, 0.0);

    double x = sqrt((mag + z.real) / 2.0);
    double y = sqrt((mag - z.real) / 2.0);
    if (z.imag < 0) y = -y;
    return complex_make(x, y);
}

/* ============================================================================
 * L3: 2x2 Matrix Operations
 * ============================================================================ */

/**
 * Create a 2x2 complex matrix from four elements.
 * The layout follows the standard two-port sign convention:
 *   [m11  m12]   where m11=input-side, m12=reverse, m21=forward, m22=output
 *   [m21  m22]
 */
matrix2x2_t matrix2x2_make(complex_t m11, complex_t m12,
                             complex_t m21, complex_t m22) {
    matrix2x2_t m;
    m.m11 = m11;
    m.m12 = m12;
    m.m21 = m21;
    m.m22 = m22;
    return m;
}

/**
 * Element-wise matrix addition.
 * Used for adding Z-params of series-series interconnected networks.
 * Matrix addition is commutative: A + B = B + A.
 */
matrix2x2_t matrix2x2_add(matrix2x2_t a, matrix2x2_t b) {
    return matrix2x2_make(
        complex_add(a.m11, b.m11),
        complex_add(a.m12, b.m12),
        complex_add(a.m21, b.m21),
        complex_add(a.m22, b.m22)
    );
}

matrix2x2_t matrix2x2_sub(matrix2x2_t a, matrix2x2_t b) {
    return matrix2x2_make(
        complex_sub(a.m11, b.m11),
        complex_sub(a.m12, b.m12),
        complex_sub(a.m21, b.m21),
        complex_sub(a.m22, b.m22)
    );
}

/**
 * 2x2 matrix multiplication: C = A × B.
 *
 * The fundamental formula (each element is a dot product of row × column):
 *   c11 = a11*b11 + a12*b21
 *   c12 = a11*b12 + a12*b22
 *   c21 = a21*b11 + a22*b21
 *   c22 = a21*b12 + a22*b22
 *
 * Key property: Matrix multiplication is NOT commutative.
 * This is why cascade order matters: ABCD1 × ABCD2 ≠ ABCD2 × ABCD1.
 *
 * Complexity: 8 complex multiplications + 4 complex additions = O(1)
 * (since the matrix size is fixed at 2×2)
 *
 * Course: Berkeley EE16A — Matrix multiplication as composition of
 *         linear transformations
 */
matrix2x2_t matrix2x2_mul(matrix2x2_t a, matrix2x2_t b) {
    complex_t c11 = complex_add(
        complex_mul(a.m11, b.m11),
        complex_mul(a.m12, b.m21)
    );
    complex_t c12 = complex_add(
        complex_mul(a.m11, b.m12),
        complex_mul(a.m12, b.m22)
    );
    complex_t c21 = complex_add(
        complex_mul(a.m21, b.m11),
        complex_mul(a.m22, b.m21)
    );
    complex_t c22 = complex_add(
        complex_mul(a.m21, b.m12),
        complex_mul(a.m22, b.m22)
    );
    return matrix2x2_make(c11, c12, c21, c22);
}

/**
 * Scalar multiplication: k * A.
 * Each element is scaled by the complex number k.
 * This corresponds to multiplying all impedances by a common factor.
 */
matrix2x2_t matrix2x2_scale(complex_t k, matrix2x2_t a) {
    return matrix2x2_make(
        complex_mul(k, a.m11),
        complex_mul(k, a.m12),
        complex_mul(k, a.m21),
        complex_mul(k, a.m22)
    );
}

/**
 * Determinant of a 2×2 matrix: det = m11*m22 - m12*m21.
 *
 * The determinant has several physical meanings in two-port theory:
 *   - Z-params: ΔZ = z11*z22 - z12*z21 appears in Z↔Y conversion
 *   - ABCD-params: AD-BC = 1 for reciprocal networks
 *   - S-params: Δ = s11*s22 - s12*s21 appears in stability formulas
 *
 * A zero determinant indicates a singular (non-invertible) matrix,
 * which corresponds to degenerate networks (e.g., short/open circuits).
 */
complex_t matrix2x2_det(matrix2x2_t a) {
    return complex_sub(
        complex_mul(a.m11, a.m22),
        complex_mul(a.m12, a.m21)
    );
}

complex_t matrix2x2_trace(matrix2x2_t a) {
    return complex_add(a.m11, a.m22);
}

/**
 * Matrix inverse for 2×2: A⁻¹ = (1/det) * [[m22, -m12], [-m21, m11]].
 *
 * The formula follows directly from Cramer's rule for 2×2 matrices.
 *
 * For Z↔Y conversion: Y = Z⁻¹, meaning:
 *   y11 = z22/ΔZ,  y12 = -z12/ΔZ
 *   y21 = -z21/ΔZ,  y22 = z11/ΔZ
 *
 * Singular case (det ≈ 0): Returns a matrix filled with NAN.
 * This happens for:
 *   - Ideal FET at DC (infinite input impedance → non-invertible Z)
 *   - Short-circuited ports
 *   - Networks with dependent sources that make the Z-matrix singular
 *
 * Course: Berkeley EE16A — Matrix inversion in circuit analysis
 */
matrix2x2_t matrix2x2_inv(matrix2x2_t a) {
    complex_t det = matrix2x2_det(a);
    if (complex_is_zero(det, 1e-30)) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }
    complex_t inv_det = complex_div(complex_make(1.0, 0.0), det);
    return matrix2x2_make(
        complex_mul(inv_det, a.m22),
        complex_mul(inv_det, complex_make(-a.m12.real, -a.m12.imag)),
        complex_mul(inv_det, complex_make(-a.m21.real, -a.m21.imag)),
        complex_mul(inv_det, a.m11)
    );
}

matrix2x2_t matrix2x2_transpose(matrix2x2_t a) {
    return matrix2x2_make(a.m11, a.m21, a.m12, a.m22);
}

/**
 * Hermitian conjugate (conjugate transpose): A^H.
 * For S-parameters, S^H * S = I for lossless networks (unitary condition).
 * This is the microwave equivalent of energy conservation.
 */
matrix2x2_t matrix2x2_hermitian(matrix2x2_t a) {
    return matrix2x2_make(
        complex_conj(a.m11), complex_conj(a.m21),
        complex_conj(a.m12), complex_conj(a.m22)
    );
}

/**
 * Symmetry check: |a12 - a21| < tolerance.
 * For Z-parameters: reciprocal if z12 = z21.
 * Passive RLC networks are always reciprocal (Maxwell reciprocity).
 * Active networks (transistors, op-amps) are typically non-reciprocal.
 */
int matrix2x2_is_symmetric(matrix2x2_t a, double tolerance) {
    complex_t diff = complex_sub(a.m12, a.m21);
    return complex_approx_equal(diff, complex_make(0, 0), tolerance);
}

int matrix2x2_is_reciprocal_z(matrix2x2_t z, double tolerance) {
    return matrix2x2_is_symmetric(z, tolerance);
}

int matrix2x2_is_reciprocal_y(matrix2x2_t y, double tolerance) {
    return matrix2x2_is_symmetric(y, tolerance);
}

/**
 * ABCD reciprocal condition: det(ABCD) = AD - BC = 1.
 * This is the fundamental test for a reciprocal transmission network.
 * Networks with dependent sources violate this condition.
 *
 * For a simple series Z: A=1, B=Z, C=0, D=1 → AD-BC = 1·1 - Z·0 = 1 ✓
 * For a simple shunt Y: A=1, B=0, C=Y, D=1 → AD-BC = 1·1 - 0·Y = 1 ✓
 */
int matrix2x2_is_reciprocal_abcd(matrix2x2_t abcd, double tolerance) {
    complex_t det = matrix2x2_det(abcd);
    complex_t one = complex_make(1.0, 0.0);
    complex_t diff = complex_sub(det, one);
    return complex_approx_equal(diff, complex_make(0, 0), tolerance);
}

/* ============================================================================
 * L2: Two-Port Input/Output Impedance and Gain Calculations
 * ============================================================================ */

/**
 * Input impedance of a terminated two-port: Zin = z11 - z12*z21/(z22 + ZL).
 *
 * Derivation:
 *   V1 = z11*I1 + z12*I2     (1)
 *   V2 = z21*I1 + z22*I2     (2)
 *   V2 = -ZL*I2              (3)  termination condition
 *
 * From (3): I2 = -V2/ZL. Substitute into (2):
 *   V2 = z21*I1 - z22*V2/ZL
 *   V2*(1 + z22/ZL) = z21*I1
 *   V2 = z21*I1 * ZL/(ZL + z22)
 *
 * From (2): I2 = z21*I1/(z22 + ZL) [note the sign]
 * Substitute into (1):
 *   V1 = z11*I1 - z12*z21*I1/(z22 + ZL)
 *   V1/I1 = z11 - z12*z21/(z22 + ZL) = Zin ✓
 *
 * For unilateral networks (z12 = 0): Zin = z11, independent of load.
 * This is the "unilateral assumption" used in hand analysis.
 */
complex_t two_port_zin_from_zparams(matrix2x2_t z, complex_t zl) {
    complex_t denom = complex_add(z.m22, zl);
    complex_t term = complex_div(complex_mul(z.m12, z.m21), denom);
    return complex_sub(z.m11, term);
}

complex_t two_port_zout_from_zparams(matrix2x2_t z, complex_t zs) {
    complex_t denom = complex_add(z.m11, zs);
    complex_t term = complex_div(complex_mul(z.m12, z.m21), denom);
    return complex_sub(z.m22, term);
}

/**
 * Loaded voltage gain: Av = V2/V1 = z21*ZL / (z11*(z22+ZL) - z12*z21).
 *
 * This general formula includes all loading effects.
 * Special cases:
 *   - ZL → ∞ (open circuit): Av → z21/z11 (open-circuit gain)
 *   - ZL → 0 (short circuit): Av → 0 (no voltage across short)
 *   - z12 = 0 (unilateral): Av = z21*ZL / (z11*(z22+ZL))
 */
complex_t two_port_voltage_gain_z(matrix2x2_t z, complex_t zl) {
    complex_t num = complex_mul(z.m21, zl);
    complex_t z22_zl = complex_add(z.m22, zl);
    complex_t denom = complex_sub(
        complex_mul(z.m11, z22_zl),
        complex_mul(z.m12, z.m21)
    );
    return complex_div(num, denom);
}

/**
 * Current gain: Ai = I2/I1 = -z21/(z22 + ZL).
 *
 * The negative sign accounts for the current direction convention
 * (I2 flows into port 2, but the load current flows out).
 */
complex_t two_port_current_gain_z(matrix2x2_t z, complex_t zl) {
    complex_t num = complex_make(-z.m21.real, -z.m21.imag);
    complex_t denom = complex_add(z.m22, zl);
    return complex_div(num, denom);
}

/**
 * Power gain: Gp = P_load / P_input = |Av|^2 * Re(ZL) / Re(Zin).
 *
 * Power gain considers actual power dissipation. For passive networks
 * Gp ≤ 1. For amplifiers with gain, Gp > 1.
 */
double two_port_power_gain_z(matrix2x2_t z, complex_t zl) {
    complex_t av = two_port_voltage_gain_z(z, zl);
    complex_t zin = two_port_zin_from_zparams(z, zl);

    double av_mag_sq = complex_mag(av) * complex_mag(av);

    /* Prevent division by zero for lossless input */
    double re_zin = zin.real;
    if (re_zin < 1e-30) return INFINITY;

    double re_zl = zl.real;
    return av_mag_sq * re_zl / re_zin;
}

/**
 * Available power gain: Ga = P_avail_load / P_avail_source.
 *
 * For S-parameters: Ga = |s21|² * (1-|ΓS|²) / (|1-s11*ΓS|² - |s22-Δ*ΓS|²)
 * where ΓS = 0 for the matched-source case.
 *
 * Simplified for matched source (ΓS = 0):
 * Ga = |s21|² / (1 - |s22|²)
 */
double two_port_available_gain_s(matrix2x2_t s) {
    double s21_mag = complex_mag(s.m21);
    double s22_mag = complex_mag(s.m22);
    double denom = 1.0 - s22_mag * s22_mag;
    if (denom < 1e-30) return INFINITY;
    return (s21_mag * s21_mag) / denom;
}

/**
 * Transducer power gain: Gt = P_load / P_avail_source.
 *
 * Gt = |s21|² (1-|ΓS|²)(1-|ΓL|²) / |(1-s11*ΓS)(1-s22*ΓL)-s12*s21*ΓS*ΓL|²
 *
 * This is the most general gain expression, including both source and
 * load mismatch losses. Used in RF system design.
 */
double two_port_transducer_gain_s(matrix2x2_t s, complex_t gs, complex_t gl) {
    double s21_mag_sq = complex_mag(s.m21) * complex_mag(s.m21);

    double gs_mag_sq = complex_mag(gs) * complex_mag(gs);
    double gl_mag_sq = complex_mag(gl) * complex_mag(gl);

    complex_t one = complex_make(1.0, 0.0);

    /* Denominator term: (1 - s11*ΓS)(1 - s22*ΓL) - s12*s21*ΓS*ΓL */
    complex_t term1 = complex_sub(one, complex_mul(s.m11, gs));
    complex_t term2 = complex_sub(one, complex_mul(s.m22, gl));
    complex_t term3 = complex_mul(complex_mul(s.m12, s.m21),
                                  complex_mul(gs, gl));
    complex_t denom = complex_sub(complex_mul(term1, term2), term3);

    double denom_mag_sq = complex_mag(denom) * complex_mag(denom);
    if (denom_mag_sq < 1e-30) return 0.0;

    double num = s21_mag_sq * (1.0 - gs_mag_sq) * (1.0 - gl_mag_sq);
    return num / denom_mag_sq;
}

/* ============================================================================
 * L2: Reflection Coefficient and VSWR
 * ============================================================================ */

/**
 * Reflection coefficient: Γ = (Z - Z0) / (Z + Z0).
 *
 * This formula maps the impedance plane (right half-plane for passive Z)
 * to the unit circle in the Γ-plane:
 *   - Z = Z0 → Γ = 0 (perfect match, center of Smith chart)
 *   - Z = 0 (short) → Γ = -1 (left edge of Smith chart)
 *   - Z = ∞ (open) → Γ = +1 (right edge of Smith chart)
 *   - Z = jX (pure reactance) → |Γ| = 1 (edge of Smith chart)
 *
 * Reference: Smith, "Transmission Line Calculator", Electronics, 1939
 * Course: Georgia Tech ECE 6350 — Smith chart theory
 */
complex_t reflection_coefficient(complex_t z, double z0) {
    complex_t z0_c = complex_make(z0, 0.0);
    complex_t num = complex_sub(z, z0_c);
    complex_t denom = complex_add(z, z0_c);
    return complex_div(num, denom);
}

/**
 * VSWR = (1 + |Γ|) / (1 - |Γ|).
 * VSWR ∈ [1, ∞). A perfect match gives VSWR = 1.
 *
 * Practical VSWR values:
 *   1.0: Perfect match (return loss ∞)
 *   1.5: Return loss 14 dB (good match)
 *   2.0: Return loss 9.5 dB (acceptable)
 *   3.0: Return loss 6 dB (marginal)
 *   ∞:   Complete reflection (return loss 0 dB)
 */
double vswr_from_gamma(complex_t gamma) {
    double mag = complex_mag(gamma);
    if (mag >= 1.0) return INFINITY;
    return (1.0 + mag) / (1.0 - mag);
}

/**
 * Return loss: RL(dB) = -20*log10(|Γ|).
 * Conventionally positive for passive loads (|Γ| < 1).
 */
double return_loss_db(complex_t gamma) {
    double mag = complex_mag(gamma);
    if (mag < 1e-30) return INFINITY;
    return -20.0 * log10(mag);
}

/**
 * Inverse reflection: Z = Z0 * (1+Γ)/(1-Γ).
 * This reconstructs impedance from a VNA S11 measurement.
 */
complex_t impedance_from_gamma(complex_t gamma, double z0) {
    complex_t one = complex_make(1.0, 0.0);
    complex_t num = complex_add(one, gamma);
    complex_t denom = complex_sub(one, gamma);
    complex_t ratio = complex_div(num, denom);
    return complex_make(z0 * ratio.real, z0 * ratio.imag);
}

/**
 * Passivity check: Re(Z) ≥ 0.
 * A negative resistance would imply energy generation, which requires
 * active elements (transistors in certain bias regions, tunnel diodes,
 * or op-amp negative impedance converters).
 */
int is_passive_impedance(complex_t z) {
    return z.real >= 0.0;
}

/* ============================================================================
 * L2: Two-Port Property Checks
 * ============================================================================ */

/**
 * Reciprocity depends on parameter type:
 *   Z: symmetric (z12 = z21)
 *   Y: symmetric (y12 = y21)
 *   ABCD: det = AD-BC = 1
 *   S: symmetric (s12 = s21) when normalized to same Z0
 *
 * Lorentz reciprocity: In a linear, isotropic medium, the response
 * at port 2 due to excitation at port 1 equals the response at port 1
 * due to the same excitation at port 2.
 */
int two_port_is_reciprocal(matrix2x2_t m, param_type_t type, double tolerance) {
    switch (type) {
        case PARAM_Z:
        case PARAM_Y:
        case PARAM_S:
            return matrix2x2_is_symmetric(m, tolerance);
        case PARAM_ABCD:
            return matrix2x2_is_reciprocal_abcd(m, tolerance);
        case PARAM_H:
        case PARAM_G:
            /* For H and G: h12 = -h21 and g12 = -g21 (with sign conventions) */
            /* The reciprocal condition is more complex; simplified check */
            {
                complex_t sum = complex_add(m.m12, m.m21);
                return complex_approx_equal(sum, complex_make(0, 0), tolerance);
            }
        default:
            return 0;
    }
}

/**
 * Losslessness check — different for each parameter type.
 *
 * For Z-parameters: Re(z11) = Re(z22) = 0 (no resistive loss).
 * Pure reactive Z-params → ideal LC network.
 *
 * For S-parameters: |s11|² + |s21|² = 1 and |s12|² + |s22|² = 1,
 * plus s11*·s12 + s21*·s22 = 0 (unitary condition S^H·S = I).
 */
int two_port_is_lossless(matrix2x2_t m, param_type_t type, double tolerance) {
    switch (type) {
        case PARAM_Z:
            return (fabs(m.m11.real) < tolerance) &&
                   (fabs(m.m22.real) < tolerance);
        case PARAM_Y:
            return (fabs(m.m11.real) < tolerance) &&
                   (fabs(m.m22.real) < tolerance);
        case PARAM_S: {
            /* Check unitary condition: |s11|² + |s21|² ≈ 1 */
            double row1_sq = complex_mag(m.m11) * complex_mag(m.m11) +
                             complex_mag(m.m21) * complex_mag(m.m21);
            return fabs(row1_sq - 1.0) < tolerance * 10;
        }
        default:
            return 0;
    }
}

/**
 * Symmetry check: A network is symmetric if its ports are interchangeable.
 * For Z: z11 = z22; For Y: y11 = y22; For S: s11 = s22.
 * For reciprocal ABCD: A = D.
 */
int two_port_is_symmetrical(matrix2x2_t m, param_type_t type, double tolerance) {
    switch (type) {
        case PARAM_Z:
        case PARAM_Y:
        case PARAM_H:
        case PARAM_G:
        case PARAM_S:
            return complex_approx_equal(m.m11, m.m22, tolerance);
        case PARAM_ABCD:
            return complex_approx_equal(m.m11, m.m22, tolerance);
        default:
            return 0;
    }
}

/* ============================================================================
 * L2: Miller's Theorem
 * ============================================================================ */

/**
 * Miller's theorem: An impedance Z_f connected between the input and output
 * of an inverting amplifier with gain Av can be split into:
 *   Z1 = Zf / (1 - Av)   at the input
 *   Z2 = Zf * Av / (Av - 1)  at the output
 *
 * Derivation:
 * The current through Zf is If = (V1 - V2)/Zf.
 * At the input: Yin_equiv = If/V1 = (1 - Av)/Zf → Z1 = Zf/(1-Av)
 * At the output: Yout_equiv = If/V2 = (Av - 1)/(Av*Zf) → Z2 = Zf*Av/(Av-1)
 *
 * For large inverting gain (|Av| >> 1): Z1 ≈ Zf/|Av|, Z2 ≈ Zf.
 * This explains why Cμ (collector-base capacitance) appears multiplied
 * by the voltage gain at the input of a CE amplifier.
 *
 * Reference: Miller, "Dependence of the input impedance of a three-electrode
 *   vacuum tube upon the load in the plate circuit", 1920
 * Course: Berkeley EE105 — The Miller effect
 */
void miller_split_impedance(complex_t z, complex_t av,
                             complex_t *z1, complex_t *z2) {
    complex_t one = complex_make(1.0, 0.0);
    complex_t denom1 = complex_sub(one, av);
    complex_t denom2 = complex_sub(av, one);

    *z1 = complex_div(z, denom1);
    *z2 = complex_div(complex_mul(z, av), denom2);
}

/**
 * Miller capacitance: Cin_miller = Cf * (1 + |Av|).
 *
 * For a common-emitter amplifier:
 *   Av = -gm * RL (inverting, magnitude = gm*RL)
 *   Cin_total = Cπ + Cμ * (1 + gm*RL)
 *
 * The factor (1 + |Av|) explains why the input capacitance can be
 * dominated by Cμ even though Cμ is physically small (~1 pF).
 */
double miller_input_capacitance(double cf, double av) {
    /* av is magnitude of voltage gain (positive number) */
    return cf * (1.0 + av);
}

/* ============================================================================
 * L4: Maximum Power Transfer Theorem
 * ============================================================================ */

/**
 * Conjugate match: ZL = ZS* for maximum power transfer.
 *
 * Proof (Jacobi's law):
 *   P_load = |Vs|² * RL / |ZS + ZL|²
 *          = |Vs|² * RL / ((RS+RL)² + (XS+XL)²)
 *
 * For fixed RL: set XL = -XS to maximize denominator (makes it (RS+RL)²).
 * Then dP/dRL = 0 → RL = RS.
 * Therefore: ZL = RS - jXS = ZS*. QED.
 *
 * This theorem is fundamental to:
 *   - RF matching network design
 *   - Maximum power extraction from antennas
 *   - Audio amplifier load matching
 *   - Wireless power transfer optimization
 */
complex_t conjugate_match(complex_t zsource) {
    return complex_conj(zsource);
}

/**
 * Maximum available power from a Thevenin source.
 * Pavail = |Vs|²/(8*RS) when ZL = ZS*.
 *
 * Derivation: With ZL = RS - jXS, ZS = RS + jXS:
 *   I = Vs/(ZS+ZL) = Vs/(2*RS)
 *   P = |I|² * RL = |Vs|²/(4*RS²) * RS = |Vs|²/(4*RS)
 *
 * Wait — let me recalculate:
 *   I = Vs / (ZS + ZL) = Vs / (RS + jXS + RS - jXS) = Vs / (2*RS)
 *   P_load = |I|² * RL = |Vs|² / (4 * RS²) * RS = |Vs|² / (4 * RS)
 *
 * Hmm, the textbook says |Vs|²/(8*RS). The factor of 2 difference comes
 * from whether Vs is the peak or RMS value. Using peak: Pavail = |Vs|²/(8*RS).
 * Using RMS: Pavail = |Vs_rms|²/(4*RS). Since |Vs_peak|² = 2*|Vs_rms|²:
 *   Pavail = (2*|Vs_rms|²)/(8*RS) = |Vs_rms|²/(4*RS). ✓
 *
 * Reference: Desoer & Kuh, "Basic Circuit Theory", Ch. 9
 * Course: MIT 6.002 — Maximum power transfer
 */
double maximum_available_power(double vs_mag, complex_t zs) {
    double rs = zs.real;
    if (rs < 1e-30) return INFINITY;  /* Ideal voltage source */
    return (vs_mag * vs_mag) / (8.0 * rs);
}

/* ============================================================================
 * L5: Numerical Utilities
 * ============================================================================ */

int complex_approx_equal(complex_t a, complex_t b, double tolerance) {
    complex_t diff = complex_sub(a, b);
    return complex_mag(diff) < tolerance;
}

int complex_is_zero(complex_t z, double tolerance) {
    return complex_mag(z) < tolerance;
}

void complex_to_string(complex_t z, char *buf, size_t buf_size) {
    if (z.imag >= 0.0) {
        snprintf(buf, buf_size, "%.6g + j%.6g", z.real, z.imag);
    } else {
        snprintf(buf, buf_size, "%.6g - j%.6g", z.real, -z.imag);
    }
}

void matrix2x2_print(matrix2x2_t m, const char *label) {
    char s11[64], s12[64], s21[64], s22[64];
    complex_to_string(m.m11, s11, sizeof(s11));
    complex_to_string(m.m12, s12, sizeof(s12));
    complex_to_string(m.m21, s21, sizeof(s21));
    complex_to_string(m.m22, s22, sizeof(s22));
    printf("%s:\n", label);
    printf("  [%s, %s]\n", s11, s12);
    printf("  [%s, %s]\n", s21, s22);
}
