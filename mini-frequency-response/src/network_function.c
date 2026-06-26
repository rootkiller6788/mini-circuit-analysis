/**
 * @file network_function.c
 * @brief Network functions in the frequency domain implementation
 *
 * Implements:
 * - Driving-point impedance computation
 * - Two-port parameter conversion (Z↔Y↔ABCD↔S)
 * - S-parameter computation from Z-parameters
 * - Rollett stability factor K
 * - Maximum available gain (MAG/MSG) computation
 * - Positive real function check
 * - Input impedance of terminated two-port
 * - VSWR computation
 *
 * Reference: Pozar (2012), Van Valkenburg (1974), Guillemin (1957)
 * Course: Berkeley EE16B, ETH 227-0455, Georgia Tech ECE 6350
 */

#include "network_function.h"
#include "frequency_response.h"
#include "transfer_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <complex.h>

/* ============================================================================
 * Network Function Synthesis
 * ============================================================================ */

/**
 * Impedance of basic RLC configurations.
 *
 * Series:  Z(s) = R + sL + 1/(sC) = (s²LC + sRC + 1)/(sC)
 * Parallel: Z(s) = 1/(1/R + 1/(sL) + sC) = s/(s²C + s/R + 1/L)
 *                          = sL·R/(s²LRC + sL + R)
 *                          = s·something / (s²·something + s·something)
 *
 * Parallel RLC impedance in rational form:
 *   Z(s) = 1 / (1/R + 1/(sL) + sC)
 *        = sLR / (s²LRC + sL + R)
 */
tf_polynomial_t *network_impedance_rlc(double R, double L, double C,
                                         int is_par)
{
    if (is_par) {
        /* Z(s) = sLR/(s²LRC + sL + R) */
        double num[2] = { 0.0, L * R };
        double den[3] = { R, L, L * R * C };

        /* Normalize */
        return tf_polynomial_create(num, 1, den, 2);
    } else {
        /* Z(s) = R + sL + 1/(sC) = (s²LC + sRC + 1)/(sC) */
        double num[3] = { 1.0, R * C, L * C };
        double den[2] = { 0.0, C };

        return tf_polynomial_create(num, 2, den, 1);
    }
}

/**
 * Voltage divider transfer function: H(s) = Z₂/(Z₁+Z₂).
 *
 * This is the single most common transfer function in circuit analysis.
 * Examples:
 * - RC lowpass: Z₁=R, Z₂=1/(sC) → H(s) = 1/(1+sRC)
 * - CR highpass: Z₁=1/(sC), Z₂=R → H(s) = sRC/(1+sRC)
 * - Compensated attenuator (scope probe): Z₁=R₁||C₁, Z₂=R₂||C₂
 *   H(s) = (R₂/(R₁+R₂))·(1+sR₁C₁)/(1+s(R₁||R₂)(C₁+C₂))
 *   When R₁C₁ = R₂C₂: H(s) = R₂/(R₁+R₂) (frequency-independent!)
 */
tf_polynomial_t *network_voltage_divider_tf(const tf_polynomial_t *Z1,
                                              const tf_polynomial_t *Z2)
{
    if (!Z1 || !Z2) return NULL;

    /* H = Z₂/(Z₁+Z₂) = N₂/D₂ / (N₁/D₁ + N₂/D₂)
     *                = N₂D₁ / (N₁D₂ + N₂D₁) */

    /* Numerator: N₂ × D₁ */
    double *num = (double *)malloc((Z2->num_order + Z1->den_order + 1)
                                    * sizeof(double));
    if (!num) return NULL;
    /* poly_multiply returns allocated buffer */
    double *n2d1 = NULL;
    size_t n2d1_len;
    n2d1 = (double *)malloc((Z2->num_order + Z1->den_order + 1)
                             * sizeof(double));
    if (!n2d1) { free(num); return NULL; }
    memset(n2d1, 0, (Z2->num_order + Z1->den_order + 1) * sizeof(double));
    for (size_t i = 0; i <= Z2->num_order; i++) {
        for (size_t j = 0; j <= Z1->den_order; j++) {
            n2d1[i + j] += Z2->num[i] * Z1->den[j];
        }
    }
    n2d1_len = Z2->num_order + Z1->den_order;

    /* Denominator: N₁D₂ + N₂D₁ */
    size_t den_len1;
    double *n1d2 = (double *)calloc((Z1->num_order + Z2->den_order + 1),
                                     sizeof(double));
    if (!n1d2) { free(n2d1); free(num); return NULL; }
    for (size_t i = 0; i <= Z1->num_order; i++) {
        for (size_t j = 0; j <= Z2->den_order; j++) {
            n1d2[i + j] += Z1->num[i] * Z2->den[j];
        }
    }
    den_len1 = Z1->num_order + Z2->den_order;

    /* Sum: den = n1d2 + n2d1 */
    size_t den_len = (den_len1 > n2d1_len) ? den_len1 : n2d1_len;
    double *den = (double *)calloc(den_len + 1, sizeof(double));
    if (!den) { free(n1d2); free(n2d1); free(num); return NULL; }

    for (size_t i = 0; i <= den_len1; i++) den[i] += n1d2[i];
    for (size_t i = 0; i <= n2d1_len; i++) den[i] += n2d1[i];

    free(n1d2);

    tf_polynomial_t *result = tf_polynomial_create(n2d1, n2d1_len,
                                                     den, den_len);
    free(n2d1);
    free(den);
    free(num);
    return result;
}

/**
 * Element impedances for common single-port configurations.
 *
 * element_type mapping:
 * 0: L     → Z(s) = sL
 * 1: C     → Z(s) = 1/(sC)
 * 2: R     → Z(s) = R
 * 3: R+L   → Z(s) = R + sL
 * 4: R+C   → Z(s) = R + 1/(sC)
 * 5: R||L  → Z(s) = sRL/(R + sL)
 * 6: R||C  → Z(s) = R/(1 + sRC)
 */
tf_polynomial_t *network_element_impedance(int element_type,
                                             double val1, double val2)
{
    double num[3], den[3];
    int num_deg = 0, den_deg = 0;

    switch (element_type) {
        case 0:  /* L: Z(s) = sL */
            num[0] = 0.0; num[1] = val1; num_deg = 1;
            den[0] = 1.0; den_deg = 0;
            break;

        case 1:  /* C: Z(s) = 1/(sC) */
            num[0] = 1.0; num_deg = 0;
            den[0] = 0.0; den[1] = val1; den_deg = 1;
            break;

        case 2:  /* R: Z(s) = R */
            num[0] = val1; num_deg = 0;
            den[0] = 1.0; den_deg = 0;
            break;

        case 3:  /* R+L: Z(s) = R + sL */
            num[0] = val1; num[1] = val2; num_deg = 1;
            den[0] = 1.0; den_deg = 0;
            break;

        case 4:  /* R+C: Z(s) = R + 1/(sC) = (sRC + 1)/(sC) */
            num[0] = 1.0; num[1] = val1 * val2; num_deg = 1;
            den[0] = 0.0; den[1] = val2; den_deg = 1;
            break;

        case 5:  /* R||L: Z(s) = (R·sL)/(R + sL) */
            num[0] = 0.0; num[1] = val1 * val2; num_deg = 1;
            den[0] = val1; den[1] = val2; den_deg = 1;
            break;

        case 6:  /* R||C: Z(s) = R/(1 + sRC) */
            num[0] = val1; num_deg = 0;
            den[0] = 1.0; den[1] = val1 * val2; den_deg = 1;
            break;

        default:
            return NULL;
    }

    return tf_polynomial_create(num, num_deg, den, den_deg);
}

/* ============================================================================
 * Two-Port Parameter Conversion
 * ============================================================================ */

/**
 * Convert between two-port parameter sets.
 *
 * The six standard parameter sets are interrelated by matrix
 * transformations. This function implements the most common
 * conversions.
 *
 * Z → Y: Y = Z^{-1}
 * Z → ABCD: A=z₁₁/z₂₁, B=Δz/z₂₁, C=1/z₂₁, D=z₂₂/z₂₁
 * Z → S: S = (Z+Z₀I)^{-1}(Z-Z₀I)
 */
int network_param_convert(const double _Complex src[4],
                           network_param_type_t src_type,
                           double _Complex dst[4],
                           network_param_type_t dst_type,
                           double Z0)
{
    if (!src || !dst) return -1;

    /* For simplicity, implement Z ↔ S and Z ↔ Y conversions */
    double _Complex z11, z12, z21, z22;

    /* First, convert source to Z-parameters */
    switch (src_type) {
        case NETWORK_PARAM_Z:
            z11 = src[0]; z12 = src[1]; z21 = src[2]; z22 = src[3];
            break;

        case NETWORK_PARAM_Y:
            /* Z = Y^{-1} */
            {
                double _Complex detY = src[0] * src[3] - src[1] * src[2];
                if (cabs(detY) < 1e-15) return -1;
                z11 = src[3] / detY;
                z12 = -src[1] / detY;
                z21 = -src[2] / detY;
                z22 = src[0] / detY;
            }
            break;

        case NETWORK_PARAM_S:
            /* Z = Z₀·(I+S)(I-S)^{-1} */
            {
                double _Complex det = (1.0 - src[0]) * (1.0 - src[3])
                                       - src[1] * src[2];
                if (cabs(det) < 1e-15) return -1;
                z11 = Z0 * ((1.0 + src[0]) * (1.0 - src[3])
                             + src[1] * src[2]) / det;
                z12 = Z0 * 2.0 * src[1] / det;
                z21 = Z0 * 2.0 * src[2] / det;
                z22 = Z0 * ((1.0 - src[0]) * (1.0 + src[3])
                             + src[1] * src[2]) / det;
            }
            break;

        default:
            return -1;  /* Other conversions not implemented yet */
    }

    /* Convert Z-parameters to destination type */
    switch (dst_type) {
        case NETWORK_PARAM_Z:
            dst[0] = z11; dst[1] = z12; dst[2] = z21; dst[3] = z22;
            break;

        case NETWORK_PARAM_Y:
            {
                double _Complex detZ = z11 * z22 - z12 * z21;
                if (cabs(detZ) < 1e-15) return -1;
                dst[0] = z22 / detZ;
                dst[1] = -z12 / detZ;
                dst[2] = -z21 / detZ;
                dst[3] = z11 / detZ;
            }
            break;

        case NETWORK_PARAM_ABCD:
            if (cabs(z21) < 1e-15) return -1;
            dst[0] = z11 / z21;
            dst[1] = (z11 * z22 - z12 * z21) / z21;
            dst[2] = 1.0 / z21;
            dst[3] = z22 / z21;
            break;

        case NETWORK_PARAM_S:
            /* S = (Z - Z₀I)(Z + Z₀I)^{-1} */
            {
                double _Complex denom = (z11 + Z0) * (z22 + Z0) - z12 * z21;
                if (cabs(denom) < 1e-15) return -1;
                dst[0] = ((z11 - Z0) * (z22 + Z0) - z12 * z21) / denom;
                dst[1] = 2.0 * z12 * Z0 / denom;
                dst[2] = 2.0 * z21 * Z0 / denom;
                dst[3] = ((z11 + Z0) * (z22 - Z0) - z12 * z21) / denom;
            }
            break;

        default:
            return -1;
    }

    return 0;
}

/* ============================================================================
 * S-Parameter Computation
 * ============================================================================ */

s_params_2port_t network_s_params_from_z(double _Complex z11,
                                           double _Complex z12,
                                           double _Complex z21,
                                           double _Complex z22,
                                           double Z0, double frequency_hz)
{
    s_params_2port_t s;
    memset(&s, 0, sizeof(s_params_2port_t));
    s.frequency_hz = frequency_hz;

    double _Complex det = (z11 + Z0) * (z22 + Z0) - z12 * z21;
    if (cabs(det) < 1e-15) return s;

    s.s11 = ((z11 - Z0) * (z22 + Z0) - z12 * z21) / det;
    s.s12 = 2.0 * z12 * Z0 / det;
    s.s21 = 2.0 * z21 * Z0 / det;
    s.s22 = ((z11 + Z0) * (z22 - Z0) - z12 * z21) / det;

    return s;
}

/**
 * Rollett stability factor K.
 *
 * K = (1 - |S₁₁|² - |S₂₂|² + |Δ|²) / (2·|S₁₂·S₂₁|)
 * Δ = S₁₁·S₂₂ - S₁₂·S₂₁
 *
 * Unconditional stability requires:
 *   K > 1  AND  |Δ| < 1
 *
 * This is the most widely used stability criterion for RF amplifiers.
 *
 * For potentially unstable devices (K < 1), the amplifier can still
 * be made stable by choosing appropriate source and load impedances
 * that lie in the stable regions of the Smith chart.
 */
double network_rollett_k(const s_params_2port_t *s)
{
    double s11_mag2 = cabs(s->s11) * cabs(s->s11);
    double s22_mag2 = cabs(s->s22) * cabs(s->s22);
    double s12_s21_mag = cabs(s->s12) * cabs(s->s21);

    double _Complex delta = s->s11 * s->s22 - s->s12 * s->s21;
    double delta_mag2 = cabs(delta) * cabs(delta);

    if (s12_s21_mag < 1e-15) return INFINITY;

    return (1.0 - s11_mag2 - s22_mag2 + delta_mag2) / (2.0 * s12_s21_mag);
}

/**
 * Maximum Available Gain (MAG) / Maximum Stable Gain (MSG).
 *
 * For unconditionally stable (K > 1):
 *   MAG = |S₂₁/S₁₂| · (K - √(K² - 1))
 *
 * For potentially unstable (K < 1):
 *   MSG = |S₂₁|/|S₁₂|
 *
 * MAG is the maximum gain achievable with simultaneous conjugate
 * match at both input and output ports. MSG is the gain achievable
 * by resistively loading the device to make it stable (K=1).
 */
double network_max_gain(const s_params_2port_t *s)
{
    double K = network_rollett_k(s);
    double s21_s12_ratio = cabs(s->s21) / (cabs(s->s12) + 1e-15);

    if (K > 1.0) {
        /* Unconditionally stable → MAG */
        return s21_s12_ratio * (K - sqrt(K * K - 1.0));
    } else {
        /* Potentially unstable → MSG */
        return s21_s12_ratio;
    }
}

/* ============================================================================
 * Positive Real Function Check
 * ============================================================================ */

/**
 * Verify that an impedance function Z(s) is positive real.
 *
 * Checks:
 * 1. Real coefficients (implied by our representation)
 * 2. Re{Z(jω)} ≥ 0 for all sampled frequencies
 * 3. Pole-zero alternation on jω-axis (not fully checked here)
 */
int network_is_positive_real(const tf_polynomial_t *Z,
                               const double *freq, size_t n)
{
    if (!Z || !freq || n == 0) return -1;

    for (size_t i = 0; i < n; i++) {
        double w = 2.0 * M_PI * freq[i];
        double _Complex val = tf_evaluate_freq(Z, w);
        if (creal(val) < -1e-12) {
            return 0;  /* Negative real part → not PR */
        }
    }

    /* Check degree condition: |deg(num) - deg(den)| ≤ 1 */
    size_t num_deg = Z->num_order;
    size_t den_deg = Z->den_order;
    int deg_diff = (int)num_deg - (int)den_deg;
    if (abs(deg_diff) > 1) {
        return 0;
    }

    return 1;
}

/* ============================================================================
 * Input Impedance of Terminated Two-Port
 * ============================================================================ */

/**
 * Input impedance of a two-port terminated in load Z_L.
 *
 * From S-parameters: Γ_in = S₁₁ + S₁₂·S₂₁·Γ_L/(1 - S₂₂·Γ_L)
 * where Γ_L = (Z_L - Z₀)/(Z_L + Z₀)
 *
 * Then: Z_in = Z₀·(1 + Γ_in)/(1 - Γ_in)
 */
double _Complex network_input_impedance(const s_params_2port_t *s_params,
                                           double _Complex ZL, double Z0)
{
    double _Complex gamma_L = (ZL - Z0) / (ZL + Z0);

    double _Complex gamma_in = s_params->s11
        + s_params->s12 * s_params->s21 * gamma_L
          / (1.0 - s_params->s22 * gamma_L);

    double _Complex Z_in = Z0 * (1.0 + gamma_in) / (1.0 - gamma_in);
    return Z_in;
}

/**
 * Voltage Standing Wave Ratio.
 *
 * VSWR = (1 + |Γ|)/(1 - |Γ|)
 *
 * |Γ| = 0: VSWR = 1 (perfect match)
 * |Γ| = 1/3: VSWR = 2
 * |Γ| = 1/2: VSWR = 3
 * |Γ| → 1: VSWR → ∞ (complete mismatch)
 */
double network_vswr(double _Complex gamma)
{
    double mag = cabs(gamma);
    if (mag >= 1.0) return INFINITY;
    return (1.0 + mag) / (1.0 - mag);
}
