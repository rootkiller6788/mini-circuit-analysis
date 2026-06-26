/**
 * @file conversion.c
 * @brief Two-Port Parameter Conversion — All 30 Directional Conversions
 *
 * Implements every mathematically meaningful conversion between the six
 * two-port parameter types. Each conversion is based on matrix algebra
 * formulas from Frickey (IEEE MTT, 1994).
 */

#include "../include/conversion.h"
#include <string.h>

/* NAN sentinel for failed conversions */
static complex_t nan_c(void) { return complex_make(NAN, NAN); }
static matrix2x2_t nan_matrix(void) {
    return matrix2x2_make(nan_c(), nan_c(), nan_c(), nan_c());
}

/* ============================================================================
 * Z ↔ Y (inversion)
 * ============================================================================ */

matrix2x2_t convert_z_to_y(matrix2x2_t z) {
    return matrix2x2_inv(z);
}

matrix2x2_t convert_y_to_z(matrix2x2_t y) {
    return matrix2x2_inv(y);
}

/* ============================================================================
 * Z ↔ H
 * ============================================================================ */

matrix2x2_t convert_z_to_h(matrix2x2_t z) {
    if (complex_is_zero(z.m22, 1e-30)) return nan_matrix();

    complex_t inv_z22 = complex_div(complex_make(1.0, 0.0), z.m22);
    complex_t det = matrix2x2_det(z);

    return matrix2x2_make(
        complex_mul(det, inv_z22),                /* h11 = ΔZ/z22 */
        complex_mul(z.m12, inv_z22),              /* h12 = z12/z22 */
        complex_mul(complex_make(-z.m21.real, -z.m21.imag), inv_z22), /* h21 = -z21/z22 */
        inv_z22                                    /* h22 = 1/z22 */
    );
}

matrix2x2_t convert_h_to_z(matrix2x2_t h) {
    if (complex_is_zero(h.m22, 1e-30)) return nan_matrix();

    complex_t inv_h22 = complex_div(complex_make(1.0, 0.0), h.m22);
    complex_t det = matrix2x2_det(h);

    return matrix2x2_make(
        complex_mul(det, inv_h22),                /* z11 = ΔH/h22 */
        complex_mul(h.m12, inv_h22),              /* z12 = h12/h22 */
        complex_mul(complex_make(-h.m21.real, -h.m21.imag), inv_h22), /* z21 = -h21/h22 */
        inv_h22                                    /* z22 = 1/h22 */
    );
}

/* ============================================================================
 * Z ↔ G
 * ============================================================================ */

matrix2x2_t convert_z_to_g(matrix2x2_t z) {
    if (complex_is_zero(z.m11, 1e-30)) return nan_matrix();

    complex_t inv_z11 = complex_div(complex_make(1.0, 0.0), z.m11);
    complex_t det = matrix2x2_det(z);

    return matrix2x2_make(
        inv_z11,                                   /* g11 = 1/z11 */
        complex_mul(complex_make(-z.m12.real, -z.m12.imag), inv_z11), /* g12 = -z12/z11 */
        complex_mul(z.m21, inv_z11),               /* g21 = z21/z11 */
        complex_mul(det, inv_z11)                  /* g22 = ΔZ/z11 */
    );
}

matrix2x2_t convert_g_to_z(matrix2x2_t g) {
    if (complex_is_zero(g.m11, 1e-30)) return nan_matrix();

    complex_t inv_g11 = complex_div(complex_make(1.0, 0.0), g.m11);
    complex_t det = matrix2x2_det(g);

    return matrix2x2_make(
        inv_g11,                                   /* z11 = 1/g11 */
        complex_mul(complex_make(-g.m12.real, -g.m12.imag), inv_g11), /* z12 = -g12/g11 */
        complex_mul(g.m21, inv_g11),               /* z21 = g21/g11 */
        complex_mul(det, inv_g11)                  /* z22 = ΔG/g11 */
    );
}

/* ============================================================================
 * Z ↔ ABCD
 * ============================================================================ */

matrix2x2_t convert_z_to_abcd(matrix2x2_t z) {
    if (complex_is_zero(z.m21, 1e-30)) return nan_matrix();

    complex_t inv_z21 = complex_div(complex_make(1.0, 0.0), z.m21);
    complex_t det = matrix2x2_det(z);

    return matrix2x2_make(
        complex_mul(z.m11, inv_z21),               /* A = z11/z21 */
        complex_mul(det, inv_z21),                 /* B = ΔZ/z21 */
        inv_z21,                                    /* C = 1/z21 */
        complex_mul(z.m22, inv_z21)                /* D = z22/z21 */
    );
}

matrix2x2_t convert_abcd_to_z(matrix2x2_t abcd) {
    if (complex_is_zero(abcd.m21, 1e-30)) return nan_matrix();

    complex_t inv_c = complex_div(complex_make(1.0, 0.0), abcd.m21);
    complex_t det = matrix2x2_det(abcd);

    return matrix2x2_make(
        complex_mul(abcd.m11, inv_c),              /* z11 = A/C */
        complex_mul(det, inv_c),                   /* z12 = ΔT/C */
        inv_c,                                      /* z21 = 1/C */
        complex_mul(abcd.m22, inv_c)               /* z22 = D/C */
    );
}

/* ============================================================================
 * Y ↔ H
 * ============================================================================ */

matrix2x2_t convert_y_to_h(matrix2x2_t y) {
    if (complex_is_zero(y.m11, 1e-30)) return nan_matrix();

    complex_t inv_y11 = complex_div(complex_make(1.0, 0.0), y.m11);
    complex_t det = matrix2x2_det(y);

    return matrix2x2_make(
        inv_y11,                                   /* h11 = 1/y11 */
        complex_mul(complex_make(-y.m12.real, -y.m12.imag), inv_y11), /* h12 = -y12/y11 */
        complex_mul(y.m21, inv_y11),               /* h21 = y21/y11 */
        complex_mul(det, inv_y11)                  /* h22 = ΔY/y11 */
    );
}

matrix2x2_t convert_h_to_y(matrix2x2_t h) {
    if (complex_is_zero(h.m11, 1e-30)) return nan_matrix();

    complex_t inv_h11 = complex_div(complex_make(1.0, 0.0), h.m11);
    complex_t det = matrix2x2_det(h);

    return matrix2x2_make(
        inv_h11,                                   /* y11 = 1/h11 */
        complex_mul(complex_make(-h.m12.real, -h.m12.imag), inv_h11), /* y12 = -h12/h11 */
        complex_mul(h.m21, inv_h11),               /* y21 = h21/h11 */
        complex_mul(det, inv_h11)                  /* y22 = ΔH/h11 */
    );
}

/* ============================================================================
 * Y ↔ ABCD
 * ============================================================================ */

/**
 * Y→ABCD conversion formula:
 *   A = -y22/y21, B = -1/y21, C = -ΔY/y21, D = -y11/y21
 */
static matrix2x2_t convert_y_to_abcd_clean(matrix2x2_t y) {
    if (complex_is_zero(y.m21, 1e-30)) return nan_matrix();

    complex_t neg_inv = complex_div(complex_make(-1.0, 0.0), y.m21);
    complex_t det_y = matrix2x2_det(y);

    return matrix2x2_make(
        complex_mul(y.m22, neg_inv),              /* A = -y22/y21 */
        neg_inv,                                   /* B = -1/y21 */
        complex_mul(det_y, neg_inv),               /* C = -ΔY/y21 */
        complex_mul(y.m11, neg_inv)               /* D = -y11/y21 */
    );
}

matrix2x2_t convert_y_to_abcd(matrix2x2_t y) {
    return convert_y_to_abcd_clean(y);
}

matrix2x2_t convert_abcd_to_y(matrix2x2_t abcd) {
    if (complex_is_zero(abcd.m12, 1e-30)) return nan_matrix();

    complex_t inv_b = complex_div(complex_make(1.0, 0.0), abcd.m12);
    complex_t det = matrix2x2_det(abcd);
    complex_t neg_det = complex_make(-det.real, -det.imag);

    /* y11 = D/B, y12 = -ΔT/B, y21 = -1/B, y22 = A/B */
    return matrix2x2_make(
        complex_mul(abcd.m22, inv_b),
        complex_mul(neg_det, inv_b),
        complex_mul(complex_make(-1.0, 0.0), inv_b),
        complex_mul(abcd.m11, inv_b)
    );
}

/* ============================================================================
 * H ↔ G (inversion)
 * ============================================================================ */

matrix2x2_t convert_h_to_g(matrix2x2_t h) {
    return matrix2x2_inv(h);
}

matrix2x2_t convert_g_to_h(matrix2x2_t g) {
    return matrix2x2_inv(g);
}

/* ============================================================================
 * S ↔ Z
 * ============================================================================ */

/**
 * S → Z: Z = Z0 * (I+S)/(I-S).
 *
 * ΔS = (1-s11)(1-s22) - s12*s21
 * z11 = Z0 * ((1+s11)(1-s22) + s12*s21) / ΔS
 * z12 = Z0 * 2*s12 / ΔS
 * z21 = Z0 * 2*s21 / ΔS
 * z22 = Z0 * ((1-s11)(1+s22) + s12*s21) / ΔS
 */
matrix2x2_t convert_s_to_z(matrix2x2_t s, double z0) {
    complex_t one = complex_make(1.0, 0.0);
    complex_t s12s21 = complex_mul(s.m12, s.m21);

    /* ΔS = (1-s11)(1-s22) - s12*s21 */
    complex_t denom = complex_sub(
        complex_mul(complex_sub(one, s.m11), complex_sub(one, s.m22)),
        s12s21
    );

    if (complex_is_zero(denom, 1e-30)) return nan_matrix();

    complex_t inv_denom = complex_div(complex_make(1.0, 0.0), denom);

    complex_t z11 = complex_mul(
        complex_make(z0, 0.0),
        complex_mul(
            complex_add(
                complex_mul(complex_add(one, s.m11), complex_sub(one, s.m22)),
                s12s21
            ),
            inv_denom
        )
    );
    complex_t z12 = complex_mul(
        complex_make(2.0 * z0, 0.0),
        complex_mul(s.m12, inv_denom)
    );
    complex_t z21 = complex_mul(
        complex_make(2.0 * z0, 0.0),
        complex_mul(s.m21, inv_denom)
    );
    complex_t z22 = complex_mul(
        complex_make(z0, 0.0),
        complex_mul(
            complex_add(
                complex_mul(complex_sub(one, s.m11), complex_add(one, s.m22)),
                s12s21
            ),
            inv_denom
        )
    );

    return matrix2x2_make(z11, z12, z21, z22);
}

/**
 * Z → S: S = (Z - Z0*I)/(Z + Z0*I).
 */
matrix2x2_t convert_z_to_s(matrix2x2_t z, double z0) {
    complex_t z0_c = complex_make(z0, 0.0);
    complex_t z11_p = complex_add(z.m11, z0_c);
    complex_t z22_p = complex_add(z.m22, z0_c);
    complex_t z11_m = complex_sub(z.m11, z0_c);
    complex_t z22_m = complex_sub(z.m22, z0_c);

    complex_t denom = complex_sub(
        complex_mul(z11_p, z22_p),
        complex_mul(z.m12, z.m21)
    );

    if (complex_is_zero(denom, 1e-30)) return nan_matrix();

    complex_t inv_denom = complex_div(complex_make(1.0, 0.0), denom);

    complex_t s11 = complex_mul(
        complex_sub(
            complex_mul(z11_m, z22_p),
            complex_mul(z.m12, z.m21)
        ),
        inv_denom
    );
    complex_t s12 = complex_mul(
        complex_make(2.0 * z0, 0.0),
        complex_mul(z.m12, inv_denom)
    );
    complex_t s21 = complex_mul(
        complex_make(2.0 * z0, 0.0),
        complex_mul(z.m21, inv_denom)
    );
    complex_t s22 = complex_mul(
        complex_sub(
            complex_mul(z11_p, z22_m),
            complex_mul(z.m12, z.m21)
        ),
        inv_denom
    );

    return matrix2x2_make(s11, s12, s21, s22);
}

/* ============================================================================
 * S ↔ Y
 * ============================================================================ */

matrix2x2_t convert_s_to_y(matrix2x2_t s, double z0) {
    matrix2x2_t z = convert_s_to_z(s, z0);
    return convert_z_to_y(z);
}

matrix2x2_t convert_y_to_s(matrix2x2_t y, double z0) {
    matrix2x2_t z = convert_y_to_z(y);
    return convert_z_to_s(z, z0);
}

/* ============================================================================
 * S ↔ ABCD
 * ============================================================================ */

/**
 * S → ABCD:
 * A = ((1+s11)(1-s22) + s12*s21) / (2*s21)
 * B = Z0 * ((1+s11)(1+s22) - s12*s21) / (2*s21)
 * C = (1/Z0) * ((1-s11)(1-s22) - s12*s21) / (2*s21)
 * D = ((1-s11)(1+s22) + s12*s21) / (2*s21)
 */
matrix2x2_t convert_s_to_abcd(matrix2x2_t s, double z0) {
    if (complex_is_zero(s.m21, 1e-30)) return nan_matrix();

    complex_t one = complex_make(1.0, 0.0);
    complex_t two_s21 = complex_make(2.0 * s.m21.real, 2.0 * s.m21.imag);
    complex_t inv_2s21 = complex_div(complex_make(1.0, 0.0), two_s21);
    complex_t s12s21 = complex_mul(s.m12, s.m21);

    complex_t a = complex_mul(
        complex_add(
            complex_mul(complex_add(one, s.m11), complex_sub(one, s.m22)),
            s12s21
        ),
        inv_2s21
    );
    complex_t b = complex_mul(
        complex_make(z0, 0.0),
        complex_mul(
            complex_sub(
                complex_mul(complex_add(one, s.m11), complex_add(one, s.m22)),
                s12s21
            ),
            inv_2s21
        )
    );
    complex_t c_val = complex_mul(
        complex_make(1.0 / z0, 0.0),
        complex_mul(
            complex_sub(
                complex_mul(complex_sub(one, s.m11), complex_sub(one, s.m22)),
                s12s21
            ),
            inv_2s21
        )
    );
    complex_t d = complex_mul(
        complex_add(
            complex_mul(complex_sub(one, s.m11), complex_add(one, s.m22)),
            s12s21
        ),
        inv_2s21
    );

    return matrix2x2_make(a, b, c_val, d);
}

/**
 * ABCD → S:
 * denom = A + B/Z0 + C*Z0 + D
 * s11 = (A + B/Z0 - C*Z0 - D) / denom
 * s12 = 2*(AD-BC) / denom
 * s21 = 2 / denom
 * s22 = (-A + B/Z0 - C*Z0 + D) / denom
 */
matrix2x2_t convert_abcd_to_s(matrix2x2_t abcd, double z0) {
    complex_t b_div_z0 = complex_make(abcd.m12.real / z0, abcd.m12.imag / z0);
    complex_t c_mul_z0 = complex_make(abcd.m21.real * z0, abcd.m21.imag * z0);

    complex_t denom = complex_add(
        complex_add(abcd.m11, b_div_z0),
        complex_add(c_mul_z0, abcd.m22)
    );

    if (complex_is_zero(denom, 1e-30)) return nan_matrix();

    complex_t inv_denom = complex_div(complex_make(1.0, 0.0), denom);

    complex_t s11 = complex_mul(
        complex_sub(
            complex_add(abcd.m11, b_div_z0),
            complex_add(c_mul_z0, abcd.m22)
        ),
        inv_denom
    );
    complex_t det = matrix2x2_det(abcd);
    complex_t s12 = complex_mul(
        complex_make(2.0 * det.real, 2.0 * det.imag),
        inv_denom
    );
    complex_t s21 = complex_mul(
        complex_make(2.0, 0.0),
        inv_denom
    );
    complex_t s22 = complex_mul(
        complex_add(
            complex_sub(b_div_z0, abcd.m11),
            complex_sub(abcd.m22, c_mul_z0)
        ),
        inv_denom
    );

    return matrix2x2_make(s11, s12, s21, s22);
}

/* ============================================================================
 * ABCD ↔ H
 * ============================================================================ */

matrix2x2_t convert_abcd_to_h(matrix2x2_t abcd) {
    if (complex_is_zero(abcd.m22, 1e-30)) return nan_matrix();

    complex_t inv_d = complex_div(complex_make(1.0, 0.0), abcd.m22);
    complex_t det = matrix2x2_det(abcd);

    return matrix2x2_make(
        complex_mul(abcd.m12, inv_d),             /* h11 = B/D */
        complex_mul(det, inv_d),                   /* h12 = ΔT/D */
        complex_mul(complex_make(-1.0, 0.0), inv_d), /* h21 = -1/D */
        complex_mul(abcd.m21, inv_d)               /* h22 = C/D */
    );
}

matrix2x2_t convert_h_to_abcd(matrix2x2_t h) {
    if (complex_is_zero(h.m21, 1e-30)) return nan_matrix();

    complex_t neg_inv_h21 = complex_div(complex_make(-1.0, 0.0), h.m21);
    complex_t det = matrix2x2_det(h);

    return matrix2x2_make(
        complex_mul(det, neg_inv_h21),             /* A = -ΔH/h21 */
        complex_mul(h.m11, neg_inv_h21),           /* B = -h11/h21 */
        complex_mul(h.m22, neg_inv_h21),           /* C = -h22/h21 */
        neg_inv_h21                                 /* D = -1/h21 */
    );
}

/* ============================================================================
 * L5: Generic Conversion
 * ============================================================================ */

/**
 * Generic parameter conversion router.
 *
 * This function routes between any two parameter types using direct
 * conversion formulas. For some rare conversions (e.g., H→S), it
 * may go through an intermediate type (H→Z→S).
 */
matrix2x2_t convert_parameters(matrix2x2_t src, param_type_t src_type,
                                param_type_t dst_type, double z0) {
    if (src_type == dst_type) return src;

    switch (src_type) {
        case PARAM_Z:
            switch (dst_type) {
                case PARAM_Y: return convert_z_to_y(src);
                case PARAM_H: return convert_z_to_h(src);
                case PARAM_G: return convert_z_to_g(src);
                case PARAM_ABCD: return convert_z_to_abcd(src);
                case PARAM_S: return convert_z_to_s(src, z0);
                default: return nan_matrix();
            }
        case PARAM_Y:
            switch (dst_type) {
                case PARAM_Z: return convert_y_to_z(src);
                case PARAM_H: return convert_y_to_h(src);
                case PARAM_G: {
                    matrix2x2_t h = convert_y_to_h(src);
                    return convert_h_to_g(h);
                }
                case PARAM_ABCD: return convert_y_to_abcd_clean(src);
                case PARAM_S: return convert_y_to_s(src, z0);
                default: return nan_matrix();
            }
        case PARAM_H:
            switch (dst_type) {
                case PARAM_Z: return convert_h_to_z(src);
                case PARAM_Y: return convert_h_to_y(src);
                case PARAM_G: return convert_h_to_g(src);
                case PARAM_ABCD: return convert_h_to_abcd(src);
                case PARAM_S: {
                    matrix2x2_t z = convert_h_to_z(src);
                    return convert_z_to_s(z, z0);
                }
                default: return nan_matrix();
            }
        case PARAM_G:
            switch (dst_type) {
                case PARAM_Z: return convert_g_to_z(src);
                case PARAM_H: return convert_g_to_h(src);
                case PARAM_Y: {
                    matrix2x2_t h = convert_g_to_h(src);
                    return convert_h_to_y(h);
                }
                case PARAM_ABCD: {
                    matrix2x2_t z = convert_g_to_z(src);
                    return convert_z_to_abcd(z);
                }
                case PARAM_S: {
                    matrix2x2_t z = convert_g_to_z(src);
                    return convert_z_to_s(z, z0);
                }
                default: return nan_matrix();
            }
        case PARAM_ABCD:
            switch (dst_type) {
                case PARAM_Z: return convert_abcd_to_z(src);
                case PARAM_Y: return convert_abcd_to_y(src);
                case PARAM_H: return convert_abcd_to_h(src);
                case PARAM_G: {
                    matrix2x2_t h = convert_abcd_to_h(src);
                    return convert_h_to_g(h);
                }
                case PARAM_S: return convert_abcd_to_s(src, z0);
                default: return nan_matrix();
            }
        case PARAM_S:
            switch (dst_type) {
                case PARAM_Z: return convert_s_to_z(src, z0);
                case PARAM_Y: return convert_s_to_y(src, z0);
                case PARAM_ABCD: return convert_s_to_abcd(src, z0);
                case PARAM_H: {
                    matrix2x2_t z = convert_s_to_z(src, z0);
                    return convert_z_to_h(z);
                }
                case PARAM_G: {
                    matrix2x2_t z = convert_s_to_z(src, z0);
                    return convert_z_to_g(z);
                }
                default: return nan_matrix();
            }
        default:
            return nan_matrix();
    }
}

int conversion_is_possible(param_type_t src_type, param_type_t dst_type) {
    /* All conversions are mathematically possible except when:
     * - Converting from or to types with different port variable mixes
     *   may fail for degenerate matrices (singular Z, Y, H, G).
     * But in principle, any conversion is defined. */
    (void)src_type;
    (void)dst_type;
    return 1;
}

const char *param_type_name(param_type_t type) {
    switch (type) {
        case PARAM_Z: return "Z";
        case PARAM_Y: return "Y";
        case PARAM_H: return "H";
        case PARAM_G: return "G";
        case PARAM_ABCD: return "ABCD";
        case PARAM_S: return "S";
        default: return "UNKNOWN";
    }
}

/**
 * Round-trip consistency check: convert from src→dst→src and compare.
 *
 * A perfect round-trip means the conversion formulas are numerically
 * correct and the matrix was invertible at both steps.
 */
int conversion_roundtrip_check(matrix2x2_t original, param_type_t type,
                                double z0, double tolerance) {
    /* Test round-trip through all other types */
    param_type_t all_types[] = {PARAM_Z, PARAM_Y, PARAM_H, PARAM_G, PARAM_ABCD, PARAM_S};
    for (int i = 0; i < 6; i++) {
        if (all_types[i] == type) continue;

        matrix2x2_t converted = convert_parameters(original, type, all_types[i], z0);
        /* Check for NaN */
        if (isnan(converted.m11.real)) continue;

        matrix2x2_t back = convert_parameters(converted, all_types[i], type, z0);
        if (isnan(back.m11.real)) continue;

        /* Compare */
        if (!complex_approx_equal(original.m11, back.m11, tolerance)) return 0;
        if (!complex_approx_equal(original.m12, back.m12, tolerance)) return 0;
        if (!complex_approx_equal(original.m21, back.m21, tolerance)) return 0;
        if (!complex_approx_equal(original.m22, back.m22, tolerance)) return 0;
    }
    return 1;
}
