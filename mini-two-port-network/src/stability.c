/**
 * @file stability.c
 * @brief Two-Port Stability Analysis Implementation
 *
 * Implements stability criteria (Rollett K, μ-factor, stability circles,
 * Nyquist analysis) for active two-port networks. These are essential
 * tools for amplifier design to prevent oscillation.
 */

#include "../include/stability.h"
#include "../include/conversion.h"
#include "../include/s_params.h"
#include <float.h>

/* ============================================================================
 * L4: Rollett Stability Analysis
 * ============================================================================ */

/**
 * Stability K-factor for any parameter type.
 *
 * For non-S types, converts to S-parameters first (using Z0 normalization),
 * then computes K from S. This provides a unified stability metric.
 *
 * The Z0 parameter is only meaningful when converting to/from S-parameters.
 * For Z/Y/H/G/ABCD types without an S-parameter reference, Z0 is a free
 * parameter that sets the normalization. Typically Z0 = 50Ω is used.
 */
double stability_rollett_k(matrix2x2_t m, param_type_t type, double z0) {
    matrix2x2_t s;
    if (type == PARAM_S) {
        s = m;
    } else {
        s = convert_parameters(m, type, PARAM_S, z0);
    }

    if (isnan(s.m11.real)) return -1.0;
    return sparams_rollett_k(s);
}

double stability_b1(matrix2x2_t s) {
    double s11_mag_sq = complex_mag(s.m11) * complex_mag(s.m11);
    double s22_mag_sq = complex_mag(s.m22) * complex_mag(s.m22);
    complex_t delta = matrix2x2_det(s);
    double delta_mag_sq = complex_mag(delta) * complex_mag(delta);

    return 1.0 + s11_mag_sq - s22_mag_sq - delta_mag_sq;
}

/**
 * Edwards-Sinsky stability classification.
 *
 * Unconditional stability (K > 1 AND |Δ| < 1):
 *   The device is stable for ALL passive source and load terminations
 *   (any Γ with |Γ| ≤ 1). This is the desired operating region for
 *   amplifiers that must work with arbitrary antennas or loads.
 *
 * Conditional stability (K < 1 OR |Δ| > 1):
 *   The device is stable only for certain terminations. A stability
 *   circle exists on the Smith chart. Terminations inside (or outside,
 *   depending on the case) the circle cause oscillation.
 *
 * Unstable: Device oscillates even with 50Ω terminations.
 *
 * Reference: Edwards & Sinsky, IEEE MTT, 1992
 */
stability_t stability_classify(matrix2x2_t m, param_type_t type, double z0) {
    matrix2x2_t s;
    if (type == PARAM_S) {
        s = m;
    } else {
        s = convert_parameters(m, type, PARAM_S, z0);
    }
    if (isnan(s.m11.real)) return UNSTABLE;

    double k = sparams_rollett_k(s);
    double delta_mag = sparams_delta_mag(s);

    if (k > 1.0 && delta_mag < 1.0) {
        return STABLE_UNCONDITIONAL;
    } else if (k > 0.0) {
        return STABLE_CONDITIONAL;
    } else {
        return UNSTABLE;
    }
}

double stability_mu_factor(matrix2x2_t s) {
    return sparams_mu_factor(s);
}

double stability_mu_prime(matrix2x2_t s) {
    /* μ' = (1 - |s22|²) / (|s11 - conj(s22)*Δ| + |s12*s21|) */
    double s22_mag_sq = complex_mag(s.m22) * complex_mag(s.m22);
    complex_t delta = matrix2x2_det(s);
    complex_t s22_conj = complex_conj(s.m22);
    complex_t term = complex_sub(s.m11, complex_mul(s22_conj, delta));
    double denom = complex_mag(term) + complex_mag(s.m12) * complex_mag(s.m21);

    if (denom < 1e-30) return INFINITY;
    return (1.0 - s22_mag_sq) / denom;
}

/* ============================================================================
 * L5: Stability Circles
 * ============================================================================ */

void stability_circles(matrix2x2_t s, complex_t *cs, double *rs,
                        complex_t *cl, double *rl) {
    sparams_source_stability_circle(s, cs, rs);
    sparams_load_stability_circle(s, cl, rl);
}

/**
 * Check if a source termination ΓS lies in the stable region.
 *
 * The Smith chart center (Γ = 0, representing Z0 = 50Ω) is typically
 * chosen to be in the stable region. A point is stable if it's on the
 * same side of the stability circle as the Smith chart center.
 *
 * If the stability circle encloses the origin:
 *   - Points INSIDE the circle are stable
 *   - Points OUTSIDE are unstable
 * If the stability circle excludes the origin:
 *   - Points INSIDE the circle are unstable
 *   - Points OUTSIDE are stable
 */
int stability_is_source_stable(matrix2x2_t s, complex_t gs) {
    complex_t cs;
    double rs;
    sparams_source_stability_circle(s, &cs, &rs);

    if (isinf(rs)) {
        /* K → ∞ or degenerate: unconditionally stable */
        return 1;
    }

    /* Distance from ΓS to circle center */
    complex_t diff = complex_sub(gs, cs);
    double d = complex_mag(diff);

    /* Distance from origin to circle center */
    double d0 = complex_mag(cs);

    /* Determine if origin is inside the circle */
    int origin_inside = (d0 < rs);

    /* Determine if ΓS is inside the circle */
    int gs_inside = (d < rs);

    /* Stable region is opposite of unstable region */
    if (origin_inside) {
        /* Origin is unstable → outside is stable */
        return !gs_inside;
    } else {
        /* Origin is stable → inside is unstable */
        return !gs_inside;
    }
}

int stability_is_load_stable(matrix2x2_t s, complex_t gl) {
    complex_t cl;
    double rl;
    sparams_load_stability_circle(s, &cl, &rl);

    if (isinf(rl)) return 1;

    complex_t diff = complex_sub(gl, cl);
    double d = complex_mag(diff);
    double d0 = complex_mag(cl);
    int origin_inside = (d0 < rl);
    int gl_inside = (d < rl);

    if (origin_inside) return !gl_inside;
    else return !gl_inside;
}

/* ============================================================================
 * L3: Nyquist Stability Criterion
 * ============================================================================ */

/**
 * Loop gain for a two-port: T = s12*s21*ΓS*ΓL / ((1-s11*ΓS)(1-s22*ΓL)).
 *
 * The Nyquist criterion: if the polar plot of T(ω) encircles the
 * point -1 + j0 in the complex plane, the closed-loop system is
 * unstable. For microwave amplifiers, stability circles provide
 * an equivalent but more convenient test.
 */
complex_t stability_loop_gain(matrix2x2_t s, complex_t gs, complex_t gl) {
    complex_t one = complex_make(1.0, 0.0);
    complex_t num = complex_mul(
        complex_mul(s.m12, s.m21),
        complex_mul(gs, gl)
    );
    complex_t denom = complex_mul(
        complex_sub(one, complex_mul(s.m11, gs)),
        complex_sub(one, complex_mul(s.m22, gl))
    );
    return complex_div(num, denom);
}

/**
 * Oscillation margin: max distance from the oscillation condition.
 *
 * Oscillation: Γin * ΓS = 1 AND Γout * ΓL = 1.
 * This function returns max(|Γin*ΓS - 1|, |Γout*ΓL - 1|).
 * A value near 0 means the network is close to oscillating.
 *
 * This is the "Ohtomo stability margin" used in microwave CAD.
 */
double stability_oscillation_margin(matrix2x2_t s, complex_t gs, complex_t gl) {
    complex_t one = complex_make(1.0, 0.0);
    complex_t gin = sparams_gamma_in(s, gl);
    complex_t gout = sparams_gamma_out(s, gs);

    double margin_in = complex_mag(complex_sub(complex_mul(gin, gs), one));
    double margin_out = complex_mag(complex_sub(complex_mul(gout, gl), one));

    return (margin_in > margin_out) ? margin_in : margin_out;
}

/* ============================================================================
 * L6: Canonical Stability Problems
 * ============================================================================ */

double stability_max_stable_gain(matrix2x2_t s) {
    return sparams_msg(s);
}

/**
 * Find the series resistance needed at the input to make K ≥ 1.
 *
 * Adding series resistance at port 1 changes s11:
 *   s11_new = s11 + (1-s11)*Rser/(Rser+2*Z0)
 *
 * But the simplest approach is to add a pad until K ≥ 1.
 * This function uses iterative search.
 */
double stability_stabilizing_series_r(matrix2x2_t s, double z0) {
    double k = sparams_rollett_k(s);
    if (k >= 1.0) return 0.0;  /* Already stable */

    /* Iterative search: increase series R until K ≥ 1 */
    double r_min = 0.0, r_max = 10.0 * z0;

    for (int iter = 0; iter < 50; iter++) {
        double r_try = (r_min + r_max) / 2.0;

        /* Add series R at input — this modifies s11 and s21 */
        /* For a series R at input in a Z0 system:
           s11' = (s11 * (Z0+R) + R - Z0) / (Z0+R + s11*(R-Z0))
           Simplified: add the series R's S-parameters */
        matrix2x2_t s_series = sparams_series_z(
            complex_make(r_try, 0.0), z0
        );
        /* Cascade: series_R → DUT */
        matrix2x2_t s_combined = convert_abcd_to_s(
            matrix2x2_mul(
                convert_s_to_abcd(s_series, z0),
                convert_s_to_abcd(s, z0)
            ),
            z0
        );

        double k_try = sparams_rollett_k(s_combined);
        if (isinf(k_try) || k_try >= 1.0) {
            r_max = r_try;
        } else {
            r_min = r_try;
        }

        if (r_max - r_min < 0.01) break;
    }

    return (r_min + r_max) / 2.0;
}

/**
 * Find shunt conductance at output needed for stabilization.
 *
 * Adding shunt conductance at port 2 reduces |s22|, increasing K.
 */
double stability_stabilizing_shunt_g(matrix2x2_t s, double z0) {
    double k = sparams_rollett_k(s);
    if (k >= 1.0) return 0.0;

    double g_min = 0.0, g_max = 1.0;  /* Up to 1 S (1 Ω shunt) */

    for (int iter = 0; iter < 50; iter++) {
        double g_try = (g_min + g_max) / 2.0;

        matrix2x2_t s_shunt = sparams_shunt_y(
            complex_make(g_try, 0.0), z0
        );
        /* Cascade: DUT → shunt */
        matrix2x2_t s_combined = convert_abcd_to_s(
            matrix2x2_mul(
                convert_s_to_abcd(s, z0),
                convert_s_to_abcd(s_shunt, z0)
            ),
            z0
        );

        double k_try = sparams_rollett_k(s_combined);
        if (isinf(k_try) || k_try >= 1.0) {
            g_max = g_try;
        } else {
            g_min = g_try;
        }

        if (g_max - g_min < 1e-6) break;
    }

    return (g_min + g_max) / 2.0;
}
