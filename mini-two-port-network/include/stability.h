/**
 * @file stability.h
 * @brief Two-Port Network Stability Analysis
 *
 * Stability analysis is critical for active two-port networks (amplifiers,
 * oscillators). An unconditionally stable amplifier remains stable for all
 * passive source and load terminations. A potentially unstable amplifier
 * may oscillate for certain termination impedances.
 *
 * Key stability metrics:
 *   - Rollett K-factor (Rollett, 1962)
 *   - μ-factor (Edwards & Sinsky, 1992)
 *   - Stability circles on the Smith chart
 *   - Nyquist stability criterion (applied to two-ports)
 *   - Normalized determinant function (NDF)
 *
 * Reference: Gonzalez, "Microwave Transistor Amplifiers: Analysis and Design"
 * Course: Stanford EE359 — RF amplifier design
 *         Georgia Tech ECE 6350 — Microwave amplifier stability
 */

#ifndef STABILITY_H
#define STABILITY_H

#include "two_port.h"

/* ============================================================================
 * L4: Rollett Stability Theory
 * ============================================================================ */

/**
 * @brief Compute the Rollett stability K-factor from any parameter type.
 *
 * General formula (from S-parameters):
 *   K = (1 - |s11|² - |s22|² + |Δ|²) / (2 * |s12 * s21|)
 *   where Δ = s11*s22 - s12*s21
 *
 * Interpretation:
 *   K > 1: Potentially unconditionally stable (need also check |Δ| < 1 or B1 > 0)
 *   K = 1: Marginal stability
 *   K < 1: Potentially unstable (conditional stability)
 *
 * @param m Parameter matrix
 * @param type Parameter type
 * @param z0 Reference impedance (required only for S-parameters; set to 0 for others)
 * @return Rollett K-factor, or -1.0 if computation fails
 *
 * Reference: Rollett, "Stability and Power-Gain Invariants", IRE Trans. CT, 1962
 */
double stability_rollett_k(matrix2x2_t m, param_type_t type, double z0);

/**
 * @brief Compute the auxiliary stability condition B1.
 *
 * For unconditional stability: K > 1 AND B1 > 0 (or equivalently |Δ| < 1)
 * B1 = 1 + |s11|² - |s22|² - |Δ|²
 *
 * If B1 < 0, the device is potentially unstable even if K > 1.
 *
 * @param s S-parameter matrix
 * @return B1 value (> 0 required)
 */
double stability_b1(matrix2x2_t s);

/**
 * @brief Determine stability classification of a two-port.
 *
 * Evaluates both K > 1 and |Δ| < 1 (the Edwards-Sinsky criteria).
 *
 * @param m Parameter matrix (any type, converted to S internally)
 * @param type Parameter type
 * @param z0 Reference impedance
 * @return STABLE_UNCONDITIONAL, STABLE_CONDITIONAL, or UNSTABLE
 */
stability_t stability_classify(matrix2x2_t m, param_type_t type, double z0);

/**
 * @brief Compute the μ stability factor (Edwards-Sinsky, 1992).
 *
 * μ = (1 - |s11|²) / (|s22 - conj(s11)*Δ| + |s12*s21|)
 *
 * μ > 1 → unconditionally stable (single-criterion test).
 * μ' = (1 - |s22|²) / (|s11 - conj(s22)*Δ| + |s12*s21|) is the dual.
 *
 * The μ-factor has the advantage over K of being a single number and
 * having a geometric interpretation: μ is the minimum distance from
 * the center of the Smith chart to the instability region.
 *
 * Larger μ = more stable (more margin against oscillation).
 *
 * @param s S-parameter matrix
 * @return μ-factor
 *
 * Reference: Edwards & Sinsky, IEEE Trans. MTT, 1992
 */
double stability_mu_factor(matrix2x2_t s);

/**
 * @brief Compute μ' (mu-prime), the dual stability factor.
 *
 * μ' = (1 - |s22|²) / (|s11 - conj(s22)*Δ| + |s12*s21|)
 *
 * For real stability, need both μ > 1 and μ' > 1, but in practice
 * μ > 1 alone is sufficient.
 *
 * @param s S-parameter matrix
 * @return μ'-factor
 */
double stability_mu_prime(matrix2x2_t s);

/* ============================================================================
 * L5: Stability Circle Analysis
 * ============================================================================ */

/**
 * @brief Compute source and load stability circles.
 *
 * Source stability circle:
 *   Center: CS = (conj(s11 - Δ*s22*)) / (|s11|² - |Δ|²)
 *   Radius: RS = |s12*s21| / | |s11|² - |Δ|² |
 *
 * Load stability circle:
 *   Center: CL = (conj(s22 - Δ*s11*)) / (|s22|² - |Δ|²)
 *   Radius: RL = |s12*s21| / | |s22|² - |Δ|² |
 *
 * Interpretation:
 *   - If the circle contains the origin (center of Smith chart), the
 *     Smith chart center is in the unstable region.
 *   - If RS + |CS| < 1: entire circle inside Smith chart, stable outside
 *   - If |CS| - RS > 1: entire circle outside Smith chart, stable everywhere
 *
 * @param s S-parameter matrix
 * @param cs Output: source stability circle center
 * @param rs Output: source stability circle radius
 * @param cl Output: load stability circle center
 * @param rl Output: load stability circle radius
 *
 * Course: Stanford EE359 — Smith chart stability design
 */
void stability_circles(matrix2x2_t s, complex_t *cs, double *rs,
                        complex_t *cl, double *rl);

/**
 * @brief Check if a given source termination ΓS is in the stable region.
 *
 * A termination is stable if it lies outside the source stability circle
 * (when the circle center + radius does not enclose the origin, i.e.,
 * the stable region is outside the circle).
 *
 * @param s S-parameter matrix
 * @param gs Source reflection coefficient to test
 * @return 1 if stable, 0 if potentially unstable
 */
int stability_is_source_stable(matrix2x2_t s, complex_t gs);

/**
 * @brief Check if a given load termination ΓL is in the stable region.
 *
 * @param s S-parameter matrix
 * @param gl Load reflection coefficient to test
 * @return 1 if stable, 0 if potentially unstable
 */
int stability_is_load_stable(matrix2x2_t s, complex_t gl);

/* ============================================================================
 * L3: Nyquist Stability Criterion for Two-Ports
 * ============================================================================ */

/**
 * @brief Calculate the open-loop gain for Nyquist stability analysis.
 *
 * For a two-port amplifier, the loop gain T = s12*s21*ΓS*ΓL / ((1-s11*ΓS)(1-s22*ΓL))
 *
 * The Nyquist criterion states: if the polar plot of T(ω) encircles -1,
 * the system is unstable when the loop is closed.
 *
 * @param s S-parameter matrix at frequency ω
 * @param gs Source reflection coefficient
 * @param gl Load reflection coefficient
 * @return Complex loop gain T
 */
complex_t stability_loop_gain(matrix2x2_t s, complex_t gs, complex_t gl);

/**
 * @brief Compute the oscillation condition for a two-port.
 *
 * Oscillation occurs when:
 *   Γin * ΓS = 1  AND  Γout * ΓL = 1
 *
 * This function checks how close a given (ΓS, ΓL) pair is to the
 * oscillation condition, returning the worst-case product magnitude.
 *
 * @param s S-parameter matrix
 * @param gs Source reflection coefficient
 * @param gl Load reflection coefficient
 * @return max(|Γin*ΓS - 1|, |Γout*ΓL - 1|), 0 means oscillation
 */
double stability_oscillation_margin(matrix2x2_t s, complex_t gs, complex_t gl);

/* ============================================================================
 * L6: Canonical Stability Problems
 * ============================================================================ */

/**
 * @brief Compute the maximum stable gain (MSG) for conditionally stable devices.
 *
 * MSG = |s21 / s12|  (for K < 1)
 *
 * This is the maximum gain achievable while maintaining stability
 * by resistive loading (adding loss to increase K to exactly 1).
 *
 * @param s S-parameter matrix
 * @return Maximum stable gain (linear, not dB)
 */
double stability_max_stable_gain(matrix2x2_t s);

/**
 * @brief Determine the stability-improving resistance to add at input.
 *
 * Adding series resistance at the input increases |s11|, which
 * increases K. This function computes the minimum series resistance
 * needed to make K = 1, assuming the device was originally K < 1.
 *
 * Uses an iterative search: increase Zser until K ≥ 1.
 *
 * @param s Original S-parameters (K < 1)
 * @param z0 Reference impedance (Ω)
 * @return Required series resistance in Ω, or -1 if cannot be stabilized
 */
double stability_stabilizing_series_r(matrix2x2_t s, double z0);

/**
 * @brief Determine the stability-improving shunt conductance to add.
 *
 * Adding shunt conductance at the output increases |s22|, which
 * also increases K.
 *
 * @param s Original S-parameters (K < 1)
 * @param z0 Reference impedance (Ω)
 * @return Required shunt conductance in S, or -1 if cannot be stabilized
 */
double stability_stabilizing_shunt_g(matrix2x2_t s, double z0);

#endif /* STABILITY_H */
