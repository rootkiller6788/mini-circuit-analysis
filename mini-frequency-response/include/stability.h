/**
 * @file stability.h
 * @brief Stability analysis for feedback circuits and systems
 *
 * Stability is the most critical property of any feedback system.
 * An unstable circuit oscillates uncontrollably and is useless.
 * This header provides the tools for stability assessment using
 * frequency-domain methods.
 *
 * Key stability criteria:
 * 1. Routh-Hurwitz (algebraic, polynomial coefficients only)
 * 2. Nyquist (graphical, open-loop frequency response)
 * 3. Bode (gain and phase margins from open-loop Bode plot)
 * 4. Root locus (pole migration with varying gain)
 *
 * Reference:
 * - Nyquist, "Regeneration Theory" (1932), Bell System Tech. J.
 * - Bode, "Network Analysis and Feedback Amplifier Design" (1945)
 * - Ogata, "Modern Control Engineering" (2010), Ch. 7-8
 *
 * Course: MIT 6.003, Stanford EE102A, Berkeley EE16B, ETH 227-0427
 */

#ifndef STABILITY_H
#define STABILITY_H

#include <stddef.h>
#include "frequency_response.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * ROUTH-HURWITZ STABILITY CRITERION
 * ============================================================================ */

/**
 * @brief Apply Routh-Hurwitz stability criterion to a polynomial.
 *
 * Constructs the Routh array and checks the first column for sign changes.
 *
 * Routh array construction (for polynomial a₀sⁿ + a₁sⁿ⁻¹ + ... + a_n):
 *
 * Row sⁿ:   a₀   a₂   a₄   ...
 * Row sⁿ⁻¹: a₁   a₃   a₅   ...
 * Row sⁿ⁻²: b₁   b₂   b₃   ...  (computed)
 * Row sⁿ⁻³: c₁   c₂   c₃   ...  (computed)
 * ...
 * Row s⁰:   ...
 *
 * Where:
 *   b₁ = (a₁a₂ - a₀a₃)/a₁
 *   b₂ = (a₁a₄ - a₀a₅)/a₁
 *   c₁ = (b₁a₃ - a₁b₂)/b₁
 *   c₂ = (b₁a₅ - a₁b₃)/b₁
 *
 * The number of RHP roots = number of sign changes in the first column.
 *
 * Special cases handled:
 * 1. Zero in first column → replace with small ε and take limit
 * 2. Entire row of zeros → use auxiliary polynomial
 *    (indicates roots symmetrically placed about origin)
 *
 * @param coeffs   Polynomial coefficients a₀, a₁, ..., a_n
 * @param order    Polynomial degree n
 * @return         Routh-Hurwitz analysis result
 *
 * L4 Fundamental Law: Routh-Hurwitz criterion is a necessary AND
 * sufficient condition for all roots of a real polynomial to have
 * negative real parts. It is purely algebraic, requiring no root-finding.
 *
 * Complexity: O(n²) for Routh array construction
 * Reference: Hurwitz, "Über die Bedingungen..." (1895), Mathematische Annalen
 */
routh_hurwitz_t *stability_routh_hurwitz(const double *coeffs, size_t order);

/**
 * @brief Free a Routh-Hurwitz result.
 */
void routh_hurwitz_free(routh_hurwitz_t *rh);

/**
 * @brief Check if a transfer function denominator is Hurwitz
 * (all poles in LHP).
 *
 * Convenience wrapper: extracts denominator coefficients and
 * applies Routh-Hurwitz criterion.
 *
 * @param tf  Transfer function
 * @return    1 = stable (all poles in LHP)
 *            0 = unstable (one or more poles in RHP or on jω-axis)
 *           -1 = error
 */
int stability_is_hurwitz(const tf_polynomial_t *tf);

/* ============================================================================
 * NYQUIST STABILITY CRITERION
 * ============================================================================ */

/**
 * @brief Compute the Nyquist contour data from a transfer function.
 *
 * The Nyquist plot is the polar plot of L(jω) where L = loop gain.
 * It maps the Nyquist contour (jω-axis + infinite semicircle in RHP)
 * to the L(s)-plane.
 *
 * Nyquist Stability Criterion (L4 Fundamental Law):
 *   Z = N + P
 * where:
 *   Z = number of closed-loop unstable poles (in RHP)
 *   N = number of clockwise encirclements of (-1, j0) by L(jω)
 *   P = number of open-loop unstable poles (in RHP)
 *
 * For a stable closed-loop system: Z = 0, so N = -P.
 * For open-loop stable systems (P=0): N must be 0 (no encirclements).
 *
 * The Nyquist criterion is more powerful than Bode because it:
 * - Handles open-loop unstable systems (P > 0)
 * - Handles non-minimum-phase systems
 * - Provides information about conditional stability
 *
 * @param loop_tf      Loop transfer function L(s)
 * @param f_start      Start frequency (Hz)
 * @param f_end        End frequency (Hz)
 * @param n_points     Number of frequency points (≥200 recommended)
 * @return             Array of complex L(jω) values (caller frees), or NULL
 *
 * Complexity: O(N·max(m,n)) for evaluation + O(N) for encirclement counting
 * Reference: Nyquist (1932)
 */
double _Complex *stability_nyquist_contour(const tf_polynomial_t *loop_tf,
                                              double f_start, double f_end,
                                              size_t n_points);

/**
 * @brief Count encirclements of the point (-1, j0) by the Nyquist plot.
 *
 * Uses the winding number algorithm (sum of angle changes).
 *
 * @param nyquist_data  Array of complex L(jω) values
 * @param n             Number of data points
 * @return              Number of clockwise encirclements (positive = CW)
 *
 * L4: The winding number (also called the Cauchy index) is a
 * topological invariant that counts how many times a closed curve
 * encircles a point. It is computed by integrating d(arg(L+1))
 * around the Nyquist contour.
 */
int stability_count_encirclements(const double _Complex *nyquist_data,
                                    size_t n);

/**
 * @brief Determine closed-loop stability from Nyquist plot.
 *
 * @param loop_tf       Loop transfer function L(s)
 * @param f_start       Start frequency (Hz)
 * @param f_end         End frequency (Hz)
 * @param n_points      Number of frequency points
 * @param is_stable     Output: 1 = stable, 0 = unstable
 * @param encirclements Output: number of CW encirclements
 * @return              0 on success, -1 on error
 */
int stability_nyquist_check(const tf_polynomial_t *loop_tf,
                             double f_start, double f_end,
                             size_t n_points,
                             int *is_stable, int *encirclements);

/* ============================================================================
 * STABILITY MARGINS (GAIN & PHASE MARGIN)
 * ============================================================================ */

/**
 * @brief Compute stability margins from open-loop frequency response.
 *
 * Finds:
 * 1. Gain crossover frequency ω_gc where |L(jω_gc)| = 1 (0 dB)
 * 2. Phase crossover frequency ω_pc where ∠L(jω_pc) = -180°
 * 3. Phase margin: PM = 180° + ∠L(jω_gc)
 * 4. Gain margin: GM = 1/|L(jω_pc)|, GM_dB = -20·log₁₀|L(jω_pc)|
 *
 * Requires open-loop frequency response data (can be computed from
 * loop_tf via bode_compute or directly from measurement).
 *
 * Interpolation is used between frequency points for more accurate
 * crossover frequency determination.
 *
 * @param resp     Open-loop frequency response (magnitude+phase)
 * @return         Stability margins, or zeroed struct on error
 *
 * L6 Canonical Problem: Determining stability margins is a fundamental
 * design verification step for any feedback amplifier or control system.
 * Typical design targets: PM ≥ 45°, GM ≥ 6 dB (factor of 2).
 *
 * Course: MIT 6.003, Berkeley EE105, Stanford EE247
 */
stability_margin_t stability_margins_from_response(
    const freq_response_t *resp);

/**
 * @brief Compute stability margins directly from loop transfer function.
 *
 * Convenience function: computes frequency response and then margins.
 *
 * @param loop_tf  Loop transfer function L(s)
 * @return         Stability margins
 */
stability_margin_t stability_margins(const tf_polynomial_t *loop_tf);

/* ============================================================================
 * ROOT LOCUS ANALYSIS (L8 Advanced)
 * ============================================================================ */

/**
 * @brief Compute root locus data for varying gain K.
 *
 * Root locus: plot of closed-loop pole locations as gain K varies
 * from 0 to ∞ for the system:
 *   H(s) = K·N(s)/(D(s) + K·N(s))
 *
 * Evans' Root Locus Rules (1948):
 * 1. Number of branches = number of poles (n)
 * 2. Branches start at open-loop poles (K=0) and end at open-loop
 *    zeros (K=∞). Excess branches go to infinity along asymptotes.
 * 3. The root locus lies on the real axis to the left of an odd
 *    number of real poles and zeros.
 * 4. Asymptote angles: θ_k = (2k+1)·180°/(n-m) for k = 0, 1, ..., n-m-1
 * 5. Asymptote centroid: σ_a = (Σp_j - Σz_i)/(n-m)
 * 6. Breakaway/break-in points occur where dK/ds = 0.
 * 7. Angle of departure from complex pole p_k:
 *    θ_dep = 180° - Σ∠(p_k - p_j) + Σ∠(p_k - z_i)
 *
 * @param loop_tf   Loop transfer function L(s) = N(s)/D(s)
 * @param K_values  Array of gain values to evaluate
 * @param n_K       Number of gain values
 * @param poles_out Output: 2D array poles_out[k][j] = j-th pole for gain K[k]
 *                 (caller allocates n_K × den_order doubles)
 * @return          0 on success, -1 on error
 *
 * L8 Advanced: Root locus was developed by Walter R. Evans at
 * North American Aviation (1950s). It remains the most intuitive
 * graphical method for understanding how feedback affects pole
 * locations and thus stability and transient response.
 *
 * Reference: Evans, "Control System Dynamics" (1954)
 * Course: MIT 6.003, Stanford EE102A
 */
int stability_root_locus(const tf_polynomial_t *loop_tf,
                          const double *K_values, size_t n_K,
                          double _Complex *poles_out);

/**
 * @brief Compute the gain margin using the root locus method.
 *
 * Finds the gain K where a pole pair crosses the jω-axis.
 * This K is the gain margin (linear).
 *
 * @param loop_tf   Loop transfer function
 * @param K_margin  Output: gain at stability boundary
 * @param freq_hz   Output: oscillation frequency at boundary (Hz)
 * @return          0 on success, -1 on error
 */
int stability_gain_margin_root_locus(const tf_polynomial_t *loop_tf,
                                       double *K_margin, double *freq_hz);

/* ============================================================================
 * STABILITY OF SWITCHED-CAPACITOR CIRCUITS (L7 Application)
 * ============================================================================ */

/**
 * @brief Check stability of a two-phase switched-capacitor filter.
 *
 * Switched-capacitor (SC) circuits use capacitors and MOSFET switches
 * to emulate resistors: R_eq = T/(C·f_clk), where T = 1 or 2 depending
 * on the integration scheme.
 *
 * Stability is checked by analyzing the z-domain transfer function
 * and verifying all poles are inside the unit circle (|z| < 1).
 *
 * This function converts the discrete-time transfer function poles
 * from s-domain (via bilinear transform) and checks the unit circle
 * condition.
 *
 * @param tf_s        Continuous-time prototype
 * @param f_clk       Switching clock frequency (Hz)
 * @param is_stable   Output: 1 = SC implementation stable, 0 = unstable
 * @return            0 on success, -1 on error
 *
 * L7 Application: Switched-capacitor filters are widely used in
 * mixed-signal ICs (e.g., anti-aliasing filters in ADC front-ends,
 * programmable gain amplifiers). Stability must be verified in the
 * z-domain after s→z mapping.
 *
 * Reference: Gregorian & Temes, "Analog MOS Integrated Circuits
 * for Signal Processing" (1986)
 * Course: Berkeley EE105, Stanford EE247
 */
int stability_sc_filter(const tf_polynomial_t *tf_s,
                          double f_clk, int *is_stable);

#ifdef __cplusplus
}
#endif

#endif /* STABILITY_H */
