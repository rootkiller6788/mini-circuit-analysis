/**
 * @file bode_plot.h
 * @brief Bode plot construction and analysis
 *
 * Bode plots (Hendrik Wade Bode, 1940) are the primary tool for
 * frequency-domain analysis and design of LTI systems. They consist
 * of two graphs:
 * 1. Magnitude: 20·log₁₀|H(jω)| dB vs. log₁₀(ω)
 * 2. Phase: ∠H(jω) degrees vs. log₁₀(ω)
 *
 * The genius of Bode's approach is the asymptotic approximation:
 * transfer function → sum of straight-line segments on log-log scale.
 *
 * Reference: Bode, H.W., "Network Analysis and Feedback Amplifier
 * Design" (1945) — the original text that introduced these plots.
 *
 * Course: MIT 6.003, Stanford EE102A, Berkeley EE16B, ETH 227-0427
 */

#ifndef BODE_PLOT_H
#define BODE_PLOT_H

#include <stddef.h>
#include "frequency_response.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * BODE PLOT CONSTRUCTION
 * ============================================================================ */

/**
 * @brief Construct exact Bode plot from transfer function.
 *
 * Evaluates H(jω) at logarithmically-spaced frequency points and
 * computes both magnitude (dB) and phase (degrees).
 *
 * The number of points per decade should be ≥ 20 for smooth plots.
 * For accurate capture of high-Q resonance peaks (>3 dB peaking),
 * use ≥ 100 points per decade.
 *
 * @param tf           Transfer function
 * @param f_start      Start frequency (Hz)
 * @param f_end        End frequency (Hz)
 * @param pts_per_dec  Points per decade for logarithmic sweep
 * @return             Bode plot data, or NULL on error
 *
 * L5 Algorithm: Direct frequency sweep evaluation.
 * For each frequency ω_k = 2π·f_start·10^{k/pts_per_dec}:
 *   1. Evaluate H(jω_k) = N(jω_k)/D(jω_k)
 *   2. magnitude_dB[k] = 20·log₁₀|H(jω_k)|
 *   3. phase_deg[k] = atan2(Im{H}, Re{H}) × 180/π
 *
 * Complexity: O(N·max(m,n)) for N frequency points
 */
bode_plot_t *bode_compute(const tf_polynomial_t *tf,
                           double f_start, double f_end,
                           int pts_per_dec);

/**
 * @brief Construct asymptotic Bode plot from pole-zero data.
 *
 * Uses the straight-line approximation:
 * - Each pole: -20 dB/decade above ω_p, -90° phase (transition over 2 decades)
 * - Each zero: +20 dB/decade above ω_z, +90° phase (transition over 2 decades)
 * - Complex pole pair: -40 dB/decade above ω₀, peaking near ω₀ for ζ < 0.707
 * - Pole at origin (1/s): constant -20 dB/decade, -90° phase at all frequencies
 *
 * L6 Canonical Problem: Asymptotic Bode plot construction by hand
 * is a fundamental skill taught in all circuits courses. The rules:
 *
 * Magnitude construction:
 * 1. Start with DC gain: 20·log₁₀(K) (where K is the DC gain)
 * 2. For each zero at ω_z, add +20 dB/decade slope starting at ω_z
 * 3. For each pole at ω_p, add -20 dB/decade slope starting at ω_p
 * 4. Sum all contributions
 *
 * Phase construction:
 * 1. Start with 0° phase (for K > 0) or -180° (for K < 0)
 * 2. For each zero ω_z: add +45°/decade from ω_z/10 to 10ω_z (total +90°)
 * 3. For each pole ω_p: add -45°/decade from ω_p/10 to 10ω_p (total -90°)
 * 4. Sum all contributions
 *
 * @param pz           Pole-zero transfer function
 * @param f_start      Start frequency (Hz)
 * @param f_end        End frequency (Hz)
 * @param pts_per_dec  Points per decade
 * @return             Bode plot with asymptotic approximation, or NULL
 *
 * Complexity: O(N·(n_zeros + n_poles))
 */
bode_plot_t *bode_asymptotic(const tf_pole_zero_t *pz,
                              double f_start, double f_end,
                              int pts_per_dec);

/**
 * @brief Free a Bode plot data structure.
 */
void bode_free(bode_plot_t *bode);

/* ============================================================================
 * BODE PLOT ANALYSIS
 * ============================================================================ */

/**
 * @brief Find the dominant pole frequency from a Bode plot.
 *
 * The dominant pole is the lowest-frequency pole, which typically
 * determines the -3 dB bandwidth of the system.
 *
 * Strategy: Find the first frequency where the magnitude drops
 * 3 dB below the DC/low-frequency value.
 *
 * @param bode        Bode plot data
 * @return            Dominant pole frequency (Hz), or -1 if not found
 */
double bode_dominant_pole_freq(const bode_plot_t *bode);

/**
 * @brief Find all corner frequencies from a Bode plot.
 *
 * A corner (break) frequency is where the slope changes by
 * ±20 dB/decade or more. Detected by computing the second
 * difference of the magnitude response.
 *
 * @param bode         Bode plot data
 * @param num_corners  Output: number of corners found
 * @return             Array of corner frequencies (Hz), caller frees
 *
 * Complexity: O(N) for N frequency points
 */
double *bode_find_corner_frequencies(const bode_plot_t *bode,
                                       size_t *num_corners);

/**
 * @brief Compute the gain-bandwidth product (GBWP).
 *
 * For an op-amp or amplifier with dominant-pole compensation:
 *   GBWP = A₀·f_p
 * where A₀ = DC gain and f_p = dominant pole frequency.
 *
 * L4: The gain-bandwidth product is constant for a dominant-pole
 * compensated amplifier. This means closed-loop bandwidth is:
 *   f_{CL} = GBWP/G_{CL}
 * where G_{CL} = 1/β is the closed-loop gain.
 *
 * For a unity-gain buffer (β=1): f_{CL} = GBWP = f_T (transition freq).
 *
 * @param bode  Bode plot data
 * @return      Gain-bandwidth product (Hz), or -1 if can't compute
 *
 * Reference: Sedra & Smith (2020), Ch. 9
 * Course: Berkeley EE105
 */
double bode_gain_bandwidth_product(const bode_plot_t *bode);

/**
 * @brief Compute the slope in dB/decade over a specified range.
 *
 * Uses linear regression of magnitude (dB) vs. log₁₀(f).
 *
 * @param bode   Bode plot
 * @param f_low  Lower bound (Hz)
 * @param f_high Upper bound (Hz)
 * @return       Slope (dB/decade): -20 = single pole, -40 = double pole, etc.
 */
double bode_slope_in_range(const bode_plot_t *bode,
                            double f_low, double f_high);

/**
 * @brief Detect resonance peaking in a Bode magnitude plot.
 *
 * Resonance peaking occurs when Q > 1/√2 ≈ 0.707. The peak magnitude
 * for a second-order lowpass:
 *   |H|_peak = Q/√(1 - 1/(4Q²))  for Q > 1/√2
 *   Peak frequency: ω_peak = ω₀·√(1 - 1/(2Q²))
 *
 * For Q >> 1: |H|_peak ≈ Q, ω_peak ≈ ω₀.
 *
 * @param bode       Bode plot
 * @param peak_freq  Output: frequency at resonance peak (Hz)
 * @param peak_db    Output: peak magnitude above baseline (dB)
 * @return           1 if peaking detected, 0 otherwise
 */
int bode_detect_peaking(const bode_plot_t *bode,
                         double *peak_freq, double *peak_db);

/**
 * @brief Verify Bode's gain-phase relationship.
 *
 * For a minimum-phase system, the phase at frequency ω₀ is related
 * to the slope of the magnitude response:
 *   φ(ω₀) ≈ (π/2) · d(log|H|)/d(log ω)|_{ω=ω₀}
 *
 * More precisely (Bode's integral):
 *   φ(ω₀) = (1/π) ∫_{-∞}^{∞} dM/du · ln(coth|u|/2) du
 * where M = ln|H| and u = ln(ω/ω₀).
 *
 * This function computes the approximate phase from the magnitude
 * slope and compares with the actual computed phase.
 *
 * @param bode         Bode plot data
 * @param approx_phase Output: phase estimated from magnitude slope (deg)
 * @param error_rms    Output: RMS error between actual and estimated phase (deg)
 * @return             0 on success, -1 on error
 *
 * L4: Bode's Gain-Phase Relation — one of the most profound results
 * in network theory. For minimum-phase systems, the magnitude response
 * determines the phase response (up to an all-pass factor).
 *
 * Reference: Bode (1945), Ch. 14
 * Course: Stanford EE102A, ETH 227-0427
 */
int bode_gain_phase_relation(const bode_plot_t *bode,
                              double *approx_phase, double *error_rms);

/* ============================================================================
 * FREQUENCY RESPONSE VISUALIZATION HELPERS
 * ============================================================================ */

/**
 * @brief Generate text-based Bode magnitude plot (ASCII art).
 *
 * Creates a simple text representation for quick terminal inspection
 * of frequency response. Each row = one frequency point.
 *
 * @param bode       Bode plot data
 * @param buffer     Output buffer
 * @param buf_size   Buffer size
 * @param plot_width Width of the plot in characters
 * @return           Number of characters written, or -1 on error
 *
 * L7 Application: Quick terminal-based visualization for embedded
 * systems where graphical displays are unavailable.
 */
int bode_ascii_plot(const bode_plot_t *bode, char *buffer,
                     size_t buf_size, int plot_width);

/**
 * @brief Export Bode plot data to CSV format.
 *
 * Columns: frequency_hz, magnitude_db, phase_deg
 *
 * @param bode    Bode plot data
 * @param buffer  Output buffer for CSV string
 * @param size    Buffer size
 * @return        Number of characters written, or -1 on error
 */
int bode_export_csv(const bode_plot_t *bode, char *buffer, size_t size);

#ifdef __cplusplus
}
#endif

#endif /* BODE_PLOT_H */
