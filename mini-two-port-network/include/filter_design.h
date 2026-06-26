/**
 * @file filter_design.h
 * @brief Filter Design Using Two-Port Network Theory
 *
 * Filters are specialized two-port networks designed to pass signals
 * in certain frequency bands and reject signals in others. The two-port
 * formalism (especially ABCD/cascade parameters) is the natural
 * framework for filter analysis and design.
 *
 * Filter types:
 *   - Low-pass (LP): Pass DC to fc, reject above
 *   - High-pass (HP): Reject DC to fc, pass above
 *   - Band-pass (BP): Pass fl to fh, reject outside
 *   - Band-stop (BS): Reject fl to fh, pass outside (notch)
 *   - All-pass (AP): Pass all frequencies, only phase shift varies
 *
 * Filter approximation types:
 *   - Butterworth: Maximally flat passband, monotonic stopband
 *   - Chebyshev I: Equiripple passband, monotonic stopband
 *   - Chebyshev II: Monotonic passband, equiripple stopband
 *   - Elliptic (Cauer): Equiripple in both passband and stopband
 *   - Bessel: Maximally flat group delay (linear phase)
 *
 * Implementation using two-port theory:
 *   - LC ladder networks (passive, low noise)
 *   - Cascaded biquad sections (active RC)
 *   - Coupled resonators (narrow bandpass)
 *   - Transmission line filters (distributed, for microwave)
 *
 * Reference: Zverev, "Handbook of Filter Synthesis" (1967)
 *            Williams & Taylor, "Electronic Filter Design Handbook"
 * Course: ETH 227-0455 — Filter design
 *         Georgia Tech ECE 6350 — Microwave filter design
 */

#ifndef FILTER_DESIGN_H
#define FILTER_DESIGN_H

#include "two_port.h"

/* ============================================================================
 * L5: Filter Approximation — Transfer Function Generation
 * ============================================================================ */

/**
 * @brief Compute the Butterworth low-pass prototype transfer function.
 *
 * N-th order Butterworth: |H(jω)|² = 1 / (1 + ω^(2N))
 *
 * Poles lie on a unit circle in the left half s-plane:
 *   s_k = -sin((2k-1)π/(2N)) + j*cos((2k-1)π/(2N)), k = 1..N
 *
 * Properties:
 *   - Maximally flat at ω = 0 (first 2N-1 derivatives zero)
 *   - -3 dB at ω = ωc (normalized to 1)
 *   - -20*N dB/decade rolloff
 *   - Monotonic in both passband and stopband
 *
 * @param order Filter order N (≥ 1)
 * @param omega Normalized frequency (1 = cutoff)
 * @return Complex transfer function H(jω)
 *
 * Course: MIT 6.003 — Butterworth approximation
 * Ref: Butterworth, "On the Theory of Filter Amplifiers", 1930
 */
complex_t filter_butterworth_lp(int order, double omega);

/**
 * @brief Compute the Chebyshev Type-I low-pass prototype transfer function.
 *
 * N-th order Chebyshev: |H(jω)|² = 1 / (1 + ε² * T_N²(ω))
 *
 * where T_N(ω) = cos(N * arccos(ω)) is the Chebyshev polynomial.
 * ε = sqrt(10^(RdB/10) - 1) where RdB is the passband ripple in dB.
 *
 * Properties:
 *   - Equiripple in passband (ripple = RdB)
 *   - Monotonic in stopband
 *   - Sharper cutoff than Butterworth (for same order)
 *   - Worse group delay variation than Butterworth
 *
 * @param order Filter order N (≥ 1)
 * @param omega Normalized frequency (1 = passband edge)
 * @param ripple_db Passband ripple in dB (typical: 0.1, 0.5, 1.0, 3.0)
 * @return Complex transfer function H(jω)
 *
 * Course: MIT 6.003 — Chebyshev filter design
 * Ref: Chebyshev polynomials applied by Cauer (1931)
 */
complex_t filter_chebyshev1_lp(int order, double omega, double ripple_db);

/**
 * @brief Compute the Chebyshev Type-II (inverse Chebyshev) low-pass prototype.
 *
 * |H(jω)|² = 1 / (1 + 1/(ε² * T_N²(1/ω)))
 *
 * Properties:
 *   - Monotonic in passband
 *   - Equiripple in stopband
 *   - Same rolloff as Chebyshev-I but different implementation
 *
 * @param order Filter order N
 * @param omega Normalized frequency
 * @param stop_db Minimum stopband attenuation (dB)
 * @return Complex transfer function H(jω)
 */
complex_t filter_chebyshev2_lp(int order, double omega, double stop_db);

/**
 * @brief Compute the elliptic (Cauer) low-pass prototype magnitude.
 *
 * |H(jω)|² = 1 / (1 + ε² * R_N²(ω, ξ))
 *
 * where R_N is the Chebyshev rational function (ratio of Jacobi
 * elliptic functions). Elliptic filters achieve the sharpest possible
 * cutoff for a given order and ripple specification.
 *
 * This implementation approximates the magnitude response.
 *
 * @param order Filter order N
 * @param omega Normalized frequency
 * @param ripple_db Passband ripple (dB)
 * @param stop_db Minimum stopband attenuation (dB)
 * @return |H(jω)| (real magnitude)
 */
double filter_elliptic_lp_mag(int order, double omega,
                               double ripple_db, double stop_db);

/**
 * @brief Compute the Bessel (Thomson) low-pass prototype transfer function.
 *
 * Bessel filters maximize the flatness of group delay at DC.
 * The transfer function is derived from Bessel polynomials.
 *
 * For N=2: H(s) = 3 / (s² + 3s + 3)
 * For N=3: H(s) = 15 / (s³ + 6s² + 15s + 15)
 *
 * Properties:
 *   - Maximally flat group delay (linear phase)
 *   - No overshoot in step response
 *   - Gentler rolloff than Butterworth
 *   - Used in: data communications, audio crossovers
 *
 * @param order Filter order N (1 ≤ N ≤ 8)
 * @param omega Normalized frequency
 * @return Complex transfer function H(jω)
 *
 * Course: Stanford EE102 — Bessel filters
 * Ref: Thomson, "Delay Networks having Maximally Flat Frequency Characteristics", 1949
 */
complex_t filter_bessel_lp(int order, double omega);

/* ============================================================================
 * L5: Frequency Transformations
 * ============================================================================ */

/**
 * @brief Transform a normalized low-pass prototype to a low-pass filter.
 *
 * s → s/ωc  (frequency scaling)
 * Component scaling: L → L/ωc, C → C/ωc
 *
 * @param omega_norm Normalized frequency in prototype domain
 * @param omega_c Actual cutoff frequency (rad/s)
 * @return Scaled frequency in prototype domain
 */
double filter_lp_to_lp(double omega_norm, double omega_c);

/**
 * @brief Transform a normalized low-pass to a high-pass filter.
 *
 * s → ωc/s  (reciprocal transformation, with frequency inversion)
 * Component: L → 1/(ωc*C), C → 1/(ωc*L)
 *
 * @param omega_norm Normalized LP frequency
 * @param omega_c HP cutoff frequency (rad/s)
 * @return Effective normalized frequency
 */
double filter_lp_to_hp(double omega_norm, double omega_c);

/**
 * @brief Transform a normalized low-pass to a band-pass filter.
 *
 * s → (s² + ω0²) / (s * BW)
 * where BW = ωh - ωl (bandwidth), ω0 = sqrt(ωl * ωh) (center)
 *
 * Each LP element becomes a series or parallel LC:
 *   Series L → series LC (resonant at ω0)
 *   Shunt C → parallel LC (resonant at ω0)
 *
 * @param omega_norm Normalized LP frequency
 * @param omega_l Lower band edge (rad/s)
 * @param omega_h Upper band edge (rad/s)
 * @return Effective normalized frequency
 */
double filter_lp_to_bp(double omega_norm, double omega_l, double omega_h);

/**
 * @brief Transform a normalized low-pass to a band-stop filter.
 *
 * s → s * BW / (s² + ω0²)  (the inverse of LP→BP)
 *
 * @param omega_norm Normalized LP frequency
 * @param omega_l Lower band edge (rad/s)
 * @param omega_h Upper band edge (rad/s)
 * @return Effective normalized frequency
 */
double filter_lp_to_bs(double omega_norm, double omega_l, double omega_h);

/* ============================================================================
 * L6: Filter Design — ABCD Cascade
 * ============================================================================ */

/**
 * @brief Design a Butterworth LC ladder low-pass filter.
 *
 * Computes the normalized component values (g-values) for a doubly-terminated
 * LC ladder network. The g-values alternate between series inductance and
 * shunt capacitance (or vice versa depending on starting element).
 *
 * g-values for N=3, Rs=Rl=1:
 *   g1 = 1.0 H  (series L — or 1.0 F if shunt C first)
 *   g2 = 2.0 F  (shunt C)
 *   g3 = 1.0 H  (series L)
 *   g4 = 1.0    (load resistance, normalized)
 *
 * General formula: g_k = 2 * sin((2k-1)π/(2N))
 *
 * @param order N
 * @param g_values Output: array of N normalized values
 * @return 0 on success
 *
 * Course: Georgia Tech ECE 6350 — Filter prototype design
 */
int filter_butterworth_g_values(int order, double *g_values);

/**
 * @brief Compute Chebyshev-I LC ladder g-values.
 *
 * The g-values for Chebyshev filters account for the passband ripple.
 * The computation requires the ripple factor ε = sqrt(10^(RdB/10)-1).
 *
 * g0 = 1 (source)
 * g1 = 2*a1 / sinh(β)  where β = ln(coth(RdB/17.37)), a1 = sin(π/(2N))
 * For k = 2..N: gk = 4*a_{k-1}*a_k / (b_{k-1} * g_{k-1})
 *   where a_k = sin((2k-1)π/(2N)), b_k = sinh²(β) + sin²(kπ/N)
 * g_{N+1} = 1 (for odd N), or coth²(β/4) (for even N)
 *
 * @param order N
 * @param ripple_db Passband ripple in dB
 * @param g_values Output: array of N normalized values
 * @return 0 on success
 */
int filter_chebyshev_g_values(int order, double ripple_db, double *g_values);

/**
 * @brief Build the cascade ABCD matrix for an N-stage LC ladder filter.
 *
 * Given g-values, construct the full filter ABCD matrix by cascading
 * alternating series-L and shunt-C sections.
 *
 * For series-L: ABCD = [[1, jω*L], [0, 1]]
 * For shunt-C: ABCD = [[1, 0], [jω*C, 1]]
 *
 * @param g_values Array of N normalized g-values
 * @param order N
 * @param omega Frequency (normalized or actual)
 * @param start_with_shunt 0 for series-L first, 1 for shunt-C first
 * @return Complete ABCD matrix of the filter
 */
matrix2x2_t filter_build_lc_ladder(const double *g_values, int order,
                                    double omega, int start_with_shunt);

/**
 * @brief Compute the filter response in dB at a given frequency.
 *
 * Converts the ABCD matrix to S21 and returns 20*log10(|s21|).
 *
 * @param abcd Filter ABCD matrix
 * @param z0 Reference impedance (Ω, usually 50 for RF, or 1 for normalized)
 * @return Insertion gain/loss in dB
 */
double filter_response_db(matrix2x2_t abcd, double z0);

/* ============================================================================
 * L6: Filter Design — Active RC Biquad
 * ============================================================================ */

/**
 * @brief Design a Sallen-Key low-pass biquad.
 *
 * A second-order active RC filter using one op-amp (Sallen-Key topology).
 * Transfer function:
 *   H(s) = K / (s²R1R2C1C2 + s(R1C1(K-1) + R1C2 + R2C2) + 1)
 *
 * With unity gain (K=1):
 *   ω0 = 1 / sqrt(R1 * R2 * C1 * C2)
 *   Q = sqrt(R1*R2*C1*C2) / (R1*C1 + R2*C2)  [simplified for K=1]
 *
 * This function computes component values given ω0 and Q.
 *
 * @param omega0 Natural frequency (rad/s)
 * @param q Quality factor (> 0, typically 0.5-10)
 * @param r1 Output: resistor R1 value (Ω)
 * @param r2 Output: resistor R2 value (Ω)
 * @param c1 Output: capacitor C1 value (F)
 * @param c2 Output: capacitor C2 value (F)
 * @return 0 on success
 *
 * Course: Stanford EE101 — Active filter design
 * Ref: Sallen & Key, "A Practical Method of Designing RC Active Filters", 1955
 */
int filter_sallen_key_lp(double omega0, double q,
                          double *r1, double *r2, double *c1, double *c2);

/**
 * @brief Design a multiple-feedback (MFB) low-pass biquad.
 *
 * MFB topology provides inverting gain with better component sensitivity
 * than Sallen-Key.
 *
 * H(s) = -R2/R1 / (s²R2R3C1C2 + s(R2R3C2/R1 + R2C2 + R3C2) + 1)
 *
 * @param omega0 Natural frequency (rad/s)
 * @param q Quality factor
 * @param gain DC gain magnitude
 * @param r1 Output: resistor R1
 * @param r2 Output: resistor R2
 * @param r3 Output: resistor R3
 * @param c1 Output: capacitor C1
 * @param c2 Output: capacitor C2
 * @return 0 on success
 */
int filter_mfb_lp(double omega0, double q, double gain,
                   double *r1, double *r2, double *r3,
                   double *c1, double *c2);

/**
 * @brief Compute the ABCD matrix of a single biquad section.
 *
 * Given a transfer function H(s) = (b2*s² + b1*s + b0) / (s² + a1*s + a0),
 * construct the equivalent ABCD matrix at frequency ω.
 *
 * This allows active filter sections to be cascaded using ABCD multiplication.
 *
 * @param b2,b1,b0 Numerator coefficients
 * @param a1,a0 Denominator coefficients (a2 = 1, normalized)
 * @param omega Angular frequency (rad/s)
 * @return ABCD matrix of the biquad
 */
matrix2x2_t filter_biquad_abcd(double b2, double b1, double b0,
                                double a1, double a0, double omega);

/* ============================================================================
 * L6: Distributed Element Filter Design
 * ============================================================================ */

/**
 * @brief Design a quarter-wave coupled-line bandpass filter.
 *
 * Uses parallel-coupled microstrip lines, each λ/4 long at center frequency.
 * This is the most common planar microwave filter topology.
 *
 * This function computes the even/odd mode impedances for each coupled
 * section given normalized g-values.
 *
 * @param g_values Normalized low-pass prototype g-values
 * @param order N (number of coupled sections = N+1)
 * @param bandwidth_fraction Fractional bandwidth Δ = BW/f0
 * @param z0 System impedance (Ω, typically 50)
 * @param z0e Output: even-mode impedances for N+1 sections
 * @param z0o Output: odd-mode impedances for N+1 sections
 * @return 0 on success
 *
 * Course: Georgia Tech ECE 6350 — Parallel-coupled line filter
 * Ref: Matthaei, Young, Jones, "Microwave Filters, Impedance-Matching
 *      Networks, and Coupling Structures", Ch. 5
 */
int filter_coupled_line_bp(const double *g_values, int order,
                            double bandwidth_fraction, double z0,
                            double *z0e, double *z0o);

/**
 * @brief Compute the response of a coupled-line bandpass filter.
 *
 * Given the even/odd impedances, compute the S-parameters at frequency ω
 * using the coupled-line theory.
 *
 * @param z0e Even-mode impedances array
 * @param z0o Odd-mode impedances array
 * @param n_sections Number of coupled sections
 * @param omega Frequency (rad/s)
 * @param omega0 Center frequency (rad/s)
 * @param z0 System impedance (Ω)
 * @return S-parameter matrix of the filter
 */
matrix2x2_t filter_coupled_line_response(const double *z0e,
                                          const double *z0o,
                                          int n_sections,
                                          double omega, double omega0,
                                          double z0);

/**
 * @brief Compute the group delay of a filter from its S-parameters.
 *
 * τg(ω) = -dφ/dω = -Im(dS21/dω / S21)
 *
 * This is the exact group delay at a single frequency, computed from
 * S-parameters at two closely spaced frequencies using finite differences.
 *
 * @param s1 S-parameters at ω
 * @param s2 S-parameters at ω + Δω
 * @param delta_omega Frequency increment (rad/s)
 * @return Group delay in seconds
 */
double filter_group_delay(matrix2x2_t s1, matrix2x2_t s2, double delta_omega);

#endif /* FILTER_DESIGN_H */
