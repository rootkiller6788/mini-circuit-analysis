/**
 * @file network_synthesis.h
 * @brief Network Synthesis from Two-Port Parameters
 *
 * Network synthesis is the inverse of network analysis: given a desired
 * two-port parameter matrix (often derived from a transfer function),
 * synthesize a physical circuit that realizes those parameters.
 *
 * Key synthesis methods:
 *   - Foster's reactance theorem (Foster canonical forms)
 *   - Cauer synthesis (continued fraction expansion)
 *   - Darlington synthesis (lossless two-port terminated in R)
 *   - Brune synthesis (general positive-real function → RLC network)
 *   - Bott-Duffin synthesis (minimum number of elements)
 *
 * This module focuses on synthesis applicable to two-port networks:
 * converting between L, T, and π equivalent circuits, and basic
 * impedance-matching network synthesis.
 *
 * Reference: Van Valkenburg, "Introduction to Modern Network Synthesis" (1960)
 * Course: Illinois ECE 451 — Network synthesis
 *         TU Munich — Netzwerksynthese
 */

#ifndef NETWORK_SYNTHESIS_H
#define NETWORK_SYNTHESIS_H

#include "two_port.h"

/* ============================================================================
 * L5: Equivalent Circuit Synthesis
 * ============================================================================ */

/**
 * @brief Synthesize Z-parameters into a T-equivalent circuit.
 *
 * Any reciprocal two-port can be represented as a T-network:
 *
 *        Za       Zb
 *   P1 ————┬———————┬——— P2
 *          │       │
 *          Zc
 *          │
 *         GND
 *
 *   Za = z11 - z12
 *   Zb = z22 - z12
 *   Zc = z12  (= z21 for reciprocal)
 *
 * This is the universal Z-parameter equivalent. The T-network is
 * physically realizable with passive components if Za, Zb, Zc satisfy
 * the positive-real conditions.
 *
 * @param z Z-parameter matrix (must be reciprocal)
 * @param za Output: first series arm impedance
 * @param zb Output: second series arm impedance
 * @param zc Output: shunt arm impedance
 * @return 0 on success, -1 if network is non-reciprocal
 */
int synthesize_t_network(matrix2x2_t z, complex_t *za, complex_t *zb,
                         complex_t *zc);

/**
 * @brief Synthesize Y-parameters into a π-equivalent circuit.
 *
 * Any reciprocal two-port can be represented as a π-network:
 *
 *        Yb
 *   P1 —┬———┬——— P2
 *       Ya  Yc
 *       │   │
 *      GND GND
 *
 *   Ya = y11 + y12  (since y12 is negative in typical π representation)
 *   Yb = -y12       (series admittance = -mutual)
 *   Yc = y22 + y12
 *
 * @param y Y-parameter matrix (must be reciprocal)
 * @param ya Output: first shunt admittance
 * @param yb Output: series admittance
 * @param yc Output: second shunt admittance
 * @return 0 on success, -1 if non-reciprocal
 */
int synthesize_pi_network(matrix2x2_t y, complex_t *ya, complex_t *yb,
                          complex_t *yc);

/**
 * @brief Convert between T-network and π-network element values.
 *
 * T → π conversion (Y-Δ transform, generalized):
 *   Ya = Zb / (Za*Zb + Zb*Zc + Zc*Za)
 *   Yb = Zc / (Za*Zb + Zb*Zc + Zc*Za)
 *   Yc = Za / (Za*Zb + Zb*Zc + Zc*Za)
 *
 * π → T conversion (Δ-Y transform):
 *   Za = Yc / (Ya*Yb + Yb*Yc + Yc*Ya)
 *   Zb = Ya / (Ya*Yb + Yb*Yc + Yc*Ya)
 *   Zc = Yb / (Ya*Yb + Yb*Yc + Yc*Ya)
 *
 * This is the electrical dual of the Δ-Y resistor transformation
 * (Kennelly's theorem), generalized to complex impedances.
 *
 * @param za_in First T-arm impedance (for T→π) or output (for π→T)
 * @param zb_in Second T-arm impedance (for T→π) or output (for π→T)
 * @param zc_in T-shunt or π-series impedance
 * @param to_pi 1 for T→π, 0 for π→T
 * @param elem1_out Output: first element of target network
 * @param elem2_out Output: second element of target network
 * @param elem3_out Output: third element of target network
 *
 * Course: Berkeley EE105 — Δ-Y transformation
 * Ref: Kennelly, "Equivalence of triangles and stars in conducting networks", 1899
 */
void synthesize_t_pi_convert(complex_t za_in, complex_t zb_in,
                              complex_t zc_in, int to_pi,
                              complex_t *elem1_out, complex_t *elem2_out,
                              complex_t *elem3_out);

/**
 * @brief Synthesize a lossless ladder network from a driving-point impedance.
 *
 * Given a rational impedance function Z(s), synthesize a lossless
 * LC ladder network terminated in a resistor using Cauer's method
 * (continued fraction expansion).
 *
 * This function performs one step of the continued fraction expansion
 * for a given impedance at a single frequency.
 *
 * For a full synthesis, this function would be called iteratively.
 *
 * @param zin Target input impedance (at frequency ω)
 * @param omega Angular frequency (rad/s)
 * @param n_stages Number of LC stages desired
 * @param series_l Output array: series inductance values (H)
 * @param shunt_c Output array: shunt capacitance values (F)
 * @return 0 on success, -1 if synthesis fails
 */
int synthesize_cauer_ladder(complex_t zin, double omega, int n_stages,
                             double *series_l, double *shunt_c);

/**
 * @brief Synthesize a matching network for conjugate match at a single frequency.
 *
 * Given source impedance ZS and load impedance ZL, design an L-network
 * that provides conjugate match at frequency ω.
 *
 * Algorithm:
 * 1. Convert ZS to admittance YS = 1/ZS
 * 2. Compute required transformation Q
 * 3. Determine L-network topology (series/shunt element types)
 * 4. Calculate component values
 *
 * Two solutions exist (low-pass and high-pass type). This function
 * returns the low-pass solution.
 *
 * @param zs Source impedance (Ω, complex)
 * @param zl Load impedance (Ω, complex)
 * @param omega Angular frequency (rad/s)
 * @param l_value Output: inductor value (H), 0 if capacitor used instead
 * @param c_value Output: capacitor value (F), 0 if inductor used instead
 * @return 0 on success, -1 if matching is impossible
 *
 * Course: Stanford EE359 — L-network matching design
 * Ref: Bowick, "RF Circuit Design", Ch. 4
 */
int synthesize_l_match(complex_t zs, complex_t zl, double omega,
                        double *l_value, double *c_value);

/**
 * @brief Synthesize a π-network matching circuit.
 *
 * For broadband or harmonic-suppression matching, π-networks provide
 * an extra degree of freedom (the loaded Q can be chosen).
 *
 * Given ZS, ZL, and desired operating Q, compute the three reactance values.
 *
 * @param zs Source impedance (Ω)
 * @param zl Load impedance (Ω)
 * @param omega Angular frequency (rad/s)
 * @param q_loaded Desired loaded quality factor (> 0)
 * @param ya Output: input shunt admittance (S)
 * @param yb Output: series admittance (S)
 * @param yc Output: output shunt admittance (S)
 * @return 0 on success, -1 if infeasible
 */
int synthesize_pi_match(complex_t zs, complex_t zl, double omega,
                         double q_loaded,
                         complex_t *ya, complex_t *yb, complex_t *yc);

/**
 * @brief Synthesize a T-network matching circuit.
 *
 * Dual of the π-network synthesis. The T-network is preferred when
 * harmonic rejection (series LC resonators) is needed.
 *
 * @param zs Source impedance (Ω)
 * @param zl Load impedance (Ω)
 * @param omega Angular frequency (rad/s)
 * @param q_loaded Desired loaded Q
 * @param za Output: first series impedance (Ω)
 * @param zb Output: second series impedance (Ω)
 * @param zc Output: shunt impedance (Ω)
 * @return 0 on success, -1 if infeasible
 */
int synthesize_t_match(complex_t zs, complex_t zl, double omega,
                        double q_loaded,
                        complex_t *za, complex_t *zb, complex_t *zc);

/* ============================================================================
 * L4: Foster's Reactance Theorem
 * ============================================================================ */

/**
 * @brief Verify Foster's reactance theorem for a given LC driving-point function.
 *
 * Foster's reactance theorem states that for a lossless one-port
 * (pure LC network), the reactance X(ω) = Im(Z(ω)) is a strictly
 * monotonically increasing function of frequency:
 *
 *   dX/dω > 0 for all ω
 *
 * Poles and zeros alternate on the jω axis.
 *
 * This function computes X(ω1) and X(ω2) and checks monotonicity.
 *
 * @param z1 Impedance at ω1
 * @param z2 Impedance at ω2 (ω2 > ω1)
 * @return 1 if X(ω2) > X(ω1) (Foster condition satisfied),
 *         0 if violated, -1 if not purely reactive
 *
 * Course: MIT 6.003 — Foster's reactance theorem
 * Ref: Foster, "A Reactance Theorem", Bell System Tech. J., 1924
 */
int foster_reactance_check(complex_t z1, complex_t z2);

/**
 * @brief Compute the poles and zeros of a two-element LC impedance.
 *
 * For Z(s) = sL || 1/(sC) = sL / (1 + s²LC):
 *   Zero at s = 0 (DC short through inductor)
 *   Pole at s = j/√(LC) (resonance)
 *
 * For Z(s) = sL + 1/(sC):
 *   Zero at s = j/√(LC) (series resonance)
 *   Pole at s = 0 (opens at DC due to capacitor)
 *
 * For a general lossless Z(s), poles and zeros are on the jω axis
 * and strictly interlace.
 *
 * @param l Inductance (H)
 * @param c Capacitance (F)
 * @param is_parallel 1 for parallel LC, 0 for series LC
 * @param pole_freq Output: pole frequency in rad/s (0 if at origin or ∞)
 * @param zero_freq Output: zero frequency in rad/s
 */
void foster_pole_zero(double l, double c, int is_parallel,
                       double *pole_freq, double *zero_freq);

/**
 * @brief Synthesize Foster-I canonical form.
 *
 * Foster-I form: series connection of parallel LC resonators.
 *   Z(s) = Σ [Ki*s / (s² + ωi²)] + K∞*s + K0/s
 *
 * The K∞ term corresponds to a series inductor.
 * The K0/s term corresponds to a series capacitor.
 * Each Ki term corresponds to a parallel LC in series.
 *
 * @param residues Array of residue values Ki
 * @param pole_freqs Array of pole frequencies ωi
 * @param n_resonators Number of LC resonators
 * @param k_inf Coeff of s term (series L)
 * @param k0 Coeff of 1/s term (series C⁻¹)
 * @param zin_out Output: synthesized Z matrix at ω (element for a single frequency)
 * @param omega Frequency at which to evaluate
 */
void synthesize_foster_one(const double *residues, const double *pole_freqs,
                            int n_resonators, double k_inf, double k0,
                            complex_t *zin_out, double omega);

/**
 * @brief Synthesize Foster-II canonical form.
 *
 * Foster-II form: parallel connection of series LC resonators.
 *   Y(s) = Σ [Ki*s / (s² + ωi²)] + K∞*s + K0/s
 *
 * The dual of Foster-I, synthesizing admittance instead of impedance.
 *
 * @param residues Array of residue values Ki
 * @param pole_freqs Array of pole frequencies ωi
 * @param n_resonators Number of LC resonators
 * @param k_inf Coeff of s term (parallel C)
 * @param k0 Coeff of 1/s term (parallel L⁻¹)
 * @param yin_out Output: synthesized Y at ω
 * @param omega Frequency at which to evaluate
 */
void synthesize_foster_two(const double *residues, const double *pole_freqs,
                            int n_resonators, double k_inf, double k0,
                            complex_t *yin_out, double omega);

/* ============================================================================
 * L6: Darlington Synthesis
 * ============================================================================ */

/**
 * @brief Perform one step of Darlington synthesis.
 *
 * Darlington's theorem states that any positive-real impedance function
 * can be realized as a lossless two-port network terminated in a 1Ω
 * resistor.
 *
 * This function extracts one reactive element (series L or shunt C)
 * from the impedance function, reducing the problem by one degree.
 *
 * @param zin Input impedance to synthesize (Ω, complex at frequency ω)
 * @param omega Angular frequency (rad/s)
 * @param extracted_type 1 for shunt C extraction, 0 for series L extraction
 * @param element_value Output: extracted L (H) or C (F)
 * @param z_remaining Output: remaining impedance after extraction
 * @return 0 on success, -1 if extraction is not possible
 *
 * Course: Illinois ECE 451 — Darlington synthesis
 * Ref: Darlington, "Synthesis of Reactance 4-Poles", 1939
 */
int synthesize_darlington_step(complex_t zin, double omega,
                                int extracted_type,
                                double *element_value,
                                complex_t *z_remaining);

/**
 * @brief Synthesize a maximally-flat (Butterworth) termination.
 *
 * Darlington synthesis applied to a normalized low-pass prototype.
 * For N = 2 (2nd order), Zin(s) = (s² + sqrt(2)*s + 1) / (s² + sqrt(2)*s + 1).
 *
 * This function computes the termination impedance for the prototype
 * at a given frequency.
 *
 * @param order Filter order
 * @param omega Normalized frequency (1 = cutoff)
 * @return Input impedance of Butterworth termination
 */
complex_t synthesize_butterworth_z(double omega, int order);

#endif /* NETWORK_SYNTHESIS_H */
