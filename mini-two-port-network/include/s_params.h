/**
 * @file s_params.h
 * @brief S-Parameter (Scattering Parameter) Analysis
 *
 * S-parameters describe a two-port network in terms of incident (a) and
 * reflected (b) power waves at each port, normalized to reference impedance Z0:
 *
 *   b1 = s11*a1 + s12*a2
 *   b2 = s21*a1 + s22*a2
 *
 * Power wave definitions:
 *   a_n = (Vn + Z0*In) / (2*sqrt(Z0))   (incident wave at port n)
 *   b_n = (Vn - Z0*In) / (2*sqrt(Z0))   (reflected wave at port n)
 *
 * The reference impedance Z0 is typically 50 Ω (RF/microwave standard),
 * 75 Ω (video/CATV standard), or system-specific.
 *
 * All S-parameters are dimensionless complex numbers. Their magnitudes are ≤1
 * for passive networks (conservation of energy).
 *
 * Physical meaning:
 *   s11 = input reflection coefficient (when port 2 is matched: a2=0)
 *   s12 = reverse transmission (isolation)
 *   s21 = forward transmission (gain/loss)
 *   s22 = output reflection coefficient (when port 1 is matched: a1=0)
 *
 * |s11| and |s22| determine return loss and VSWR.
 * |s21| determines gain (amplifier) or loss (attenuator, filter).
 * |s12| determines reverse isolation.
 *
 * S-parameters are the ONLY parameters measurable at microwave frequencies,
 * where open/short circuit conditions are impractical (parasitic effects).
 * Vector Network Analyzers (VNAs) directly measure S-parameters.
 *
 * Reference: Kurokawa, "Power Waves and the Scattering Matrix", IEEE MTT, 1965
 * Course: ETH 227-0455 — Scattering parameters
 *         Georgia Tech ECE 6350 — Microwave measurements
 *         TU Munich HF Engineering — Streuparameter
 */

#ifndef S_PARAMS_H
#define S_PARAMS_H

#include "two_port.h"

/* ============================================================================
 * L1: S-parameter Creation Functions
 * ============================================================================ */

/**
 * @brief Create S-parameter matrix from four complex elements.
 *
 * @param s11 Input reflection coefficient (port 2 matched)
 * @param s12 Reverse transmission coefficient
 * @param s21 Forward transmission coefficient
 * @param s22 Output reflection coefficient (port 1 matched)
 * @return S-parameter matrix
 */
matrix2x2_t sparams_create(complex_t s11, complex_t s12,
                            complex_t s21, complex_t s22);

/**
 * @brief Create S-parameters for an ideal series impedance in a Z0 system.
 *
 * For series impedance Z: s11 = Z/(2*Z0 + Z), s21 = 2*Z0/(2*Z0 + Z)
 *
 *   s11 = s22 = Z / (Z + 2*Z0)
 *   s12 = s21 = 2*Z0 / (Z + 2*Z0)
 *
 * For Z = 0 (short): s11 = 0 (no reflection), s21 = 1 (full transmission)
 * Wait — for a series short across the line, actually:
 *   Z = 0 → s11 = 0/(2Z0) = 0, s21 = 2Z0/(2Z0) = 1 (no effect, just a wire)
 * For Z = ∞ (open): s11 = 1 (total reflection), s21 = 0 (no transmission)
 *
 * @param z Series impedance (Ω)
 * @param z0 Reference impedance (Ω)
 * @return S-parameter matrix
 */
matrix2x2_t sparams_series_z(complex_t z, double z0);

/**
 * @brief Create S-parameters for an ideal shunt admittance in a Z0 system.
 *
 * For shunt admittance Y: Y_total = Y + 2/Z0 (at the node, two Z0 in parallel)
 *
 *   s11 = s22 = -Z0*Y / (2 + Z0*Y) = -Y / (2/Z0 + Y)
 *   s12 = s21 = 2 / (2 + Z0*Y) = 2*Z0/(2*Z0 + Z0²Y)
 *
 * Corrected formula:
 *   s11 = s22 = -Y*Z0 / (2 + Y*Z0)
 *   s12 = s21 = 2 / (2 + Y*Z0)
 *
 * @param y Shunt admittance (S)
 * @param z0 Reference impedance (Ω)
 * @return S-parameter matrix
 */
matrix2x2_t sparams_shunt_y(complex_t y, double z0);

/**
 * @brief Create S-parameters for a transmission line section.
 *
 * For a lossless TL of length l, phase constant β:
 *
 *   s11 = s22 = 0                  (matched — no reflection)
 *   s12 = s21 = e^(-jβl)           (pure phase shift)
 *
 * For a lossy line with attenuation α and phase β:
 *
 *   s11 = s22 = 0                  (still matched)
 *   s12 = s21 = e^(-αl) * e^(-jβl)  (attenuation + phase shift)
 *
 * @param z0 Line characteristic impedance (must equal system Z0 for match)
 * @param alpha Attenuation constant (Np/m)
 * @param beta Phase constant (rad/m)
 * @param length Line length (m)
 * @param z0_sys System reference impedance (Ω)
 * @return S-parameter matrix
 *
 * Course: Georgia Tech ECE 6350 — Transmission line S-parameters
 */
matrix2x2_t sparams_transmission_line(double z0, double alpha,
                                       double beta, double length,
                                       double z0_sys);

/**
 * @brief Create S-parameters for a discrete attenuator in a Z0 system.
 *
 * For an attenuator with attenuation ATT (dB):
 *   |s21| = 10^(-ATT/20)            (voltage ratio)
 *   |s11| = 0 (ideally matched)     (well-matched attenuator)
 *
 * @param att_db Attenuation in dB (positive number)
 * @return S-parameter matrix for ideal attenuator
 */
matrix2x2_t sparams_attenuator(double att_db);

/**
 * @brief Create S-parameters for an ideal amplifier.
 *
 * Ideal model:
 *   s11 = 0         (perfect input match)
 *   s12 = 0         (infinite reverse isolation)
 *   s21 = gain      (forward gain, magnitude > 1)
 *   s22 = 0         (perfect output match)
 *
 * Realistic model includes finite reflection coefficients.
 *
 * @param gain Forward voltage gain (complex, |gain| > 1 for amplifier)
 * @param s11_input Input reflection (e.g., 0.1 for -20 dB return loss)
 * @param s22_input Output reflection
 * @param s12_input Reverse transmission (e.g., 0.01 for -40 dB isolation)
 * @return S-parameter matrix
 */
matrix2x2_t sparams_amplifier(complex_t gain, complex_t s11_input,
                               complex_t s22_input, complex_t s12_input);

/* ============================================================================
 * L3: S-parameter Analysis Functions
 * ============================================================================ */

/**
 * @brief Compute the forward gain in dB from S-parameters.
 *
 * G_db = 20 * log10(|s21|)
 *
 * Standard convention in RF engineering. Negative dB = attenuation.
 * |s21|² = power gain into matched load.
 *
 * @param s S-parameter matrix
 * @return Forward gain in dB
 */
double sparams_s21_db(matrix2x2_t s);

/**
 * @brief Compute the input return loss in dB.
 *
 * RL_in = -20 * log10(|s11|)
 * Typical specs: > 10 dB (acceptable), > 15 dB (good), > 20 dB (excellent)
 *
 * A return loss of ∞ means s11 = 0 (perfect match).
 *
 * @param s S-parameter matrix
 * @return Input return loss in dB
 */
double sparams_input_return_loss_db(matrix2x2_t s);

/**
 * @brief Compute the output return loss in dB.
 *
 * RL_out = -20 * log10(|s22|)
 *
 * @param s S-parameter matrix
 * @return Output return loss in dB
 */
double sparams_output_return_loss_db(matrix2x2_t s);

/**
 * @brief Compute the reverse isolation in dB.
 *
 * ISO = -20 * log10(|s12|)
 *
 * Reverse isolation quantifies how much signal leaks from output to input.
 * Important for amplifier stability and preventing local oscillator leakage.
 *
 * @param s S-parameter matrix
 * @return Reverse isolation in dB (higher = better)
 */
double sparams_isolation_db(matrix2x2_t s);

/**
 * @brief Compute the maximum stable gain (MSG) from S-parameters.
 *
 * MSG = |s21| / |s12|
 *
 * This is the maximum gain achievable while maintaining stability when
 * K < 1 (conditionally stable). MSG is a widely-used figure of merit
 * for active devices.
 *
 * @param s S-parameter matrix
 * @return MSG (linear, not dB)
 */
double sparams_msg(matrix2x2_t s);

/**
 * @brief Compute the maximum available gain (MAG) for K > 1.
 *
 * MAG = |s21/s12| * (K - sqrt(K² - 1))
 *
 * This is the maximum gain achievable with simultaneous conjugate matching
 * at both ports. Valid only when the device is unconditionally stable (K > 1).
 *
 * @param s S-parameter matrix
 * @return MAG (linear, not dB), or -1.0 if K <= 1
 */
double sparams_mag(matrix2x2_t s);

/**
 * @brief Compute the input reflection coefficient with arbitrary load.
 *
 * Γin = s11 + (s12 * s21 * ΓL) / (1 - s22 * ΓL)
 *
 * This is the fundamental formula relating input reflection to load.
 * For a well-matched device (small s11, s22): Γin ≈ s11 + s12*s21*ΓL.
 * The second term represents the "load-pull" effect.
 *
 * @param s S-parameter matrix
 * @param gl Load reflection coefficient
 * @return Input reflection coefficient
 *
 * Course: ETH 227-0455 — Reflection coefficient transformation
 */
complex_t sparams_gamma_in(matrix2x2_t s, complex_t gl);

/**
 * @brief Compute the output reflection coefficient with arbitrary source.
 *
 * Γout = s22 + (s12 * s21 * ΓS) / (1 - s11 * ΓS)
 *
 * @param s S-parameter matrix
 * @param gs Source reflection coefficient
 * @return Output reflection coefficient
 */
complex_t sparams_gamma_out(matrix2x2_t s, complex_t gs);

/**
 * @brief Compute the Rollett stability factor K.
 *
 * K = (1 - |s11|² - |s22|² + |Δ|²) / (2 * |s12 * s21|)
 * where Δ = s11*s22 - s12*s21
 *
 * Unconditional stability requires K > 1 AND |Δ| < 1.
 *
 * This is the most important stability metric in RF amplifier design.
 *
 * @param s S-parameter matrix
 * @return Rollett K-factor (dimensionless)
 *
 * Reference: Rollett, IRE Trans. Circuit Theory, 1962
 * Course: Stanford EE359 — RF amplifier stability
 */
double sparams_rollett_k(matrix2x2_t s);

/**
 * @brief Compute the auxiliary stability condition |Δ|.
 *
 * |Δ| = |s11*s22 - s12*s21|
 *
 * For unconditional stability: K > 1 AND |Δ| < 1 (Edwards-Sinsky criteria).
 *
 * @param s S-parameter matrix
 * @return |Δ| (determinant magnitude)
 */
double sparams_delta_mag(matrix2x2_t s);

/**
 * @brief Compute the μ stability factor (Edwards-Sinsky, 1992).
 *
 * μ = (1 - |s11|²) / (|s22 - Δ*s11*| + |s12*s21|)
 *
 * A single number: μ > 1 means unconditional stability.
 * μ has the advantage over K of being a single criterion and having
 * physical meaning (distance from center of Smith chart to instability).
 *
 * @param s S-parameter matrix
 * @return μ-factor (> 1 = unconditionally stable)
 *
 * Reference: Edwards & Sinsky, IEEE MTT, 1992
 * Course: Georgia Tech ECE 6350 — μ stability
 */
double sparams_mu_factor(matrix2x2_t s);

/**
 * @brief Compute the source stability circle parameters.
 *
 * Center: Cs = (s11 - Δ*s22*)* / (|s11|² - |Δ|²)
 * Radius: Rs = |s12*s21| / | |s11|² - |Δ|² |
 *
 * The source stability circle divides the ΓS Smith chart into stable
 * and potentially unstable regions. Used in amplifier design to choose
 * safe source terminations.
 *
 * @param s S-parameter matrix
 * @param cs Output: center of source stability circle
 * @param rs Output: radius of source stability circle
 *
 * Course: Stanford EE359 — Stability circles
 * Ref: Gonzalez, "Microwave Transistor Amplifiers", Ch. 3
 */
void sparams_source_stability_circle(matrix2x2_t s, complex_t *cs, double *rs);

/**
 * @brief Compute the load stability circle parameters.
 *
 * Center: CL = (s22 - Δ*s11*)* / (|s22|² - |Δ|²)
 * Radius: RL = |s12*s21| / | |s22|² - |Δ|² |
 *
 * @param s S-parameter matrix
 * @param cl Output: center of load stability circle
 * @param rl Output: radius of load stability circle
 */
void sparams_load_stability_circle(matrix2x2_t s, complex_t *cl, double *rl);

/**
 * @brief Compute the bilateral conjugate match at input and output.
 *
 * For an unconditionally stable two-port, the simultaneous conjugate
 * match gives maximum gain. The optimal source and load reflection
 * coefficients are:
 *
 * ΓMS = (B1 ± sqrt(B1² - 4*|C1|²)) / (2*C1)
 * where B1 = 1 + |s11|² - |s22|² - |Δ|²
 *       C1 = s11 - Δ*s22*
 *
 * ΓML = (B2 ± sqrt(B2² - 4*|C2|²)) / (2*C2)
 * where B2 = 1 + |s22|² - |s11|² - |Δ|²
 *       C2 = s22 - Δ*s11*
 *
 * The minus sign is used for the physically realizable solution.
 *
 * @param s S-parameter matrix (must be unconditionally stable)
 * @param gms Output: optimal source reflection coefficient
 * @param gml Output: optimal load reflection coefficient
 * @return 0 on success, -1 if not unconditionally stable
 *
 * Course: Stanford EE359 — Simultaneous conjugate match
 */
int sparams_conjugate_match(matrix2x2_t s, complex_t *gms, complex_t *gml);

/**
 * @brief Compute noise figure of a two-port given source reflection ΓS.
 *
 * NF = NFmin + 4*Rn/Z0 * |ΓS - Γopt|² / ((1-|ΓS|²) * |1+Γopt|²)
 *
 * where NFmin = minimum noise figure, Rn = noise resistance, Γopt = optimum source reflection for noise.
 *
 * This is the fundamental noise parameter equation for two-ports.
 *
 * @param nfmin Minimum noise figure (linear, not dB), e.g., 1.58 for 2 dB
 * @param rn Normalized noise resistance Rn/Z0
 * @param gopt Optimum source reflection coefficient for minimum noise
 * @param gs Actual source reflection coefficient
 * @return Noise figure (linear, not dB)
 *
 * Course: Stanford EE359 — Low-noise amplifier design
 * Ref: Haus, "IRE Standards on Methods of Measuring Noise", 1960
 */
double sparams_noise_figure(double nfmin, double rn,
                             complex_t gopt, complex_t gs);

/* ============================================================================
 * L2: S-parameter to other parameter conversion helpers
 * ============================================================================ */

/**
 * @brief Convert S-parameters to impedance at port 1 with port 2 terminated.
 *
 * Zin = Z0 * (1 + Γin) / (1 - Γin)
 *
 * @param s S-parameter matrix
 * @param z0 Reference impedance (Ω)
 * @param gl Load reflection coefficient at port 2
 * @return Input impedance (Ω)
 */
complex_t sparams_input_impedance(matrix2x2_t s, double z0, complex_t gl);

/**
 * @brief Compute the transmission phase through a two-port.
 *
 * φ = arg(s21) in degrees. Used for phase shifter and delay line analysis.
 *
 * @param s S-parameter matrix
 * @return Transmission phase in degrees [-180, 180]
 */
double sparams_transmission_phase_deg(matrix2x2_t s);

/**
 * @brief Compute group delay from phase slope.
 *
 * τg = -dφ/dω ≈ -(Δφ)/(Δω)
 *
 * Group delay measures the time delay of the signal envelope.
 * Constant group delay = linear phase = no dispersion.
 *
 * @param s1 S-parameters at frequency ω1
 * @param s2 S-parameters at frequency ω2
 * @param omega1 Angular frequency 1 (rad/s)
 * @param omega2 Angular frequency 2 (rad/s)
 * @return Group delay in seconds
 *
 * Course: Michigan EECS 411 — Group delay in filters
 */
double sparams_group_delay(matrix2x2_t s1, matrix2x2_t s2,
                            double omega1, double omega2);

#endif /* S_PARAMS_H */
