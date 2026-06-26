/**
 * @file network_function.h
 * @brief Network functions in the frequency domain
 *
 * Network functions relate voltages and currents at the ports of
 * a linear network. In the frequency domain, these are rational
 * functions of the complex frequency variable s = σ + jω.
 *
 * Types of network functions:
 * - Driving-point impedance Z(s) = V(s)/I(s) at same port
 * - Driving-point admittance Y(s) = I(s)/V(s) = 1/Z(s)
 * - Transfer impedance Z₂₁(s) = V₂(s)/I₁(s) (I₂ = 0)
 * - Transfer admittance Y₂₁(s) = I₂(s)/V₁(s) (V₂ = 0)
 * - Voltage transfer ratio H(s) = V₂(s)/V₁(s)
 * - Current transfer ratio = I₂(s)/I₁(s)
 *
 * Reference:
 * - Van Valkenburg, "Network Analysis" (1974), Ch. 12-14
 * - Guillemin, "Synthesis of Passive Networks" (1957)
 * - Hayt, Kemmerly & Durbin (2019), Ch. 16
 *
 * Course: Berkeley EE16B, MIT 6.003, Illinois ECE 451
 */

#ifndef NETWORK_FUNCTION_H
#define NETWORK_FUNCTION_H

#include <stddef.h>
#include <complex.h>
#include "frequency_response.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1 DEFINITIONS: Network function types
 * ============================================================================ */

/**
 * @brief Two-port network parameter sets.
 *
 * Six common parameter sets characterize linear two-port networks:
 *
 * Z-parameters (impedance, open-circuit):
 *   V₁ = Z₁₁I₁ + Z₁₂I₂
 *   V₂ = Z₂₁I₁ + Z₂₂I₂
 *   Measured with open-circuit terminations.
 *
 * Y-parameters (admittance, short-circuit):
 *   I₁ = Y₁₁V₁ + Y₁₂V₂
 *   I₂ = Y₂₁V₁ + Y₂₂V₂
 *   Measured with short-circuit terminations.
 *
 * H-parameters (hybrid):
 *   V₁ = h₁₁I₁ + h₁₂V₂
 *   I₂ = h₂₁I₁ + h₂₂V₂
 *   Popular for BJT small-signal modeling.
 *
 * G-parameters (inverse hybrid):
 *   I₁ = g₁₁V₁ + g₁₂I₂
 *   V₂ = g₂₁V₁ + g₂₂I₂
 *
 * ABCD-parameters (transmission, chain):
 *   V₁ = A·V₂ - B·I₂
 *   I₁ = C·V₂ - D·I₂
 *   Convenient for cascading: ABCD_total = ABCD₁ × ABCD₂.
 *
 * S-parameters (scattering):
 *   b₁ = S₁₁a₁ + S₁₂a₂
 *   b₂ = S₂₁a₁ + S₂₂a₂
 *   Preferred at microwave frequencies (traveling waves).
 *
 * L2 Concept: The choice of parameter set depends on the application:
 * - Z/Y for circuit analysis and synthesis
 * - h for BJTs (manufacturer datasheets)
 * - ABCD for cascaded networks
 * - S for RF/microwave (network analyzer measurements)
 */
typedef enum {
    NETWORK_PARAM_Z,     /**< Open-circuit impedance parameters */
    NETWORK_PARAM_Y,     /**< Short-circuit admittance parameters */
    NETWORK_PARAM_H,     /**< Hybrid parameters */
    NETWORK_PARAM_G,     /**< Inverse hybrid parameters */
    NETWORK_PARAM_ABCD,  /**< Transmission (chain) parameters */
    NETWORK_PARAM_S      /**< Scattering parameters */
} network_param_type_t;

/**
 * @brief Frequency-dependent two-port network parameters.
 *
 * Each parameter is a complex function of frequency.
 */
typedef struct {
    double _Complex z11, z12, z21, z22;  /**< Z-parameters (Ω) */
    double _Complex y11, y12, y21, y22;  /**< Y-parameters (S) */
    double _Complex h11, h12, h21, h22;  /**< H-parameters (mixed) */
    double _Complex s11, s12, s21, s22;  /**< S-parameters (unitless) */
    double frequency_hz;                 /**< Frequency (Hz) */
} network_params_t;

/* ============================================================================
 * NETWORK FUNCTION SYNTHESIS
 * ============================================================================ */

/**
 * @brief Compute the driving-point impedance of common RLC networks.
 *
 * For a series RLC: Z(s) = R + sL + 1/(sC)
 * For a parallel RLC: Z(s) = 1/(1/R + 1/(sL) + sC)
 *
 * @param R       Resistance (Ω)
 * @param L       Inductance (H)
 * @param C       Capacitance (F)
 * @param is_par  0=series, 1=parallel
 * @return        Impedance transfer function Z(s)
 */
tf_polynomial_t *network_impedance_rlc(double R, double L, double C,
                                         int is_par);

/**
 * @brief Compute the voltage transfer function of a voltage divider.
 *
 * H(s) = Z₂(s)/(Z₁(s) + Z₂(s))
 *
 * where Z₁ and Z₂ are arbitrary impedances.
 *
 * Used for:
 * - Compensated attenuator (oscilloscope probe): Z₁ = R₁||C₁, Z₂ = R₂||C₂
 *   10× probe: R₁=9MΩ, R₂=1MΩ, compensation: R₁C₁ = R₂C₂
 *
 * - Lead-lag compensator: Z₁ = R₁, Z₂ = R₂ + 1/(sC)
 * - Lag-lead compensator: Z₁ = R₁||(1/(sC₁)), Z₂ = R₂
 *
 * @param Z1  Impedance Z₁(s)
 * @param Z2  Impedance Z₂(s)
 * @return    Voltage divider transfer function H(s) = Z₂/(Z₁+Z₂)
 */
tf_polynomial_t *network_voltage_divider_tf(const tf_polynomial_t *Z1,
                                              const tf_polynomial_t *Z2);

/**
 * @brief Compute impedance of common single-port elements.
 *
 * L: Z(s) = sL
 * C: Z(s) = 1/(sC)
 * R: Z(s) = R
 * R-L series: Z(s) = R + sL
 * R-C series: Z(s) = R + 1/(sC)
 * R||L: Z(s) = (R·sL)/(R + sL)
 * R||C: Z(s) = R/(1 + sRC)
 *
 * @param element_type 0=L, 1=C, 2=R, 3=R+L, 4=R+C, 5=R||L, 6=R||C
 * @param val1         First element value (R for composite, otherwise L/C)
 * @param val2         Second element value (0 if single element)
 * @return             Impedance Z(s), or NULL
 */
tf_polynomial_t *network_element_impedance(int element_type,
                                             double val1, double val2);

/* ============================================================================
 * PARAMETER CONVERSION
 * ============================================================================ */

/**
 * @brief Convert between two-port parameter sets.
 *
 * Standard conversion formulas (for reciprocal networks, Z₁₂ = Z₂₁):
 *
 * Z → Y: Y = Z^{-1}
 *   y₁₁ = z₂₂/Δz, y₁₂ = -z₁₂/Δz
 *   y₂₁ = -z₂₁/Δz, y₂₂ = z₁₁/Δz
 *   where Δz = z₁₁z₂₂ - z₁₂z₂₁
 *
 * Z → ABCD:
 *   A = z₁₁/z₂₁, B = Δz/z₂₁
 *   C = 1/z₂₁, D = z₂₂/z₂₁
 *
 * Z → S (with reference impedance Z₀):
 *   S = (Z + Z₀I)^{-1}(Z - Z₀I)
 *
 * Y → Z: Z = Y^{-1}
 * H → Z, Z → H: derived from defining equations
 * ABCD → Z, ABCD → Y, ABCD → S
 *
 * @param src        Source parameters
 * @param src_type   Source parameter type
 * @param dst        Destination parameters (output)
 * @param dst_type   Destination parameter type
 * @param Z0         Reference impedance for S-parameters (Ω, typically 50)
 * @return           0 on success, -1 on error
 *
 * L5 Algorithm: Two-port parameter conversion is essential for
 * circuit analysis because different parameter sets are convenient
 * for different interconnections:
 * - Series connection: Z_total = Z₁ + Z₂
 * - Parallel connection: Y_total = Y₁ + Y₂
 * - Cascade connection: ABCD_total = ABCD₁ × ABCD₂
 * - Series-parallel: H_total = H₁ + H₂
 * - Parallel-series: G_total = G₁ + G₂
 *
 * Reference: Pozar (2012), Appendix D
 */
int network_param_convert(const double _Complex src[4],
                           network_param_type_t src_type,
                           double _Complex dst[4],
                           network_param_type_t dst_type,
                           double Z0);

/* ============================================================================
 * S-PARAMETER COMPUTATION (L6 Canonical Problem)
 * ============================================================================ */

/**
 * @brief Compute S-parameters from a two-port impedance description.
 *
 * Given a two-port network with impedance matrix Z at frequency f:
 *   S = (Z + Z₀I)^{-1} · (Z - Z₀I)
 *
 * where Z₀ is the reference impedance (typically 50 Ω for RF,
 * 75 Ω for video/CATV).
 *
 * For a simple series impedance Z_s:
 *   S₁₁ = Z_s/(Z_s + 2Z₀)
 *   S₂₁ = 2Z₀/(Z_s + 2Z₀)
 *   S₁₂ = S₂₁ (reciprocal)
 *   S₂₂ = S₁₁
 *
 * For a simple shunt admittance Y_p:
 *   S₁₁ = -Y_p/(Y_p + 2Y₀) where Y₀ = 1/Z₀
 *   S₂₁ = 2Y₀/(Y_p + 2Y₀)
 *
 * @param z11, z12, z21, z22  Z-parameters at frequency f (Ω)
 * @param Z0                  Reference impedance (Ω)
 * @return                    S-parameters at frequency f
 *
 * L6 Canonical Problem: Computing S-parameters from circuit theory
 * bridges the gap between lumped-element circuit analysis and
 * microwave network analysis. This is essential for:
 * - Matching network design
 * - Filter design verification
 * - Amplifier stability analysis (K-factor, μ-factor)
 *
 * Course: ETH 227-0455, Michigan EECS 411, Georgia Tech ECE 6350
 */
s_params_2port_t network_s_params_from_z(double _Complex z11,
                                           double _Complex z12,
                                           double _Complex z21,
                                           double _Complex z22,
                                           double Z0, double frequency_hz);

/**
 * @brief Compute the Rollett stability factor K from S-parameters.
 *
 * K = (1 - |S₁₁|² - |S₂₂|² + |Δ|²) / (2|S₁₂S₂₁|)
 * where Δ = S₁₁S₂₂ - S₁₂S₂₁
 *
 * Unconditional stability requires:
 *   K > 1 AND |Δ| < 1 (Rollett, 1962)
 *
 * Alternative: μ-factor (Edwards & Sinsky, 1992):
 *   μ = (1 - |S₁₁|²)/(|S₂₂ - Δ·S₁₁*| + |S₁₂S₂₁|)
 *   μ > 1 = unconditional stability (single test)
 *
 * @param s      S-parameters at a single frequency
 * @return       K-factor. K>1 + |Δ|<1 → unconditionally stable
 *
 * L8 Advanced: Amplifier stability analysis using S-parameters
 * is critical for RF/microwave design. The K-factor must be >1
 * across the entire frequency range (not just at the operating
 * frequency) for unconditional stability.
 *
 * Reference: Rollett, "Stability and Power-Gain Invariants of
 * Linear Twoports" (1962), IRE Trans. Circuit Theory
 * Course: ETH 227-0455, Georgia Tech ECE 6350
 */
double network_rollett_k(const s_params_2port_t *s);

/**
 * @brief Compute maximum available gain (MAG) from S-parameters.
 *
 * For unconditionally stable device (K > 1):
 *   MAG = |S₂₁/S₁₂| · (K - √(K² - 1))  (linear)
 *   MAG_dB = 10·log₁₀(MAG)
 *
 * For potentially unstable device (K < 1):
 *   MSG = |S₂₁|/|S₁₂|  (maximum stable gain)
 *
 * @param s   S-parameters at a single frequency
 * @return    MAG if K>1, MSG if K<1
 *
 * L6 Canonical Problem: Computing the maximum gain available from
 * a transistor at a given frequency, assuming optimal simultaneous
 * conjugate matching at input and output.
 */
double network_max_gain(const s_params_2port_t *s);

/* ============================================================================
 * POSITIVE REAL FUNCTION CHECK
 * ============================================================================ */

/**
 * @brief Check if a rational function Z(s) is positive real (PR).
 *
 * Necessary and sufficient conditions for Z(s) = N(s)/D(s) to be PR
 * (Brune, 1931):
 *
 * 1. Z(s) is real when s is real (coefficients are real — always
 *    satisfied for physical circuits)
 *
 * 2. Re{Z(jω)} ≥ 0 for all ω (real part of impedance is non-negative
 *    on the imaginary axis)
 *
 * 3. Poles on the jω-axis must be simple with positive real residues
 *
 * 4. The difference between numerator and denominator degrees ≤ 1
 *
 * A function is PR iff it can be realized as a driving-point impedance
 * of a network containing only passive R, L, C elements (and ideal
 * transformers).
 *
 * @param Z     Impedance transfer function
 * @param freq  Frequency array for checking Re{Z(jω)} ≥ 0 (Hz)
 * @param n     Number of frequency points
 * @return      1 = positive real, 0 = not positive real, -1 = error
 *
 * L4: Brune's theorem is a fundamental result in network synthesis.
 * It establishes the exact class of functions realizable as passive
 * impedances, providing both a test and a synthesis procedure.
 *
 * Reference: Brune, "Synthesis of a Finite Two-terminal Network
 * whose Driving-point Impedance is a Prescribed Function of
 * Frequency" (1931), J. Math. Phys.
 */
int network_is_positive_real(const tf_polynomial_t *Z,
                               const double *freq, size_t n);

/* ============================================================================
 * INPUT IMPEDANCE OF TERMINATED TWO-PORT (L6 Canonical Problem)
 * ============================================================================ */

/**
 * @brief Compute input impedance of a two-port terminated in Z_L.
 *
 * Z_in = V₁/I₁ = Z₁₁ - (Z₁₂Z₂₁)/(Z₂₂ + Z_L)
 *
 * This is one of the most frequently used formulas in RF circuit
 * design. For a transistor amplifier:
 * - Z_in determines the input matching network design
 * - Z_out = Z₂₂ - (Z₁₂Z₂₁)/(Z₁₁ + Z_s) determines output matching
 *
 * Conjugate match condition for maximum power transfer:
 *   Z_s = Z_in* (source matched to input)
 *   Z_L = Z_out* (load matched to output)
 *
 * @param s_params  S-parameters of the two-port
 * @param ZL        Load impedance (Ω, complex)
 * @param Z0        Reference impedance (Ω)
 * @return          Input impedance Z_in (Ω, complex)
 *
 * L6 Canonical Problem: Determining the input impedance of a
 * loaded two-port is a fundamental building block for:
 * - Amplifier design (matching networks)
 * - Filter design (termination effects)
 * - Oscillator design (negative resistance condition)
 */
double _Complex network_input_impedance(const s_params_2port_t *s_params,
                                           double _Complex ZL, double Z0);

/**
 * @brief Compute the voltage standing wave ratio (VSWR).
 *
 * VSWR = (1 + |Γ|)/(1 - |Γ|)
 *
 * where Γ = reflection coefficient = (Z_L - Z₀)/(Z_L + Z₀)
 * or Γ = S₁₁ for a one-port measurement.
 *
 * VSWR ranges from 1 (perfect match, Γ=0) to ∞ (complete mismatch,
 * |Γ|=1, open or short).
 *
 * Return loss: RL = -20·log₁₀|Γ| dB
 * Good match: VSWR < 1.5 (RL > 14 dB)
 * Acceptable: VSWR < 2.0 (RL > 10 dB)
 *
 * @param gamma  Reflection coefficient Γ (complex)
 * @return       VSWR (≥1), INFINITY if |Γ| ≥ 1
 *
 * L1 Definition: VSWR is the fundamental measure of impedance
 * matching quality in RF systems.
 */
double network_vswr(double _Complex gamma);

#ifdef __cplusplus
}
#endif

#endif /* NETWORK_FUNCTION_H */
