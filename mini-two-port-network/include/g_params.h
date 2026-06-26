/**
 * @file g_params.h
 * @brief G-Parameter (Inverse Hybrid Parameter) Analysis
 *
 * G-parameters are the dual of H-parameters, mixing input admittance with
 * output impedance:
 *
 *   I1 = g11*V1 + g12*I2
 *   V2 = g21*V1 + g22*I2
 *
 * Matrix form: [[I1], [V2]] = [[g11, g12], [g21, g22]] * [[V1], [I2]]
 *
 * Physical interpretation:
 *   g11 = I1/V1 | I2=0  (input admittance, output open) — S
 *   g12 = I1/I2 | V1=0  (reverse current ratio, input shorted) — dimensionless
 *   g21 = V2/V1 | I2=0  (forward voltage gain, output open) — dimensionless
 *   g22 = V2/I2 | V1=0  (output impedance, input shorted) — Ω
 *
 * G-parameters are natural for:
 *   - FET common-source amplifiers (gate voltage in, drain current out)
 *   - Parallel-series connected networks (G_total = G1 + G2)
 *   - Feedback networks with shunt mixing and series sampling
 *
 * Relationship to other parameters: G = H^(-1) for invertible H.
 *
 * Reference: Carson & DeBuda, "Circuit Theory" (1990), Ch. 11
 * Course: Georgia Tech ECE 6350 — Inverse hybrid analysis
 */

#ifndef G_PARAMS_H
#define G_PARAMS_H

#include "two_port.h"

/* ============================================================================
 * L1: G-parameter Creation Functions
 * ============================================================================ */

/**
 * @brief Create G-parameter matrix from four mixed parameters.
 *
 * @param g11 Input admittance with output open (S)
 * @param g12 Reverse current ratio with input shorted (dimensionless)
 * @param g21 Forward voltage gain with output open (dimensionless)
 * @param g22 Output impedance with input shorted (Ω)
 * @return G-parameter matrix
 */
matrix2x2_t gparams_create(complex_t g11, complex_t g12,
                            complex_t g21, complex_t g22);

/**
 * @brief Create G-parameters for a common-source FET.
 *
 * In the G-parameter representation:
 *   g11 = 1/RG + jω*CGD_total (input admittance)
 *   g12 = 0 (reverse current ratio, negligible in simplified model)
 *   g21 = -gm * ro = intrinsic gain (forward voltage gain)
 *   g22 = ro (output resistance)
 *
 * At DC: g11 = 0 (infinite input impedance), g21 = -gm*ro (intrinsic gain).
 * The intrinsic gain gm*ro is a technology parameter, typically 20-100.
 *
 * @param gm Transconductance (S)
 * @param ro Output resistance (Ω)
 * @param rg Gate bias resistance (Ω), typically 1MΩ for discrete
 * @return G-parameter matrix for CS FET
 *
 * Course: Berkeley EE140 — Analog IC Design (common-source analysis)
 * Ref: Razavi, "Design of Analog CMOS Integrated Circuits", Ch. 3
 */
matrix2x2_t gparams_fet_cs(double gm, double ro, double rg);

/**
 * @brief Create G-parameters for a common-gate FET.
 *
 * CG stage (gate common, input at source, output at drain):
 *   g11 = gm (input admittance ≈ gm, low Zin)
 *   g12 = 0
 *   g21 = gm * ro (voltage gain, positive — non-inverting)
 *   g22 = ro + (1+gm*ro)*RS (output impedance, boosted by source degeneration)
 *
 * @param gm Transconductance (S)
 * @param ro Output resistance (Ω)
 * @param rs Source degeneration resistance (Ω, typically 0 for ideal CG)
 * @return G-parameter matrix for CG FET
 */
matrix2x2_t gparams_fet_cg(double gm, double ro, double rs);

/**
 * @brief Create G-parameters for a common-drain (source follower) FET.
 *
 * CD stage (drain common, input at gate, output at source):
 *   g11 = 1/RG + jω*(Cgs+Cgd) (input admittance, very low at DC)
 *   g12 = 0
 *   g21 = gm*rs/(1+gm*rs) (voltage gain, < 1, near unity for large gm*rs)
 *   g22 = 1/gm || rs (output impedance, ~1/gm for large gm*rs)
 *
 * where rs is the source load resistance (bias resistor or current source).
 *
 * @param gm Transconductance (S)
 * @param rg Input bias resistance (Ω)
 * @param rs Source load resistance (Ω)
 * @return G-parameter matrix for CD FET
 *
 * Course: Stanford EE214 — Source follower analysis
 */
matrix2x2_t gparams_fet_cd(double gm, double rg, double rs);

/* ============================================================================
 * L3: G-parameter Analysis Functions
 * ============================================================================ */

/**
 * @brief Compute voltage gain from G-parameters.
 *
 * Av = g21 * ZL / (g22 + ZL)
 *
 * For CS FET with ZL << ro: Av ≈ g21 * ZL / ro = -gm * ZL
 * which is the classic gain formula.
 *
 * @param g G-parameter matrix
 * @param zl Load impedance (Ω)
 * @return Complex voltage gain
 */
complex_t gparams_voltage_gain(matrix2x2_t g, complex_t zl);

/**
 * @brief Compute input admittance from G-parameters.
 *
 * Yin = g11 - (g12 * g21) / (g22 + ZL)
 *
 * For CS FET (g12 ≈ 0): Yin ≈ g11 (very high input impedance).
 * The g12*g21 term represents feedback through the Miller capacitance.
 *
 * @param g G-parameter matrix
 * @param zl Load impedance (Ω)
 * @return Input admittance (S)
 */
complex_t gparams_input_admittance(matrix2x2_t g, complex_t zl);

/**
 * @brief Compute output impedance from G-parameters.
 *
 * Zout = g22 - (g12 * g21) / (g11 + YS)
 *
 * For CS FET with g12=0: Zout = g22 = ro.
 *
 * @param g G-parameter matrix
 * @param ys Source admittance (S) = 1/ZS
 * @return Output impedance (Ω)
 */
complex_t gparams_output_impedance(matrix2x2_t g, complex_t ys);

/**
 * @brief Compute source follower output impedance.
 *
 * For a CD stage: Rout = 1/gm (approximately, for large gm*rs).
 * This is the key advantage of the source follower — low output impedance
 * for driving capacitive or low-impedance loads.
 *
 * The exact formula includes body effect: Rout = 1/(gm + gmb) || rs
 * where gmb = body transconductance (typically 0.1-0.3 × gm).
 *
 * @param gm Transconductance (S)
 * @param gmb Body transconductance (S), 0 if not applicable
 * @param rs Source load resistance (Ω)
 * @return Output impedance (Ω)
 *
 * Course: Berkeley EE140 — Output impedance of source followers
 */
double gparams_source_follower_rout(double gm, double gmb, double rs);

/**
 * @brief Compute the intrinsic gain of a MOSFET from G-parameters.
 *
 * A_int = |g21(open circuit)| = gm * ro
 *
 * The intrinsic gain gm*ro is a fundamental technology parameter:
 *   Long-channel: gm*ro ≈ 2/λ*Vdsat ≈ 50-200
 *   Short-channel: gm*ro ≈ 10-50
 *
 * It sets the maximum single-stage voltage gain achievable.
 *
 * @param g G-parameter matrix
 * @return Intrinsic gain (linear, dimensionless)
 */
double gparams_intrinsic_gain(matrix2x2_t g);

#endif /* G_PARAMS_H */
