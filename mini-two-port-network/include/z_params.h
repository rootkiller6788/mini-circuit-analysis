/**
 * @file z_params.h
 * @brief Z-Parameter (Impedance Parameter) Analysis
 *
 * Z-parameters (open-circuit impedance parameters) are the most fundamental
 * two-port representation, defined by port voltages in terms of port currents:
 *
 *   V1 = z11*I1 + z12*I2
 *   V2 = z21*I1 + z22*I2
 *
 * Matrix form: [V] = [Z] * [I]
 *
 * Physical interpretation:
 *   z11 = V1/I1 | I2=0  (input impedance, output open)
 *   z12 = V1/I2 | I1=0  (reverse transfer impedance, input open)
 *   z21 = V2/I1 | I2=0  (forward transfer impedance, output open)
 *   z22 = V2/I2 | I1=0  (output impedance, input open)
 *
 * All Z-parameters have units of ohms (Ω).
 *
 * Z-parameters are natural for:
 *   - Series-connected networks (Z_total = Z1 + Z2)
 *   - T-equivalent circuit models
 *   - Transformer coupling analysis
 *   - Feedback network analysis (series-series)
 *
 * Reference: Kuo, "Network Analysis and Synthesis" (2006), Ch. 13
 * Course: MIT 6.002 — Equivalent circuits and two-port networks
 */

#ifndef Z_PARAMS_H
#define Z_PARAMS_H

#include "two_port.h"

/* ============================================================================
 * L1: Z-parameter Creation Functions
 * ============================================================================ */

/**
 * @brief Create Z-parameter matrix from four complex impedances.
 *
 * @param z11 Input impedance with output open-circuit (Ω)
 * @param z12 Reverse transfer impedance with input open-circuit (Ω)
 * @param z21 Forward transfer impedance with output open-circuit (Ω)
 * @param z22 Output impedance with input open-circuit (Ω)
 * @return Z-parameter matrix
 */
matrix2x2_t zparams_create(complex_t z11, complex_t z12,
                            complex_t z21, complex_t z22);

/**
 * @brief Create Z-parameters for a series impedance Zs at port 1.
 *
 * Circuit: port1 — Zs — port2 (with common ground)
 * Z = [[Zs, Zs], [Zs, Zs]]
 *
 * This represents a floating impedance between two ports with shared ground.
 * Used for: modeling PCB trace impedance, bond wire inductance.
 *
 * Course: Berkeley EE117 — Transmission line modeling
 */
matrix2x2_t zparams_series_impedance(complex_t zs);

/**
 * @brief Create Z-parameters for a shunt impedance Zp from port to ground.
 *
 * Circuit: port1 — + — port2; shunt Zp from line to ground
 * Z = [[Zp, Zp], [Zp, Zp]]
 *
 * Course: Michigan EECS 411 — Shunt element in microwave circuits
 */
matrix2x2_t zparams_shunt_impedance(complex_t zp);

/**
 * @brief Create Z-parameters for a T-network with three impedances.
 *
 *        Za       Zb
 *   P1 ————┬———————┬——— P2
 *          │       │
 *          Zc
 *          │
 *         GND
 *
 * Z-matrix:
 *   z11 = Za + Zc
 *   z12 = Zc = z21 (reciprocal)
 *   z22 = Zb + Zc
 *
 * The T-network is the natural Z-parameter equivalent circuit.
 * Used for: attenuator pads, impedance matching networks.
 *
 * Course: TU Munich HF Engineering — T-equivalent circuit
 */
matrix2x2_t zparams_t_network(complex_t za, complex_t zb, complex_t zc);

/**
 * @brief Create Z-parameters for a π-network (using Z-parameters).
 *
 *        Yb
 *   P1 —┬———┬——— P2
 *       Ya  Yc
 *       │   │
 *      GND GND
 *
 * The π-network is more naturally described by Y-parameters,
 * but Z-parameters exist:
 *   z11 = (Yb+Yc) / (Ya*Yb + Yb*Yc + Yc*Ya)
 *   z12 = Yb / (Ya*Yb + Yb*Yc + Yc*Ya)
 *   z21 = Yb / (Ya*Yb + Yb*Yc + Yc*Ya)  (= z12 for reciprocal)
 *   z22 = (Ya+Yb) / (Ya*Yb + Yb*Yc + Yc*Ya)
 *
 * @param ya Admittance Ya (first shunt, Siemens)
 * @param yb Admittance Yb (series, Siemens)
 * @param yc Admittance Yc (second shunt, Siemens)
 * @return Z-parameter matrix
 */
matrix2x2_t zparams_pi_network(complex_t ya, complex_t yb, complex_t yc);

/**
 * @brief Create Z-parameters for an ideal transformer.
 *
 * An ideal transformer with turns ratio n:1 has:
 *   V1 = (1/n) * V2
 *   I1 = -n * I2
 *
 * However, ideal transformers cannot be represented by Z-parameters alone
 * (they are not impedance-admittance describable). This function returns
 * Z-parameters for a practical transformer model with magnetizing inductance Lm,
 * leakage inductance Ll, and winding resistance Rw.
 *
 * Model: Primary: Rw1 + jωLl1, Secondary: Rw2 + jωLl2, Mutual: jωM
 *   z11 = Rw1 + jω(Ll1 + Lm)
 *   z12 = jωM = jωk*sqrt(Lp*Ls)  (mutual coupling)
 *   z21 = z12 (reciprocal)
 *   z22 = Rw2 + jω(Ll2 + n²*Lm)
 *
 * @param rw1 Primary winding resistance (Ω)
 * @param rw2 Secondary winding resistance (Ω)
 * @param l_lk1 Primary leakage inductance (H)
 * @param l_lk2 Secondary leakage inductance (H)
 * @param lm Magnetizing inductance (H)
 * @param k Coupling coefficient (0 ≤ k ≤ 1)
 * @param n Turns ratio (Ns/Np)
 * @param omega Angular frequency (rad/s)
 * @return Z-parameter matrix at frequency ω
 *
 * Course: MIT 6.002 — Transformer equivalent circuits
 * Ref: Erickson & Maksimovic, "Fundamentals of Power Electronics", Ch. 13
 */
matrix2x2_t zparams_transformer(double rw1, double rw2,
                                 double l_lk1, double l_lk2,
                                 double lm, double k, double n,
                                 double omega);

/**
 * @brief Create Z-parameters for a common-emitter BJT at low frequencies.
 *
 * Simplified hybrid-π model mapped to Z-parameters:
 *   z11 = rπ  (base-emitter input resistance)
 *   z12 = 0   (negligible reverse transmission at low freq)
 *   z21 = -gm * rπ * ro  (forward transimpedance, ~ -β*ro)
 *   z22 = ro   (collector-emitter output resistance)
 *
 * @param rpi Base-emitter small-signal resistance rπ (Ω)
 * @param gm Transconductance (S)
 * @param ro Early-effect output resistance (Ω)
 * @return Z-parameter matrix for CE BJT
 *
 * Course: Berkeley EE105 — BJT small-signal models
 * Ref: Sedra & Smith §7.2
 */
matrix2x2_t zparams_bjt_ce(double rpi, double gm, double ro);

/**
 * @brief Create Z-parameters for a common-source FET at low frequencies.
 *
 * Simplified model:
 *   z11 = ∞ (gate is open-circuit at DC/LF, represented as very large)
 *   z12 = 0 (no reverse transmission in simplified model)
 *   z21 = -gm * ro (forward transimpedance)
 *   z22 = ro (drain-source output resistance)
 *
 * For practical computation, z11 is set to a very large value (e.g., 1e12).
 *
 * @param gm Transconductance (S)
 * @param ro Output resistance (Ω) — from channel-length modulation
 * @return Z-parameter matrix for CS FET
 *
 * Course: Berkeley EE105 — MOSFET small-signal models
 * Ref: Sedra & Smith §6.3
 */
matrix2x2_t zparams_fet_cs(double gm, double ro);

/* ============================================================================
 * L3: Z-parameter Analysis Functions
 * ============================================================================ */

/**
 * @brief Compute the Z-parameter representation of a resistor.
 *
 * For a resistor R, Z = [[R, R], [R, R]] (series model).
 * This is the simplest passive two-port.
 *
 * @param r Resistance in ohms (Ω)
 * @return Z-parameter matrix
 */
matrix2x2_t zparams_resistor(double r);

/**
 * @brief Compute the Z-parameter representation of an inductor at frequency ω.
 *
 * Z = [[jωL, jωL], [jωL, jωL]] (series inductor)
 * Inductive impedance: ZL = jωL, where ω = 2πf.
 *
 * @param l Inductance in henries (H)
 * @param omega Angular frequency in rad/s
 * @return Z-parameter matrix
 *
 * Course: MIT 6.002 — Frequency-domain impedance
 */
matrix2x2_t zparams_inductor(double l, double omega);

/**
 * @brief Compute the Z-parameter representation of a capacitor at frequency ω.
 *
 * Z = [[1/(jωC), 1/(jωC)], [1/(jωC), 1/(jωC)]] (series capacitor)
 * Capacitive impedance: ZC = 1/(jωC) = -j/(ωC).
 *
 * @param c Capacitance in farads (F)
 * @param omega Angular frequency in rad/s
 * @return Z-parameter matrix
 */
matrix2x2_t zparams_capacitor(double c, double omega);

/**
 * @brief Compute Z-parameters for a series RLC branch at frequency ω.
 *
 * Zseries = R + j(ωL - 1/(ωC))
 *
 * Resonance occurs at ω0 = 1/sqrt(LC), where Z is purely real (R).
 * Below resonance: capacitive (negative imag), above: inductive (positive imag).
 *
 * @param r Resistance (Ω)
 * @param l Inductance (H)
 * @param c Capacitance (F)
 * @param omega Angular frequency (rad/s)
 * @return Z-parameter matrix for series RLC
 */
matrix2x2_t zparams_series_rlc(double r, double l, double c, double omega);

/**
 * @brief Compute Z-parameters for a parallel RLC tank at frequency ω.
 *
 * For a parallel RLC: Ztank = 1 / (1/R + 1/(jωL) + jωC)
 *
 * At resonance (ω0 = 1/sqrt(LC)): Ztank = R (maximum impedance, pure real).
 * The quality factor Q = R*sqrt(C/L) determines bandwidth and selectivity.
 *
 * @param r Resistance (Ω)
 * @param l Inductance (H)
 * @param c Capacitance (F)
 * @param omega Angular frequency (rad/s)
 * @return Z-parameter matrix
 */
matrix2x2_t zparams_parallel_rlc(double r, double l, double c, double omega);

/**
 * @brief Compute the open-circuit voltage gain from Z-parameters.
 *
 * Av_oc = V2/V1 | I2=0 = z21 / z11
 *
 * This is the unloaded (open-circuit) voltage gain. For a loaded amplifier,
 * use two_port_voltage_gain_z() which includes the load effect.
 *
 * @param z Z-parameter matrix
 * @return Open-circuit voltage gain (complex)
 */
complex_t zparams_open_circuit_gain(matrix2x2_t z);

/**
 * @brief Compute the short-circuit current gain from Z-parameters.
 *
 * Ai_sc = I2/I1 | V2=0
 *
 * With V2 = 0 (short circuit at output): I2 = (-z21/z22) * I1
 * (derived from the Z-parameter equations with V2 = 0)
 *
 * @param z Z-parameter matrix
 * @return Short-circuit current gain (complex)
 */
complex_t zparams_short_circuit_gain(matrix2x2_t z);

/**
 * @brief Extract equivalent T-network elements from Z-parameters.
 *
 * Given Z-parameters of a reciprocal network, extract Za, Zb, Zc:
 *   Za = z11 - z12
 *   Zb = z22 - z12
 *   Zc = z12
 *
 * This is the inverse of zparams_t_network() and demonstrates the
 * equivalence between any reciprocal two-port and a T-network.
 *
 * @param z Z-parameter matrix (must be reciprocal: z12 ≈ z21)
 * @param za Output: first series arm
 * @param zb Output: second series arm
 * @param zc Output: shunt arm
 * @return 0 on success, -1 if non-reciprocal
 *
 * Course: Illinois ECE 451 — Equivalent network synthesis
 */
int zparams_to_t_network(matrix2x2_t z, complex_t *za, complex_t *zb,
                         complex_t *zc);

/**
 * @brief Compute the impedance looking into port 1 with port 2 terminated.
 *
 * Zin = z11 - (z12 * z21) / (z22 + ZL)
 *
 * This is the fundamental formula for impedance transformation in a
 * feedback network. When z12 = 0 (unilateral), Zin = z11, independent of load.
 *
 * @param z Z-parameter matrix
 * @param zl Load impedance at port 2
 * @return Input impedance at port 1
 */
complex_t zparams_input_impedance(matrix2x2_t z, complex_t zl);

/**
 * @brief Compute the Thevenin output impedance at port 2.
 *
 * Zout = z22 - (z12 * z21) / (z11 + ZS)
 *
 * This includes the effect of source impedance ZS reflected to the output.
 * For unilateral networks (z12=0): Zout = z22.
 *
 * @param z Z-parameter matrix
 * @param zs Source impedance at port 1
 * @return Output impedance at port 2
 */
complex_t zparams_output_impedance(matrix2x2_t z, complex_t zs);

/**
 * @brief Compute the loaded voltage gain V2/V1.
 *
 * Av = z21 * ZL / (z11 * z22 - z12 * z21 + z11 * ZL)
 *
 * This general formula accounts for finite source and load impedance effects.
 *
 * @param z Z-parameter matrix
 * @param zl Load impedance
 * @return Loaded voltage gain (complex)
 */
complex_t zparams_voltage_gain(matrix2x2_t z, complex_t zl);

/**
 * @brief Compute the reverse voltage gain V1/V2 when port 2 is driven.
 *
 * Av_rev = z12 * ZS / (z11 * z22 - z12 * z21 + z22 * ZS)
 *
 * For a unilateral amplifier, this should be zero. Non-zero reverse gain
 * indicates feedback effects in the two-port.
 *
 * @param z Z-parameter matrix
 * @param zs Source impedance at port 1
 * @return Reverse voltage gain
 */
complex_t zparams_reverse_voltage_gain(matrix2x2_t z, complex_t zs);

/**
 * @brief Compute the condition for unconditional stability from Z-parameters.
 *
 * Using the Rollett stability factor adapted for Z-parameters:
 * K = (2*Re(z11)*Re(z22) - Re(z12*z21)) / |z12*z21|
 *
 * Unconditional stability requires K > 1.
 *
 * @param z Z-parameter matrix
 * @return Rollett K-factor for Z-parameters (dimensionless)
 */
double zparams_rollett_k(matrix2x2_t z);

#endif /* Z_PARAMS_H */
