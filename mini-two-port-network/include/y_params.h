/**
 * @file y_params.h
 * @brief Y-Parameter (Admittance Parameter) Analysis
 *
 * Y-parameters (short-circuit admittance parameters) are defined by port
 * currents in terms of port voltages:
 *
 *   I1 = y11*V1 + y12*V2
 *   I2 = y21*V1 + y22*V2
 *
 * Matrix form: [I] = [Y] * [V]
 *
 * Physical interpretation (with port short-circuited):
 *   y11 = I1/V1 | V2=0  (input admittance, output shorted)
 *   y12 = I1/V2 | V1=0  (reverse transfer admittance, input shorted)
 *   y21 = I2/V1 | V2=0  (forward transfer admittance, output shorted)
 *   y22 = I2/V2 | V1=0  (output admittance, input shorted)
 *
 * All Y-parameters have units of siemens (S).
 *
 * Y-parameters are natural for:
 *   - Parallel-connected networks (Y_total = Y1 + Y2)
 *   - π-equivalent circuit models
 *   - Nodal analysis in SPICE
 *   - FET small-signal models (yfs = forward transadmittance)
 *
 * Y = Z^(-1) for invertible Z-matrices.
 *
 * Reference: Desoer & Kuh, "Basic Circuit Theory" (1969), Ch. 10
 * Course: Berkeley EE16B — Admittance and nodal analysis
 */

#ifndef Y_PARAMS_H
#define Y_PARAMS_H

#include "two_port.h"

/* ============================================================================
 * L1: Y-parameter Creation Functions
 * ============================================================================ */

/**
 * @brief Create Y-parameter matrix from four complex admittances.
 *
 * @param y11 Input admittance with output shorted (S)
 * @param y12 Reverse transfer admittance with input shorted (S)
 * @param y21 Forward transfer admittance with output shorted (S)
 * @param y22 Output admittance with input shorted (S)
 * @return Y-parameter matrix
 */
matrix2x2_t yparams_create(complex_t y11, complex_t y12,
                            complex_t y21, complex_t y22);

/**
 * @brief Create Y-parameters for a shunt admittance Yp from port to ground.
 *
 * Circuit: port1 — + — port2; shunt Yp from line to ground
 * Y = [[Yp, -Yp], [-Yp, Yp]]
 *
 * The negative off-diagonal signs arise from the current divider effect.
 * This is the simplest Y-parameter network (single shunt element).
 *
 * @param yp Shunt admittance (S)
 * @return Y-parameter matrix
 *
 * Course: Georgia Tech ECE 6350 — Shunt admittance in microwave NWA
 */
matrix2x2_t yparams_shunt_admittance(complex_t yp);

/**
 * @brief Create Y-parameters for a series admittance Ys between ports.
 *
 * Circuit: port1 — Ys — port2; Ys is a series admittance.
 * Y = [[Ys, -Ys], [-Ys, Ys]]
 *
 * Used for: modeling series coupling capacitors, DC blocks in RF.
 *
 * @param ys Series admittance (S)
 * @return Y-parameter matrix
 */
matrix2x2_t yparams_series_admittance(complex_t ys);

/**
 * @brief Create Y-parameters for a π-network with three admittances.
 *
 *        Yb
 *   P1 —┬———┬——— P2
 *       Ya  Yc
 *       │   │
 *      GND GND
 *
 * Y-matrix:
 *   y11 = Ya + Yb         (input admittance, output shorted)
 *   y12 = -Yb = y21       (negative mutual admittance, reciprocal)
 *   y22 = Yb + Yc         (output admittance, input shorted)
 *
 * The π-network is the natural Y-parameter equivalent circuit.
 * Used for: matching networks, attenuators, filter sections.
 *
 * @param ya First shunt admittance (S)
 * @param yb Series admittance (S)
 * @param yc Second shunt admittance (S)
 * @return Y-parameter matrix
 *
 * Course: TU Munich HF Engineering — π-equivalent circuit
 */
matrix2x2_t yparams_pi_network(complex_t ya, complex_t yb, complex_t yc);

/**
 * @brief Create Y-parameters for a T-network (using Y-parameters).
 *
 *        Za       Zb
 *   P1 ————┬———————┬——— P2
 *          │       │
 *          Zc
 *          │
 *         GND
 *
 * The T-network is more naturally described by Z-parameters,
 * but Y-parameters exist:
 *   y11 = (Zb+Zc) / (Za*Zb + Zb*Zc + Zc*Za)
 *   y12 = -Zc / (Za*Zb + Zb*Zc + Zc*Za)
 *   y21 = y12 (reciprocal)
 *   y22 = (Za+Zc) / (Za*Zb + Zb*Zc + Zc*Za)
 *
 * @param za Impedance Za (first series arm, Ω)
 * @param zb Impedance Zb (second series arm, Ω)
 * @param zc Impedance Zc (shunt arm, Ω)
 * @return Y-parameter matrix
 */
matrix2x2_t yparams_t_network(complex_t za, complex_t zb, complex_t zc);

/**
 * @brief Create Y-parameters for a MOSFET small-signal model.
 *
 * Simplified high-frequency MOSFET model using Y-parameters:
 *   y11 = jω*Cgs (gate-source capacitance)
 *   y12 = jω*Cgd (gate-drain capacitance, Miller capacitance)
 *   y21 = gm     (transconductance, dominant term)
 *   y22 = go + jω*Cds (output conductance + drain-source capacitance)
 *
 * where go = 1/ro ≈ λ*Id (channel-length modulation).
 *
 * @param gm Transconductance (S)
 * @param go Output conductance (S) = 1/ro
 * @param cgs Gate-source capacitance (F)
 * @param cgd Gate-drain capacitance (F)
 * @param cds Drain-source capacitance (F)
 * @param omega Angular frequency (rad/s)
 * @return Y-parameter matrix for MOSFET
 *
 * Course: Stanford EE214 — High-frequency MOSFET models
 * Ref: Razavi, "Design of Analog CMOS Integrated Circuits", Ch. 6
 */
matrix2x2_t yparams_mosfet_hf(double gm, double go,
                               double cgs, double cgd, double cds,
                               double omega);

/**
 * @brief Create Y-parameters for a BJT hybrid-π model at high frequencies.
 *
 *   y11 = gπ + jω*(Cπ + Cμ)    (input admittance)
 *   y12 = -jω*Cμ              (reverse transfer via Cμ)
 *   y21 = gm - jω*Cμ           (forward transfer)
 *   y22 = go + jω*Cμ           (output admittance)
 *
 * where gπ = 1/rπ, go = 1/ro, Cπ = Cbe + Cdiff, Cμ = Cbc.
 *
 * @param rpi Base-emitter resistance rπ (Ω)
 * @param gm Transconductance (S) = Ic/VT
 * @param ro Output resistance (Ω) from Early effect
 * @param cpi Base-emitter capacitance Cπ (F)
 * @param cmu Collector-base capacitance Cμ (F)
 * @param omega Angular frequency (rad/s)
 * @return Y-parameter matrix for BJT
 *
 * Course: Berkeley EE105 — BJT frequency response
 * Ref: Sedra & Smith §10.2 — High-frequency BJT model
 */
matrix2x2_t yparams_bjt_hf(double rpi, double gm, double ro,
                            double cpi, double cmu, double omega);

/* ============================================================================
 * L3: Y-parameter Analysis Functions
 * ============================================================================ */

/**
 * @brief Compute the short-circuit current gain from Y-parameters.
 *
 * Ai_sc = I2/I1 | V2=0 = y21 / y11
 *
 * For a BJT: Ai_sc = gm / (gπ + jωCπ) ≈ gm*rπ = β at low frequencies.
 * The magnitude rolls off at the β cutoff frequency fβ.
 *
 * @param y Y-parameter matrix
 * @return Short-circuit current gain (complex)
 */
complex_t yparams_short_circuit_current_gain(matrix2x2_t y);

/**
 * @brief Compute the open-circuit voltage gain from Y-parameters.
 *
 * Av_oc = V2/V1 | I2=0 = -y21 / y22
 *
 * For a CS MOSFET: Av_oc = -gm/go = -gm*ro (intrinsic gain).
 * This represents the maximum voltage gain achievable from a single transistor.
 *
 * @param y Y-parameter matrix
 * @return Open-circuit voltage gain (complex)
 */
complex_t yparams_open_circuit_voltage_gain(matrix2x2_t y);

/**
 * @brief Compute the input admittance with load YL.
 *
 * Yin = y11 - (y12 * y21) / (y22 + YL)
 *
 * This is the admittance analog of Zin = z11 - z12*z21/(z22 + ZL).
 * Miller effect appears through y12*y21 term.
 *
 * @param y Y-parameter matrix
 * @param yl Load admittance (S)
 * @return Input admittance (S)
 */
complex_t yparams_input_admittance(matrix2x2_t y, complex_t yl);

/**
 * @brief Compute the output admittance with source YS.
 *
 * Yout = y22 - (y12 * y21) / (y11 + YS)
 *
 * @param y Y-parameter matrix
 * @param ys Source admittance (S)
 * @return Output admittance (S)
 */
complex_t yparams_output_admittance(matrix2x2_t y, complex_t ys);

/**
 * @brief Compute unity-gain frequency fT from Y-parameters.
 *
 * fT is the frequency where |Ai_sc(fT)| = 1 (0 dB current gain).
 *
 * For a MOSFET: fT ≈ gm / (2π * (Cgs + Cgd))
 * For a BJT: fT ≈ gm / (2π * (Cπ + Cμ))
 *
 * Using Y-parameters: |y21/y11| = 1 at fT.
 * This function solves fT numerically by binary search over frequency.
 *
 * @param gm Transconductance (S)
 * @param cgs Gate-source capacitance (F)
 * @param cgd Gate-drain capacitance (F)
 * @return fT in Hz
 *
 * Course: Stanford EE314 — RF IC Design
 * Ref: Lee, "The Design of CMOS Radio-Frequency Integrated Circuits", Ch. 3
 */
double yparams_f_t_mosfet(double gm, double cgs, double cgd);

/**
 * @brief Compute the maximum available gain (MAG) from Y-parameters.
 *
 * MAG = |y21 / y12|
 *
 * This is the maximum gain achievable when the device is simultaneously
 * conjugately matched at input and output. Requires stability (K > 1).
 *
 * For stable amplifiers (K > 1): MAG = |y21/y12| * (K - sqrt(K²-1)⁻¹)
 * For conditionally stable: MSG = |y21/y12| (maximum stable gain)
 *
 * @param y Y-parameter matrix
 * @return MAG (linear, not dB)
 */
double yparams_max_available_gain(matrix2x2_t y);

/* ============================================================================
 * L4: Admittance Analysis Theorems
 * ============================================================================ */

/**
 * @brief Compute the driving-point admittance at port 1 for arbitrary termination.
 *
 * Yin(s) = y11(s) - y12(s)*y21(s)*YL(s) / (1 + y22(s)*YL(s))
 *
 * This is the general formula for admittance transformation through a two-port.
 * At resonance: Im(Yin) = 0, giving the resonant frequency condition.
 *
 * @param y Y-parameter matrix (at single frequency)
 * @param yl Load admittance
 * @return Driving-point input admittance
 */
complex_t yparams_driving_point_admittance(matrix2x2_t y, complex_t yl);

/**
 * @brief Check if a Y-parameter matrix represents a passive network.
 *
 * A network is passive if Re(Y) is positive semi-definite:
 *   Re(y11) ≥ 0, Re(y22) ≥ 0
 *   |y12 + y21*|² ≤ 4 * Re(y11) * Re(y22)
 *
 * This condition ensures the network does not generate power.
 *
 * @param y Y-parameter matrix
 * @param tolerance Numerical tolerance
 * @return 1 if passive, 0 if active
 *
 * Course: Illinois ECE 451 — Passivity criteria
 */
int yparams_is_passive(matrix2x2_t y, double tolerance);

#endif /* Y_PARAMS_H */
