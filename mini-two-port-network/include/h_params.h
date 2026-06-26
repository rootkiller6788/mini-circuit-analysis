/**
 * @file h_params.h
 * @brief H-Parameter (Hybrid Parameter) Analysis
 *
 * H-parameters (hybrid parameters) are the classic BJT small-signal model,
 * mixing an input impedance with output admittance:
 *
 *   V1 = h11*I1 + h12*V2
 *   I2 = h21*I1 + h22*V2
 *
 * Matrix form: [[V1], [I2]] = [[h11, h12], [h21, h22]] * [[I1], [V2]]
 *
 * Physical interpretation:
 *   h11 = V1/I1 | V2=0  (input impedance, output shorted) — hie/Ω
 *   h12 = V1/V2 | I1=0  (reverse voltage ratio, input open) — hre/dimensionless
 *   h21 = I2/I1 | V2=0  (forward current gain, output shorted) — hfe/dimensionless
 *   h22 = I2/V2 | I1=0  (output admittance, input open) — hoe/S
 *
 * BJT notation:
 *   hie: input impedance (common-emitter), ~rπ
 *   hre: reverse voltage ratio, ~0 (typically 10⁻⁴)
 *   hfe: forward current gain, β (typically 50-500)
 *   hoe: output admittance, ~1/ro (typically 10⁻⁵ S)
 *
 * H-parameters are natural for:
 *   - BJT CE amplifiers (the original application)
 *   - Series-parallel connected networks (H_total = H1 + H2)
 *   - Feedback analysis with series mixing and shunt sampling
 *
 * Reference: Sedra & Smith, "Microelectronic Circuits" (2020), Ch. 7
 * Course: Berkeley EE105 — BJT amplifier analysis
 */

#ifndef H_PARAMS_H
#define H_PARAMS_H

#include "two_port.h"

/* ============================================================================
 * L1: H-parameter Creation Functions
 * ============================================================================ */

/**
 * @brief Create H-parameter matrix from four mixed parameters.
 *
 * @param h11 Input impedance with output shorted (Ω)
 * @param h12 Reverse voltage ratio with input open (dimensionless)
 * @param h21 Forward current gain with output shorted (dimensionless)
 * @param h22 Output admittance with input open (S)
 * @return H-parameter matrix
 */
matrix2x2_t hparams_create(complex_t h11, complex_t h12,
                            complex_t h21, complex_t h22);

/**
 * @brief Create H-parameters for a common-emitter BJT.
 *
 * h_ie = rbb' + rπ        (input impedance, ≈ rπ when rbb' is small)
 * h_re ≈ 0                (reverse voltage transfer, typically ~10⁻⁴)
 * h_fe = β = gm*rπ        (forward current gain)
 * h_oe = 1/ro ≈ Ic/VA    (output admittance)
 *
 * These are the classic "h-parameters" found in transistor datasheets.
 * h_re is often neglected in hand analysis but can be computed from
 * the Early effect as h_re ≈ rπ / (ro + rπ).
 *
 * @param rbb Base spreading resistance (Ω), typically 10-100 Ω
 * @param rpi Base-emitter small-signal resistance rπ (Ω)
 * @param beta DC current gain hFE ≈ ic/ib
 * @param ro Early-effect output resistance (Ω)
 * @return H-parameter matrix for CE BJT
 *
 * Course: Berkeley EE105 — H-parameter BJT model
 * Ref: Gray, Hurst, Lewis, Meyer, "Analysis and Design of Analog ICs", Ch. 1
 */
matrix2x2_t hparams_bjt_ce(double rbb, double rpi, double beta, double ro);

/**
 * @brief Create H-parameters for a common-base BJT.
 *
 * Converted from CE h-parameters via standard transistor configuration
 * transformation formulas:
 *
 * h_ib = h_ie / (1 + h_fe)          (input impedance, very low ~rπ/β)
 * h_rb = h_ie*h_oe/(1+h_fe) - h_re  (≈ 0)
 * h_fb = -h_fe / (1 + h_fe)         (current gain ≈ -α, just under -1)
 * h_ob = h_oe / (1 + h_fe)          (output admittance, very low)
 *
 * where α = β/(1+β) ≈ 0.99 for typical β=100.
 * CB configuration features: low Zin, high Zout, current gain < 1.
 *
 * @param h_ce CE H-parameter matrix
 * @return H-parameter matrix for CB
 *
 * Course: MIT 6.002 — Transistor configurations
 */
matrix2x2_t hparams_ce_to_cb(matrix2x2_t h_ce);

/**
 * @brief Create H-parameters for a common-collector (emitter follower) BJT.
 *
 * h_ic = h_ie                      (input impedance, same as CE, very high)
 * h_rc = 1 - h_re ≈ 1             (voltage gain ≈ 1)
 * h_fc = -(1 + h_fe)              (current gain ≈ -β)
 * h_oc = h_oe                     (output admittance, same as CE)
 *
 * CC configuration features: high Zin, low Zout, voltage gain ≈ 1.
 * Used for: impedance buffering, voltage regulators.
 *
 * @param h_ce CE H-parameter matrix
 * @return H-parameter matrix for CC
 *
 * Course: Stanford EE101 — Emitter follower analysis
 */
matrix2x2_t hparams_ce_to_cc(matrix2x2_t h_ce);

/**
 * @brief Create H-parameters for a BJT including high-frequency effects.
 *
 * At high frequencies, the CE H-parameters become complex:
 *   h11 = rbb' + rπ / (1 + jω/ωβ)   (input impedance with frequency rolloff)
 *   h12 = jω*Cμ / (gπ + jω*Cπ)      (reverse via Cμ, ≈0 at low f)
 *   h21 = hfe / (1 + jω/ωβ)         (current gain rolloff, single-pole)
 *   h22 = hoe + jω*(Cμ+Ccs)         (output admittance increases with f)
 *
 * @param rbb Base spreading resistance (Ω)
 * @param rpi Base-emitter rπ (Ω)
 * @param beta DC current gain
 * @param ro Output resistance (Ω)
 * @param cpi Cπ capacitance (F)
 * @param cmu Cμ capacitance (F)
 * @param omega Angular frequency (rad/s)
 * @return Frequency-dependent H-parameter matrix
 *
 * Course: Stanford EE214 — Frequency response of BJT amplifiers
 */
matrix2x2_t hparams_bjt_ce_hf(double rbb, double rpi, double beta, double ro,
                               double cpi, double cmu, double omega);

/* ============================================================================
 * L3: H-parameter Analysis Functions
 * ============================================================================ */

/**
 * @brief Compute the voltage gain of a CE amplifier using H-parameters.
 *
 * Av = V2/V1 = -h21 / (h11 * h22 - h12 * h21 + h11 * YL)
 *
 * Simplified for typical BJT (h12 ≈ 0, h22 small):
 * Av ≈ -h21 / (h11 * YL) = -hfe * RL / hie = -gm * RL
 *
 * @param h H-parameter matrix
 * @param yl Load admittance (S) = 1/RL
 * @return Complex voltage gain
 */
complex_t hparams_voltage_gain(matrix2x2_t h, complex_t yl);

/**
 * @brief Compute the current gain of a CE amplifier using H-parameters.
 *
 * Ai = I2/I1 = h21 / (1 + h22 * ZL)
 *
 * For ZL → 0 (short circuit): Ai → h21 = hfe
 * For typical loads (ZL << 1/h22): Ai ≈ hfe
 *
 * @param h H-parameter matrix
 * @param zl Load impedance (Ω)
 * @return Complex current gain
 */
complex_t hparams_current_gain(matrix2x2_t h, complex_t zl);

/**
 * @brief Compute the input impedance using H-parameters.
 *
 * Zin = h11 - (h12 * h21) / (h22 + YL)
 *
 * For typical BJT with negligible h12: Zin ≈ h11 = hie
 *
 * @param h H-parameter matrix
 * @param yl Load admittance (S) = 1/ZL
 * @return Input impedance (Ω)
 */
complex_t hparams_input_impedance(matrix2x2_t h, complex_t yl);

/**
 * @brief Compute the output admittance using H-parameters.
 *
 * Yout = h22 - (h12 * h21) / (h11 + ZS)
 *
 * For typical BJT with negligible h12: Yout ≈ h22 = hoe
 *
 * @param h H-parameter matrix
 * @param zs Source impedance (Ω)
 * @return Output admittance (S)
 */
complex_t hparams_output_admittance(matrix2x2_t h, complex_t zs);

/**
 * @brief Compute the power gain using H-parameters.
 *
 * Gp = |h21|^2 * Re(YL) / (Re(Yin) * |1 + h22*ZL|^2)
 *
 * where Yin = h11 - h12*h21*YL/(1+h22*ZL)
 *
 * @param h H-parameter matrix
 * @param yl Load admittance (S)
 * @return Power gain (linear)
 */
double hparams_power_gain(matrix2x2_t h, complex_t yl);

/**
 * @brief Determine β cutoff frequency from H-parameters.
 *
 * At f = fβ, |hfe(fβ)| = hfe_DC / √2 (3 dB drop).
 *
 * Using the single-pole model: hfe(s) = hfe_DC / (1 + s/ωβ)
 * fβ = ωβ / (2π) where ωβ = 1 / (rπ * Cπ)
 *
 * @param rpi rπ (Ω)
 * @param cpi Cπ total (Cdiff + Cje depletion), in F
 * @return fβ in Hz
 *
 * Course: Stanford EE314 — BJT frequency limits
 */
double hparams_beta_cutoff_freq(double rpi, double cpi);

/**
 * @brief Compute the hfe vs frequency magnitude at a given frequency.
 *
 * |hfe(ω)| = hfe_DC / sqrt(1 + (ω/ωβ)²)
 *
 * Single-pole magnitude response. The phase is -atan(ω/ωβ).
 *
 * @param beta_dc DC current gain
 * @param f_beta Beta cutoff frequency (Hz)
 * @param freq Frequency of interest (Hz)
 * @return |hfe| at the given frequency
 */
double hparams_hfe_magnitude(double beta_dc, double f_beta, double freq);

#endif /* H_PARAMS_H */
