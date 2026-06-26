/**
 * @file resonance.h
 * @brief Resonance analysis for RLC circuits
 *
 * Resonance is the phenomenon where the inductive and capacitive
 * reactances cancel, resulting in purely resistive impedance.
 * It is fundamental to:
 * - Tuned amplifiers (radio, TV receivers)
 * - Oscillators (crystal, LC, RC)
 * - Filters (bandpass, bandstop)
 * - Impedance matching networks
 * - Wireless power transfer
 *
 * Resonance frequency (undamped natural frequency):
 *   ω₀ = 1/√(LC)    f₀ = 1/(2π√(LC))
 *
 * Reference:
 * - Hayt, Kemmerly & Durbin, "Engineering Circuit Analysis" (2019), Ch. 15
 * - Sedra & Smith, "Microelectronic Circuits" (2020), Ch. 15
 * - Terman, "Radio Engineers' Handbook" (1943)
 *
 * Course: Berkeley EE16B, MIT 6.003, Stanford EE102A
 */

#ifndef RESONANCE_H
#define RESONANCE_H

#include <stddef.h>
#include "frequency_response.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1 DEFINITIONS: Resonance circuit configurations
 * ============================================================================ */

/**
 * @brief RLC circuit parameters for resonance analysis.
 *
 * All RLC resonance circuits can be parameterized by:
 * - R: resistance (Ω) — dissipative element
 * - L: inductance (H) — stores magnetic energy
 * - C: capacitance (F) — stores electric energy
 *
 * The resonance frequency ω₀ = 1/√(LC) is independent of R.
 * However, R determines:
 * - Quality factor Q (sharpness of resonance)
 * - Bandwidth BW = ω₀/Q
 * - Damping factor ζ = R/(2L) for series, ζ = 1/(2RC) for parallel
 */
typedef struct {
    double R;          /**< Resistance (Ω) */
    double L;          /**< Inductance (H) */
    double C;          /**< Capacitance (F) */
    double source_V;   /**< Source voltage amplitude (V), AC analysis */
} rlc_params_t;

/* ============================================================================
 * RESONANCE COMPUTATION
 * ============================================================================ */

/**
 * @brief Analyze series RLC resonance.
 *
 * Series RLC circuit: R, L, C connected in series with voltage source.
 *
 * Impedance: Z(s) = R + sL + 1/(sC)
 *
 * At resonance (ω = ω₀):
 *   Z(jω₀) = R (purely resistive, minimum impedance)
 *   I_max = V_s/R (maximum current)
 *   V_L = V_C = Q·V_s (voltage magnification, may exceed source!)
 *   |V_L| = |V_C| = ω₀L·I = I/(ω₀C) = Q·V_s
 *
 * Phase of Z(jω):
 *   φ_Z(ω) = atan((ωL - 1/(ωC))/R)
 *   φ_Z < 0 (capacitive) for ω < ω₀
 *   φ_Z = 0 (resistive) at ω = ω₀
 *   φ_Z > 0 (inductive) for ω > ω₀
 *
 * Quality factor: Q = ω₀L/R = 1/(ω₀RC)
 *
 * L6 Canonical Problem: Series RLC analysis is the first resonance
 * problem encountered in circuits courses. It demonstrates:
 * - Cancellation of reactances at resonance
 * - Voltage magnification (Q times source voltage across L and C)
 * - Bandwidth-BPF response of the current
 *
 * @param params  RLC component values
 * @return        Resonance analysis result
 *
 * Reference: Hayt et al. (2019), Example 15.1-15.4
 */
resonance_result_t resonance_series(const rlc_params_t *params);

/**
 * @brief Analyze parallel RLC resonance.
 *
 * Parallel RLC circuit: R, L, C connected in parallel with current source
 * (or voltage source with series resistance R).
 *
 * Admittance: Y(s) = 1/R + 1/(sL) + sC
 *
 * At resonance (ω = ω₀):
 *   Y(jω₀) = 1/R (purely conductive, minimum admittance = maximum impedance)
 *   Z_max = R (maximum impedance)
 *   I_L = I_C = Q·I_source (current magnification in LC tank)
 *
 * Phase of Y(jω):
 *   φ_Y(ω) = atan((ωC - 1/(ωL))/(1/R))
 *   φ_Y < 0 (inductive) for ω < ω₀
 *   φ_Y = 0 (resistive) at ω = ω₀
 *   φ_Y > 0 (capacitive) for ω > ω₀
 *
 * Quality factor: Q = R/(ω₀L) = ω₀RC
 *
 * Note: Q formula is the reciprocal of the series case.
 * For high-Q parallel resonance, the circulating (tank) current
 * is Q times the source current.
 *
 * L6 Canonical Problem: Parallel RLC is the dual of series RLC.
 * Used in tuned amplifiers where high impedance at resonance
 * provides voltage gain.
 *
 * @param params  RLC component values
 * @return        Resonance analysis result
 */
resonance_result_t resonance_parallel(const rlc_params_t *params);

/**
 * @brief Compute the transfer function of a series RLC circuit.
 *
 * Five possible transfer functions for series RLC, each with
 * different filtering characteristics:
 *
 * 1. V_R/V_in (across R): Bandpass
 *    H(s) = R/(R + sL + 1/(sC)) = (R/L)s/(s² + (R/L)s + 1/(LC))
 *
 * 2. V_L/V_in (across L): Highpass with resonance
 *    H(s) = sL/(R + sL + 1/(sC)) = s²/(s² + (R/L)s + 1/(LC))·L
 *
 * 3. V_C/V_in (across C): Lowpass with resonance
 *    H(s) = (1/(sC))/(R + sL + 1/(sC)) = 1/(LC)/(s² + (R/L)s + 1/(LC))
 *
 * 4. V_L+V_C/V_in: Bandstop (notch)
 *    H(s) = (sL + 1/(sC))/(R + sL + 1/(sC))
 *
 * 5. I/V_in (current): Bandpass
 *    Y(s) = 1/(R + sL + 1/(sC))
 *
 * @param params    RLC parameters
 * @param output    0=VR, 1=VL, 2=VC, 3=VL+VC, 4=I
 * @return          Transfer function H(s), or NULL
 */
tf_polynomial_t *resonance_series_tf(const rlc_params_t *params, int output);

/**
 * @brief Compute the transfer function of a parallel RLC circuit.
 *
 * Parallel RLC with input current source (or Thevenin equivalent):
 *
 * 1. V_out (voltage across tank): Bandpass
 *    Z(s) = 1/(1/R + 1/(sL) + sC) = s/(C)/(s² + s/(RC) + 1/(LC))
 *
 * 2. I_R (current through R): Same shape as V_out (bandpass)
 * 3. I_L (current through L): Lowpass
 * 4. I_C (current through C): Highpass
 *
 * @param params    RLC parameters
 * @param output    0=V_out, 1=I_R, 2=I_L, 3=I_C
 * @return          Transfer function H(s), or NULL
 */
tf_polynomial_t *resonance_parallel_tf(const rlc_params_t *params, int output);

/* ============================================================================
 * UNIVERSAL RESONANCE CURVE
 * ============================================================================ */

/**
 * @brief Compute the universal resonance curve.
 *
 * The magnitude response of any second-order bandpass filter
 * (series or parallel RLC) can be expressed in normalized form:
 *
 *   |H(jω)|/|H_max| = 1/√(1 + Q²(ω/ω₀ - ω₀/ω)²)
 *
 * This is the "universal resonance curve" because it depends only
 * on Q and the normalized frequency Ω = ω/ω₀.
 *
 * At the half-power points:
 *   Ω_{H,L} = √(1 + 1/(4Q²)) ± 1/(2Q)
 *   For Q >> 1: Ω_{H,L} ≈ 1 ± 1/(2Q)
 *
 * The sharpness of resonance:
 *   BW/f₀ = 1/Q
 *   For Q = 10: BW = f₀/10 (10% bandwidth)
 *   For Q = 100: BW = f₀/100 (1% bandwidth)
 *
 * @param norm_freq  Array of normalized frequencies Ω = f/f₀
 * @param Q          Quality factor
 * @param n          Number of frequency points
 * @return           Normalized magnitude |H|/|H_max| array (caller frees)
 *
 * L2 Concept: The universal resonance curve encapsulates the
 * frequency-domain behavior of ALL second-order resonant systems,
 * from electrical RLC circuits to mechanical mass-spring-damper
 * systems to optical Fabry-Perot cavities. This is an example of
 * the unifying power of transfer function analysis.
 */
double *resonance_universal_curve(const double *norm_freq,
                                    double Q, size_t n);

/* ============================================================================
 * QUALITY FACTOR ANALYSIS
 * ============================================================================ */

/**
 * @brief Compute quality factor Q from -3 dB bandwidth measurement.
 *
 * Q = f₀/(f_H - f_L)
 *
 * This is the most common method for measuring Q experimentally.
 *
 * @param f0    Resonance frequency (Hz)
 * @param f_low Lower -3 dB frequency (Hz)
 * @param f_high Upper -3 dB frequency (Hz)
 * @return      Quality factor Q
 */
double resonance_q_from_bandwidth(double f0, double f_low, double f_high);

/**
 * @brief Compute quality factor from energy storage perspective.
 *
 * Q = 2π·(maximum energy stored)/(energy dissipated per cycle)
 *
 * For series RLC: Q = ω₀L/R (energy stored in L and C / energy lost in R)
 * For parallel RLC: Q = ω₀RC
 *
 * This energy definition is the most fundamental — it applies to
 * any resonant system, not just electrical.
 *
 * L1 Definition: The quality factor Q fundamentally represents the
 * ratio of stored to dissipated energy. High Q means low loss,
 * sharp resonance, and long ring-down time.
 *
 * Ring-down (decay) time constant: τ = 2Q/ω₀ = Q/(πf₀)
 * Number of cycles to decay to 1/e: N_cycles = Q/π
 *
 * @param params       RLC parameters
 * @param topology     RESONANCE_SERIES_RLC or RESONANCE_PARALLEL_RLC
 * @return             Q from energy definition
 */
double resonance_q_from_energy(const rlc_params_t *params,
                                 resonance_topology_t topology);

/**
 * @brief Compute the damping factor ζ from Q.
 *
 * ζ = 1/(2Q)
 *
 * ζ < 1:  underdamped (oscillatory, peaking in frequency response)
 * ζ = 1:  critically damped (fastest settling without overshoot)
 * ζ > 1:  overdamped (sluggish, no oscillation)
 * ζ = 0:  undamped (sustained oscillation, poles on jω-axis)
 * ζ = 1/√2 ≈ 0.707: maximally flat magnitude (Butterworth)
 * ζ = 1/2 = 0.5: Bessel (maximally flat group delay)
 *
 * @param Q  Quality factor
 * @return   Damping factor ζ
 */
double resonance_damping_from_q(double Q);

/* ============================================================================
 * STEP RESPONSE OF RESONANT CIRCUITS
 * ============================================================================ */

/**
 * @brief Compute the step response of a series RLC circuit.
 *
 * The step response of a second-order system with damping ζ:
 *
 * For ζ < 1 (underdamped):
 *   v(t) = V_final·[1 - e^{-ζω₀t}·(cos(ω_d·t) + ζ/√(1-ζ²)·sin(ω_d·t))]
 *   where ω_d = ω₀·√(1-ζ²) is the damped natural frequency.
 *
 * Peak overshoot: M_p = exp(-πζ/√(1-ζ²))
 * Peak time: t_p = π/ω_d
 * Settling time (2%): t_s ≈ 4/(ζω₀)
 *
 * For ζ = 1 (critically damped):
 *   v(t) = V_final·[1 - (1 + ω₀t)·e^{-ω₀t}]
 *
 * For ζ > 1 (overdamped):
 *   v(t) = V_final·[1 - (τ₂/(τ₂-τ₁))·e^{-t/τ₁} - (τ₁/(τ₁-τ₂))·e^{-t/τ₂}]
 *   where τ₁,₂ = 1/(ω₀(ζ ± √(ζ²-1)))
 *
 * @param params    RLC parameters
 * @param t         Time array (seconds)
 * @param v_out     Output voltage array (allocated by caller)
 * @param n         Number of time points
 * @param topology  Series or parallel topology
 *
 * L6 Canonical Problem: The step response of a resonant circuit
 * reveals the complete dynamic behavior: overshoot, ringing frequency,
 * settling time. This bridges frequency-domain (Q, ζ) and time-domain
 * (overshoot, settling) characterizations.
 *
 * Course: Berkeley EE16B, MIT 6.003
 */
void resonance_step_response(const rlc_params_t *params,
                              const double *t, double *v_out,
                              size_t n, resonance_topology_t topology);

/* ============================================================================
 * COUPLED RESONATORS (L8 Advanced)
 * ============================================================================ */

/**
 * @brief Analyze two magnetically coupled RLC resonators.
 *
 * Coupled resonators with mutual inductance M:
 * - Coupling coefficient: k = M/√(L₁L₂)
 * - Critical coupling: k_c = 1/√(Q₁Q₂)  (for identical Q: k_c = 1/Q)
 * - Over-coupled (k > k_c): double-peaked response
 * - Critically coupled (k = k_c): maximally flat (transitional)
 * - Under-coupled (k < k_c): single peak, reduced bandwidth
 *
 * The split resonant frequencies for over-coupled identical resonators:
 *   ω_{1,2} = ω₀/√(1 ± k)
 *
 * L8 Advanced: Coupled resonator theory is fundamental to:
 * - Bandpass filter design (coupled resonator filters)
 * - Wireless power transfer (magnetic resonance coupling)
 * - Impedance matching networks
 * - Oscillator injection locking
 *
 * @param primary    Primary RLC parameters
 * @param secondary  Secondary RLC parameters
 * @param M          Mutual inductance (H)
 * @param freq       Frequency array (Hz)
 * @param n_freq     Number of frequency points
 * @return           Frequency response of secondary voltage, or NULL
 *
 * Reference: Terman (1943), Dishal (1951)
 * Course: ETH 227-0455, Michigan EECS 411
 */
freq_response_t *resonance_coupled(const rlc_params_t *primary,
                                     const rlc_params_t *secondary,
                                     double M,
                                     const double *freq, size_t n_freq);

/* ============================================================================
 * CRYSTAL RESONATOR EQUIVALENT CIRCUIT (L7 Application)
 * ============================================================================ */

/**
 * @brief Model a quartz crystal resonator.
 *
 * A quartz crystal is electrically modeled as:
 * - C₀: static (holder) capacitance (~1-10 pF)
 * - L₁: motional inductance (very large, ~mH-H equivalent)
 * - C₁: motional capacitance (very small, ~fF)
 * - R₁: motional resistance (loss, ~10-100 Ω)
 *
 * The crystal has TWO resonant frequencies:
 *   Series resonance: f_s = 1/(2π√(L₁C₁))
 *   Parallel resonance: f_p = f_s·√(1 + C₁/C₀) ≈ f_s·(1 + C₁/(2C₀))
 *
 * Typically f_p - f_s is very small (~100-1000 ppm), making crystals
 * excellent for precision frequency references (Q ~ 10⁴ to 10⁶).
 *
 * The figure shows the reactance changes from inductive (between
 * f_s and f_p) to capacitive (outside this range). The crystal is
 * used as an inductive element in Pierce and Colpitts oscillators.
 *
 * @param C0    Static capacitance (F)
 * @param L1    Motional inductance (H)
 * @param C1    Motional capacitance (F)
 * @param R1    Motional resistance (Ω)
 * @param fs    Output: series resonant frequency (Hz)
 * @param fp    Output: parallel resonant frequency (Hz)
 * @param Q     Output: quality factor
 * @return      Impedance transfer function Z(s), or NULL
 *
 * L7 Application: Quartz crystal resonators are ubiquitous in
 * electronic systems for frequency control:
 * - Microcontroller clock oscillators (8-40 MHz fundamental)
 * - Real-time clocks (32.768 kHz tuning fork crystals)
 * - RF frequency synthesis (TCXO, OCXO for GPS, cellular)
 * - Filtering (crystal ladder filters for IF stages)
 *
 * Reference: Bottom, "Introduction to Quartz Crystal Unit Design" (1982)
 * Course: Berkeley EE105 oscillator design, Stanford EE247
 */
tf_polynomial_t *resonance_crystal_model(double C0, double L1, double C1,
                                           double R1,
                                           double *fs, double *fp, double *Q);

#ifdef __cplusplus
}
#endif

#endif /* RESONANCE_H */
