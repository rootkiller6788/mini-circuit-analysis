/**
 * @file frequency_response.h
 * @brief Core frequency response definitions and types
 *
 * This header defines the fundamental data structures for frequency
 * response analysis in linear time-invariant (LTI) circuits and systems.
 *
 * Key concepts covered:
 * - Frequency response H(jω) = |H(jω)|∠H(jω)
 * - Magnitude response in dB: 20·log₁₀|H(jω)|
 * - Phase response in degrees or radians
 * - Linear and logarithmic frequency scales
 * - Cutoff frequency, bandwidth, roll-off rate
 *
 * Reference textbooks:
 * - Oppenheim & Willsky, "Signals and Systems" (1997), Ch. 6
 * - Sedra & Smith, "Microelectronic Circuits" (2020), Ch. 8, 9
 * - Hayt, Kemmerly & Durbin, "Engineering Circuit Analysis" (2019), Ch. 15
 *
 * Course alignment:
 * - MIT 6.003 Signal Processing: Frequency response fundamentals
 * - Berkeley EE16B: Frequency domain analysis
 * - Stanford EE102A: LTI system frequency response
 * - ETH 227-0427: Signal Processing frequency analysis
 */

#ifndef FREQUENCY_RESPONSE_H
#define FREQUENCY_RESPONSE_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1 DEFINITIONS: Core frequency response data types
 * ============================================================================ */

/**
 * @brief Frequency scale type enumeration.
 *
 * Two fundamental frequency representations exist in circuit analysis:
 *
 * 1. LINEAR: f (Hz) — Direct frequency in hertz, useful for narrow-band
 *    analysis and time-domain ↔ frequency-domain conversion via ω = 2πf.
 *
 * 2. LOGARITHMIC: log₁₀(f) — Bode plots use logarithmic frequency axis
 *    to compress wide frequency ranges. One decade = factor of 10 in
 *    frequency. One octave = factor of 2 in frequency.
 *
 * In Bode plot asymptotes, the frequency response is approximated by
 * straight lines on log-log (magnitude) and log-linear (phase) scales.
 *
 * L2 Concept: Bode plot construction relies on logarithmic frequency
 * representation to turn multiplicative pole/zero effects into additive
 * asymptotic contributions.
 */
typedef enum {
    FREQ_SCALE_LINEAR,       /**< Linear frequency in Hz */
    FREQ_SCALE_LOGARITHMIC,  /**< Logarithmic frequency, log₁₀(f) */
    FREQ_SCALE_RAD_PER_SEC   /**< Angular frequency ω (rad/s) */
} freq_scale_t;

/**
 * @brief Magnitude representation format.
 *
 * L1 Definition: The magnitude of a frequency response can be expressed
 * in several equivalent forms:
 *
 * - ABSOLUTE: |H(jω)| — raw voltage/current gain
 * - DB: 20·log₁₀|H(jω)| — decibel scale, standard in Bode plots
 * - NEPER: ln|H(jω)| — natural log scale, used in filter theory
 *
 * The decibel (dB) is defined by the International System of Units (SI)
 * as a logarithmic unit for ratios. In circuit analysis:
 *   Gain_dB = 20·log₁₀(V_out/V_in) for voltage
 *   Gain_dB = 10·log₁₀(P_out/P_in) for power
 *
 * Roll-off rate is measured in dB/decade or dB/octave:
 *   - Single pole: -20 dB/decade (-6 dB/octave)
 *   - Single zero: +20 dB/decade (+6 dB/octave)
 *   - Double pole: -40 dB/decade (-12 dB/octave)
 */
typedef enum {
    MAG_FORMAT_ABSOLUTE,  /**< |H(jω)| — linear magnitude */
    MAG_FORMAT_DB,        /**< 20·log₁₀|H(jω)| — decibel */
    MAG_FORMAT_NEPER      /**< ln|H(jω)| — neper */
} mag_format_t;

/**
 * @brief Phase representation format.
 *
 * L1 Definition: Phase response φ(ω) = ∠H(jω) represents the phase shift
 * introduced by the circuit at angular frequency ω. For a transfer function
 * H(jω) = |H(jω)|·e^{jφ(ω)}, the phase is:
 *   φ(ω) = atan2(Im{H(jω)}, Re{H(jω)})
 *
 * Phase can be expressed in:
 * - Degrees: [0°, 360°) or [-180°, 180°)
 * - Radians: [0, 2π) or [-π, π)
 *
 * Group delay τ_g(ω) = -dφ(ω)/dω is the derivative of phase with respect
 * to angular frequency, measuring the time delay of the envelope of a
 * narrow-band signal passing through the circuit.
 */
typedef enum {
    PHASE_FORMAT_DEGREES,       /**< Phase in degrees */
    PHASE_FORMAT_RADIANS,       /**< Phase in radians */
    PHASE_FORMAT_GROUP_DELAY    /**< Group delay τ_g = -dφ/dω (seconds) */
} phase_format_t;

/**
 * @brief Filter classification by passband characteristic.
 *
 * L1 Definition: Filters are classified by which frequency bands they
 * pass (low attenuation) and which they reject (high attenuation):
 *
 * - LOWPASS: Passes DC to f_c, rejects above f_c
 * - HIGHPASS: Rejects DC to f_c, passes above f_c
 * - BANDPASS: Passes f_L to f_H, rejects outside
 * - BANDSTOP: Rejects f_L to f_H, passes outside (notch filter)
 * - ALLPASS: Passes all frequencies equally in magnitude,
 *            but introduces frequency-dependent phase shift.
 *            Used for phase equalization / delay equalization.
 *
 * The cutoff frequency f_c (or -3dB frequency) is defined as the frequency
 * where |H(jω_c)| = |H_max|/√2, equivalent to -3.01 dB below the maximum.
 * For bandpass/bandstop, two cutoff frequencies exist: f_L (lower) and
 * f_H (upper). The bandwidth BW = f_H - f_L.
 *
 * L2 Concept: Filter order N determines the asymptotic roll-off rate:
 * roll-off = ±N × 20 dB/decade = ±N × 6 dB/octave.
 */
typedef enum {
    FILTER_TYPE_LOWPASS,    /**< Low-pass filter: passes f < f_c */
    FILTER_TYPE_HIGHPASS,   /**< High-pass filter: passes f > f_c */
    FILTER_TYPE_BANDPASS,   /**< Band-pass filter: passes f_L < f < f_H */
    FILTER_TYPE_BANDSTOP,   /**< Band-stop filter: rejects f_L < f < f_H */
    FILTER_TYPE_ALLPASS     /**< All-pass filter: |H|=1, phase varies */
} filter_type_t;

/**
 * @brief Filter approximation type.
 *
 * L5 Algorithm: Different polynomial approximations are used to design
 * analog filters with specific magnitude response characteristics:
 *
 * - BUTTERWORTH: Maximally flat magnitude in passband.
 *   |H(jω)|² = 1 / (1 + (ω/ω_c)^{2N})
 *   Poles equally spaced on a circle in the s-plane.
 *
 * - CHEBYSHEV_I: Equiripple in passband, monotonic in stopband.
 *   |H(jω)|² = 1 / (1 + ε²·T_N²(ω/ω_c))
 *   where T_N(x) is the Chebyshev polynomial of order N.
 *   Sharper cutoff than Butterworth but with passband ripple.
 *
 * - CHEBYSHEV_II: Monotonic in passband, equiripple in stopband.
 *   Also called Inverse Chebyshev.
 *
 * - ELLIPTIC: Equiripple in both passband and stopband.
 *   Also called Cauer filter. Uses Jacobi elliptic functions.
 *   Provides the sharpest transition band for a given order.
 *   |H(jω)|² = 1 / (1 + ε²·R_N²(ξ, ω/ω_c))
 *   where R_N is the Chebyshev rational function.
 *
 * - BESSEL: Maximally flat group delay (linear phase response).
 *   Uses Bessel polynomials. Best for preserving waveform shape
 *   (minimal overshoot, constant group delay).
 *
 * - GAUSSIAN: Gaussian magnitude response.
 *   |H(jω)| = exp(-α·ω²). No overshoot in step response, but slow roll-off.
 *
 * Reference: Zverev, "Handbook of Filter Synthesis" (1967)
 * Course: Stanford EE247, Berkeley EE105, ETH 227-0455
 */
typedef enum {
    APPROX_BUTTERWORTH,   /**< Maximally flat magnitude */
    APPROX_CHEBYSHEV_I,   /**< Equiripple passband */
    APPROX_CHEBYSHEV_II,  /**< Equiripple stopband */
    APPROX_ELLIPTIC,      /**< Equiripple passband & stopband (Cauer) */
    APPROX_BESSEL,        /**< Maximally flat group delay */
    APPROX_GAUSSIAN       /**< Gaussian magnitude response */
} filter_approx_t;

/* ============================================================================
 * L1 DEFINITIONS: Complex frequency response data structure
 * ============================================================================ */

/**
 * @brief Single frequency point in a frequency response measurement.
 *
 * Represents the complex frequency response H(jω) = |H|·e^{jφ} at
 * a single frequency point. This is the atomic unit of frequency
 * response data.
 *
 * The complex value is stored in both rectangular (real, imag) and
 * polar (magnitude, phase) forms for computational convenience.
 * When computing from a transfer function H(s) with s = jω:
 *   real = Re{H(jω)}, imag = Im{H(jω)}
 *   magnitude = sqrt(real² + imag²) = |H(jω)|
 *   phase = atan2(imag, real) = ∠H(jω)
 *
 * L3 Math Structure: Complex number representation — the fundamental
 * mathematical structure underlying all frequency response analysis.
 */
typedef struct {
    double frequency;       /**< Frequency in Hz */
    double angular_freq;    /**< Angular frequency ω = 2πf (rad/s) */
    double real;            /**< Real part: Re{H(jω)} */
    double imag;            /**< Imaginary part: Im{H(jω)} */
    double magnitude;       /**< Magnitude: |H(jω)| */
    double phase_rad;       /**< Phase in radians: ∠H(jω) */
    double magnitude_db;    /**< Magnitude in dB: 20·log₁₀|H(jω)| */
} freq_point_t;

/**
 * @brief Complete frequency response data set.
 *
 * A frequency response is a set of complex-valued measurements
 * at discrete frequency points spanning the range of interest.
 * This structure holds the full sweep data along with metadata.
 *
 * The frequency sweep can be linear (equally spaced in Hz) or
 * logarithmic (equally spaced in log₁₀(f)), depending on the
 * application. Bode plots typically use logarithmic spacing.
 *
 * L2 Concept: The frequency response completely characterizes
 * an LTI system's steady-state response to sinusoidal inputs.
 * For input x(t) = A·cos(ωt + θ), the output is:
 *   y(t) = A·|H(jω)|·cos(ωt + θ + ∠H(jω))
 */
typedef struct {
    freq_point_t *points;   /**< Array of frequency response points */
    size_t num_points;      /**< Number of points in the sweep */
    double freq_start;      /**< Start frequency (Hz) */
    double freq_end;        /**< End frequency (Hz) */
    freq_scale_t scale;     /**< Frequency scale type */
    char label[128];        /**< Descriptive label for this response */
} freq_response_t;

/* ============================================================================
 * L1 DEFINITIONS: Bode plot data
 * ============================================================================ */

/**
 * @brief Asymptotic Bode plot corner (break) point.
 *
 * L5 Algorithm: In asymptotic Bode plot construction, each pole and zero
 * contributes a corner frequency where the slope changes:
 *
 * - Pole at ω₀: slope changes by -20 dB/decade starting at ω₀
 * - Zero at ω₀: slope changes by +20 dB/decade starting at ω₀
 * - Complex conjugate pole pair at ω₀: -40 dB/decade at ω₀
 *   (with peaking near ω₀ depending on damping factor ζ)
 * - Complex conjugate zero pair: +40 dB/decade at ω₀
 *
 * Phase contribution of a single pole:
 *   φ(ω) ≈ 0°       for ω << ω₀
 *   φ(ω) ≈ -45°     for ω = ω₀
 *   φ(ω) ≈ -90°     for ω >> ω₀
 *   Transition: ±45° per decade around ω₀
 *
 * Phase contribution of a single zero:
 *   φ(ω) ≈ 0°       for ω << ω₀
 *   φ(ω) ≈ +45°     for ω = ω₀
 *   φ(ω) ≈ +90°     for ω >> ω₀
 */
typedef struct {
    double corner_freq;          /**< Corner (break) frequency ω₀ (rad/s) */
    int is_pole;                 /**< 1 = pole, 0 = zero */
    int is_complex_pair;         /**< 1 = complex conjugate pair */
    double damping_factor;       /**< Damping factor ζ (for complex pairs) */
    double slope_change_db_dec;  /**< Slope change in dB/decade */
    double phase_shift_deg;      /**< Total phase shift contribution (deg) */
    double quality_factor;       /**< Q = 1/(2ζ) for complex poles */
} bode_corner_t;

/**
 * @brief Bode plot data structure.
 *
 * A Bode plot consists of two graphs:
 * 1. Magnitude plot: 20·log₁₀|H(jω)| vs. log₁₀(ω)
 * 2. Phase plot: ∠H(jω) vs. log₁₀(ω)
 *
 * L6 Canonical Problem: Bode plot construction from a transfer function
 * is a fundamental skill in circuit analysis. The asymptotic approximation
 * allows quick hand-drawn sketches by summing the contributions of
 * individual poles and zeros.
 */
typedef struct {
    freq_response_t magnitude_response;  /**< Magnitude frequency response */
    freq_response_t phase_response;      /**< Phase frequency response */
    bode_corner_t *corners;              /**< Array of corner points */
    size_t num_corners;                  /**< Number of corner points */
    double dc_gain_db;                   /**< DC gain in dB (ω→0) */
    double hf_slope_db_dec;              /**< High-frequency asymptotic slope */
} bode_plot_t;

/* ============================================================================
 * L1 DEFINITIONS: Resonance parameters
 * ============================================================================ */

/**
 * @brief Resonance analysis result.
 *
 * L1 Definition: Resonance occurs when the capacitive and inductive
 * reactances cancel: X_C = X_L, i.e., 1/(ω₀C) = ω₀L.
 *
 * Resonance frequency (undamped natural frequency):
 *   ω₀ = 1/√(LC)    (rad/s)
 *   f₀ = 1/(2π√(LC)) (Hz)
 *
 * Quality factor Q — measures sharpness of resonance:
 *   Q = ω₀/(BW) = f₀/(f_H - f_L)
 *
 * For series RLC:  Q = ω₀L/R = 1/(ω₀RC)
 * For parallel RLC: Q = R/(ω₀L) = ω₀RC
 *
 * Bandwidth:  BW = ω₀/Q = f₀/Q
 * Half-power (-3 dB) frequencies: f_L, f_H
 *   f_L = f₀·(√(1+1/(4Q²)) - 1/(2Q))
 *   f_H = f₀·(√(1+1/(4Q²)) + 1/(2Q))
 * For Q >> 1:  f_L ≈ f₀ - BW/2,  f_H ≈ f₀ + BW/2
 *
 * Damping factor ζ:  ζ = 1/(2Q)
 * For ζ < 1: underdamped (oscillatory, peaking at resonance)
 * For ζ = 1: critically damped (no overshoot, fastest settling)
 * For ζ > 1: overdamped (sluggish, no oscillation)
 *
 * L4 Fundamental Law: The resonance condition ω₀ = 1/√(LC) follows
 * directly from the impedance of L and C being equal in magnitude
 * and opposite in sign at resonance.
 */
typedef struct {
    double resonant_freq_hz;     /**< Resonance frequency f₀ (Hz) */
    double resonant_freq_rad;    /**< Resonance frequency ω₀ (rad/s) */
    double quality_factor;       /**< Quality factor Q */
    double bandwidth_hz;         /**< 3-dB bandwidth BW (Hz) */
    double half_power_low_hz;    /**< Lower half-power frequency f_L (Hz) */
    double half_power_high_hz;   /**< Upper half-power frequency f_H (Hz) */
    double damping_factor;       /**< Damping factor ζ = 1/(2Q) */
    double peak_magnitude;       /**< Magnitude at resonance |H(jω₀)| */
    double peak_magnitude_db;    /**< Peak magnitude in dB */
    double impedance_at_res;     /**< Impedance magnitude at resonance (Ω) */
} resonance_result_t;

/**
 * @brief Resonance circuit topology.
 *
 * L2 Concept: The behavior at resonance depends critically on whether
 * L and C are in series or parallel:
 *
 * - SERIES RLC: At resonance, Z = R (minimum impedance).
 *   Current is maximum. Voltage across L and C are equal in magnitude
 *   but 180° out of phase, and Q times the source voltage.
 *   Used in: series-tuned RF amplifiers, notch filters.
 *
 * - PARALLEL RLC: At resonance, Z = R (maximum impedance).
 *   Current is minimum. Circulating current between L and C is Q
 *   times the line current.
 *   Used in: parallel-tuned RF amplifiers, bandpass filters.
 */
typedef enum {
    RESONANCE_SERIES_RLC,    /**< Series RLC resonance */
    RESONANCE_PARALLEL_RLC,  /**< Parallel RLC resonance */
    RESONANCE_GENERAL        /**< General second-order resonance */
} resonance_topology_t;

/* ============================================================================
 * L2 DEFINITIONS: Network characterization in frequency domain
 * ============================================================================ */

/**
 * @brief Driving-point impedance function Z(s) parameters.
 *
 * L2 Concept: A driving-point impedance Z(s) = V(s)/I(s) at a port
 * is a positive real function — it maps the right half s-plane to
 * the right half Z-plane. This property is both necessary and
 * sufficient for a rational function to be realizable as a passive
 * network of R, L, C components (Brune's theorem, 1931).
 *
 * Properties of positive real functions:
 * 1. Z(s) is real when s is real
 * 2. Re{Z(s)} ≥ 0 when Re{s} ≥ 0
 * 3. Poles on the jω-axis must be simple with positive real residues
 * 4. The degree difference between numerator and denominator ≤ 1
 */
typedef struct {
    size_t num_poles;            /**< Number of poles */
    size_t num_zeros;            /**< Number of zeros */
    double *poles_real;          /**< Real parts of poles */
    double *poles_imag;          /**< Imaginary parts of poles */
    double *zeros_real;          /**< Real parts of zeros */
    double *zeros_imag;          /**< Imaginary parts of zeros */
    double dc_resistance;        /**< DC resistance Z(0) (Ω) */
    double hf_behavior;          /**< High-frequency behavior type */
    int is_positive_real;        /**< 1 if Z(s) is positive real */
} impedance_function_t;

/* ============================================================================
 * L1 DEFINITIONS: Stability analysis types
 * ============================================================================ */

/**
 * @brief Stability margin data.
 *
 * L1 Definition: Stability margins quantify how far a system is from
 * instability. They are measured from the open-loop frequency response.
 *
 * Gain Margin (GM): The amount of gain increase required to make the
 * system marginally stable. Measured at the phase crossover frequency
 * ω_pc where ∠L(jω_pc) = -180° (or -π rad):
 *   GM = 1/|L(jω_pc)| (linear)
 *   GM_dB = -20·log₁₀|L(jω_pc)| (dB)
 *   For stability: GM_dB > 0 (|L(jω_pc)| < 1 at -180°)
 *
 * Phase Margin (PM): The additional phase lag required to bring the
 * system to marginal stability. Measured at the gain crossover
 * frequency ω_gc where |L(jω_gc)| = 1 (0 dB):
 *   PM = 180° + ∠L(jω_gc)
 *   For stability: PM > 0°
 *   Typical design target: PM = 45° to 60°
 *
 * L4: Nyquist Stability Criterion — A closed-loop system is stable
 * iff the Nyquist plot of L(jω) encircles the point (-1, j0) exactly
 * P times counterclockwise, where P is the number of open-loop
 * unstable poles.
 *
 * Course: MIT 6.003, Stanford EE102A, Berkeley EE16B
 * Reference: Ogata, "Modern Control Engineering" (2010), Ch. 8
 */
typedef struct {
    double gain_margin_db;       /**< Gain margin (dB) */
    double gain_margin_linear;   /**< Gain margin (linear) */
    double phase_margin_deg;     /**< Phase margin (degrees) */
    double phase_margin_rad;     /**< Phase margin (radians) */
    double gain_crossover_hz;    /**< Gain crossover frequency ω_gc (Hz) */
    double phase_crossover_hz;   /**< Phase crossover frequency ω_pc (Hz) */
    int is_stable;               /**< 1 = stable, 0 = unstable */
} stability_margin_t;

/**
 * @brief Routh-Hurwitz table for polynomial stability analysis.
 *
 * L4 Fundamental Law: Routh-Hurwitz Stability Criterion (Routh, 1874;
 * Hurwitz, 1895) — A polynomial a₀sⁿ + a₁sⁿ⁻¹ + ... + aₙ has all
 * roots with negative real parts iff:
 * 1. All coefficients aᵢ have the same sign (necessary condition)
 * 2. All elements in the first column of the Routh array are positive
 *    (necessary and sufficient condition)
 *
 * The number of sign changes in the first column equals the number
 * of roots with positive real parts (unstable poles).
 *
 * For a second-order polynomial a₀s² + a₁s + a₂:
 *   Stability requires: a₀>0, a₁>0, a₂>0 (all same sign)
 *
 * For a third-order polynomial a₀s³ + a₁s² + a₂s + a₃:
 *   Stability requires: a₀>0, a₁>0, a₂>0, a₃>0, and a₁a₂ > a₀a₃
 *
 * This criterion is algebraic (no root-finding needed), making it
 * essential for symbolic stability analysis in circuit design.
 *
 * Course: MIT 6.003, Stanford EE102A
 * Reference: Dorf & Bishop, "Modern Control Systems" (2016), Ch. 6
 */
typedef struct {
    size_t order;                /**< Polynomial order */
    double *first_column;        /**< First column of Routh array */
    size_t sign_changes;         /**< Number of sign changes (unstable poles) */
    int is_stable;               /**< 1 = all roots in LHP, 0 = unstable */
    int is_marginally_stable;    /**< 1 = poles on jω-axis, 0 = no */
} routh_hurwitz_t;

/* ============================================================================
 * CORE API: Frequency response computation
 * ============================================================================ */

/**
 * @brief Allocate and initialize a frequency response data set.
 *
 * @param num_points  Number of frequency points
 * @param f_start     Start frequency (Hz)
 * @param f_end       End frequency (Hz)
 * @param scale       Frequency scale type
 * @return            Pointer to allocated freq_response_t, or NULL on error
 *
 * Complexity: O(n) where n = num_points (memory allocation)
 * Theorem: Nyquist-Shannon sampling in reverse — the frequency spacing
 * must be fine enough to capture the narrowest feature of |H(jω)|.
 */
freq_response_t *freq_response_alloc(size_t num_points, double f_start,
                                      double f_end, freq_scale_t scale);

/**
 * @brief Free a frequency response data set.
 * @param resp  Pointer to freq_response_t to free (safe to pass NULL)
 */
void freq_response_free(freq_response_t *resp);

/**
 * @brief Generate logarithmically-spaced frequency points.
 *
 * Creates an array of frequencies uniformly spaced on a log₁₀ scale,
 * which is the standard for Bode plots. For N points per decade:
 *   f[i] = f_start × 10^{i/N_points_per_decade}
 *
 * @param f_start       Start frequency (Hz), must be > 0
 * @param f_end         End frequency (Hz), must be > f_start
 * @param pts_per_dec   Points per decade (typically 10-100)
 * @param num_points    Output: number of points generated
 * @return              Array of frequencies (Hz), caller must free()
 *
 * L5 Algorithm: Logarithmic frequency grid generation for Bode plot
 * construction. The number of decades is log₁₀(f_end/f_start).
 *
 * Complexity: O(n) where n is the number of frequency points.
 * Usage: Standard for Bode plot frequency axis generation.
 */
double *freq_logspace(double f_start, double f_end, int pts_per_dec,
                       size_t *num_points);

/**
 * @brief Convert between linear frequency f (Hz) and angular frequency ω (rad/s).
 *
 * ω = 2πf
 * f = ω/(2π)
 *
 * @param freq_hz  Frequency in Hz, or pass 0 to convert from rad/s
 * @param omega    Angular frequency in rad/s
 * @return         Converted value based on which input is nonzero
 *
 * L3: The relationship ω = 2πf is fundamental to AC circuit analysis.
 * It bridges the time-domain (ω used in phasors and Laplace transforms)
 * and practical frequency measurements (Hz used in instrumentation).
 */
double freq_hz_to_rad(double freq_hz);
double freq_rad_to_hz(double omega);

/**
 * @brief Compute the magnitude in dB from a linear magnitude.
 *
 * dB = 20·log₁₀(|H|)
 *
 * For |H| ≤ 0, returns -INFINITY (dB of zero magnitude is -∞ dB).
 * For power ratios: dB = 10·log₁₀(P₂/P₁)
 *
 * L1 Definition: The decibel is a logarithmic unit expressing the
 * ratio of two quantities. In electronics, 0 dB = unity gain,
 * 20 dB = gain of 10, -20 dB = attenuation of 10.
 *
 * Common dB values to remember:
 *   +3 dB   ≈ ×1.414 (√2, half-power point)
 *   +6 dB   ≈ ×2
 *   +20 dB  = ×10
 *   +40 dB  = ×100
 *   -3 dB   ≈ ×0.707 (1/√2, -3dB cutoff)
 */
double magnitude_to_db(double magnitude);

/**
 * @brief Compute the linear magnitude from dB.
 *
 * |H| = 10^{dB/20}
 */
double db_to_magnitude(double db);

/**
 * @brief Compute phase in degrees from a complex value.
 *
 * ∠H = atan2(imag, real) × 180/π
 *
 * Uses atan2 for correct quadrant determination.
 * Range: (-180°, 180°] or [-180°, 180°)
 *
 * L3: Phase unwrapping may be needed for continuous phase plots,
 * as atan2 produces discontinuities at ±180°. Phase unwrapping
 * adds/subtracts multiples of 360° to ensure continuity.
 */
double phase_degrees(double real, double imag);

/**
 * @brief Unwrap phase to remove ±180° discontinuities.
 *
 * L5 Algorithm: Phase unwrapping is essential for computing
 * group delay τ_g = -dφ/dω. The raw atan2 output has 360° jumps
 * that must be removed before differentiation.
 *
 * @param phase_deg  Array of wrapped phase values (degrees)
 * @param n          Number of points
 * @param tolerance  Jump threshold (typically 180°)
 * @return           Array of unwrapped phase (caller must free)
 *
 * Complexity: O(n)
 * Reference: Oppenheim & Schafer, "Discrete-Time Signal Processing" (2010) Ch. 10
 */
double *phase_unwrap(const double *phase_deg, size_t n, double tolerance);

/**
 * @brief Compute group delay from unwrapped phase data.
 *
 * τ_g(ω_k) = -Δφ/Δω = -(φ_{k+1} - φ_{k-1})/(ω_{k+1} - ω_{k-1})
 *
 * Uses central difference for interior points, forward/backward
 * difference at endpoints.
 *
 * L1 Definition: Group delay measures the time delay experienced by
 * the envelope of a narrow-band signal. Constant group delay (linear
 * phase) preserves the waveform shape — crucial for pulse transmission.
 *
 * @param phase_rad  Unwrapped phase in radians
 * @param omega      Angular frequencies (rad/s)
 * @param n          Number of points
 * @return           Group delay array (seconds), caller must free
 *
 * Complexity: O(n)
 */
double *group_delay(const double *phase_rad, const double *omega, size_t n);

/* ============================================================================
 * FREQUENCY SWEEP OPERATIONS
 * ============================================================================ */

/**
 * @brief Find the -3 dB cutoff frequency from a frequency response.
 *
 * Searches for the frequency where |H(jω)| drops to 1/√2 of the
 * passband value (or rises to 1/√2 of maximum for highpass).
 *
 * Uses linear interpolation between frequency points for accuracy.
 *
 * L6 Canonical Problem: Determining cutoff frequency is essential for
 * characterizing filters and amplifiers. The -3 dB point corresponds
 * to half-power: P_out = P_in/2.
 *
 * @param resp        Frequency response data
 * @param is_highpass 0 for lowpass, 1 for highpass
 * @return            Cutoff frequency in Hz, or -1 if not found
 */
double find_cutoff_freq(const freq_response_t *resp, int is_highpass);

/**
 * @brief Find the -3 dB bandwidth from a frequency response.
 *
 * Identifies f_L and f_H where gain is 3 dB below maximum for
 * bandpass/bandstop responses.
 *
 * @param resp   Frequency response data
 * @param f_low  Output: lower -3 dB frequency (Hz)
 * @param f_high Output: upper -3 dB frequency (Hz)
 * @return       Bandwidth f_H - f_L in Hz, or -1 if not found
 */
double find_bandwidth(const freq_response_t *resp,
                       double *f_low, double *f_high);

/**
 * @brief Find peak magnitude and the frequency at which it occurs.
 *
 * @param resp       Frequency response data
 * @param peak_freq  Output: frequency at peak (Hz)
 * @return           Peak magnitude (linear), or -1 if no peak
 */
double find_peak_response(const freq_response_t *resp, double *peak_freq);

/**
 * @brief Compute the roll-off rate in dB/decade over a specified range.
 *
 * Uses linear regression on the dB magnitude in the stopband.
 *
 * @param resp   Frequency response data
 * @param f_low  Lower bound of stopband (Hz)
 * @param f_high Upper bound of stopband (Hz)
 * @return       Roll-off rate in dB/decade (negative for lowpass)
 */
double compute_rolloff(const freq_response_t *resp,
                        double f_low, double f_high);

/**
 * @brief Interpolate frequency response at an arbitrary frequency.
 *
 * Uses linear interpolation in log-frequency and dB-magnitude space
 * for accurate Bode plot interpolation.
 *
 * @param resp  Frequency response data
 * @param freq  Target frequency (Hz)
 * @return      Interpolated freq_point_t, or zeroed struct if out of range
 */
freq_point_t freq_response_interpolate(const freq_response_t *resp,
                                         double freq);

#ifdef __cplusplus
}
#endif

#endif /* FREQUENCY_RESPONSE_H */
