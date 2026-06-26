/**
 * @file transfer_function.h
 * @brief Transfer function representation and analysis
 *
 * A transfer function H(s) = N(s)/D(s) completely characterizes an
 * LTI system in the s-domain (Laplace domain). For frequency response,
 * evaluate at s = jω: H(jω) = H(s)|_{s=jω}.
 *
 * Representation forms:
 * 1. Polynomial ratio: H(s) = (bₘsᵐ + ... + b₁s + b₀)/(aₙsⁿ + ... + a₁s + a₀)
 * 2. Pole-zero-gain: H(s) = K·Π(s - zᵢ)/Π(s - pⱼ)
 * 3. Partial fraction: H(s) = Σ rᵢ/(s - pᵢ) + direct term
 * 4. Biquad cascade: H(s) = Π (b₂ᵢs² + b₁ᵢs + b₀ᵢ)/(a₂ᵢs² + a₁ᵢs + a₀ᵢ)
 *
 * Reference:
 * - Ogata, "Modern Control Engineering" (2010), Ch. 3
 * - Oppenheim & Willsky, "Signals and Systems" (1997), Ch. 9
 * - Sedra & Smith, "Microelectronic Circuits" (2020), Ch. 8
 *
 * Course: MIT 6.003, Stanford EE102A, Berkeley EE16B, ETH 227-0427
 */

#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

/* Need freq_response_t from frequency_response.h */
#include "frequency_response.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1 DEFINITIONS: Transfer function representation forms
 * ============================================================================ */

/**
 * @brief Transfer function representation form.
 *
 * L2 Concept: A transfer function can be represented in several
 * mathematically equivalent forms, each useful for different
 * analysis tasks:
 *
 * - POLYNOMIAL: Direct coefficient form. Best for evaluation.
 *   H(s) = (bₘsᵐ + ... + b₁s + b₀)/(aₙsⁿ + ... + a₁s + a₀)
 *
 * - POLE_ZERO_GAIN: Factored form reveals poles and zeros directly.
 *   Essential for stability analysis and Bode plot construction.
 *   H(s) = K·Π(s - zᵢ)/Π(s - pⱼ)
 *
 * - PARTIAL_FRACTION: Sum of first-order terms.
 *   Essential for inverse Laplace transform and time-domain response.
 *   H(s) = Σ rᵢ/(s - pᵢ) (for simple poles)
 *
 * - BIQUAD_CASCADE: Product of second-order sections.
 *   Standard form for digital filter implementation (reduces
 *   numerical sensitivity compared to direct form).
 */
typedef enum {
    TF_FORM_POLYNOMIAL,       /**< H(s) = N(s)/D(s) polynomial ratio */
    TF_FORM_POLE_ZERO_GAIN,   /**< H(s) = K·Π(s-zᵢ)/Π(s-pⱼ) */
    TF_FORM_PARTIAL_FRACTION, /**< H(s) = Σ rᵢ/(s-pᵢ) */
    TF_FORM_BIQUAD_CASCADE    /**< Cascade of 2nd-order sections */
} tf_form_t;

/* ============================================================================
 * L1 DEFINITIONS: Polynomial ratio representation
 * ============================================================================ */

/**
 * @brief Transfer function in polynomial ratio form.
 *
 * H(s) = N(s)/D(s) where:
 *   N(s) = b₀ + b₁s + b₂s² + ... + b_ms^m  (numerator)
 *   D(s) = a₀ + a₁s + a₂s² + ... + a_ns^n  (denominator)
 *
 * For a physically realizable (causal) system: m ≤ n.
 * Coefficients are stored in ascending power order:
 *   b[0] = b₀, b[1] = b₁, ..., b[m] = b_m
 *   a[0] = a₀, a[1] = a₁, ..., a[n] = a_n
 *
 * L4: For a stable system, all poles (roots of D(s)) must have
 * negative real parts, i.e., lie in the left half s-plane.
 */
typedef struct {
    double *num;          /**< Numerator coefficients b₀, b₁, ..., b_m */
    double *den;          /**< Denominator coefficients a₀, a₁, ..., a_n */
    size_t num_order;     /**< Numerator order m */
    size_t den_order;     /**< Denominator order n */
} tf_polynomial_t;

/* ============================================================================
 * L1 DEFINITIONS: Pole-zero-gain representation
 * ============================================================================ */

/**
 * @brief Transfer function in pole-zero-gain form.
 *
 * H(s) = K · Π_{i=1}^m (s - z_i) / Π_{j=1}^n (s - p_j)
 *
 * where:
 *   K   = gain factor (leading coefficient ratio: b_m/a_n)
 *   z_i = zeros (numerator roots), may be complex
 *   p_j = poles (denominator roots), may be complex
 *
 * L2 Concept: Poles and zeros are the most intuitive representation
 * for understanding frequency response behavior:
 * - Each pole contributes -20 dB/decade and -90° phase at high frequencies
 * - Each zero contributes +20 dB/decade and +90° phase at high frequencies
 * - Pole locations determine stability (must be in LHP)
 * - Zero locations determine minimum-phase vs. non-minimum-phase behavior
 *
 * Complex poles/zeros always appear in conjugate pairs for real-coefficient
 * transfer functions (real systems).
 */
typedef struct {
    double _Complex *zeros;   /**< Array of zeros z_i */
    double _Complex *poles;   /**< Array of poles p_j */
    size_t num_zeros;         /**< Number of zeros m */
    size_t num_poles;         /**< Number of poles n */
    double gain;              /**< Gain factor K */
} tf_pole_zero_t;

/* ============================================================================
 * L1 DEFINITIONS: Partial fraction expansion
 * ============================================================================ */

/**
 * @brief Residue term in partial fraction expansion.
 *
 * For a simple pole at p_k with residue r_k:
 *   r_k/(s - p_k)
 *
 * For a pole of multiplicity M:
 *   r_{k,1}/(s - p_k) + r_{k,2}/(s - p_k)^2 + ... + r_{k,M}/(s - p_k)^M
 *
 * L3 Math Structure: Partial fraction expansion decomposes a rational
 * function into a sum of simpler fractions. This is the key to:
 * 1. Inverse Laplace transform (each term ↦ exponential in time domain)
 * 2. Parallel-form implementation of digital filters
 * 3. Modal analysis in state-space representations
 *
 * The residue r_k for a simple pole p_k is computed via:
 *   r_k = lim_{s→p_k} (s - p_k)·H(s)
 * or equivalently (Heaviside cover-up method):
 *   r_k = N(p_k) / D'(p_k)  where D'(s) = dD/ds
 */
typedef struct {
    double _Complex pole;       /**< Pole location p_k */
    double _Complex residue;    /**< Residue r_k */
    int multiplicity;           /**< Pole multiplicity (1 = simple) */
} tf_residue_t;

/**
 * @brief Transfer function in partial fraction form.
 *
 * H(s) = D₀ + Σ_{k=1}^P Σ_{m=1}^{M_k} r_{k,m}/(s - p_k)^m
 *
 * where D₀ is the direct feedthrough term (non-zero when m = n).
 */
typedef struct {
    tf_residue_t *residues;     /**< Array of residue terms */
    size_t num_residues;        /**< Number of residue terms */
    double _Complex direct_term; /**< Direct feedthrough D₀ */
} tf_partial_fraction_t;

/* ============================================================================
 * L1 DEFINITIONS: Biquad section
 * ============================================================================ */

/**
 * @brief Second-order (biquad) filter section.
 *
 * H(s) = (b₂s² + b₁s + b₀)/(a₂s² + a₁s + a₀)
 *
 * L5 Algorithm: Biquad sections are the fundamental building block
 * of high-order filters. Any even-order transfer function can be
 * decomposed into a cascade of biquad sections. Odd-order filters
 * include one first-order section.
 *
 * Biquad parameters can be expressed in terms of:
 * - Natural frequency ω₀ and quality factor Q:
 *   D(s) = s² + (ω₀/Q)s + ω₀²  (normalized with a₂=1)
 *
 * - Pole frequency f₀ and pole Q:
 *   f₀ = ω₀/(2π),  Q = ω₀/(2α) where poles are at -α ± jω_d
 *
 * Standard biquad topologies for active filter implementation:
 * - Sallen-Key (positive feedback, low component sensitivity)
 * - Multiple Feedback (MFB) (inverting, good for high Q)
 * - Tow-Thomas (state-variable, simultaneous LP/BP/HP outputs)
 * - KHN (Kerwin-Huelsman-Newcomb, state-variable)
 *
 * Course: Berkeley EE105, Stanford EE247
 * Reference: Sedra & Smith (2020) Ch. 16, Zverev (1967)
 */
typedef struct {
    double b2, b1, b0;         /**< Numerator coefficients */
    double a2, a1, a0;         /**< Denominator coefficients */
    double natural_freq_rad;   /**< Natural frequency ω₀ (rad/s) */
    double quality_factor;     /**< Quality factor Q */
    double dc_gain;            /**< DC gain b₀/a₀ */
    int filter_type;           /**< Resulting filter type (LP/HP/BP/BS) */
} biquad_section_t;

/**
 * @brief Cascade of biquad sections forming a high-order filter.
 *
 * H(s) = gain · Π_{i=1}^{N_sections} H_i(s)
 *
 * where each H_i(s) is a biquad_section_t.
 */
typedef struct {
    biquad_section_t *sections;  /**< Array of biquad sections */
    size_t num_sections;         /**< Number of sections */
    double overall_gain;         /**< Overall gain factor */
} biquad_cascade_t;

/* ============================================================================
 * L1 DEFINITIONS: S-parameters in frequency domain
 * ============================================================================ */

/**
 * @brief S-parameter matrix at a single frequency.
 *
 * L1 Definition: Scattering parameters (S-parameters) describe the
 * reflection and transmission behavior of a network at its ports.
 * For a 2-port network:
 *   S₁₁ = input reflection coefficient (return loss = -20log|S₁₁|)
 *   S₂₁ = forward transmission coefficient (gain/insertion loss)
 *   S₁₂ = reverse transmission coefficient (isolation)
 *   S₂₂ = output reflection coefficient
 *
 * In frequency response analysis, S-parameters characterize:
 * - Filter insertion loss: IL = -20·log₁₀|S₂₁| dB
 * - Return loss: RL = -20·log₁₀|S₁₁| dB
 * - Voltage Standing Wave Ratio: VSWR = (1+|S₁₁|)/(1-|S₁₁|)
 *
 * L2 Concept: S-parameters are preferred at high frequencies (RF/microwave)
 * because they are based on traveling waves (power), whereas Z, Y, H
 * parameters require open/short circuit terminations that are difficult
 * to realize at high frequencies.
 *
 * Reference: Pozar, "Microwave Engineering" (2012), Ch. 4
 * Course: ETH 227-0455, Georgia Tech ECE 6350
 */
typedef struct {
    double _Complex s11;       /**< Input reflection coefficient */
    double _Complex s21;       /**< Forward transmission */
    double _Complex s12;       /**< Reverse transmission */
    double _Complex s22;       /**< Output reflection coefficient */
    double frequency_hz;       /**< Frequency at which S-params are measured */
} s_params_2port_t;

/* ============================================================================
 * CORE API: Transfer function creation and evaluation
 * ============================================================================ */

/**
 * @brief Create a transfer function in polynomial form.
 *
 * @param num         Numerator coefficients (ascending powers)
 * @param num_order   Numerator order m (degree = m, length = m+1)
 * @param den         Denominator coefficients (ascending powers)
 * @param den_order   Denominator order n (degree = n, length = n+1)
 * @return            Allocated tf_polynomial_t, or NULL on error
 *
 * L4: For a physically realizable system, den_order ≥ num_order.
 * The function validates this constraint.
 *
 * Complexity: O(n) for coefficient copying
 */
tf_polynomial_t *tf_polynomial_create(const double *num, size_t num_order,
                                       const double *den, size_t den_order);

/**
 * @brief Free a polynomial transfer function.
 */
void tf_polynomial_free(tf_polynomial_t *tf);

/**
 * @brief Evaluate H(s) at a complex frequency s = σ + jω.
 *
 * Evaluates both numerator and denominator polynomials using
 * Horner's method for numerical stability:
 *   N(s) = b₀ + s·(b₁ + s·(b₂ + ... + s·b_m))
 *
 * @param tf  Transfer function
 * @param s   Complex frequency s = σ + jω
 * @return    Complex value H(s)
 *
 * L3: Horner's method — evaluates a polynomial of degree n with
 * only n multiplications and n additions, reducing both operations
 * count and floating-point error accumulation.
 *
 * Complexity: O(max(m,n))
 */
double _Complex tf_evaluate(const tf_polynomial_t *tf, double _Complex s);

/**
 * @brief Evaluate frequency response H(jω) at angular frequency ω.
 *
 * Convenience wrapper: tf_evaluate(tf, I * omega)
 *
 * @param tf     Transfer function
 * @param omega  Angular frequency (rad/s)
 * @return       Complex frequency response H(jω)
 */
double _Complex tf_evaluate_freq(const tf_polynomial_t *tf, double omega);

/**
 * @brief Compute magnitude response at frequency f.
 *
 * @param tf     Transfer function
 * @param freq   Frequency (Hz)
 * @return       |H(j2πf)| (linear magnitude)
 */
double tf_magnitude_at(const tf_polynomial_t *tf, double freq);

/**
 * @brief Compute phase response at frequency f.
 *
 * @param tf     Transfer function
 * @param freq   Frequency (Hz)
 * @return       ∠H(j2πf) (degrees)
 */
double tf_phase_at(const tf_polynomial_t *tf, double freq);

/**
 * @brief Generate complete frequency response from transfer function.
 *
 * Evaluates H(jω) at all frequency points in a logarithmic sweep.
 *
 * @param tf           Transfer function
 * @param f_start      Start frequency (Hz)
 * @param f_end        End frequency (Hz)
 * @param pts_per_dec  Points per decade
 * @return             freq_response_t with computed data, or NULL
 *
 * Complexity: O(N·max(m,n)) where N = number of frequency points
 */
freq_response_t *tf_compute_freq_response(const tf_polynomial_t *tf,
                                            double f_start, double f_end,
                                            int pts_per_dec);

/**
 * @brief Normalize transfer function so aₙ = 1 (monic denominator).
 *
 * Divides all coefficients by a_n. This is standard form for
 * pole-zero analysis.
 *
 * @param tf  Transfer function to normalize (modified in place)
 *
 * L3: Monic polynomial form simplifies root-finding and pole-zero
 * analysis by removing one degree of freedom.
 */
void tf_normalize(tf_polynomial_t *tf);

/**
 * @brief Compute the DC gain H(0) = b₀/a₀.
 *
 * @param tf  Transfer function
 * @return    DC gain (linear), or INFINITY if a₀ = 0 (pole at origin)
 */
double tf_dc_gain(const tf_polynomial_t *tf);

/**
 * @brief Compute the high-frequency asymptotic gain.
 *
 * For m < n: lim_{ω→∞} |H(jω)| = 0 (roll-off)
 * For m = n: lim_{ω→∞} H(jω) = b_m/a_n (constant)
 * For m > n: diverges (non-causal, physically unrealizable)
 *
 * @param tf  Transfer function
 * @return    High-frequency gain, INFINITY if m > n
 */
double tf_hf_gain(const tf_polynomial_t *tf);

/* ============================================================================
 * POLE-ZERO ANALYSIS
 * ============================================================================ */

/**
 * @brief Convert polynomial form to pole-zero-gain form.
 *
 * Finds roots of numerator and denominator polynomials using
 * the Jenkins-Traub algorithm (via companion matrix eigenvalue method).
 *
 * @param tf  Polynomial transfer function
 * @return    Pole-zero-gain representation, or NULL on error
 *
 * L5 Algorithm: Root-finding for polynomials. For degree ≤ 4,
 * analytical formulas exist. For higher degrees, numerical methods
 * (eigenvalue of companion matrix, Jenkins-Traub, Laguerre) are used.
 * The companion matrix method is robust and leverages LAPACK-quality
 * eigenvalue algorithms.
 *
 * Complexity: O(n³) for eigenvalue computation (companion matrix)
 * Reference: Press et al., "Numerical Recipes" (2007), Ch. 9.5
 */
tf_pole_zero_t *tf_to_pole_zero(const tf_polynomial_t *tf);

/**
 * @brief Convert pole-zero-gain form to polynomial form.
 *
 * Expands Π(s - pⱼ) and Π(s - zᵢ) into polynomial coefficients.
 *
 * @param pz  Pole-zero-gain transfer function
 * @return    Polynomial transfer function, or NULL
 *
 * Complexity: O(n²) for polynomial expansion
 */
tf_polynomial_t *tf_from_pole_zero(const tf_pole_zero_t *pz);

/**
 * @brief Free a pole-zero-gain transfer function.
 */
void tf_pole_zero_free(tf_pole_zero_t *pz);

/**
 * @brief Determine if a transfer function is minimum-phase.
 *
 * A minimum-phase system has all poles AND all zeros in the left
 * half s-plane (LHP). Minimum-phase systems have the property that
 * for a given magnitude response, the phase lag is minimized.
 *
 * Conversely, non-minimum-phase (NMP) systems have one or more zeros
 * in the right half-plane (RHP), resulting in excess phase lag.
 *
 * L2 Concept: Minimum-phase property is important because:
 * 1. For minimum-phase systems, magnitude uniquely determines phase
 *    (via Bode's gain-phase relation / Hilbert transform)
 * 2. NMP zeros impose fundamental limitations on feedback control
 *    (waterbed effect — reducing sensitivity at one frequency
 *     increases it at another)
 *
 * @param pz  Pole-zero transfer function
 * @return    1 = minimum-phase, 0 = non-minimum-phase
 *
 * L4: Bode's Integral Theorem — For a stable minimum-phase LTI system,
 * the integral of ln|S(jω)| over all frequencies is zero, where S is
 * the sensitivity function. This imposes a conservation law on feedback.
 */
int tf_is_minimum_phase(const tf_pole_zero_t *pz);

/**
 * @brief Compute pole and zero frequencies in Hz.
 *
 * For a pole at s = σ + jω: the natural frequency is |s| = √(σ²+ω²).
 * For real poles: freq = |σ|/(2π) Hz.
 * For complex conjugate poles: freq = √(σ²+ω²)/(2π) Hz.
 *
 * @param pz         Pole-zero transfer function
 * @param pole_freqs Output: array of pole frequencies (Hz), caller frees
 * @param zero_freqs Output: array of zero frequencies (Hz), caller frees
 */
void tf_pole_zero_frequencies(const tf_pole_zero_t *pz,
                               double **pole_freqs, double **zero_freqs);

/* ============================================================================
 * PARTIAL FRACTION EXPANSION
 * ============================================================================ */

/**
 * @brief Compute partial fraction expansion of a transfer function.
 *
 * Decomposes H(s) = N(s)/D(s) into:
 *   H(s) = D₀ + Σ r_k/(s - p_k)
 *
 * Uses the residue formula for simple poles:
 *   r_k = N(p_k) / D'(p_k)
 * where D'(s) is the derivative of D(s).
 *
 * @param tf  Polynomial transfer function
 * @return    Partial fraction expansion, or NULL on error
 *
 * L3 Math Structure: The residue theorem from complex analysis.
 * For a rational function with simple poles:
 *   f(s) = Σ Res(f, p_k)/(s - p_k)
 * where Res(f, p_k) = lim_{s→p_k} (s-p_k)·f(s)
 *
 * L5 Algorithm: For repeated poles, the residues are computed via:
 *   r_{k,m} = 1/(M-m)! · d^{M-m}/ds^{M-m} [(s-p_k)^M·H(s)]|_{s=p_k}
 *
 * Complexity: O(n²) due to polynomial evaluation at each pole
 * Reference: Proakis & Manolakis, "Digital Signal Processing" (2007), Ch. 6
 */
tf_partial_fraction_t *tf_partial_fraction(const tf_polynomial_t *tf);

/**
 * @brief Free a partial fraction expansion.
 */
void tf_partial_fraction_free(tf_partial_fraction_t *pf);

/* ============================================================================
 * BIQUAD DECOMPOSITION
 * ============================================================================ */

/**
 * @brief Decompose a transfer function into biquad cascade.
 *
 * Pairs complex conjugate poles into second-order sections.
 * Real poles are paired to minimize peak gain in each section.
 *
 * Standard pairing strategy:
 * 1. Pair the pole with highest Q with the zero closest to it
 * 2. Pair the pole with next highest Q with the next closest zero
 * 3. Continue until all poles are paired
 * 4. Remaining real poles form first-order sections (or are paired)
 *
 * @param tf  Polynomial transfer function
 * @return    Biquad cascade, or NULL on error
 *
 * L5 Algorithm: Pole-zero pairing for optimal dynamic range.
 * Improper pairing can cause internal overflow even when overall
 * transfer function has unity gain. The L₂ scaling rule ensures
 * no internal node saturates.
 *
 * Reference: Jackson, "Digital Filters and Signal Processing" (1996), Ch. 9
 * Course: MIT 6.003, Stanford EE102A
 */
biquad_cascade_t *tf_to_biquad_cascade(const tf_polynomial_t *tf);

/**
 * @brief Free a biquad cascade.
 */
void biquad_cascade_free(biquad_cascade_t *cascade);

/**
 * @brief Evaluate a biquad cascade at complex frequency s.
 *
 * @param cascade  Biquad cascade
 * @param s        Complex frequency
 * @return         H(s) = gain · Π H_i(s)
 *
 * Complexity: O(N) where N = number of sections
 */
double _Complex biquad_cascade_evaluate(const biquad_cascade_t *cascade,
                                          double _Complex s);

/**
 * @brief Create a single biquad section from specifications.
 *
 * Given filter type, natural frequency ω₀, quality factor Q, and DC gain,
 * computes the coefficients for a normalized biquad (a₂ = 1):
 *
 * Lowpass:   H(s) = (b₀)/(s² + (ω₀/Q)s + ω₀²)
 * Highpass:  H(s) = (b₂s²)/(s² + (ω₀/Q)s + ω₀²)
 * Bandpass:  H(s) = (b₁s)/(s² + (ω₀/Q)s + ω₀²)
 * Bandstop:  H(s) = b₂(s² + ω_z²)/(s² + (ω₀/Q)s + ω₀²)
 *
 * @param type       Filter type
 * @param omega0     Natural frequency (rad/s)
 * @param q_factor   Quality factor Q
 * @param dc_gain    DC gain for LP, HF gain for HP, peak gain for BP
 * @return           Biquad section
 */
biquad_section_t biquad_create(int type, double omega0,
                                double q_factor, double dc_gain);

/* ============================================================================
 * TRANSFER FUNCTION ARITHMETIC
 * ============================================================================ */

/**
 * @brief Multiply two transfer functions: H(s) = H₁(s)·H₂(s).
 *
 * Series/cascade connection: H_total = H₁ × H₂.
 * Polynomial multiplication (convolution of coefficients).
 *
 * @param tf1  First transfer function
 * @param tf2  Second transfer function
 * @return     Product transfer function, or NULL
 *
 * Complexity: O(m₁m₂ + n₁n₂) where mᵢ,nᵢ are numerator/denominator orders
 */
tf_polynomial_t *tf_multiply(const tf_polynomial_t *tf1,
                              const tf_polynomial_t *tf2);

/**
 * @brief Add two transfer functions: H(s) = H₁(s) + H₂(s).
 *
 * Parallel connection: H_total = H₁ + H₂.
 * Requires common denominator computation.
 *
 * @param tf1  First transfer function
 * @param tf2  Second transfer function
 * @return     Sum transfer function, or NULL
 *
 * Complexity: O((n₁+n₂)²) due to polynomial multiplication for
 * common denominator
 */
tf_polynomial_t *tf_add(const tf_polynomial_t *tf1,
                         const tf_polynomial_t *tf2);

/**
 * @brief Compute feedback connection: H(s) = H₁/(1 + H₁·H₂).
 *
 * Negative feedback with H₁ in forward path, H₂ in feedback path.
 * For unity feedback, H₂(s) = 1.
 * For positive feedback, H(s) = H₁/(1 - H₁·H₂).
 *
 * @param forward     Forward path H₁(s)
 * @param feedback    Feedback path H₂(s)
 * @param is_positive 0 = negative feedback, 1 = positive feedback
 * @return            Closed-loop transfer function, or NULL
 *
 * L2 Concept: Negative feedback is the most important concept in
 * analog circuit design. It:
 * - Reduces sensitivity to parameter variations
 * - Extends bandwidth (gain-bandwidth product is constant)
 * - Reduces nonlinear distortion
 * - Modifies input/output impedance
 * - BUT may cause instability if not properly compensated
 *
 * L4: Black's Formula (Harold Black, 1927):
 *   H_closed = A/(1 + Aβ) for negative feedback
 * where A = forward gain, β = feedback factor.
 * Loop gain L = Aβ determines stability.
 *
 * Course: Berkeley EE105, Stanford EE247, MIT 6.003
 * Reference: Sedra & Smith (2020), Ch. 10
 */
tf_polynomial_t *tf_feedback(const tf_polynomial_t *forward,
                              const tf_polynomial_t *feedback,
                              int is_positive);

#ifdef __cplusplus
}
#endif

#endif /* TRANSFER_FUNCTION_H */
