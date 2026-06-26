/**
 * @file filter_design.h
 * @brief Analog filter design and approximation
 *
 * This header provides the complete analog filter design API:
 * 1. Lowpass prototype generation (normalized ω_c = 1 rad/s)
 * 2. Frequency transformation (LP→HP, LP→BP, LP→BS)
 * 3. Impedance denormalization to real component values
 * 4. Active filter realization (Sallen-Key, MFB, Tow-Thomas)
 *
 * Theory: The design of analog filters begins with a normalized
 * lowpass prototype (cutoff at 1 rad/s, terminated in 1 Ω). This
 * prototype is then transformed to the desired filter type, cutoff
 * frequency, and impedance level.
 *
 * Reference:
 * - Zverev, "Handbook of Filter Synthesis" (1967)
 * - Williams & Taylor, "Electronic Filter Design Handbook" (2006)
 * - Van Valkenburg, "Analog Filter Design" (1982)
 * - Sedra & Smith, "Microelectronic Circuits" (2020), Ch. 16
 *
 * Course: Berkeley EE105, Stanford EE247, ETH 227-0455, MIT 6.003
 */

#ifndef FILTER_DESIGN_H
#define FILTER_DESIGN_H

#include <stddef.h>
#include "frequency_response.h"
#include "transfer_function.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * L1 DEFINITIONS: Filter specification
 * ============================================================================ */

/**
 * @brief Filter design specification.
 *
 * Complete specification for designing an analog lowpass filter
 * using the insertion loss method.
 *
 * Parameters follow the standard filter specification template:
 * - Passband: 0 ≤ f ≤ f_pass, with attenuation ≤ A_pass dB
 *   (A_pass is typically 0.1-3 dB, "passband ripple")
 * - Stopband: f ≥ f_stop, with attenuation ≥ A_stop dB
 *   (A_stop is typically 20-80 dB, "stopband attenuation")
 * - Transition band: f_pass < f < f_stop
 *
 * The filter order N is chosen to meet both constraints with the
 * selected approximation type.
 *
 * L4 Fundamental Law: The sharper the cutoff (smaller transition band),
 * the higher the required filter order. There is a fundamental trade-off
 * between filter complexity (order) and selectivity.
 */
typedef struct {
    filter_approx_t approx;   /**< Approximation type */
    double f_pass;            /**< Passband edge frequency (Hz) */
    double f_stop;            /**< Stopband edge frequency (Hz) */
    double a_pass;            /**< Maximum passband attenuation (dB) */
    double a_stop;            /**< Minimum stopband attenuation (dB) */
    double dc_gain;           /**< Desired DC (passband) gain (linear) */
    double source_impedance;  /**< Source impedance (Ω), for passive filters */
    double load_impedance;    /**< Load impedance (Ω), for passive filters */
} filter_spec_t;

/**
 * @brief Filter design result with component values.
 *
 * For passive LC ladder filters, provides normalized element values.
 * For active RC filters, provides resistor and capacitor values.
 *
 * Passive LC ladder (doubly-terminated):
 *   Elements are alternating series L and shunt C, normalized to
 *   1 Ω termination and 1 rad/s cutoff. g₀ = source resistance = 1,
 *   g_{N+1} = load resistance = 1 (or calculated).
 *
 *   Denormalization:
 *     L_actual = (g_k · R₀)/ω_c    (for series inductor g_k)
 *     C_actual = g_k/(R₀ · ω_c)    (for shunt capacitor g_k)
 *     R_actual = g_k · R₀          (for termination resistor)
 *
 *   where R₀ is the impedance scaling factor and ω_c is the
 *   frequency scaling factor (2πf_c).
 *
 * L5 Algorithm: Filter synthesis — the process of converting a
 * transfer function into a physical circuit. For LC ladder filters,
 * element values are found by expanding the input impedance Z_in(s)
 * as a continued fraction (Cauer synthesis).
 */
typedef struct {
    double *g_values;        /**< Normalized g-values (prototype elements) */
    size_t n;                /**< Filter order (= number of reactive elements) */
    double *L_values;        /**< Actual inductor values (H) */
    double *C_values;        /**< Actual capacitor values (H) */
    double *R_values;        /**< Actual resistor values (Ω) */
    filter_type_t type;      /**< Resulting filter type */
    double cutoff_hz;        /**< Actual cutoff frequency (Hz) */
} filter_design_t;

/* ============================================================================
 * CORE API: Filter order computation
 * ============================================================================ */

/**
 * @brief Compute required Butterworth filter order.
 *
 * For Butterworth (maximally flat) lowpass:
 *   |H(jω)|² = 1/(1 + (ω/ω_c)^{2N})
 *
 * Required order:
 *   N ≥ log₁₀((10^{A_stop/10} - 1)/(10^{A_pass/10} - 1)) /
 *        (2·log₁₀(ω_stop/ω_pass))
 *
 * @param spec  Filter specification
 * @return      Required filter order N (≥1)
 *
 * L5 Algorithm: Butterworth approximation — poles are equally spaced
 * on a circle of radius ω_c in the s-plane. The denominator polynomial
 * is the Butterworth polynomial B_N(s).
 *
 * B₁(s) = s + 1
 * B₂(s) = s² + √2·s + 1
 * B₃(s) = s³ + 2s² + 2s + 1 = (s+1)(s²+s+1)
 * B₄(s) = s⁴ + 2.613s³ + 3.414s² + 2.613s + 1
 *
 * Reference: Butterworth, "On the Theory of Filter Amplifiers" (1930)
 * Course: Berkeley EE105, Stanford EE247
 */
int filter_order_butterworth(const filter_spec_t *spec);

/**
 * @brief Compute required Chebyshev Type I filter order.
 *
 * For Chebyshev I (equiripple passband):
 *   |H(jω)|² = 1/(1 + ε²·T_N²(ω/ω_c))
 *
 * where ε = √(10^{A_pass/10} - 1) and T_N(x) is the Chebyshev
 * polynomial of the first kind: T_N(cos θ) = cos(Nθ).
 *
 * Required order:
 *   N ≥ acosh(√((10^{A_stop/10}-1)/(10^{A_pass/10}-1))) /
 *        acosh(ω_stop/ω_pass)
 *
 * @param spec  Filter specification
 * @return      Required filter order N (≥1)
 *
 * L5 Algorithm: Chebyshev approximation provides sharper cutoff than
 * Butterworth for the same order, at the cost of passband ripple.
 * Poles lie on an ellipse in the s-plane.
 *
 * Reference: Chebyshev, P.L. (1854)
 * Course: Berkeley EE105, Stanford EE247, ETH 227-0455
 */
int filter_order_chebyshev1(const filter_spec_t *spec);

/**
 * @brief Compute required Chebyshev Type II filter order.
 *
 * Chebyshev II (Inverse Chebyshev): monotonic passband, equiripple stopband.
 * Same order formula as Chebyshev I (the types are mathematically dual).
 *
 * @param spec  Filter specification
 * @return      Required filter order N (≥1)
 */
int filter_order_chebyshev2(const filter_spec_t *spec);

/**
 * @brief Compute required Elliptic (Cauer) filter order.
 *
 * For Elliptic (equiripple passband and stopband):
 *   |H(jω)|² = 1/(1 + ε²·R_N²(ξ, ω/ω_c))
 *
 * where R_N is the Chebyshev rational function.
 *
 * Required order (approximate):
 *   N ≥ K(1/k₁)·K'(1/k)/(K(k₁)·K'(1/k₁))
 * where k = ω_pass/ω_stop, k₁ = ε/√(10^{A_stop/10}-1),
 * K(k) is the complete elliptic integral of the first kind.
 *
 * @param spec  Filter specification
 * @return      Required filter order N (≥1)
 *
 * L5 Algorithm: Elliptic (Cauer) filters provide the sharpest
 * transition band for a given order. They use Jacobi elliptic
 * functions to place poles and zeros optimally.
 *
 * Reference: Cauer, "Theorie der linearen Wechselstromschaltungen" (1941)
 * Course: Stanford EE247, ETH 227-0455
 */
int filter_order_elliptic(const filter_spec_t *spec);

/**
 * @brief Compute required Bessel filter order.
 *
 * Bessel filters approximate a constant group delay (linear phase)
 * in the passband. The order is determined by the required group
 * delay flatness.
 *
 * For Bessel filters, the -3 dB frequency is approximately:
 *   ω_{-3dB} ≈ √((2N-1)·ln 2) / (normalized delay)
 *
 * @param spec        Filter specification
 * @param delay_flat  Allowed group delay variation (% in passband)
 * @return            Required filter order N (≥1)
 *
 * L5 Algorithm: Bessel approximation uses Bessel polynomials
 * derived from the Bessel function of the first kind. The poles
 * lie on a set of contours in the s-plane, not on a circle.
 *
 * The transfer function uses Bessel polynomials y_N(s):
 *   H(s) = y_N(0)/y_N(s)
 * where y_N(s) = Σ_{k=0}^N (N+k)!/(2^{N-k}·k!·(N-k)!) · s^k
 *
 * Reference: Thomson, "Delay Networks having Maximally Flat
 * Frequency Characteristics" (1949)
 * Course: Stanford EE247
 */
int filter_order_bessel(const filter_spec_t *spec, double delay_flat);

/* ============================================================================
 * PROTOTYPE GENERATION
 * ============================================================================ */

/**
 * @brief Generate Butterworth lowpass prototype transfer function.
 *
 * Normalized cutoff ω_c = 1 rad/s.
 * Denominator = Butterworth polynomial B_N(s).
 * Numerator = 1 (all-pole filter).
 *
 * @param order  Filter order N
 * @return       Normalized transfer function H(s), or NULL
 *
 * L5 Algorithm: Butterworth poles in the s-plane are located at:
 *   p_k = -sin((2k-1)π/(2N)) + j·cos((2k-1)π/(2N))
 * for k = 1, 2, ..., N. These lie on the unit circle in the left
 * half-plane.
 */
tf_polynomial_t *filter_butterworth_prototype(int order);

/**
 * @brief Generate Chebyshev Type I lowpass prototype.
 *
 * Passband ripple ε determines the pole locations on an ellipse.
 * The ellipse major axis (real) = sinh(φ), minor axis (imag) = cosh(φ),
 * where φ = (1/N)·asinh(1/ε).
 *
 * @param order     Filter order N
 * @param ripple_db Passband ripple in dB (typically 0.1, 0.5, 1, 3)
 * @return          Normalized transfer function H(s), or NULL
 */
tf_polynomial_t *filter_chebyshev1_prototype(int order, double ripple_db);

/**
 * @brief Generate Chebyshev Type II lowpass prototype.
 *
 * Inverse Chebyshev: has finite zeros on the jω-axis (stopband).
 * The stopband attenuation determines zero/pole locations.
 *
 * @param order        Filter order N
 * @param stopband_db  Minimum stopband attenuation (dB)
 * @return             Normalized transfer function H(s), or NULL
 */
tf_polynomial_t *filter_chebyshev2_prototype(int order, double stopband_db);

/**
 * @brief Generate Elliptic (Cauer) lowpass prototype.
 *
 * Provides equiripple behavior in both passband and stopband.
 * Uses Jacobi elliptic functions for pole-zero placement.
 *
 * The selectivity factor k = ω_c/ω_s determines the transition ratio.
 * The discrimination factor k₁ = ε/√(10^{A_stop/10}-1).
 *
 * @param order        Filter order N
 * @param ripple_db    Passband ripple (dB)
 * @param stopband_db  Minimum stopband attenuation (dB)
 * @return             Normalized transfer function H(s), or NULL
 *
 * L8 Advanced: The computation of elliptic filter poles and zeros
 * requires Jacobi elliptic functions sn(u,k), cn(u,k), dn(u,k),
 * which are computed via the arithmetic-geometric mean (AGM) algorithm.
 *
 * Reference: Orfanidis, "Lecture Notes on Elliptic Filter Design" (2006)
 */
tf_polynomial_t *filter_elliptic_prototype(int order, double ripple_db,
                                             double stopband_db);

/**
 * @brief Generate Bessel lowpass prototype.
 *
 * Maximally flat group delay response.
 * Normalized such that group delay = 1 second at DC.
 *
 * Bessel polynomials y_N(s) satisfy the recurrence:
 *   y_N(s) = (2N-1)·y_{N-1}(s) + s²·y_{N-2}(s)
 * with y₀(s) = 1, y₁(s) = s + 1.
 *
 * @param order  Filter order N
 * @return       Normalized transfer function H(s), or NULL
 */
tf_polynomial_t *filter_bessel_prototype(int order);

/* ============================================================================
 * FREQUENCY TRANSFORMATIONS
 * ============================================================================ */

/**
 * @brief Transform a normalized lowpass to a lowpass with cutoff ω_c.
 *
 * Substitution: s → s/ω_c
 *
 * This scales the frequency axis by ω_c. A normalized LP with
 * cutoff at 1 rad/s becomes a LP with cutoff at ω_c rad/s.
 *
 * @param lp_proto    Normalized lowpass prototype
 * @param cutoff_rad  Desired cutoff ω_c (rad/s)
 * @return            Denormalized lowpass transfer function
 *
 * L5 Algorithm: Frequency scaling is the simplest transformation.
 * Each coefficient a_k in the denominator is multiplied by ω_c^{N-k}.
 */
tf_polynomial_t *filter_lp_to_lp(const tf_polynomial_t *lp_proto,
                                   double cutoff_rad);

/**
 * @brief Transform normalized lowpass to highpass with cutoff ω_c.
 *
 * Substitution: s → ω_c/s
 *
 * This maps ω=0 (DC) to ω=∞ and ω=∞ to ω=0.
 * A normalized LP becomes a HP with cutoff at ω_c.
 *
 * The transformation:
 *   H_HP(s) = H_LP(ω_c/s)
 *
 * For an Nth-order LP, the resulting HP has N zeros at s=0
 * (capacitive coupling at low frequencies).
 *
 * @param lp_proto    Normalized lowpass prototype
 * @param cutoff_rad  Desired cutoff ω_c (rad/s)
 * @return            Highpass transfer function
 *
 * L5 Algorithm: LP→HP transformation is a conformal mapping
 * that interchanges s=0 and s=∞ while preserving the imaginary axis.
 */
tf_polynomial_t *filter_lp_to_hp(const tf_polynomial_t *lp_proto,
                                   double cutoff_rad);

/**
 * @brief Transform normalized lowpass to bandpass.
 *
 * Substitution: s → (s² + ω₀²)/(BW·s)
 *
 * where ω₀ = √(ω_L·ω_H) is the geometric center frequency,
 * BW = ω_H - ω_L is the bandwidth.
 *
 * An Nth-order LP becomes a 2Nth-order BP.
 * Each LP pole at p_k becomes two BP poles.
 *
 * @param lp_proto    Normalized lowpass prototype
 * @param center_rad  Center frequency ω₀ (rad/s)
 * @param bw_rad      Bandwidth BW (rad/s)
 * @return            Bandpass transfer function (order 2N)
 *
 * L5 Algorithm: LP→BP transformation doubles the filter order.
 * It maps the LP frequency range [0, ∞) to the BP range [0, ∞)
 * with Ω=0↔ω=ω₀, Ω=1↔ω=(BW±√(BW²+4ω₀²))/2.
 *
 * For narrow-band filters (BW << ω₀), the BP response is
 * approximately geometrically symmetric about ω₀.
 */
tf_polynomial_t *filter_lp_to_bp(const tf_polynomial_t *lp_proto,
                                   double center_rad, double bw_rad);

/**
 * @brief Transform normalized lowpass to bandstop.
 *
 * Substitution: s → BW·s/(s² + ω₀²)
 *
 * An Nth-order LP becomes a 2Nth-order BS (notch filter).
 *
 * @param lp_proto    Normalized lowpass prototype
 * @param center_rad  Center frequency ω₀ (rad/s)
 * @param bw_rad      Bandwidth BW (rad/s)
 * @return            Bandstop transfer function (order 2N)
 */
tf_polynomial_t *filter_lp_to_bs(const tf_polynomial_t *lp_proto,
                                   double center_rad, double bw_rad);

/* ============================================================================
 * COMPLETE FILTER DESIGN
 * ============================================================================ */

/**
 * @brief Design a complete filter from specification.
 *
 * This is the top-level design function that:
 * 1. Computes the required order N
 * 2. Generates the prototype transfer function
 * 3. Applies frequency transformation
 * 4. Computes normalized g-values (for LC ladder)
 * 5. Denormalizes to actual component values
 *
 * @param spec   Filter specification
 * @param type   Desired filter type
 * @return       Complete filter design, or NULL on error
 *
 * L6 Canonical Problem: End-to-end filter design from specification
 * to component values. This is the central problem of analog filter
 * synthesis, combining approximation theory, network synthesis,
 * and practical implementation constraints.
 */
filter_design_t *filter_design(const filter_spec_t *spec, filter_type_t type);

/**
 * @brief Free a filter design result.
 */
void filter_design_free(filter_design_t *design);

/* ============================================================================
 * G-VALUE COMPUTATION (LC LADDER ELEMENTS)
 * ============================================================================ */

/**
 * @brief Compute normalized g-values for Butterworth LC ladder.
 *
 * For a doubly-terminated LC ladder:
 *   g₀ = 1 (source resistance)
 *   g_k = 2·sin((2k-1)π/(2N)) for k = 1, 2, ..., N
 *   g_{N+1} = 1 (load resistance for N odd)
 *   g_{N+1} = may differ for N even
 *
 * Alternating: g₁, g₃, g₅, ... are series L or shunt C depending
 * on the first element choice.
 *
 * @param order  Filter order N
 * @param n_elements  Output: number of g-values (N+2 for doubly-terminated)
 * @return            Array of g-values (caller frees), or NULL
 *
 * L5 Algorithm: The g-values for Butterworth are analytically known
 * (they don't require numerical synthesis). This is a special property
 * of Butterworth filters.
 *
 * Reference: Matthaei, Young, Jones, "Microwave Filters, Impedance-
 * Matching Networks, and Coupling Structures" (1980)
 */
double *filter_g_values_butterworth(int order, size_t *n_elements);

/**
 * @brief Compute normalized g-values for Chebyshev LC ladder.
 *
 * Requires iterative computation using Darlington synthesis.
 *
 * @param order     Filter order N
 * @param ripple_db Passband ripple (dB)
 * @param n_elements Output: number of g-values
 * @return          Array of g-values (caller frees), or NULL
 */
double *filter_g_values_chebyshev(int order, double ripple_db,
                                    size_t *n_elements);

/**
 * @brief Denormalize g-values to actual component values.
 *
 * Applies impedance and frequency scaling:
 *   For series inductor g_k:  L = g_k·R₀/ω_c
 *   For shunt capacitor g_k:  C = g_k/(R₀·ω_c)
 *   For termination resistor:  R = g_k·R₀
 *
 * @param g_values       Normalized g-value array
 * @param n              Number of g-values
 * @param R0             Impedance scaling factor (Ω, typically 50)
 * @param cutoff_rad     Cutoff frequency ω_c (rad/s)
 * @param L_values       Output: inductor values (H), length = n/2
 * @param C_values       Output: capacitor values (F), length = n/2
 * @param R_load         Output: load resistance (Ω)
 * @return               0 on success, -1 on error
 */
int filter_denormalize(const double *g_values, size_t n,
                        double R0, double cutoff_rad,
                        double **L_values, double **C_values,
                        double *R_load);

/* ============================================================================
 * ACTIVE FILTER COMPONENT CALCULATION
 * ============================================================================ */

/**
 * @brief Compute Sallen-Key lowpass component values.
 *
 * Sallen-Key topology (1955): single op-amp, positive feedback
 * through a capacitor. Provides a non-inverting 2nd-order section.
 *
 * For a normalized biquad s² + (ω₀/Q)s + ω₀²:
 * Choose C₁, C₂ (arbitrary), then:
 *   R₁ = R₂ = 1/(ω₀·√(C₁·C₂))  (equal-R design)
 * or:
 *   Choose m = C₂/C₁ ratio, then compute R₁, R₂.
 *
 * @param biquad        Biquad specification
 * @param C1            Chosen C₁ value (F)
 * @param C2            Chosen C₂ value (F)
 * @param gain           Desired passband gain (≥1)
 * @param R1            Output: resistor R₁ (Ω)
 * @param R2            Output: resistor R₂ (Ω)
 * @param R3            Output: resistor R₃ for gain (Ω)
 * @param R4            Output: resistor R₄ for gain (Ω)
 * @return              0 on success, -1 if values infeasible
 *
 * L6 Canonical Problem: Sallen-Key active filter design.
 * Component sensitivity to tolerance is a key practical concern.
 *
 * Reference: Sallen & Key, "A Practical Method of Designing
 * RC Active Filters" (1955), IRE Trans. Circuit Theory
 * Course: Berkeley EE105
 */
int filter_sallen_key_lp(const biquad_section_t *biquad,
                          double C1, double C2, double gain,
                          double *R1, double *R2, double *R3, double *R4);

/**
 * @brief Compute Multiple Feedback (MFB) lowpass component values.
 *
 * MFB topology: single op-amp, negative feedback through two paths.
 * Provides an inverting 2nd-order section with lower sensitivity
 * than Sallen-Key for high Q values.
 *
 * For normalized biquad with DC gain = -R₂/R₁ (inverting):
 *   C₁, C₂ = chosen
 *   R₂ = 1/(ω₀·C₁·√(...)) — solved from design equations
 *   R₁ = R₂/|gain|
 *   R₃ = 1/(ω₀²·R₂·C₁·C₂)
 *
 * @param biquad  Biquad specification
 * @param C1      Chosen C₁ (F)
 * @param C2      Chosen C₂ (F)
 * @param R1      Output: input resistor (Ω)
 * @param R2      Output: feedback resistor (Ω)
 * @param R3      Output: ground resistor (Ω)
 * @return        0 on success, -1 if infeasible
 *
 * Reference: Sedra & Smith (2020), Ch. 16.4
 */
int filter_mfb_lp(const biquad_section_t *biquad,
                   double C1, double C2,
                   double *R1, double *R2, double *R3);

/**
 * @brief Compute Tow-Thomas biquad component values.
 *
 * Tow-Thomas topology (1968): three op-amp state-variable biquad
 * providing simultaneous LP and BP outputs. Excellent for high-Q
 * filters due to low component sensitivity.
 *
 * Uses integrator-based implementation of the state equations:
 *   s²·V_out + (ω₀/Q)·s·V_out + ω₀²·V_out = ω₀²·V_in (for LP)
 *
 * @param biquad      Biquad specification
 * @param C           Integration capacitor (F, same for both integrators)
 * @param R_freq      Output: frequency-setting resistor (Ω)
 * @param R_q         Output: Q-setting resistor (Ω)
 * @param R_input     Output: input resistor (Ω)
 * @param R_feedback  Output: feedback resistor (Ω)
 * @return            0 on success, -1 if infeasible
 *
 * Reference: Tow, "A Step-by-Step Active Filter Design" (1969)
 * Course: Berkeley EE105, Stanford EE247
 */
int filter_tow_thomas(const biquad_section_t *biquad, double C,
                       double *R_freq, double *R_q,
                       double *R_input, double *R_feedback);

/* ============================================================================
 * FILTER SENSITIVITY ANALYSIS
 * ============================================================================ */

/**
 * @brief Compute component sensitivity for a filter parameter.
 *
 * Sensitivity S_x^P = (∂P/∂x)·(x/P) = fractional change in parameter P
 * per fractional change in component x.
 *
 * Low sensitivity |S| < 1 indicates robustness to component tolerances.
 * High sensitivity |S| > 1 may require precision components or trimming.
 *
 * For Sallen-Key LP:
 *   S_R1^ω₀ = S_R2^ω₀ = S_C1^ω₀ = S_C2^ω₀ = -1/2
 *   S_R1^Q and S_R2^Q depend on the gain setting
 *
 * @param biquad      Biquad specification
 * @param topology    0=Sallen-Key, 1=MFB, 2=Tow-Thomas
 * @param component   Component index (0=R1, 1=R2, 2=C1, 3=C2)
 * @param S_omega0    Output: sensitivity to ω₀
 * @param S_Q         Output: sensitivity to Q
 * @return            0 on success, -1 on error
 *
 * L8 Advanced: Sensitivity analysis quantifies the effect of component
 * variations on filter performance. This is critical for production
 * yield prediction and design for manufacturability (DFM).
 *
 * Reference: Sedra & Brackett, "Filter Theory and Design: Active
 * and Passive" (1978)
 */
int filter_sensitivity(const biquad_section_t *biquad, int topology,
                        int component, double *S_omega0, double *S_Q);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_DESIGN_H */
