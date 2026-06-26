/**
 * @file harmonic_power.h
 * @brief Harmonic analysis, THD computation, and distortion power theory
 *
 * In non-sinusoidal AC systems, the power factor decomposes into:
 *   PF_true = PF_displacement × PF_distortion
 *
 * Where PF_displacement = cos φ₁ (fundamental phase shift) and
 *       PF_distortion = 1 / √(1 + THD_v²) × 1 / √(1 + THD_i²)
 *
 * Harmonic power flow follows from Parseval's theorem:
 *   P_total = Σ P_h  (real power at each harmonic)
 *   Cross-frequency products integrate to zero over a fundamental period.
 *
 * References:
 *   - IEEE Std 519-2014, "Recommended Practice for Harmonic Control"
 *   - IEEE Std 1459-2010, "Definitions for Power Quantities"
 *   - Arrillaga & Watson, "Power System Harmonics" (2nd ed, 2003)
 *   - Grady, "Understanding Power System Harmonics" (UTexas course notes)
 *
 * Knowledge coverage:
 *   L1 (Definitions): THD, TDD, harmonic order, interharmonic
 *   L2 (Concepts):   Distortion PF, harmonic power flow, K-factor
 *   L3 (Math):       Fourier series, Parseval, Goertzel algorithm
 *   L4 (Laws):       Parseval's theorem for power, harmonic orthogonality
 *   L5 (Algorithms): FFT for harmonics, Goertzel detector, THD computation
 */

#ifndef HARMONIC_POWER_H
#define HARMONIC_POWER_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HP_MAX_HARMONICS   63
#define HP_FFT_MAX_SIZE    4096

/* ==========================================================================
 * L1: Harmonic Distortion Definitions
 * ========================================================================== */

/**
 * @brief Harmonic measurement standard enumeration.
 */
typedef enum {
    HP_STD_IEEE519    = 0,   /**< IEEE 519-2014 (current harmonics at PCC) */
    HP_STD_IEC61000_32 = 1,  /**< IEC 61000-3-2 (equipment < 16A) */
    HP_STD_IEC61000_34 = 2,  /**< IEC 61000-3-4 (equipment > 16A) */
    HP_STD_IEC61000_312 = 3, /**< IEC 61000-3-12 (equipment 16-75A) */
    HP_STD_EN50160    = 4,   /**< EN 50160 (voltage quality) */
} hp_standard_t;

/**
 * @brief Harmonic limits per IEEE 519-2014 for general distribution systems.
 *
 * Current distortion limits (I_sc / I_L ratio dependent):
 *   I_sc / I_L < 20:  TDD ≤ 5.0%
 *   20 ≤ I_sc/I_L < 50:  TDD ≤ 8.0%
 *   50 ≤ I_sc/I_L < 100: TDD ≤ 12.0%
 *   100 ≤ I_sc/I_L < 1000: TDD ≤ 15.0%
 *   I_sc/I_L ≥ 1000: TDD ≤ 20.0%
 *
 * Individual harmonic limits vary by odd/even and order.
 */
typedef struct {
    double tdd_limit_percent;           /**< Total demand distortion limit [%] */
    double individual_limit_odd[50];    /**< Limit per odd harmonic order [%] */
    double individual_limit_even[50];   /**< Limit per even harmonic order [%] */
    hp_standard_t standard;             /**< Which standard applies */
    double i_sc;                        /**< Short-circuit current at PCC [A] */
    double i_l_max;                     /**< Maximum demand load current [A] */
} hp_limits_t;

/**
 * @brief Complete harmonic spectrum analysis result.
 *
 * Contains the Fourier decomposition of voltage and current waveforms,
 * along with derived power quality metrics.
 */
typedef struct {
    /* Fundamental */
    double   f_fundamental;             /**< Fundamental frequency [Hz] */
    double   v1_mag;                    /**< Fundamental voltage magnitude [V] */
    double   v1_phase_deg;              /**< Fundamental voltage phase [deg] */
    double   i1_mag;                    /**< Fundamental current magnitude [A] */
    double   i1_phase_deg;              /**< Fundamental current phase [deg] */

    /* THD metrics */
    double   thd_v_percent;             /**< Voltage THD = V_H / V₁ × 100 [%] */
    double   thd_i_percent;             /**< Current THD = I_H / I₁ × 100 [%] */
    double   tdd_i_percent;             /**< Total demand distortion [%] (I_H / I_L_max) */
    double   pwrhd_v_percent;           /**< Power-related HD for voltage */
    double   crest_factor_v;            /**< Voltage crest factor */
    double   crest_factor_i;            /**< Current crest factor */

    /* Distortion power factor */
    double   pf_displacement;           /**< cos(φ₁) — displacement PF */
    double   pf_distortion;             /**< 1/√((1+THDv²)(1+THDi²)) */
    double   pf_true;                   /**< PF_displacement × PF_distortion */

    /* K-factor (transformer derating) */
    double   k_factor;                  /**< K-factor = Σ h² × (I_h/I₁)² */

    /* Individual harmonic components */
    uint32_t num_harmonics;
    struct {
        uint32_t order;                 /**< Harmonic order h */
        double   v_mag;                 /**< Voltage magnitude [V] */
        double   v_phase_deg;           /**< Voltage phase [deg] */
        double   i_mag;                 /**< Current magnitude [A] */
        double   i_phase_deg;           /**< Current phase [deg] */
        double   p_h;                   /**< Real power at harmonic h [W] */
        double   i_h_percent;           /**< Current as % of fundamental */
        uint8_t  within_limit;          /**< Flag: within IEEE 519 limit */
    } harmonics[HP_MAX_HARMONICS];

    /* Interharmonic detection */
    uint32_t num_interharmonics;
    double   interharmonic_thd_percent;
} hp_spectrum_t;

/* ==========================================================================
 * L2: Total Harmonic Distortion (THD) Computation
 * ========================================================================== */

/**
 * @brief Compute THD from harmonic magnitudes.
 *
 * THD = sqrt( Σ_{h=2}^{N} M_h² ) / M₁ × 100%
 *
 * where M_h is the magnitude (voltage or current) of harmonic h.
 *
 * @param harmonics   Array of harmonic magnitudes (index 0 = fundamental)
 * @param n_harmonics Number of harmonics (including fundamental)
 * @return THD as percentage [0-∞%), -1 on error
 */
double hp_compute_thd(const double *harmonics, size_t n_harmonics);

/**
 * @brief Compute Total Demand Distortion (TDD).
 *
 * TDD = sqrt( Σ_{h=2}^{N} I_h² ) / I_L_max × 100%
 *
 * Unlike THD which uses I₁ as reference, TDD uses the maximum demand
 * load current I_L_max as the denominator per IEEE 519.
 *
 * @param i_harmonics  Current harmonic magnitudes [A]
 * @param n_harmonics  Number of harmonics
 * @param i_l_max      Maximum demand load current [A]
 * @return TDD as percentage, -1 on error
 */
double hp_compute_tdd(const double *i_harmonics, size_t n_harmonics,
                       double i_l_max);

/**
 * @brief Compute partial weighted harmonic distortion (PWHD).
 *
 * PWHD = sqrt( Σ_{h=14}^{40} h × V_h² ) / V₁ × 100%
 *
 * Gives higher weighting to higher-order harmonics (14th-40th),
 * which are more likely to cause interference in communication circuits.
 *
 * @param v_harmonics  Voltage harmonic magnitudes [V]
 * @param n_harmonics  Total number of harmonics (must include up to h=40)
 * @return PWHD as percentage, -1 on error
 */
double hp_compute_pwhd(const double *v_harmonics, size_t n_harmonics);

/* ==========================================================================
 * L3: Fourier-Based Harmonic Power Analysis
 * ========================================================================== */

/**
 * @brief Compute real power at each harmonic from sampled waveforms.
 *
 * Uses DFT (discrete Fourier transform) to decompose voltage and current
 * into harmonic components and computes P_h = V_h × I_h × cos(φ_h) for each.
 *
 * Assumes an integer number of fundamental periods in the sample window
 * (synchronous sampling) to avoid spectral leakage.
 *
 * @param v_samples    Voltage samples [V]
 * @param i_samples    Current samples [A]
 * @param n_samples    Number of samples (must be power of 2 for FFT speed)
 * @param f_sample_hz  Sampling frequency [Hz]
 * @param f_fund_hz    Fundamental frequency [Hz]
 * @param spectrum     [out] Harmonic spectrum analysis
 * @return 0 on success, -1 on error
 */
int hp_analyze_spectrum(const double *v_samples, const double *i_samples,
                         size_t n_samples, double f_sample_hz,
                         double f_fund_hz, hp_spectrum_t *spectrum);

/**
 * @brief Goertzel algorithm for single-harmonic power detection.
 *
 * The Goertzel algorithm (1958) efficiently computes a single DFT bin,
 * requiring only 2N+2 real multiplications per bin versus N·log₂N for FFT.
 * Optimal for real-time monitoring of specific harmonics (3rd, 5th, 7th).
 *
 * Transfer function: H_k(z) = (1 - W_N^k z^{-1}) / (1 - 2cos(2πk/N)z^{-1} + z^{-2})
 *
 * @param samples      Input sample array
 * @param n_samples    Number of samples
 * @param target_hz    Target frequency to detect [Hz]
 * @param sample_hz    Sampling frequency [Hz]
 * @param magnitude    [out] Detected magnitude
 * @param phase_deg    [out] Detected phase [degrees]
 * @return 0 on success, -1 on error
 */
int hp_goertzel_detect(const double *samples, size_t n_samples,
                        double target_hz, double sample_hz,
                        double *magnitude, double *phase_deg);

/**
 * @brief Compute harmonic power using FFT.
 *
 * Performs radix-2 DIT FFT on both voltage and current arrays,
 * then computes complex power at each harmonic bin.
 *
 * Uses the Cooley-Tukey (1965) in-place algorithm.
 *
 * @param v_complex    Complex voltage samples (in-place modified)
 * @param i_complex    Complex current samples (in-place modified)
 * @param n_fft        FFT size (must be power of 2)
 * @param s_harmonic   [out] Complex power at each harmonic bin
 * @return 0 on success, -1 on error
 */
int hp_fft_power(double complex *v_complex, double complex *i_complex,
                  size_t n_fft, double complex *s_harmonic);

/**
 * @brief In-place radix-2 DIT FFT (Cooley-Tukey).
 *
 * Core DFT computation using the butterfly algorithm.
 * Complexity: O(N log₂ N)
 *
 * @param data     Complex array (modified in place)
 * @param n        FFT size (must be power of 2)
 * @param inverse  Non-zero for IFFT
 * @return 0 on success
 */
int hp_fft_radix2(double complex *data, size_t n, int inverse);

/* ==========================================================================
 * L4: Distortion Power and Parseval's Theorem
 * ========================================================================== */

/**
 * @brief Compute distortion power factor components.
 *
 * From IEEE 1459-2010:
 *   PF = P_total / S_total = (P₁ / S₁) × (S₁ / S_total) × (P_total / P₁)
 *
 * For sinusoidal voltage: PF = cos φ₁ / √(1 + THD_i²)
 *   where cos φ₁ is the displacement PF
 *   and 1/√(1 + THD_i²) is the distortion PF due to current harmonics.
 *
 * For general case (both V and I distorted):
 *   PF_distortion = 1 / √((1 + THD_v²)(1 + THD_i²))
 *
 * @param thd_v       Voltage THD (as decimal, e.g., 0.05 for 5%)
 * @param thd_i       Current THD (as decimal)
 * @param cos_phi1    Displacement power factor (cos of fundamental angle)
 * @param pf_distortion [out] Distortion power factor
 * @param pf_true     [out] True power factor
 * @return 0 on success
 */
int hp_decompose_pf(double thd_v, double thd_i, double cos_phi1,
                     double *pf_distortion, double *pf_true);

/**
 * @brief Verify Parseval's theorem for power.
 *
 * Parseval: (1/N) Σ |x[n]|² = Σ |X[k]|²
 *
 * This means the total power computed in the time domain equals
 * the sum of powers in the frequency domain (each harmonic contributes
 * independently — cross-frequency products integrate to zero).
 *
 * @param time_domain_power  Average power from time samples [W]
 * @param harmonic_powers    Array of per-harmonic real powers [W]
 * @param n_harmonics        Number of harmonics
 * @param tolerance          Allowed numerical tolerance
 * @return 0 if Parseval holds, 1 if violation
 */
int hp_verify_parseval_power(double time_domain_power,
                              const double *harmonic_powers,
                              size_t n_harmonics, double tolerance);

/**
 * @brief Compute K-factor for transformer harmonic derating.
 *
 * K-factor = Σ_{h=1}^{N} h² × (I_h / I₁)²
 *
 * Transformers supplying non-linear loads must be derated because
 * eddy-current losses increase with the square of both the harmonic
 * order and the current magnitude. A K-4 transformer handles loads
 * with K-factor ≤ 4; K-13 handles ≤ 13, etc.
 *
 * Per UL 1561 and IEEE C57.110.
 *
 * @param i_harmonics  Current magnitudes per harmonic [A]
 * @param n_harmonics  Number of harmonics (index 0 = fundamental)
 * @return K-factor, -1 on error
 */
double hp_k_factor(const double *i_harmonics, size_t n_harmonics);

/**
 * @brief Compute transformer derating factor from K-factor.
 *
 * Derating = 1 / (1 + α × (K_factor - 1))
 *
 * where α is the ratio of eddy-current loss to I²R loss at fundamental
 * (typically 0.05-0.10 for dry-type transformers per IEEE C57.110).
 *
 * @param k_factor     Computed K-factor
 * @param alpha        Eddy/R loss ratio (typically 0.08)
 * @return Derating factor (0-1), multiply rated kVA by this
 */
double hp_transformer_derating(double k_factor, double alpha);

/* ==========================================================================
 * L5: IEEE 519 Compliance Checking
 * ========================================================================== */

/**
 * @brief Initialize IEEE 519 harmonic limits for a given system.
 *
 * Sets up per-harmonic limits based on I_sc/I_L ratio and bus voltage.
 *
 * @param limits       [out] Initialized limits structure
 * @param standard     Standard to apply
 * @param i_sc         Short-circuit current at PCC [A]
 * @param i_l_max      Maximum demand load current [A]
 * @param v_nominal    Nominal bus voltage [V]
 * @return 0 on success
 */
int hp_init_ieee519_limits(hp_limits_t *limits, hp_standard_t standard,
                            double i_sc, double i_l_max, double v_nominal);

/**
 * @brief Check if measured harmonics comply with IEEE 519 limits.
 *
 * Evaluates both TDD and individual harmonic limits.
 *
 * @param spectrum  Measured harmonic spectrum
 * @param limits    Applicable limits
 * @param violations [out] Number of violations found
 * @return 0 if compliant, >0 for TDD violation, <0 for individual violation
 */
int hp_check_ieee519(const hp_spectrum_t *spectrum,
                      const hp_limits_t *limits, int *violations);

/**
 * @brief Identify dominant harmonic sources from spectrum.
 *
 * Returns the orders of the top-N most significant harmonic
 * contributors (by current magnitude or power).
 *
 * @param spectrum     Measured spectrum
 * @param top_n        Number of top contributors to identify
 * @param orders_out   [out] Harmonic orders, sorted by contribution
 * @return Number of contributors found
 */
int hp_dominant_harmonics(const hp_spectrum_t *spectrum, int top_n,
                           int *orders_out);

/**
 * @brief Generate harmonic mitigation recommendations.
 *
 * Based on harmonic profile, suggests:
 * - Passive filters (tuned to 5th, 7th, 11th)
 * - Active harmonic filters
 * - Phase-shifting transformers for 12-pulse operation
 * - Line reactors to reduce THD
 *
 * @param spectrum       Measured spectrum
 * @param recommendation [out] String buffer for recommendation
 * @param buf_len        Buffer length
 * @return 0 on success
 */
int hp_recommend_mitigation(const hp_spectrum_t *spectrum,
                             char *recommendation, size_t buf_len);

#ifdef __cplusplus
}
#endif

#endif /* HARMONIC_POWER_H */
