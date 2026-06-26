/**
 * @file harmonic_power.c
 * @brief Harmonic analysis, THD/TDD computation, distortion power, and FFT
 *
 * In modern power systems, non-linear loads (rectifiers, VFDs, LED drivers,
 * SMPS) generate harmonic currents that distort the voltage waveform.
 * IEEE 519-2014 governs harmonic limits at the point of common coupling (PCC).
 *
 * The power factor in a non-sinusoidal system decomposes into:
 *   PF_true = PF_displacement × PF_distortion
 *
 * where PF_distortion accounts for the reduction in power factor caused
 * by harmonic currents that do not contribute to real power transfer.
 *
 * Knowledge points:
 *   L1: THD, TDD, PWHD definitions and computation
 *   L2: Distortion power factor decomposition
 *   L3: DFT/FFT for harmonic power analysis
 *   L3: Goertzel algorithm for single-bin detection
 *   L3: Radix-2 DIT FFT (Cooley-Tukey 1965)
 *   L4: Parseval's theorem verification for power
 *   L5: K-factor for transformer harmonic derating
 *   L5: IEEE 519 compliance checking
 *   L6: Harmonic source identification and mitigation
 *
 * References:
 *   - IEEE Std 519-2014
 *   - IEEE Std 1459-2010
 *   - Arrillaga & Watson (2003)
 *   - Cooley & Tukey, "An Algorithm for the Machine Calculation of
 *     Complex Fourier Series" (Math. Comp., 1965)
 *   - Goertzel, "An Algorithm for the Evaluation of Finite Trigonometric
 *     Series" (American Math. Monthly, 1958)
 */

#include "harmonic_power.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define HP_EPS 1e-12

/* ==========================================================================
 * L1: Total Harmonic Distortion (THD)
 * ========================================================================== */

double hp_compute_thd(const double *harmonics, size_t n_harmonics)
{
    if (!harmonics || n_harmonics < 2) return -1.0;

    double fundamental = harmonics[0];
    if (fundamental < HP_EPS) return -1.0;

    double sum_h2 = 0.0;
    for (size_t h = 1; h < n_harmonics; h++) {
        sum_h2 += harmonics[h] * harmonics[h];
    }

    return 100.0 * sqrt(sum_h2) / fundamental;
}

/* ==========================================================================
 * L1: Total Demand Distortion (TDD)
 * ========================================================================== */

double hp_compute_tdd(const double *i_harmonics, size_t n_harmonics,
                       double i_l_max)
{
    if (!i_harmonics || n_harmonics < 2) return -1.0;
    if (i_l_max < HP_EPS) return -1.0;

    /* TDD = sqrt(Σ_{h=2}^{N} I_h²) / I_L_max × 100% */
    double sum_h2 = 0.0;
    for (size_t h = 1; h < n_harmonics; h++) {
        sum_h2 += i_harmonics[h] * i_harmonics[h];
    }

    return 100.0 * sqrt(sum_h2) / i_l_max;
}

/* ==========================================================================
 * L1: Partial Weighted Harmonic Distortion (PWHD)
 * ========================================================================== */

double hp_compute_pwhd(const double *v_harmonics, size_t n_harmonics)
{
    if (!v_harmonics || n_harmonics < 2) return -1.0;

    double fundamental = v_harmonics[0];
    if (fundamental < HP_EPS) return -1.0;

    /* PWHD = sqrt(Σ_{h=14}^{40} h × V_h²) / V₁ × 100% */
    double sum_weighted = 0.0;
    size_t h_start = 14;
    size_t h_end = (n_harmonics < 41) ? n_harmonics : 41;

    for (size_t h = h_start; h < h_end && h < n_harmonics; h++) {
        /* Index h corresponds to harmonic order h (0-indexed, so h=1 is idx 0) */
        int order = (int)(h);
        sum_weighted += order * v_harmonics[h] * v_harmonics[h];
    }

    return 100.0 * sqrt(sum_weighted) / fundamental;
}

/* ==========================================================================
 * L3: Bit-Reversal Permutation for FFT
 * ========================================================================== */

static void hp_bit_reverse(double complex *data, size_t n)
{
    size_t j = 0;
    size_t n_over_2 = n >> 1;
    for (size_t i = 1; i < n - 1; i++) {
        size_t bit = n_over_2;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            double complex tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
        }
    }
}

/* ==========================================================================
 * L3: Radix-2 DIT FFT (Cooley-Tukey 1965)
 * ========================================================================== */

int hp_fft_radix2(double complex *data, size_t n, int inverse)
{
    if (!data || n == 0) return -1;

    /* Check n is power of 2 */
    size_t tmp = n;
    while ((tmp & 1) == 0) tmp >>= 1;
    if (tmp != 1 && n > 1) return -1;

    /* Bit-reversal permutation */
    hp_bit_reverse(data, n);

    /* Butterfly loops */
    double sign = inverse ? 1.0 : -1.0;

    for (size_t len = 2; len <= n; len <<= 1) {
        double angle = sign * 2.0 * M_PI / (double)len;
        double complex wlen = cos(angle) + I * sin(angle);

        for (size_t i = 0; i < n; i += len) {
            double complex w = 1.0 + I * 0.0;
            size_t half = len >> 1;

            for (size_t j = 0; j < half; j++) {
                double complex u = data[i + j];
                double complex v = data[i + j + half] * w;

                data[i + j]         = u + v;
                data[i + j + half]  = u - v;

                w *= wlen;
            }
        }
    }

    /* Scale for IFFT */
    if (inverse) {
        double inv_n = 1.0 / (double)n;
        for (size_t i = 0; i < n; i++) {
            data[i] *= inv_n;
        }
    }

    return 0;
}

/* ==========================================================================
 * L3: FFT-Based Harmonic Power Computation
 * ========================================================================== */

int hp_fft_power(double complex *v_complex, double complex *i_complex,
                  size_t n_fft, double complex *s_harmonic)
{
    if (!v_complex || !i_complex || !s_harmonic) return -1;
    if (n_fft == 0) return -1;

    /* Forward FFT on both */
    int rv = hp_fft_radix2(v_complex, n_fft, 0);
    if (rv < 0) return rv;

    int ri = hp_fft_radix2(i_complex, n_fft, 0);
    if (ri < 0) return ri;

    /* Complex power at each bin: S[k] = V[k] × conj(I[k]) / N */
    /* Note: We want the RMS phasor, so divide by N (not N²) */
    double inv_n = 1.0 / (double)n_fft;

    for (size_t k = 0; k < n_fft; k++) {
        /* DC and Nyquist are real; conjugate symmetric bins are redundent */
        s_harmonic[k] = v_complex[k] * conj(i_complex[k]) * inv_n;
    }

    /* DC bin (k=0): average real power */
    /* Fundamental (k=1 for integer-period window): positive-frequency
     * real power = 2 × Re(S[1]) for single-sided spectrum */
    /* For negative frequencies (k > n/2), phase information is conjugate */

    return 0;
}

/* ==========================================================================
 * L5: Harmonic Spectrum Analysis from Sampled Waveforms
 * ========================================================================== */

int hp_analyze_spectrum(const double *v_samples, const double *i_samples,
                         size_t n_samples, double f_sample_hz,
                         double f_fund_hz, hp_spectrum_t *spectrum)
{
    if (!v_samples || !i_samples || !spectrum) return -1;
    if (n_samples < 2 || f_sample_hz < HP_EPS || f_fund_hz < HP_EPS) return -1;

    memset(spectrum, 0, sizeof(*spectrum));
    spectrum->f_fundamental = f_fund_hz;

    /* Determine max harmonic: Nyquist / fundamental */
    double nyquist = f_sample_hz / 2.0;
    int max_harmonic = (int)(nyquist / f_fund_hz);
    if (max_harmonic > (int)HP_MAX_HARMONICS) max_harmonic = (int)HP_MAX_HARMONICS;

    /* Use DFT for each harmonic individually.
     * For a window with exactly k cycles of the fundamental:
     *   V_h = (2/N) Σ v[n] × e^{-j 2π h n / N}
     *
     * where N = n_samples, h = harmonic order.
     *
     * We compute using the quadrature method (sine/cosine correlation).
     */
    for (int h = 1; h <= max_harmonic && h <= (int)HP_MAX_HARMONICS; h++) {
        int idx = h - 1;
        double re_v = 0.0, im_v = 0.0;
        double re_i = 0.0, im_i = 0.0;

        for (size_t n = 0; n < n_samples; n++) {
            double theta = 2.0 * M_PI * (double)h * (double)n / (double)n_samples;
            double cos_th = cos(theta);
            double sin_th = sin(theta);

            re_v += v_samples[n] * cos_th;
            im_v -= v_samples[n] * sin_th;
            re_i += i_samples[n] * cos_th;
            im_i -= i_samples[n] * sin_th;
        }

        /* Normalize: (2/N) for single-sided spectrum */
        double norm = 2.0 / (double)n_samples;
        re_v *= norm;
        im_v *= norm;
        re_i *= norm;
        im_i *= norm;

        double v_mag = sqrt(re_v * re_v + im_v * im_v);
        double v_phase = atan2(im_v, re_v);
        double i_mag = sqrt(re_i * re_i + im_i * im_i);
        double i_phase = atan2(im_i, re_i);

        spectrum->harmonics[idx].order       = h;
        spectrum->harmonics[idx].v_mag       = v_mag;
        spectrum->harmonics[idx].v_phase_deg = v_phase * 180.0 / M_PI;
        spectrum->harmonics[idx].i_mag       = i_mag;
        spectrum->harmonics[idx].i_phase_deg = i_phase * 180.0 / M_PI;

        /* Real power at harmonic h: P_h = V_h × I_h × cos(φv_h - φi_h) */
        double phi_diff = v_phase - i_phase;
        spectrum->harmonics[idx].p_h = v_mag * i_mag * cos(phi_diff);

        spectrum->num_harmonics++;
    }

    /* Extract fundamental and THD */
    if (spectrum->num_harmonics > 0) {
        spectrum->v1_mag       = spectrum->harmonics[0].v_mag;
        spectrum->v1_phase_deg = spectrum->harmonics[0].v_phase_deg;
        spectrum->i1_mag       = spectrum->harmonics[0].i_mag;
        spectrum->i1_phase_deg = spectrum->harmonics[0].i_phase_deg;

        /* THD_i = sqrt(Σ I_h² for h>1) / I₁ */
        double sum_i2_harm = 0.0;
        double sum_v2_harm = 0.0;
        for (uint32_t h = 1; h < spectrum->num_harmonics; h++) {
            sum_i2_harm += spectrum->harmonics[h].i_mag * spectrum->harmonics[h].i_mag;
            sum_v2_harm += spectrum->harmonics[h].v_mag * spectrum->harmonics[h].v_mag;
        }

        if (spectrum->i1_mag > HP_EPS) {
            spectrum->thd_i_percent = 100.0 * sqrt(sum_i2_harm) / spectrum->i1_mag;
        }
        if (spectrum->v1_mag > HP_EPS) {
            spectrum->thd_v_percent = 100.0 * sqrt(sum_v2_harm) / spectrum->v1_mag;
        }

        /* Displacement PF = cos(φ₁) */
        double phi1 = spectrum->harmonics[0].v_phase_deg * M_PI / 180.0
                       - spectrum->harmonics[0].i_phase_deg * M_PI / 180.0;
        spectrum->pf_displacement = cos(phi1);

        /* Distortion PF */
        double thdv = spectrum->thd_v_percent / 100.0;
        double thdi = spectrum->thd_i_percent / 100.0;
        spectrum->pf_distortion = 1.0 / sqrt((1.0 + thdv * thdv) * (1.0 + thdi * thdi));
        spectrum->pf_true = spectrum->pf_displacement * spectrum->pf_distortion;
    }

    return 0;
}

/* ==========================================================================
 * L3: Goertzel Algorithm for Single-Harmonic Detection
 * ========================================================================== */

int hp_goertzel_detect(const double *samples, size_t n_samples,
                        double target_hz, double sample_hz,
                        double *magnitude, double *phase_deg)
{
    if (!samples || !magnitude || !phase_deg) return -1;
    if (n_samples < 3 || target_hz <= 0.0 || sample_hz <= 0.0) return -1;
    if (target_hz >= sample_hz / 2.0) return -1; /* above Nyquist */

    /* Goertzel algorithm (1958):
     *
     * Computes a single DFT bin for frequency f = k * Fs / N.
     *
     * Coefficient: c = 2 × cos(2πk/N)
     * Recursion:   s[n] = x[n] + c × s[n-1] - s[n-2]
     *              for n = 0...N-1, with s[-1] = s[-2] = 0
     *
     * Output: |X[k]|² = s[N-1]² + s[N-2]² - c × s[N-1] × s[N-2]
     *         arg(X[k]) = atan2(sin_term, cos_term)
     */

    double k_frac = target_hz * (double)n_samples / sample_hz;
    size_t k = (size_t)(k_frac + 0.5); /* nearest bin */
    if (k >= n_samples / 2) k = n_samples / 2 - 1;

    double omega = 2.0 * M_PI * (double)k / (double)n_samples;
    double coeff = 2.0 * cos(omega);

    double s0 = 0.0, s1 = 0.0, s2 = 0.0;

    for (size_t n = 0; n < n_samples; n++) {
        s0 = samples[n] + coeff * s1 - s2;
        s2 = s1;
        s1 = s0;
    }

    /* |X[k]|² = s1² + s2² - coeff × s1 × s2 */
    double mag_sq = s1 * s1 + s2 * s2 - coeff * s1 * s2;
    if (mag_sq < 0.0) mag_sq = 0.0;
    *magnitude = 2.0 * sqrt(mag_sq) / (double)n_samples;

    /* Phase: arg(X[k]) = atan2(sin(ω) × s2, s1 - cos(ω) × s2) */
    double real_part = s1 - cos(omega) * s2;
    double imag_part = sin(omega) * s2;
    *phase_deg = atan2(imag_part, real_part) * 180.0 / M_PI;

    return 0;
}

/* ==========================================================================
 * L2: Decompose Power Factor into Displacement and Distortion
 * ========================================================================== */

int hp_decompose_pf(double thd_v, double thd_i, double cos_phi1,
                     double *pf_distortion, double *pf_true)
{
    if (!pf_distortion || !pf_true) return -1;

    /* PF_distortion = 1 / √((1 + THD_v²)(1 + THD_i²)) */
    double factor = (1.0 + thd_v * thd_v) * (1.0 + thd_i * thd_i);
    if (factor < HP_EPS) {
        *pf_distortion = 1.0;
    } else {
        *pf_distortion = 1.0 / sqrt(factor);
        if (*pf_distortion > 1.0) *pf_distortion = 1.0;
    }

    /* PF_true = PF_displacement × PF_distortion */
    *pf_true = cos_phi1 * (*pf_distortion);
    if (*pf_true > 1.0) *pf_true = 1.0;
    if (*pf_true < -1.0) *pf_true = -1.0;

    return 0;
}

/* ==========================================================================
 * L4: Verify Parseval's Theorem for Power
 * ========================================================================== */

int hp_verify_parseval_power(double time_domain_power,
                              const double *harmonic_powers,
                              size_t n_harmonics, double tolerance)
{
    if (!harmonic_powers || n_harmonics == 0) return -1;

    double freq_domain_sum = 0.0;
    for (size_t h = 0; h < n_harmonics; h++) {
        freq_domain_sum += harmonic_powers[h];
    }

    /* Parseval: total time-domain average power = sum of per-harmonic powers */
    if (fabs(time_domain_power - freq_domain_sum) > tolerance) return 1;
    return 0;
}

/* ==========================================================================
 * L5: K-Factor (Transformer Harmonic Derating)
 * ========================================================================== */

double hp_k_factor(const double *i_harmonics, size_t n_harmonics)
{
    if (!i_harmonics || n_harmonics == 0) return -1.0;

    double i1 = i_harmonics[0];
    if (i1 < HP_EPS) return -1.0;

    /* K-factor = Σ h² × (I_h/I₁)²
     *
     * The factor h² accounts for the fact that eddy-current losses
     * increase with the square of frequency. A pure sinusoidal load
     * has K-factor = 1.
     */
    double k_sum = 1.0; /* fundamental always contributes 1 × (I₁/I₁)² = 1 */
    for (size_t h = 1; h < n_harmonics; h++) {
        double order = (double)(h + 1); /* h=0 is fundamental (order 1) */
        double ratio = i_harmonics[h] / i1;
        k_sum += order * order * ratio * ratio;
    }

    return k_sum;
}

/* ==========================================================================
 * L5: Transformer Derating Factor
 * ========================================================================== */

double hp_transformer_derating(double k_factor, double alpha)
{
    if (k_factor < 0.0 || alpha < 0.0) return -1.0;

    /* Derating = P_rated_with_harmonics / P_rated_sinusoidal
     *
     * = sqrt( (1 + α)/(1 + α·K) )
     *
     * where α is eddy-current loss / I²R loss ratio at fundamental.
     * Typical α = 0.05-0.10 for dry-type transformers (IEEE C57.110).
     *
     * For K = 1 (pure sine): derating = 1.0 (no derating needed)
     * For K = 9 (heavy harmonics): derating ≈ 0.7-0.8
     */
    if (k_factor < 1.0) k_factor = 1.0;

    double num = 1.0 + alpha;
    double den = 1.0 + alpha * k_factor;
    if (den < HP_EPS) return 0.0;

    return sqrt(num / den);
}

/* ==========================================================================
 * L5: Initialize IEEE 519 Limits
 * ========================================================================== */

int hp_init_ieee519_limits(hp_limits_t *limits, hp_standard_t standard,
                            double i_sc, double i_l_max, double v_nominal)
{
    if (!limits) return -1;

    memset(limits, 0, sizeof(*limits));
    limits->standard = standard;
    limits->i_sc     = i_sc;
    limits->i_l_max  = i_l_max;

    /* IEEE 519-2014 Table 1: Current distortion limits
     *
     * Voltage level affects which table applies:
     *   V ≤ 1.0 kV (low voltage): Table 1 — most stringent
     *   1.0 kV < V ≤ 69 kV: Table 2
     *   69 kV < V ≤ 161 kV: Table 3
     *   V > 161 kV: Table 4
     *
     * This implementation uses Table 1 (≤1kV) as default.
     * For higher voltage systems, limits are relaxed.
     */
    double voltage_factor = 1.0;
    if (v_nominal > 69000.0) voltage_factor = 2.5;
    else if (v_nominal > 1000.0) voltage_factor = 1.5;

    /* I_sc / I_L ratio determines the allowable TDD:
     */
    double ratio = (i_l_max > 0.0) ? (i_sc / i_l_max) : 1000.0;

    if (ratio < 20.0) {
        limits->tdd_limit_percent = 5.0 * voltage_factor;
    } else if (ratio < 50.0) {
        limits->tdd_limit_percent = 8.0 * voltage_factor;
    } else if (ratio < 100.0) {
        limits->tdd_limit_percent = 12.0 * voltage_factor;
    } else if (ratio < 1000.0) {
        limits->tdd_limit_percent = 15.0 * voltage_factor;
    } else {
        limits->tdd_limit_percent = 20.0 * voltage_factor;
    }

    /* Individual harmonic limits (odd, non-triplen, for <11th):
     * For ratio < 20: 4.0%
     * For 20-50: 7.0%
     * etc.
     */
    for (int h = 3; h < 50; h += 2) {
        if (h % 3 == 0) continue; /* triplen harmonics: separate limits */

        if (ratio < 20.0) {
            limits->individual_limit_odd[h] = (h < 11) ? 4.0 :
                                              (h < 17) ? 2.0 :
                                              (h < 23) ? 1.5 :
                                              (h < 35) ? 0.6 : 0.3;
        } else if (ratio < 50.0) {
            limits->individual_limit_odd[h] = (h < 11) ? 7.0 :
                                              (h < 17) ? 3.5 :
                                              (h < 23) ? 2.5 :
                                              (h < 35) ? 1.0 : 0.5;
        } else if (ratio < 100.0) {
            limits->individual_limit_odd[h] = (h < 11) ? 10.0 :
                                              (h < 17) ? 4.5 :
                                              (h < 23) ? 4.0 :
                                              (h < 35) ? 1.5 : 0.7;
        } else {
            limits->individual_limit_odd[h] = (h < 11) ? 12.0 :
                                              (h < 17) ? 5.5 :
                                              (h < 23) ? 5.0 :
                                              (h < 35) ? 2.0 : 1.0;
        }
    }

    /* Even harmonics: 25% of odd limit */
    for (int h = 2; h < 50; h += 2) {
        limits->individual_limit_even[h] = limits->individual_limit_odd[h+1] * 0.25;
    }

    return 0;
}

/* ==========================================================================
 * L5: Check IEEE 519 Compliance
 * ========================================================================== */

int hp_check_ieee519(const hp_spectrum_t *spectrum,
                      const hp_limits_t *limits, int *violations)
{
    if (!spectrum || !limits || !violations) return -1;

    *violations = 0;

    /* Check TDD */
    if (spectrum->tdd_i_percent > limits->tdd_limit_percent) {
        (*violations)++;
    }

    /* Check individual harmonics */
    for (uint32_t i = 0; i < spectrum->num_harmonics; i++) {
        int order = (int)(spectrum->harmonics[i].order);
        if (order >= 50 || order <= 1) continue;

        double limit;
        if (order % 2 == 0) {
            limit = limits->individual_limit_even[order];
        } else {
            limit = limits->individual_limit_odd[order];
        }

        if (limit > 0.0) {
            double i_percent = spectrum->harmonics[i].i_mag
                                / spectrum->i1_mag * 100.0;
            if (i_percent > limit) {
                (*violations)++;
            }
        }
    }

    if (*violations > 0) return 1;
    return 0;
}

/* ==========================================================================
 * L6: Identify Dominant Harmonics
 * ========================================================================== */

int hp_dominant_harmonics(const hp_spectrum_t *spectrum, int top_n,
                           int *orders_out)
{
    if (!spectrum || !orders_out || top_n <= 0) return -1;

    /* Simple selection sort to find top-N by current magnitude */
    typedef struct { int order; double mag; } hp_pair_t;
    hp_pair_t ranked[HP_MAX_HARMONICS];
    uint32_t count = 0;

    for (uint32_t i = 0; i < spectrum->num_harmonics; i++) {
        ranked[count].order = (int)(spectrum->harmonics[i].order);
        ranked[count].mag   = spectrum->harmonics[i].i_mag;
        count++;
    }

    /* Sort descending */
    for (uint32_t i = 0; i < count; i++) {
        uint32_t max_idx = i;
        for (uint32_t j = i + 1; j < count; j++) {
            if (ranked[j].mag > ranked[max_idx].mag) {
                max_idx = j;
            }
        }
        hp_pair_t tmp = ranked[i];
        ranked[i] = ranked[max_idx];
        ranked[max_idx] = tmp;
    }

    int n_out = (int)(top_n < (int)count ? top_n : (int)count);
    for (int i = 0; i < n_out; i++) {
        orders_out[i] = ranked[i].order;
    }

    return n_out;
}

/* ==========================================================================
 * L6: Harmonic Mitigation Recommendations
 * ========================================================================== */

int hp_recommend_mitigation(const hp_spectrum_t *spectrum,
                             char *recommendation, size_t buf_len)
{
    if (!spectrum || !recommendation || buf_len == 0) return -1;

    int dominant[5];
    int n_dom = hp_dominant_harmonics(spectrum, 5, dominant);

    char temp[512];
    int written = snprintf(temp, sizeof(temp),
        "Harmonic Mitigation Analysis:\n"
        "  THDv = %.1f%%,  THDi = %.1f%%\n"
        "  PF_displacement = %.3f,  PF_true = %.3f\n",
        spectrum->thd_v_percent, spectrum->thd_i_percent,
        spectrum->pf_displacement, spectrum->pf_true);

    if (n_dom > 0) {
        written += snprintf(temp + written, sizeof(temp) - written,
            "  Dominant harmonics: ");
        for (int i = 0; i < n_dom; i++) {
            written += snprintf(temp + written, sizeof(temp) - written,
                "%d%s", dominant[i], (i < n_dom - 1) ? ", " : "\n");
        }
    }

    /* Recommendations based on dominant harmonics */
    int has_5th = 0, has_7th = 0, has_11th = 0;
    for (int i = 0; i < n_dom; i++) {
        if (dominant[i] == 5) has_5th = 1;
        if (dominant[i] == 7) has_7th = 1;
        if (dominant[i] == 11) has_11th = 1;
    }

    if (has_5th || has_7th) {
        written += snprintf(temp + written, sizeof(temp) - written,
            "  Recommend: Passive tuned filter at 5th/7th harmonic.\n");
    }
    if (has_11th) {
        written += snprintf(temp + written, sizeof(temp) - written,
            "  Recommend: Passive filter or 12-pulse rectifier.\n");
    }
    if (spectrum->pf_displacement < 0.90) {
        written += snprintf(temp + written, sizeof(temp) - written,
            "  Recommend: PF correction capacitors with detuning reactors.\n");
    }
    if (spectrum->thd_i_percent > 8.0) {
        written += snprintf(temp + written, sizeof(temp) - written,
            "  Recommend: Active harmonic filter (AHF) for dynamic loads.\n");
    }

    snprintf(recommendation, buf_len, "%s", temp);
    return 0;
}
