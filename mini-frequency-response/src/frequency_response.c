/**
 * @file frequency_response.c
 * @brief Core frequency response computation implementation
 *
 * Implements the fundamental frequency response operations:
 * - Data structure allocation and management
 * - Frequency grid generation (linear and logarithmic)
 * - Basic conversions (Hz↔rad/s, magnitude↔dB, phase computation)
 * - Frequency response analysis (cutoff, bandwidth, peaking, roll-off)
 *
 * All implementations are numerically robust and handle edge cases.
 */

#include "frequency_response.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Memory management
 * ============================================================================ */

freq_response_t *freq_response_alloc(size_t num_points, double f_start,
                                      double f_end, freq_scale_t scale)
{
    if (num_points == 0 || f_start <= 0.0 || f_end <= f_start) {
        return NULL;
    }

    freq_response_t *resp = (freq_response_t *)calloc(1, sizeof(freq_response_t));
    if (!resp) return NULL;

    resp->points = (freq_point_t *)calloc(num_points, sizeof(freq_point_t));
    if (!resp->points) {
        free(resp);
        return NULL;
    }

    resp->num_points = num_points;
    resp->freq_start = f_start;
    resp->freq_end = f_end;
    resp->scale = scale;
    snprintf(resp->label, sizeof(resp->label), "FreqResp_%.1f-%.1fHz",
             f_start, f_end);

    return resp;
}

void freq_response_free(freq_response_t *resp)
{
    if (resp) {
        free(resp->points);
        free(resp);
    }
}

/* ============================================================================
 * Frequency grid generation
 * ============================================================================ */

/**
 * Logarithmic frequency spacing algorithm:
 *
 * For N points spanning D = log₁₀(f_end/f_start) decades:
 *   f[k] = f_start × 10^{k·D/(N-1)}    for k = 0, 1, ..., N-1
 *
 * Alternative specification via points-per-decade:
 *   N_total = pts_per_dec × ceil(D) + 1
 *   f[k] = f_start × 10^{k/pts_per_dec}
 *
 * This function uses the points-per-decade form.
 */
double *freq_logspace(double f_start, double f_end, int pts_per_dec,
                       size_t *num_points)
{
    if (f_start <= 0.0 || f_end <= f_start || pts_per_dec < 1) {
        if (num_points) *num_points = 0;
        return NULL;
    }

    double decades = log10(f_end / f_start);
    /* Number of complete decades, at least 1 */
    size_t n_dec = (size_t)ceil(decades);
    if (n_dec < 1) n_dec = 1;

    /* Total points: pts_per_dec points per decade, plus the endpoint */
    size_t n = (size_t)(n_dec * pts_per_dec) + 1;
    double *freqs = (double *)malloc(n * sizeof(double));
    if (!freqs) {
        if (num_points) *num_points = 0;
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        double exponent = (double)i / (double)pts_per_dec;
        freqs[i] = f_start * pow(10.0, exponent);
        /* Clamp the last point to f_end to avoid overshoot */
        if (i == n - 1 && freqs[i] > f_end) {
            freqs[i] = f_end;
        }
    }

    if (num_points) *num_points = n;
    return freqs;
}

/* ============================================================================
 * Basic frequency/angular frequency conversion
 * ============================================================================ */

double freq_hz_to_rad(double freq_hz)
{
    return 2.0 * M_PI * freq_hz;
}

double freq_rad_to_hz(double omega)
{
    return omega / (2.0 * M_PI);
}

/* ============================================================================
 * Magnitude/dB conversion
 * ============================================================================ */

double magnitude_to_db(double magnitude)
{
    if (magnitude <= 0.0) {
        return -INFINITY;
    }
    return 20.0 * log10(magnitude);
}

double db_to_magnitude(double db)
{
    /* Handle -infinity dB gracefully */
    if (db <= -300.0) {  /* -300 dB is practically zero */
        return 0.0;
    }
    return pow(10.0, db / 20.0);
}

/* ============================================================================
 * Phase computation
 * ============================================================================ */

double phase_degrees(double real, double imag)
{
    return atan2(imag, real) * (180.0 / M_PI);
}

/**
 * Phase unwrapping algorithm:
 *
 * The raw phase from atan2 has discontinuities of ±360° at the
 * [-180°, 180°] boundary. Unwrapping removes these by detecting
 * jumps exceeding a tolerance (typically 180° = π rad) and adding
 * or subtracting multiples of 360° to ensure continuity.
 *
 * Algorithm (Itoh, 1982):
 * For k = 1, 2, ..., N-1:
 *   Δ = φ[k] - φ[k-1]
 *   if Δ > threshold: add correction -= 360°, repeat if needed
 *   if Δ < -threshold: add correction += 360°, repeat if needed
 *   φ_unwrapped[k] = φ[k] + correction
 *
 * This is a 1D phase unwrapping — 2D unwrapping (for images/surfaces)
 * is significantly more complex due to path dependence.
 *
 * Reference: Itoh, "Analysis of the phase unwrapping algorithm"
 * (1982), Applied Optics
 */
double *phase_unwrap(const double *phase_deg, size_t n, double tolerance)
{
    if (!phase_deg || n == 0) return NULL;

    double *unwrapped = (double *)malloc(n * sizeof(double));
    if (!unwrapped) return NULL;

    unwrapped[0] = phase_deg[0];
    double cumulative_correction = 0.0;

    for (size_t i = 1; i < n; i++) {
        double diff = phase_deg[i] - phase_deg[i - 1];

        /* Check for negative jump (e.g., +179° → -179°, diff = -358°) */
        while (diff < -tolerance) {
            cumulative_correction += 360.0;
            diff += 360.0;
        }
        /* Check for positive jump (e.g., -179° → +179°, diff = +358°) */
        while (diff > tolerance) {
            cumulative_correction -= 360.0;
            diff -= 360.0;
        }

        unwrapped[i] = phase_deg[i] + cumulative_correction;
    }

    return unwrapped;
}

/**
 * Group delay computation using central difference:
 *
 * τ_g(ω_k) = -dφ/dω|_{ω=ω_k}
 *
 * Forward difference at k=0:     -(φ₁ - φ₀)/(ω₁ - ω₀)
 * Central difference at 0<k<N-1: -(φ_{k+1} - φ_{k-1})/(ω_{k+1} - ω_{k-1})
 * Backward difference at k=N-1:  -(φ_{N-1} - φ_{N-2})/(ω_{N-1} - ω_{N-2})
 *
 * Input phase must be unwrapped (continuous) in radians.
 * Input omega must be in rad/s.
 *
 * Constant group delay (linear phase) is the defining characteristic
 * of Bessel filters and is essential for preserving pulse shapes in
 * digital communications.
 */
double *group_delay(const double *phase_rad, const double *omega, size_t n)
{
    if (!phase_rad || !omega || n < 2) return NULL;

    double *tau = (double *)malloc(n * sizeof(double));
    if (!tau) return NULL;

    /* Forward difference at first point */
    tau[0] = -(phase_rad[1] - phase_rad[0]) / (omega[1] - omega[0]);

    /* Central difference for interior points */
    for (size_t i = 1; i < n - 1; i++) {
        tau[i] = -(phase_rad[i + 1] - phase_rad[i - 1])
                  / (omega[i + 1] - omega[i - 1]);
    }

    /* Backward difference at last point */
    tau[n - 1] = -(phase_rad[n - 1] - phase_rad[n - 2])
                  / (omega[n - 1] - omega[n - 2]);

    return tau;
}

/* ============================================================================
 * Frequency sweep analysis
 * ============================================================================ */

/**
 * Find cutoff frequency using linear interpolation.
 *
 * For a lowpass response, the -3 dB cutoff is where |H| = |H_max|/√2,
 * i.e., magnitude_dB = max_dB - 3.01 dB.
 *
 * For a highpass response, the cutoff is where the response rises
 * to within 3 dB of the high-frequency (passband) value.
 *
 * Searches for the first crossing point and interpolates between
 * the two adjacent frequency points for sub-bin accuracy.
 */
double find_cutoff_freq(const freq_response_t *resp, int is_highpass)
{
    if (!resp || resp->num_points < 2) return -1.0;

    /* Find the reference level (passband magnitude) */
    double ref_db;
    if (is_highpass) {
        /* Highpass: reference is the high-frequency asymptote */
        ref_db = resp->points[resp->num_points - 1].magnitude_db;
    } else {
        /* Lowpass: reference is the DC / low-frequency value */
        ref_db = resp->points[0].magnitude_db;
    }

    double cutoff_db = ref_db - 3.01029995664;  /* -10·log₁₀(2) */

    for (size_t i = 0; i < resp->num_points - 1; i++) {
        double db_a = resp->points[i].magnitude_db;
        double db_b = resp->points[i + 1].magnitude_db;

        /* Check if cutoff lies between these two points */
        if ((db_a >= cutoff_db && db_b <= cutoff_db) ||
            (db_a <= cutoff_db && db_b >= cutoff_db)) {

            /* Linear interpolation in dB vs log-frequency */
            double log_f_a = log10(resp->points[i].frequency);
            double log_f_b = log10(resp->points[i + 1].frequency);

            if (fabs(db_b - db_a) < 1e-12) {
                return resp->points[i].frequency;  /* Flat segment */
            }

            double alpha = (cutoff_db - db_a) / (db_b - db_a);
            double log_f_cutoff = log_f_a + alpha * (log_f_b - log_f_a);
            return pow(10.0, log_f_cutoff);
        }
    }

    return -1.0;  /* Cutoff not found in sweep range */
}

/**
 * Find bandwidth of bandpass/bandstop response.
 *
 * For a bandpass response, the -3 dB bandwidth is the frequency
 * range between the two half-power points.
 *
 * For a bandstop response, this finds the -3 dB rejection bandwidth.
 *
 * Algorithm:
 * 1. Find the peak (for BP) or dip (for BS) magnitude
 * 2. Find where magnitude drops to -3 dB of peak (BP) or rises
 *    to -3 dB of passband (BS)
 */
double find_bandwidth(const freq_response_t *resp,
                       double *f_low, double *f_high)
{
    if (!resp || resp->num_points < 3) return -1.0;

    /* Find maximum and minimum magnitude in dB */
    double max_db = -INFINITY;
    double min_db = INFINITY;
    size_t max_idx = 0;

    for (size_t i = 0; i < resp->num_points; i++) {
        if (resp->points[i].magnitude_db > max_db) {
            max_db = resp->points[i].magnitude_db;
            max_idx = i;
        }
        if (resp->points[i].magnitude_db < min_db) {
            min_db = resp->points[i].magnitude_db;
        }
    }

    double ref_db = max_db;
    double threshold_db = ref_db - 3.01029995664;

    /* Search left of peak for lower -3 dB frequency */
    double fl = -1.0, fh = -1.0;

    for (size_t i = max_idx; i > 0; i--) {
        double db_a = resp->points[i - 1].magnitude_db;
        double db_b = resp->points[i].magnitude_db;
        if ((db_a <= threshold_db && db_b >= threshold_db)) {
            double log_fa = log10(resp->points[i - 1].frequency);
            double log_fb = log10(resp->points[i].frequency);
            double alpha = (threshold_db - db_a) / (db_b - db_a);
            double log_fl = log_fa + alpha * (log_fb - log_fa);
            fl = pow(10.0, log_fl);
            break;
        }
    }

    /* Search right of peak for upper -3 dB frequency */
    for (size_t i = max_idx; i < resp->num_points - 1; i++) {
        double db_a = resp->points[i].magnitude_db;
        double db_b = resp->points[i + 1].magnitude_db;
        if ((db_a >= threshold_db && db_b <= threshold_db)) {
            double log_fa = log10(resp->points[i].frequency);
            double log_fb = log10(resp->points[i + 1].frequency);
            double alpha = (threshold_db - db_a) / (db_b - db_a);
            double log_fh = log_fa + alpha * (log_fb - log_fa);
            fh = pow(10.0, log_fh);
            break;
        }
    }

    if (f_low) *f_low = fl;
    if (f_high) *f_high = fh;

    if (fl > 0.0 && fh > 0.0) {
        return fh - fl;
    }
    return -1.0;
}

/**
 * Find peak magnitude in frequency response.
 *
 * Simply scans all points and returns the maximum magnitude.
 * For more accurate peak location (sub-bin), quadratic interpolation
 * around the maximum point can be used (parabolic peak fitting).
 *
 * Quadratic peak interpolation (for dB magnitude):
 *   f_peak = f_k + (Δf/2)·(M_{k-1} - M_{k+1})/(M_{k-1} - 2M_k + M_{k+1})
 *   M_peak = M_k - (M_{k-1} - M_{k+1})²/(8(M_{k-1} - 2M_k + M_{k+1}))
 * where f_k is the frequency at the discrete maximum M_k, and Δf
 * is the frequency spacing (assumed uniform on log scale).
 */
double find_peak_response(const freq_response_t *resp, double *peak_freq)
{
    if (!resp || resp->num_points == 0) return -1.0;

    double max_mag = resp->points[0].magnitude;
    size_t max_idx = 0;

    for (size_t i = 1; i < resp->num_points; i++) {
        if (resp->points[i].magnitude > max_mag) {
            max_mag = resp->points[i].magnitude;
            max_idx = i;
        }
    }

    /* Quadratic interpolation for sub-bin peak accuracy */
    if (max_idx > 0 && max_idx < resp->num_points - 1) {
        double M_minus = resp->points[max_idx - 1].magnitude;
        double M_0     = resp->points[max_idx].magnitude;
        double M_plus  = resp->points[max_idx + 1].magnitude;

        double denom = M_minus - 2.0 * M_0 + M_plus;
        if (fabs(denom) > 1e-15) {
            double log_f_minus = log10(resp->points[max_idx - 1].frequency);
            double log_f_0     = log10(resp->points[max_idx].frequency);
            double delta_log_f = log_f_0 - log_f_minus;

            double offset = 0.5 * delta_log_f * (M_minus - M_plus) / denom;
            double log_f_peak = log_f_0 + offset;

            if (peak_freq) {
                *peak_freq = pow(10.0, log_f_peak);
            }
            /* Interpolated peak magnitude */
            double peak_mag = M_0 - (M_minus - M_plus) * (M_minus - M_plus)
                                     / (8.0 * denom);
            return peak_mag;
        }
    }

    if (peak_freq) *peak_freq = resp->points[max_idx].frequency;
    return max_mag;
}

/**
 * Compute roll-off rate using linear regression on dB vs log-freq.
 *
 * For a lowpass stopband:
 *   dB(f) = m·log₁₀(f) + b
 *   roll-off = m dB/decade
 *
 * Linear regression formula:
 *   m = (N·Σxᵢyᵢ - Σxᵢ·Σyᵢ) / (N·Σxᵢ² - (Σxᵢ)²)
 *
 * where xᵢ = log₁₀(fᵢ), yᵢ = magnitude_dBᵢ.
 */
double compute_rolloff(const freq_response_t *resp,
                        double f_low, double f_high)
{
    if (!resp || resp->num_points < 2) return 0.0;

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < resp->num_points; i++) {
        double f = resp->points[i].frequency;
        if (f >= f_low && f <= f_high) {
            double x = log10(f);
            double y = resp->points[i].magnitude_db;
            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_xx += x * x;
            count++;
        }
    }

    if (count < 2) return 0.0;

    double denom = (double)count * sum_xx - sum_x * sum_x;
    if (fabs(denom) < 1e-15) return 0.0;

    double slope = ((double)count * sum_xy - sum_x * sum_y) / denom;
    return slope;
}

/**
 * Interpolate frequency response at arbitrary frequency using
 * linear interpolation in log-frequency, dB-magnitude space.
 *
 * This is the most natural interpolation for Bode plots because
 * both axes (dB and log₁₀(f)) are linear on a Bode plot grid.
 */
freq_point_t freq_response_interpolate(const freq_response_t *resp,
                                         double freq)
{
    freq_point_t result;
    memset(&result, 0, sizeof(freq_point_t));

    if (!resp || resp->num_points < 2 || freq <= 0.0) return result;

    /* Find bracketing points */
    if (freq < resp->points[0].frequency || 
        freq > resp->points[resp->num_points - 1].frequency) {
        return result;  /* Out of range */
    }

    size_t lo = 0, hi = resp->num_points - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (resp->points[mid].frequency <= freq) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    /* Interpolate in log-frequency space */
    double log_flo = log10(resp->points[lo].frequency);
    double log_fhi = log10(resp->points[hi].frequency);
    double log_f   = log10(freq);

    double alpha = (log_f - log_flo) / (log_fhi - log_flo);

    result.frequency = freq;
    result.angular_freq = 2.0 * M_PI * freq;
    result.magnitude_db = resp->points[lo].magnitude_db
                          + alpha * (resp->points[hi].magnitude_db
                                     - resp->points[lo].magnitude_db);
    result.magnitude = pow(10.0, result.magnitude_db / 20.0);

    /* For phase, interpolate unwrapped values if available, otherwise raw */
    double phase_lo = resp->points[lo].phase_rad;
    double phase_hi = resp->points[hi].phase_rad;

    /* Handle phase wrap-around in interpolation */
    double phase_diff = phase_hi - phase_lo;
    if (phase_diff > M_PI) phase_diff -= 2.0 * M_PI;
    if (phase_diff < -M_PI) phase_diff += 2.0 * M_PI;

    result.phase_rad = phase_lo + alpha * phase_diff;
    result.real = result.magnitude * cos(result.phase_rad);
    result.imag = result.magnitude * sin(result.phase_rad);

    return result;
}
