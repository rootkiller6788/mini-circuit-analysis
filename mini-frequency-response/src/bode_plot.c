/**
 * @file bode_plot.c
 * @brief Bode plot construction and analysis implementation
 *
 * Implements exact and asymptotic Bode plot construction,
 * corner frequency detection, gain-bandwidth product computation,
 * resonance peaking detection, Bode's gain-phase relationship
 * verification, and visualization helpers.
 *
 * Reference: Bode, "Network Analysis and Feedback Amplifier Design" (1945)
 * Course: MIT 6.003, Stanford EE102A, Berkeley EE16B
 */

#include "bode_plot.h"
#include "frequency_response.h"
#include "transfer_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Bode Plot Construction
 * ============================================================================ */

bode_plot_t *bode_compute(const tf_polynomial_t *tf,
                           double f_start, double f_end,
                           int pts_per_dec)
{
    if (!tf || f_start <= 0.0 || f_end <= f_start || pts_per_dec < 1) {
        return NULL;
    }

    /* Generate frequency grid */
    size_t n_freq = 0;
    double *freqs = freq_logspace(f_start, f_end, pts_per_dec, &n_freq);
    if (!freqs || n_freq == 0) return NULL;

    bode_plot_t *bode = (bode_plot_t *)calloc(1, sizeof(bode_plot_t));
    if (!bode) { free(freqs); return NULL; }

    /* Allocate magnitude response data */
    bode->magnitude_response.points = (freq_point_t *)calloc(n_freq,
        sizeof(freq_point_t));
    bode->phase_response.points = (freq_point_t *)calloc(n_freq,
        sizeof(freq_point_t));

    if (!bode->magnitude_response.points || !bode->phase_response.points) {
        bode_free(bode);
        free(freqs);
        return NULL;
    }

    bode->magnitude_response.num_points = n_freq;
    bode->magnitude_response.freq_start = f_start;
    bode->magnitude_response.freq_end = f_end;
    bode->magnitude_response.scale = FREQ_SCALE_LOGARITHMIC;

    bode->phase_response.num_points = n_freq;
    bode->phase_response.freq_start = f_start;
    bode->phase_response.freq_end = f_end;
    bode->phase_response.scale = FREQ_SCALE_LOGARITHMIC;

    /* Evaluate H(jω) at each frequency */
    for (size_t i = 0; i < n_freq; i++) {
        double f = freqs[i];
        double omega = 2.0 * M_PI * f;
        double _Complex H = tf_evaluate_freq(tf, omega);

        double mag = cabs(H);
        double mag_db = magnitude_to_db(mag);
        double phase = atan2(cimag(H), creal(H));

        /* Magnitude data */
        bode->magnitude_response.points[i].frequency = f;
        bode->magnitude_response.points[i].angular_freq = omega;
        bode->magnitude_response.points[i].magnitude = mag;
        bode->magnitude_response.points[i].magnitude_db = mag_db;
        bode->magnitude_response.points[i].phase_rad = phase;
        bode->magnitude_response.points[i].real = creal(H);
        bode->magnitude_response.points[i].imag = cimag(H);

        /* Phase data */
        bode->phase_response.points[i].frequency = f;
        bode->phase_response.points[i].angular_freq = omega;
        bode->phase_response.points[i].magnitude = mag;
        bode->phase_response.points[i].magnitude_db = mag_db;
        bode->phase_response.points[i].phase_rad = phase;
        bode->phase_response.points[i].real = creal(H);
        bode->phase_response.points[i].imag = cimag(H);
    }

    /* Compute DC gain */
    bode->dc_gain_db = magnitude_to_db(tf_dc_gain(tf));

    /* Compute high-frequency slope */
    double hf_slope = 0.0;
    if (n_freq >= 2) {
        size_t last = n_freq - 1;
        double df_log = log10(freqs[last]) - log10(freqs[last - 1]);
        double dmag = bode->magnitude_response.points[last].magnitude_db
                       - bode->magnitude_response.points[last - 1].magnitude_db;
        if (fabs(df_log) > 1e-15) {
            hf_slope = dmag / df_log;
        }
    }
    bode->hf_slope_db_dec = hf_slope;

    free(freqs);
    return bode;
}

/**
 * Construct asymptotic Bode plot from pole-zero data.
 *
 * The asymptotic (straight-line) approximation is:
 *
 * Magnitude:
 * - Start with flat line at 20·log₁₀|K| dB.
 * - At each zero frequency ω_z, slope increases by +20 dB/dec.
 * - At each pole frequency ω_p, slope decreases by -20 dB/dec.
 * - Pole/zero at origin (s=0): slope of ±20 dB/dec from ω=0.
 *
 * Phase:
 * - Each zero: phase rises from 0° to +90° over ω_z/10 to 10ω_z.
 *   At ω_z, phase = +45°.
 * - Each pole: phase falls from 0° to -90° over ω_p/10 to 10ω_p.
 *   At ω_p, phase = -45°.
 * - Multiple poles/zeros: phases add.
 * - Pole at origin: constant -90° phase.
 *
 * For complex poles with damping ζ:
 * - Magnitude has peaking: |H|_peak ≈ 1/(2ζ√(1-ζ²)) for ζ < 1/√2.
 * - Phase transition is steeper for lower ζ.
 */
bode_plot_t *bode_asymptotic(const tf_pole_zero_t *pz,
                              double f_start, double f_end,
                              int pts_per_dec)
{
    if (!pz || f_start <= 0.0 || f_end <= f_start || pts_per_dec < 1) {
        return NULL;
    }

    size_t n_freq = 0;
    double *freqs = freq_logspace(f_start, f_end, pts_per_dec, &n_freq);
    if (!freqs || n_freq == 0) return NULL;

    bode_plot_t *bode = (bode_plot_t *)calloc(1, sizeof(bode_plot_t));
    if (!bode) { free(freqs); return NULL; }

    bode->magnitude_response.points = (freq_point_t *)calloc(n_freq,
        sizeof(freq_point_t));
    bode->phase_response.points = (freq_point_t *)calloc(n_freq,
        sizeof(freq_point_t));

    if (!bode->magnitude_response.points || !bode->phase_response.points) {
        bode_free(bode);
        free(freqs);
        return NULL;
    }

    bode->magnitude_response.num_points = n_freq;
    bode->magnitude_response.freq_start = f_start;
    bode->magnitude_response.freq_end = f_end;
    bode->magnitude_response.scale = FREQ_SCALE_LOGARITHMIC;
    bode->phase_response.num_points = n_freq;
    bode->phase_response.freq_start = f_start;
    bode->phase_response.freq_end = f_end;
    bode->phase_response.scale = FREQ_SCALE_LOGARITHMIC;

    /* Count poles and zeros at origin */
    int poles_at_origin = 0;
    int zeros_at_origin = 0;

    for (size_t i = 0; i < pz->num_poles; i++) {
        if (cabs(pz->poles[i]) < 1e-12) poles_at_origin++;
    }
    for (size_t i = 0; i < pz->num_zeros; i++) {
        if (cabs(pz->zeros[i]) < 1e-12) zeros_at_origin++;
    }

    /* Initial slope from poles/zeros at origin:
     * n·pole at origin → -20n dB/dec, m·zero at origin → +20m dB/dec */
    double origin_slope = -20.0 * poles_at_origin + 20.0 * zeros_at_origin;
    double origin_phase = -90.0 * poles_at_origin + 90.0 * zeros_at_origin;

    /* DC gain in dB */
    double K_db = 20.0 * log10(fabs(pz->gain));
    bode->dc_gain_db = K_db;

    /* Compute magnitude and phase at each frequency */
    double ref_freq = 1.0;  /* Reference frequency for origin slope */
    for (size_t i = 0; i < n_freq; i++) {
        double f = freqs[i];
        double w = 2.0 * M_PI * f;

        /* Start with DC gain + origin contributions */
        double mag_db = K_db + origin_slope * log10(f / ref_freq);
        double phase_deg = origin_phase;

        /* Add non-origin pole contributions */
        for (size_t j = 0; j < pz->num_poles; j++) {
            double wp = cabs(pz->poles[j]);
            if (wp < 1e-12) continue;  /* Skip poles at origin */

            /* Magnitude: -20 dB/dec above wp */
            if (w > wp) {
                mag_db -= 20.0 * log10(w / wp);
            }

            /* Phase: -45°/dec from wp/10 to 10*wp, total -90° */
            if (w >= wp / 10.0 && w <= wp * 10.0) {
                phase_deg -= 45.0 * log10(w / (wp / 10.0));
            } else if (w > wp * 10.0) {
                phase_deg -= 90.0;
            }
        }

        /* Add zero contributions */
        for (size_t j = 0; j < pz->num_zeros; j++) {
            double wz = cabs(pz->zeros[j]);
            if (wz < 1e-12) continue;

            if (w > wz) {
                mag_db += 20.0 * log10(w / wz);
            }

            if (w >= wz / 10.0 && w <= wz * 10.0) {
                phase_deg += 45.0 * log10(w / (wz / 10.0));
            } else if (w > wz * 10.0) {
                phase_deg += 90.0;
            }
        }

        double mag_lin = pow(10.0, mag_db / 20.0);
        double phase_rad = phase_deg * M_PI / 180.0;

        bode->magnitude_response.points[i].frequency = f;
        bode->magnitude_response.points[i].angular_freq = w;
        bode->magnitude_response.points[i].magnitude = mag_lin;
        bode->magnitude_response.points[i].magnitude_db = mag_db;
        bode->magnitude_response.points[i].phase_rad = phase_rad;
        bode->magnitude_response.points[i].real = mag_lin * cos(phase_rad);
        bode->magnitude_response.points[i].imag = mag_lin * sin(phase_rad);

        bode->phase_response.points[i].frequency = f;
        bode->phase_response.points[i].angular_freq = w;
        bode->phase_response.points[i].magnitude = mag_lin;
        bode->phase_response.points[i].magnitude_db = mag_db;
        bode->phase_response.points[i].phase_rad = phase_rad;
        bode->phase_response.points[i].real = mag_lin * cos(phase_rad);
        bode->phase_response.points[i].imag = mag_lin * sin(phase_rad);
    }

    bode->hf_slope_db_dec = origin_slope
                            - 20.0 * (pz->num_poles - poles_at_origin)
                            + 20.0 * (pz->num_zeros - zeros_at_origin);

    free(freqs);
    return bode;
}

void bode_free(bode_plot_t *bode)
{
    if (bode) {
        free(bode->magnitude_response.points);
        free(bode->phase_response.points);
        free(bode->corners);
        free(bode);
    }
}

/* ============================================================================
 * Bode Plot Analysis
 * ============================================================================ */

double bode_dominant_pole_freq(const bode_plot_t *bode)
{
    if (!bode || bode->magnitude_response.num_points < 2) return -1.0;

    const freq_response_t *resp = &bode->magnitude_response;
    double dc_db = resp->points[0].magnitude_db;
    double cutoff_db = dc_db - 3.0103;

    for (size_t i = 1; i < resp->num_points; i++) {
        if (resp->points[i].magnitude_db <= cutoff_db) {
            /* Interpolate */
            double log_f1 = log10(resp->points[i - 1].frequency);
            double log_f2 = log10(resp->points[i].frequency);
            double db1 = resp->points[i - 1].magnitude_db;
            double db2 = resp->points[i].magnitude_db;
            double alpha = (cutoff_db - db1) / (db2 - db1);
            return pow(10.0, log_f1 + alpha * (log_f2 - log_f1));
        }
    }
    return -1.0;
}

double *bode_find_corner_frequencies(const bode_plot_t *bode,
                                       size_t *num_corners)
{
    if (!bode || bode->magnitude_response.num_points < 3) {
        if (num_corners) *num_corners = 0;
        return NULL;
    }

    const freq_response_t *resp = &bode->magnitude_response;
    size_t n = resp->num_points;

    /* Detect slope changes by computing second difference of dB magnitude.
     * A corner frequency occurs where the slope changes significantly. */

    /* First compute slopes between adjacent points */
    double *slopes = (double *)malloc((n - 1) * sizeof(double));
    if (!slopes) { if (num_corners) *num_corners = 0; return NULL; }

    for (size_t i = 0; i < n - 1; i++) {
        double dlogf = log10(resp->points[i + 1].frequency)
                        - log10(resp->points[i].frequency);
        double ddb = resp->points[i + 1].magnitude_db
                      - resp->points[i].magnitude_db;
        slopes[i] = (fabs(dlogf) > 1e-15) ? (ddb / dlogf) : 0.0;
    }

    /* Detect changes in slope indicating corner frequencies */
    size_t max_corners = n;
    double *corners = (double *)malloc(max_corners * sizeof(double));
    if (!corners) { free(slopes); if (num_corners) *num_corners = 0; return NULL; }

    size_t count = 0;
    for (size_t i = 1; i < n - 1 && count < max_corners; i++) {
        double slope_change = slopes[i] - slopes[i - 1];
        /* Significant slope change (> 5 dB/dec) indicates a corner */
        if (fabs(slope_change) > 5.0) {
            corners[count++] = resp->points[i].frequency;
            i++;  /* Skip next point to avoid double-detection */
        }
    }

    free(slopes);
    if (num_corners) *num_corners = count;

    if (count == 0) {
        free(corners);
        return NULL;
    }
    return corners;
}

double bode_gain_bandwidth_product(const bode_plot_t *bode)
{
    if (!bode) return -1.0;

    /* GBWP = DC gain × dominant pole frequency */
    double dc_gain_lin = pow(10.0, bode->dc_gain_db / 20.0);
    double dom_pole = bode_dominant_pole_freq(bode);
    if (dom_pole <= 0.0) return -1.0;

    return dc_gain_lin * dom_pole;
}

double bode_slope_in_range(const bode_plot_t *bode,
                            double f_low, double f_high)
{
    if (!bode) return 0.0;
    return compute_rolloff(&bode->magnitude_response, f_low, f_high);
}

/**
 * Detect resonance peaking in Bode magnitude plot.
 *
 * Peaking occurs when a complex pole pair has ζ < 1/√2 ≈ 0.707.
 * The peak magnitude above the low-frequency asymptote is detected
 * by finding a local maximum that exceeds neighboring values.
 *
 * For a second-order lowpass:
 *   |H|_peak_dB = 10·log₁₀(Q²/(1 - 1/(4Q²))) for Q > 1/√2
 *   ≈ 20·log₁₀(Q) for Q >> 1
 */
int bode_detect_peaking(const bode_plot_t *bode,
                         double *peak_freq, double *peak_db)
{
    if (!bode || bode->magnitude_response.num_points < 5) return 0;

    const freq_response_t *resp = &bode->magnitude_response;
    size_t n = resp->num_points;

    /* Baseline: DC gain (first point) */
    double baseline_db = resp->points[0].magnitude_db;

    int found = 0;
    double max_above = 0.0;
    double peak_f = 0.0;

    for (size_t i = 2; i < n - 2; i++) {
        /* Check if this point is a local maximum */
        double db_m1 = resp->points[i - 1].magnitude_db;
        double db_0  = resp->points[i].magnitude_db;
        double db_p1 = resp->points[i + 1].magnitude_db;

        if (db_0 > db_m1 && db_0 > db_p1 && db_0 > baseline_db) {
            /* Local peak detected */
            double excursion = db_0 - baseline_db;
            if (excursion > 0.5 && excursion > max_above) {  /* > 0.5 dB threshold */
                max_above = excursion;
                peak_f = resp->points[i].frequency;
                found = 1;
            }
        }
    }

    if (found) {
        if (peak_freq) *peak_freq = peak_f;
        if (peak_db) *peak_db = max_above;
    }
    return found;
}

/**
 * Bode's gain-phase relationship (Bode, 1945, Ch. 14):
 *
 * For a minimum-phase transfer function, the phase at frequency ω₀
 * is determined by the magnitude slope:
 *
 *   φ(ω₀) ≈ (π/2) · n  where n = d(log|H|)/d(log ω) at ω=ω₀
 *
 * More precisely, if |H(jω)| ∝ ωⁿ over a wide frequency range
 * centered at ω₀, then φ(ω₀) ≈ n·90°.
 *
 * This function computes the approximate phase from the local
 * magnitude slope and compares with the actual phase, returning
 * the RMS error which should be small for minimum-phase systems.
 */
int bode_gain_phase_relation(const bode_plot_t *bode,
                              double *approx_phase, double *error_rms)
{
    if (!bode || bode->magnitude_response.num_points < 3) return -1;

    const freq_response_t *mag_resp = &bode->magnitude_response;
    const freq_response_t *phase_resp = &bode->phase_response;
    size_t n = mag_resp->num_points;

    double *approx = (double *)malloc(n * sizeof(double));
    if (!approx) return -1;

    double sum_sq_error = 0.0;
    size_t count = 0;

    for (size_t i = 1; i < n - 1; i++) {
        /* Compute local magnitude slope n = d(log|H|)/d(log ω) */
        double dlog_mag = mag_resp->points[i + 1].magnitude_db
                           - mag_resp->points[i - 1].magnitude_db;
        double dlog_w = log10(mag_resp->points[i + 1].frequency)
                         - log10(mag_resp->points[i - 1].frequency);

        if (fabs(dlog_w) < 1e-15) continue;

        /* dB slope / 20 = d(log10|H|)/d(log10 ω) = n */
        double n_slope = dlog_mag / (20.0 * dlog_w);

        /* Approximate phase: n × 90° */
        approx[i] = n_slope * 90.0;

        /* Actual unwrapped phase at this point */
        double actual_phase = phase_resp->points[i].phase_rad * 180.0 / M_PI;

        double error = approx[i] - actual_phase;
        sum_sq_error += error * error;
        count++;
    }

    if (approx_phase && count > 0) {
        /* Return the approximation at the midpoint */
        *approx_phase = approx[n / 2];
    }
    if (error_rms && count > 0) {
        *error_rms = sqrt(sum_sq_error / count);
    }

    free(approx);
    return 0;
}

/* ============================================================================
 * Visualization Helpers
 * ============================================================================ */

int bode_ascii_plot(const bode_plot_t *bode, char *buffer,
                     size_t buf_size, int plot_width)
{
    if (!bode || !buffer || buf_size < 128 || plot_width < 20) return -1;

    const freq_response_t *resp = &bode->magnitude_response;
    size_t n = resp->num_points;

    /* Find dB range */
    double db_min = INFINITY, db_max = -INFINITY;
    for (size_t i = 0; i < n; i++) {
        if (resp->points[i].magnitude_db > db_max
            && isfinite(resp->points[i].magnitude_db))
            db_max = resp->points[i].magnitude_db;
        if (resp->points[i].magnitude_db < db_min
            && isfinite(resp->points[i].magnitude_db))
            db_min = resp->points[i].magnitude_db;
    }

    if (!isfinite(db_min) || !isfinite(db_max)) return -1;

    /* Build ASCII plot */
    size_t offset = 0;
    int step = (n > 50) ? (int)(n / 30) : 1;

    for (size_t i = 0; i < n; i += step) {
        if (offset + 80 >= buf_size) break;

        double db = resp->points[i].magnitude_db;
        if (!isfinite(db)) continue;

        int pos = (int)((db - db_min) / (db_max - db_min) * (plot_width - 1));
        if (pos < 0) pos = 0;
        if (pos >= plot_width) pos = plot_width - 1;

        offset += snprintf(buffer + offset, buf_size - offset,
                           "%8.1f Hz |", resp->points[i].frequency);
        for (int j = 0; j < pos; j++) {
            if (offset + 1 < buf_size) buffer[offset++] = ' ';
        }
        if (offset < buf_size) buffer[offset++] = '*';
        if (offset < buf_size) buffer[offset++] = '\n';
    }

    if (offset < buf_size) buffer[offset] = '\0';
    return (int)offset;
}

int bode_export_csv(const bode_plot_t *bode, char *buffer, size_t size)
{
    if (!bode || !buffer || size < 64) return -1;

    const freq_response_t *mag = &bode->magnitude_response;
    const freq_response_t *phase = &bode->phase_response;
    size_t n = mag->num_points;

    size_t offset = 0;
    offset += snprintf(buffer + offset, size - offset,
                       "frequency_hz,magnitude_db,phase_deg\n");

    for (size_t i = 0; i < n; i++) {
        if (offset + 80 >= size) break;
        offset += snprintf(buffer + offset, size - offset,
                           "%.6e,%.6f,%.6f\n",
                           mag->points[i].frequency,
                           mag->points[i].magnitude_db,
                           phase->points[i].phase_rad * 180.0 / M_PI);
    }

    return (int)offset;
}
