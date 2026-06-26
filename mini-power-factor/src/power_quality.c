/**
 * @file power_quality.c
 * @brief Power quality assessment — ITIC/CBEMA curves, flicker, reliability indices
 *
 * Power quality encompasses all deviations from the ideal sinusoidal
 * voltage waveform: sags, swells, interruptions, harmonics, flicker,
 * unbalance, and transients. This module implements the standard
 * assessment methods from IEEE 1159-2019 and IEC 61000-4-30.
 *
 * Knowledge points:
 *   L1: IEEE 1159 event types and classification
 *   L2: ITIC (CBEMA) voltage tolerance curve evaluation
 *   L4: Flicker severity Pst/Plt per IEC 61000-4-15
 *   L4: Reliability indices SAIFI/SAIDI/MAIFI per IEEE 1366
 *   L5: IEC 61000-4-30 Class A PQ aggregator
 *   L7: Data center PQ assessment (80 PLUS, Energy Star)
 *   L7: EV charging station PQ compliance (SAE J2894)
 *   L7: Cost of poor power quality (EPRI methodology)
 *
 * References:
 *   - IEEE Std 1159-2019
 *   - IEEE Std 1459-2010
 *   - IEEE Std 1366-2012
 *   - IEC 61000-4-15 (Flickermeter)
 *   - IEC 61000-4-30 (PQ Measurement)
 *   - ITIC (CBEMA) Curve, rev. 2000
 *   - EPRI, "Cost of Power Quality" studies
 */

#include "power_quality.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define PQ_EPS 1e-12

/* ==========================================================================
 * L2: ITIC (CBEMA) Curve Upper Envelope
 * ========================================================================== */

double pq_itic_upper_envelope(double duration_sec)
{
    if (duration_sec < 0.0) return 0.0;

    /* ITIC upper envelope (revised 2000):
     *   t ≤ 1 ms:         500%
     *   1 ms < t ≤ 10 ms: 200% (ramp down from 500 to 200)
     *   10 ms < t ≤ 0.5 s: 200% to 120%
     *   0.5 s < t ≤ 10 s: 120% to 110%
     *   t > 10 s:         110%
     *
     * Using log-linear interpolation between breakpoints.
     */

    if (duration_sec <= 0.001) {
        return 500.0;
    } else if (duration_sec <= 0.01) {
        /* Log interpolation: 500% @ 1ms → 200% @ 10ms */
        double log_t    = log10(duration_sec);
        double log_t1   = log10(0.001);
        double log_t2   = log10(0.01);
        double frac = (log_t - log_t1) / (log_t2 - log_t1);
        return 500.0 + frac * (200.0 - 500.0);
    } else if (duration_sec <= 0.5) {
        double log_t    = log10(duration_sec);
        double log_t1   = log10(0.01);
        double log_t2   = log10(0.5);
        double frac = (log_t - log_t1) / (log_t2 - log_t1);
        return 200.0 + frac * (120.0 - 200.0);
    } else if (duration_sec <= 10.0) {
        double log_t    = log10(duration_sec);
        double log_t1   = log10(0.5);
        double log_t2   = log10(10.0);
        double frac = (log_t - log_t1) / (log_t2 - log_t1);
        return 120.0 + frac * (110.0 - 120.0);
    } else {
        return 110.0;
    }
}

/* ==========================================================================
 * L2: ITIC Lower Envelope
 * ========================================================================== */

double pq_itic_lower_envelope(double duration_sec)
{
    if (duration_sec < 0.0) return 0.0;

    /* ITIC lower envelope (revised 2000):
     *   t ≤ 0.02 s (1 cycle @ 60Hz): 0%
     *   0.02 s < t ≤ 0.5 s: 0% to 70%
     *   0.5 s < t ≤ 10 s:   70% to 80%
     *   t > 10 s:           90%
     *
     * Events below the lower envelope cause equipment malfunction.
     * Events that are 0% for > 0.5 s are considered sustained interruptions.
     */

    if (duration_sec <= 0.02) {
        return 0.0;
    } else if (duration_sec <= 0.5) {
        double log_t    = log10(duration_sec);
        double log_t1   = log10(0.02);
        double log_t2   = log10(0.5);
        double frac = (log_t - log_t1) / (log_t2 - log_t1);
        return 0.0 + frac * 70.0;
    } else if (duration_sec <= 10.0) {
        double log_t    = log10(duration_sec);
        double log_t1   = log10(0.5);
        double log_t2   = log10(10.0);
        double frac = (log_t - log_t1) / (log_t2 - log_t1);
        return 70.0 + frac * (80.0 - 70.0);
    } else {
        return 90.0;
    }
}

/* ==========================================================================
 * L2: ITIC Curve Evaluation
 * ========================================================================== */

pq_itic_result_t pq_check_itic(double duration_sec, double v_percent)
{
    if (duration_sec < 0.0 || v_percent < 0.0) return PQ_ITIC_BOUNDARY;

    double upper = pq_itic_upper_envelope(duration_sec);
    double lower = pq_itic_lower_envelope(duration_sec);

    if (v_percent > upper) return PQ_ITIC_ABOVE_UPPER;
    if (v_percent < lower && duration_sec > 0.02) return PQ_ITIC_BELOW_LOWER;

    /* Zero voltage for < 1 cycle is acceptable per ITIC */
    if (v_percent < lower && duration_sec <= 0.02) return PQ_ITIC_ACCEPTABLE;

    return PQ_ITIC_ACCEPTABLE;
}

/* ==========================================================================
 * L4: Short-Term Flicker Severity Pst (IEC 61000-4-15)
 * ========================================================================== */

double pq_compute_pst(const double *v_rms_samples, size_t n_samples,
                       double nominal_v)
{
    if (!v_rms_samples || n_samples < 10 || nominal_v < PQ_EPS) return -1.0;

    /* Standard Pst computation involves:
     * 1. Demodulate voltage to extract flicker envelope
     * 2. Weight envelope by lamp-eye-brain transfer function
     * 3. Statistical analysis of instantaneous flicker sensation
     *
     * Simplified method using voltage variance as proxy:
     * Pst ≈ sqrt(0.0314·P0.1 + 0.0525·P1 + 0.0657·P3 + 0.28·P10 + 0.08·P50)
     *
     * where Px is the x-th percentile of the flicker sensation over 10 min.
     *
     * For coarse estimation: use normalized voltage variance.
     */

    double sum = 0.0, sum_sq = 0.0;
    double max_dev = 0.0;

    for (size_t i = 0; i < n_samples; i++) {
        double dev_pu = (v_rms_samples[i] - nominal_v) / nominal_v;
        sum += dev_pu;
        sum_sq += dev_pu * dev_pu;
        if (fabs(dev_pu) > max_dev) max_dev = fabs(dev_pu);
    }

    double mean_dev = sum / (double)n_samples;
    double var = sum_sq / (double)n_samples - mean_dev * mean_dev;
    if (var < 0.0) var = 0.0;

    /* IEC 61000-4-15: Pst ∝ voltage fluctuation magnitude.
     * A ΔV/V = 0.5% at 8.8 Hz modulation → Pst ≈ 1.0 (irritability threshold)
     *
     * Simple model: Pst ≈ 200 × σ (where σ = std dev of V/V_nom)
     */
    double pst = 200.0 * sqrt(var);

    /* Clamp to reasonable range */
    if (pst < 0.0) pst = 0.0;
    if (pst > 20.0) pst = 20.0;

    return pst;
}

/* ==========================================================================
 * L4: Long-Term Flicker Severity Plt
 * ========================================================================== */

double pq_compute_plt(const double *pst_values)
{
    if (!pst_values) return -1.0;

    /* Plt = ³√( (1/12) Σ Pst_i³ )
     *
     * Computed over 2 hours from 12 consecutive 10-min Pst values.
     */
    double sum_cube = 0.0;
    int valid = 0;

    for (int i = 0; i < 12; i++) {
        if (pst_values[i] >= 0.0) {
            sum_cube += pst_values[i] * pst_values[i] * pst_values[i];
            valid++;
        }
    }

    if (valid == 0) return -1.0;

    double avg_cube = sum_cube / (double)valid;
    return cbrt(avg_cube);
}

/* ==========================================================================
 * L4: SAIFI — System Average Interruption Frequency Index
 * ========================================================================== */

double pq_compute_saifi(const uint32_t *customers_per_event,
                         size_t n_events, uint32_t total_customers)
{
    if (!customers_per_event || total_customers == 0) return -1.0;

    double sum = 0.0;
    for (size_t i = 0; i < n_events; i++) {
        sum += (double)customers_per_event[i];
    }

    return sum / (double)total_customers;
}

/* ==========================================================================
 * L4: SAIDI — System Average Interruption Duration Index
 * ========================================================================== */

double pq_compute_saidi(const double *customer_min_per_event,
                         size_t n_events, uint32_t total_customers)
{
    if (!customer_min_per_event || total_customers == 0) return -1.0;

    double sum = 0.0;
    for (size_t i = 0; i < n_events; i++) {
        sum += customer_min_per_event[i];
    }

    return sum / (double)total_customers;
}

/* ==========================================================================
 * L4: MAIFI — Momentary Average Interruption Frequency Index
 * ========================================================================== */

double pq_compute_maifi(const uint32_t *momentary_events,
                         size_t n_groups, uint32_t total_customers)
{
    if (!momentary_events || total_customers == 0) return -1.0;

    double sum = 0.0;
    for (size_t i = 0; i < n_groups; i++) {
        sum += (double)momentary_events[i];
    }

    return sum / (double)total_customers;
}

/* ==========================================================================
 * L7: Data Center Power Quality Assessment
 * ========================================================================== */

int pq_assess_datacenter(double pf_measured, double thd_i_percent,
                          double v_deviation_pct, double efficiency_pct,
                          int *tier, double *score)
{
    if (!tier || !score) return -1;

    int issues = 0;
    *score = 100.0;

    /* Check PF */
    if (pf_measured < 0.95) { issues |= 1; *score -= 10.0; }
    if (pf_measured < 0.90) { issues |= 1; *score -= 15.0; }

    /* Check THD_i */
    if (thd_i_percent > 5.0) { issues |= 2; *score -= 15.0; }
    if (thd_i_percent > 8.0) { issues |= 2; *score -= 15.0; }

    /* Check voltage deviation */
    if (v_deviation_pct > 10.0) { issues |= 4; *score -= 15.0; }
    if (v_deviation_pct > 15.0) { issues |= 4; *score -= 10.0; }

    /* 80 PLUS efficiency tiers:
     *   Titanium (4): ≥ 96% at 50% load
     *   Platinum (3): ≥ 94%
     *   Gold     (2): ≥ 92%
     *   Silver   (1): ≥ 88%
     *   Bronze   (0): ≥ 85%
     *   None     (6): < 85%
     */
    if (efficiency_pct >= 96.0) *tier = 4;
    else if (efficiency_pct >= 94.0) *tier = 3;
    else if (efficiency_pct >= 92.0) *tier = 2;
    else if (efficiency_pct >= 88.0) *tier = 1;
    else if (efficiency_pct >= 85.0) *tier = 0;
    else *tier = 6;

    if (*tier < 2) { issues |= 8; *score -= 10.0; }

    if (*score < 0.0) *score = 0.0;
    return issues;
}

/* ==========================================================================
 * L7: EV Charger Power Quality Assessment (SAE J2894)
 * ========================================================================== */

int pq_assess_ev_charger(double pf_rated, double thd_i_percent,
                          double v_sag_remaining, double sag_duration_ms,
                          int *compliance)
{
    if (!compliance) return -1;

    *compliance = 1; /* start compliant */

    /* SAE J2894 / IEC 61851 requirements:
     * - PF ≥ 0.95 at rated power
     * - THD_i ≤ 5%
     * - Ride through: 70% voltage for 500 ms
     * - Ride through: 0% voltage for 10 ms (1/2 cycle at 50Hz)
     */

    if (pf_rated < 0.95) *compliance = 0;
    if (thd_i_percent > 5.0) *compliance = 0;

    /* Voltage sag ride-through check */
    if (v_sag_remaining < 70.0 && sag_duration_ms > 10.0) {
        /* Below 70% for more than half-cycle → must ride through? */
        /* Actually the requirement is the charger must RIDE THROUGH these,
         * meaning it should keep operating. We check if the sag is severe. */
        if (v_sag_remaining < 10.0 && sag_duration_ms > 10.0) {
            *compliance = 0; /* interruption → allowed to disconnect */
        }
    }

    /* If voltage goes to zero for > 500ms, it's not a sag but an interruption */
    if (v_sag_remaining < 10.0 && sag_duration_ms > 500.0) {
        *compliance = 1; /* interruption: allowed to trip */
    }

    return 0;
}

/* ==========================================================================
 * L7: Cost of Poor Power Quality (EPRI Methodology)
 * ========================================================================== */

double pq_cost_of_poor_quality(double avg_pf, double annual_kwh,
                                double cost_per_kwh, double pf_penalty_rate,
                                int sag_events_yr, double cost_per_sag)
{
    if (annual_kwh <= 0.0) return 0.0;

    double total_cost = 0.0;

    /* 1. I²R loss penalty from low PF
     * Loss increase factor = (1/PF² - 1)
     * Typical I²R loss ≈ 2% of total energy at PF=1.0
     */
    double pf_loss_factor = 0.0;
    if (avg_pf > 0.0) {
        pf_loss_factor = (1.0 / (avg_pf * avg_pf)) - 1.0;
    }
    double base_loss_kwh = 0.02 * annual_kwh;
    double extra_loss_kwh = base_loss_kwh * pf_loss_factor;
    total_cost += extra_loss_kwh * cost_per_kwh;

    /* 2. Utility PF penalty charges
     * Typical penalty: $X per kVAR per month for PF < 0.85
     */
    if (avg_pf < 0.85) {
        double avg_kw = annual_kwh / 8760.0;
        double avg_kva = (avg_pf > 0.0) ? (avg_kw / avg_pf) : avg_kw;
        double avg_kvar = sqrt(avg_kva * avg_kva - avg_kw * avg_kw);
        if (avg_kvar < 0.0) avg_kvar = 0.0;
        total_cost += avg_kvar * pf_penalty_rate * 12.0;
    }

    /* 3. Production downtime from voltage sags */
    total_cost += (double)sag_events_yr * cost_per_sag;

    /* 4. Harmonic losses (simplified: 0.5% of energy for THD > 5%) */
    /* For now, embedded in I²R factor */

    return total_cost;
}

/* ==========================================================================
 * L5: Power Quality Report Generation
 * ========================================================================== */

int pq_generate_report(const pq_monitor_t *monitor, char *report_buf,
                        size_t buf_len)
{
    if (!monitor || !report_buf || buf_len == 0) return -1;

    int written = snprintf(report_buf, buf_len,
        "========================================\n"
        "  POWER QUALITY MONITORING REPORT\n"
        "  Standard: IEC 61000-4-30 Class A\n"
        "========================================\n\n");

    written += snprintf(report_buf + written,
        (buf_len > (size_t)written) ? buf_len - written : 0,
        "Steady-State Metrics (10-min avg):\n"
        "  Vrms        = %.1f V\n"
        "  THDv        = %.2f %%\n"
        "  Frequency   = %.3f Hz\n"
        "  Power Factor = %.3f\n"
        "  Unbalance   = %.2f %%\n"
        "  Pst (flicker) = %.3f\n\n",
        monitor->v_rms_10min_avg,
        monitor->thd_v_10min_avg,
        monitor->freq_10min_avg,
        monitor->pf_10min_avg,
        monitor->unbalance_10min_avg,
        monitor->flicker_pst);

    if (monitor->flicker_plt > 0.8) {
        written += snprintf(report_buf + written,
            (buf_len > (size_t)written) ? buf_len - written : 0,
            "  *** WARNING: Plt = %.3f exceeds limit (0.8) ***\n\n",
            monitor->flicker_plt);
    }

    written += snprintf(report_buf + written,
        (buf_len > (size_t)written) ? buf_len - written : 0,
        "Event Summary (24 hours):\n"
        "  Sags:           %u\n"
        "  Swells:         %u\n"
        "  Interruptions:  %u\n"
        "  Total events:   %u\n\n",
        monitor->sag_count_24h,
        monitor->swell_count_24h,
        monitor->interruption_count_24h,
        monitor->event_count);

    written += snprintf(report_buf + written,
        (buf_len > (size_t)written) ? buf_len - written : 0,
        "========================================\n");

    return 0;
}
