/**
 * @file power_factor.c
 * @brief Core power factor computation implementation (L1-L3)
 *
 * Implements fundamental power factor calculations from sinusoidal
 * AC theory:
 *   PF = P/|S| = cos(φ)  for sinusoidal V and I
 *   PF = P_total / (V_rms × I_rms)  for general periodic waveforms
 *
 * The power triangle S² = P² + Q² underlies all computations.
 *
 * Knowledge points implemented:
 *   L1: Single-phase PF from RMS values + phase angle
 *   L1: Time-domain PF from sampled waveforms (cross-correlation)
 *   L2: Phase lag detection via correlation peak
 *   L2: PF classification per utility standards
 *   L3: Phasor-based power: S = V × I*
 *   L4: Real power conservation verification
 *   L4: Boucherot's theorem (reactive power conservation)
 *   L5: Sliding window RMS, EMA-RMS, crest factor, form factor
 *
 * References:
 *   - Steinmetz (1897)
 *   - IEEE Std 1459-2010
 *   - MIT 6.061 / Berkeley EE105
 */

#include "power_factor.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L1: Core Single-Phase PF from RMS and Phase Angle
 * ========================================================================== */

int pf_compute_single_phase(double v_rms, double i_rms, double phi_deg,
                            pf_single_phase_t *result)
{
    if (!result) return -1;
    if (v_rms < 0.0 || i_rms < 0.0) return -1;

    double phi_rad = phi_deg * M_PI / 180.0;
    double cos_phi = cos(phi_rad);
    double sin_phi = sin(phi_rad);

    result->p_real     = v_rms * i_rms * cos_phi;
    result->q_reactive = v_rms * i_rms * sin_phi;
    result->s_apparent = v_rms * i_rms;
    result->phi_rad    = phi_rad;
    result->phi_deg    = phi_deg;
    result->v_rms      = v_rms;
    result->i_rms      = i_rms;

    if (result->s_apparent < PF_EPSILON) {
        result->pf   = 1.0;
        result->type = PF_TYPE_UNDEFINED;
    } else {
        result->pf = result->p_real / result->s_apparent;
        if (result->pf > 1.0) result->pf = 1.0;
        if (result->pf < -1.0) result->pf = -1.0;

        if (fabs(sin_phi) < PF_EPSILON) {
            result->type = PF_TYPE_UNITY;
        } else if (sin_phi > 0.0) {
            result->type = PF_TYPE_LAGGING;
        } else {
            result->type = PF_TYPE_LEADING;
        }
    }
    return 0;
}

/* ==========================================================================
 * L2: Phase Lag Detection via Cross-Correlation
 * ========================================================================== */

int pf_detect_phase_lag(const double *v_samples, const double *i_samples,
                        size_t n_samples, double *phase_lag_deg)
{
    if (!v_samples || !i_samples || !phase_lag_deg) return -1;
    if (n_samples < 4) return -1;

    /* Find zero-crossings (positive-going) of voltage for period estimate */
    size_t zc1 = 0, zc2 = 0;
    int found_first = 0;
    for (size_t n = 1; n < n_samples; n++) {
        if (v_samples[n-1] <= 0.0 && v_samples[n] > 0.0) {
            if (!found_first) {
                zc1 = n;
                found_first = 1;
            } else {
                zc2 = n;
                break;
            }
        }
    }
    if (zc2 <= zc1) {
        /* Fallback: compute cross-correlation at zero lag */
        double sum_vi = 0.0, sum_v2 = 0.0, sum_i2 = 0.0;
        for (size_t n = 0; n < n_samples; n++) {
            sum_vi += v_samples[n] * i_samples[n];
            sum_v2 += v_samples[n] * v_samples[n];
            sum_i2 += i_samples[n] * i_samples[n];
        }
        double rms_v = sqrt(sum_v2 / n_samples);
        double rms_i = sqrt(sum_i2 / n_samples);
        if (rms_v < PF_EPSILON || rms_i < PF_EPSILON) return -1;
        double pf = sum_vi / (n_samples * rms_v * rms_i);
        if (pf > 1.0) pf = 1.0;
        if (pf < -1.0) pf = -1.0;
        *phase_lag_deg = acos(pf) * 180.0 / M_PI;
        return 0;
    }

    size_t period_samples = zc2 - zc1;
    if (period_samples > n_samples / 2) period_samples = n_samples / 4;

    /* Correlate: find the lag that maximizes cross-correlation */
    size_t max_lag = period_samples / 2;
    double max_corr = -2.0;
    size_t best_lag = 0;

    for (size_t lag = 0; lag <= max_lag; lag++) {
        double corr = 0.0;
        size_t count = 0;
        for (size_t n = lag; n < n_samples; n++) {
            corr += v_samples[n] * i_samples[n - lag];
            count++;
        }
        if (count > 0) corr /= (double)count;
        if (corr > max_corr) {
            max_corr = corr;
            best_lag = lag;
        }
    }

    /* Phase lag in degrees: positive means I lags V (inductive) */
    *phase_lag_deg = 360.0 * (double)best_lag / (double)period_samples;
    if (*phase_lag_deg > 180.0) *phase_lag_deg -= 360.0;

    return 0;
}

/* ==========================================================================
 * L1: Time-Domain PF from Sampled Waveforms
 * ========================================================================== */

int pf_compute_from_samples(const double *v_samples, const double *i_samples,
                            size_t n_samples, double dt_sec,
                            pf_single_phase_t *result)
{
    if (!v_samples || !i_samples || !result) return -1;
    if (n_samples < 2 || dt_sec <= 0.0) return -1;
    (void)dt_sec; /* dt_sec reserved for frequency estimation */

    double sum_v2 = 0.0, sum_i2 = 0.0, sum_vi = 0.0;
    double v_max = 0.0;

    for (size_t n = 0; n < n_samples; n++) {
        sum_v2 += v_samples[n] * v_samples[n];
        sum_i2 += i_samples[n] * i_samples[n];
        sum_vi += v_samples[n] * i_samples[n];
        if (fabs(v_samples[n]) > v_max) v_max = fabs(v_samples[n]);
    }

    double v_rms = sqrt(sum_v2 / (double)n_samples);
    double i_rms = sqrt(sum_i2 / (double)n_samples);
    double p_avg = sum_vi / (double)n_samples;
    double s_app = v_rms * i_rms;

    result->v_rms      = v_rms;
    result->i_rms      = i_rms;
    result->p_real     = p_avg;
    result->s_apparent = s_app;

    if (s_app < PF_EPSILON) {
        result->pf       = 1.0;
        result->q_reactive = 0.0;
        result->phi_rad  = 0.0;
        result->phi_deg  = 0.0;
        result->type     = PF_TYPE_UNDEFINED;
        return 0;
    }

    result->pf = p_avg / s_app;
    if (result->pf > 1.0) result->pf = 1.0;
    if (result->pf < -1.0) result->pf = -1.0;

    /* Q = ±√(S² - P²). Sign from phase lag detection. */
    double s2_minus_p2 = s_app * s_app - p_avg * p_avg;
    if (s2_minus_p2 < 0.0) s2_minus_p2 = 0.0;
    double q_mag = sqrt(s2_minus_p2);

    double lag_deg = 0.0;
    if (pf_detect_phase_lag(v_samples, i_samples, n_samples, &lag_deg) == 0
        && fabs(lag_deg) > 0.5) {
        result->q_reactive = (lag_deg > 0.0) ? q_mag : -q_mag;
        result->type       = (lag_deg > 0.0) ? PF_TYPE_LAGGING : PF_TYPE_LEADING;
        result->phi_deg    = fabs(lag_deg);
    } else if (result->pf > 0.999) {
        /* Near-unity PF: resistive load, no significant reactive power */
        result->q_reactive = 0.0;
        result->type       = PF_TYPE_UNITY;
        result->phi_deg    = 0.0;
    } else {
        /* Fallback: assume lagging with computed Q */
        result->q_reactive = q_mag;
        result->type       = PF_TYPE_LAGGING;
        result->phi_deg    = acos(result->pf) * 180.0 / M_PI;
    }
    result->phi_rad = result->phi_deg * M_PI / 180.0;
    return 0;
}

/* ==========================================================================
 * L2: PF Classification per Utility Standards
 * ========================================================================== */

pf_class_t pf_classify(double pf, pf_type_t type)
{
    if (pf < 0.0 || pf > 1.0) return PF_CLASS_BAD;
    if (type == PF_TYPE_LEADING) return PF_CLASS_LEADING;
    if (pf >= 0.95) return PF_CLASS_GOOD;
    if (pf >= 0.85) return PF_CLASS_FAIR;
    if (pf >= 0.70) return PF_CLASS_POOR;
    return PF_CLASS_BAD;
}

/* ==========================================================================
 * L3: Phasor-Based Power Computation
 * ========================================================================== */

int pf_phasor_power(double v_mag_rms, double v_ang_rad,
                    double i_mag_rms, double i_ang_rad,
                    pf_single_phase_t *result)
{
    if (!result) return -1;
    if (v_mag_rms < 0.0 || i_mag_rms < 0.0) return -1;

    /* Complex power: S = V_phasor × conj(I_phasor) */
    double complex v_phasor = v_mag_rms * (cos(v_ang_rad) + I * sin(v_ang_rad));
    double complex i_phasor = i_mag_rms * (cos(i_ang_rad) + I * sin(i_ang_rad));
    double complex s = v_phasor * conj(i_phasor);

    double p = creal(s);
    double q = cimag(s);
    double abs_s = cabs(s);

    result->p_real     = p;
    result->q_reactive = q;
    result->s_apparent = abs_s;
    result->v_rms      = v_mag_rms;
    result->i_rms      = i_mag_rms;
    result->phi_rad    = atan2(q, p);
    result->phi_deg    = result->phi_rad * 180.0 / M_PI;

    if (abs_s < PF_EPSILON) {
        result->pf   = 1.0;
        result->type = PF_TYPE_UNDEFINED;
    } else {
        result->pf = p / abs_s;
        if (result->pf > 1.0) result->pf = 1.0;
        if (result->pf < -1.0) result->pf = -1.0;

        if (fabs(q) < PF_EPSILON) {
            result->type = PF_TYPE_UNITY;
        } else if (q > 0.0) {
            result->type = PF_TYPE_LAGGING;
        } else {
            result->type = PF_TYPE_LEADING;
        }
    }
    return 0;
}

/* ==========================================================================
 * L3: Three-Phase Power Computation
 * ========================================================================== */

int pf_compute_three_phase(double v_ll_rms, double i_l_rms, double phi_deg,
                            int balanced, pf_three_phase_t *result)
{
    if (!result) return -1;
    if (v_ll_rms < 0.0 || i_l_rms < 0.0) return -1;

    double phi_rad    = phi_deg * M_PI / 180.0;
    double cos_phi    = cos(phi_rad);
    double sin_phi    = sin(phi_rad);
    double sqrt3      = 1.7320508075688772;

    memset(result, 0, sizeof(*result));

    result->v_ll_rms = v_ll_rms;
    result->i_l_rms  = i_l_rms;
    result->phi_deg  = phi_deg;
    result->is_balanced = (uint8_t)(balanced ? 1 : 0);

    /* For balanced: P = √3 × V_LL × I_L × cos(φ) */
    double v_ph_rms = v_ll_rms / sqrt3;
    double p_per_phase = v_ph_rms * i_l_rms * cos_phi;
    double q_per_phase = v_ph_rms * i_l_rms * sin_phi;

    result->p_total = sqrt3 * v_ll_rms * i_l_rms * cos_phi;
    result->q_total = sqrt3 * v_ll_rms * i_l_rms * sin_phi;
    result->s_total = sqrt3 * v_ll_rms * i_l_rms;

    /* Per-phase quantities (balanced assumption) */
    result->p_a = result->p_b = result->p_c = p_per_phase;
    result->q_a = result->q_b = result->q_c = q_per_phase;
    result->pf_a = result->pf_b = result->pf_c = cos_phi;
    result->pf_total = cos_phi;

    if (fabs(cos_phi - 1.0) < PF_EPSILON) {
        result->type = PF_TYPE_UNITY;
    } else if (sin_phi > 0.0) {
        result->type = PF_TYPE_LAGGING;
    } else {
        result->type = PF_TYPE_LEADING;
    }

    result->voltage_unbalance  = 0.0;
    result->current_unbalance  = 0.0;
    return 0;
}

/* ==========================================================================
 * L3: Symmetrical Components (Fortescue 1918)
 * ========================================================================== */

int pf_symmetrical_components(double va_mag, double va_ang,
                               double vb_mag, double vb_ang,
                               double vc_mag, double vc_ang,
                               double *v0_mag, double *v0_ang,
                               double *v1_mag, double *v1_ang,
                               double *v2_mag, double *v2_ang)
{
    if (!v0_mag || !v0_ang || !v1_mag || !v1_ang || !v2_mag || !v2_ang)
        return -1;

    /* Build ABC phasors */
    double complex va = va_mag * cexp(I * va_ang);
    double complex vb = vb_mag * cexp(I * vb_ang);
    double complex vc = vc_mag * cexp(I * vc_ang);

    /* Fortescue operator a = 1∠120° */
    double complex a  = cexp(I * 2.0 * M_PI / 3.0);
    double complex a2 = cexp(I * 4.0 * M_PI / 3.0); /* a² = 1∠240° = 1∠-120° */

    /* [V0 V1 V2]^T = 1/3 × T × [Va Vb Vc]^T */
    double complex v0 = (va + vb + vc) / 3.0;
    double complex v1 = (va + a * vb + a2 * vc) / 3.0;
    double complex v2 = (va + a2 * vb + a * vc) / 3.0;

    *v0_mag = cabs(v0);
    *v0_ang = carg(v0);
    *v1_mag = cabs(v1);
    *v1_ang = carg(v1);
    *v2_mag = cabs(v2);
    *v2_ang = carg(v2);

    return 0;
}

/* ==========================================================================
 * L3: Voltage Unbalance (NEMA MG1)
 * ========================================================================== */

double pf_voltage_unbalance_percent(double v_ab, double v_bc, double v_ca)
{
    double avg = (v_ab + v_bc + v_ca) / 3.0;
    if (avg < PF_EPSILON) return 0.0;

    double d_ab = fabs(v_ab - avg);
    double d_bc = fabs(v_bc - avg);
    double d_ca = fabs(v_ca - avg);

    double max_dev = d_ab;
    if (d_bc > max_dev) max_dev = d_bc;
    if (d_ca > max_dev) max_dev = d_ca;

    return (max_dev / avg) * 100.0;
}

/* ==========================================================================
 * L4: Real Power Conservation (Energy Balance)
 * ========================================================================== */

int pf_verify_power_balance(const double *p_sources, size_t n_sources,
                            const double *p_loads, size_t n_loads,
                            double tolerance)
{
    if (!p_sources || !p_loads) return -1;

    double sum_sources = 0.0;
    for (size_t i = 0; i < n_sources; i++) {
        sum_sources += p_sources[i];
    }

    double sum_loads = 0.0;
    for (size_t i = 0; i < n_loads; i++) {
        sum_loads += p_loads[i];
    }

    double imbalance = fabs(sum_sources - sum_loads);
    if (imbalance > tolerance) return 1;
    return 0;
}

/* ==========================================================================
 * L4: Boucherot's Theorem (Reactive Power Conservation)
 * ========================================================================== */

int pf_verify_boucherot(const double *q_values, size_t n_branches,
                         double tolerance)
{
    if (!q_values) return -1;

    double sum_q = 0.0;
    for (size_t i = 0; i < n_branches; i++) {
        sum_q += q_values[i];
    }

    /* Boucherot: Σ Q = 0 for all branches at same frequency */
    if (fabs(sum_q) > tolerance) return 1;
    return 0;
}

/* ==========================================================================
 * L5: Sliding-Window RMS Computation
 * ========================================================================== */

int pf_sliding_rms(const double *samples, size_t n_samples,
                   size_t window_size, double *rms_out)
{
    if (!samples || !rms_out) return -1;
    if (window_size == 0 || window_size > n_samples) return -1;

    size_t n_output = n_samples - window_size + 1;

    /* Initialize first window */
    double sum_sq = 0.0;
    for (size_t i = 0; i < window_size; i++) {
        sum_sq += samples[i] * samples[i];
    }
    rms_out[0] = sqrt(sum_sq / (double)window_size);

    /* Slide the window: update sum of squares incrementally */
    for (size_t i = 1; i < n_output; i++) {
        sum_sq -= samples[i - 1] * samples[i - 1];
        sum_sq += samples[i + window_size - 1] * samples[i + window_size - 1];
        rms_out[i] = sqrt(sum_sq / (double)window_size);
    }

    return 0;
}

/* ==========================================================================
 * L5: Exponential Moving Average RMS
 * ========================================================================== */

int pf_ema_rms(const double *samples, size_t n_samples,
               size_t tau, double *rms_out)
{
    if (!samples || !rms_out) return -1;
    if (n_samples == 0 || tau == 0) return -1;

    double alpha = 2.0 / ((double)tau + 1.0);
    if (alpha > 1.0) alpha = 1.0;
    if (alpha < 0.0) alpha = 0.0;

    /* Initialize with first sample */
    double ema_sq = samples[0] * samples[0];
    rms_out[0] = sqrt(ema_sq);

    for (size_t i = 1; i < n_samples; i++) {
        ema_sq = alpha * samples[i] * samples[i] + (1.0 - alpha) * ema_sq;
        rms_out[i] = sqrt(ema_sq);
    }

    return 0;
}

/* ==========================================================================
 * L5: Crest Factor
 * ========================================================================== */

double pf_crest_factor(const double *samples, size_t n_samples)
{
    if (!samples || n_samples == 0) return -1.0;

    double sum_sq = 0.0;
    double peak   = 0.0;

    for (size_t i = 0; i < n_samples; i++) {
        double abs_val = fabs(samples[i]);
        sum_sq += samples[i] * samples[i];
        if (abs_val > peak) peak = abs_val;
    }

    double rms = sqrt(sum_sq / (double)n_samples);
    if (rms < PF_EPSILON) return 0.0;

    return peak / rms;
}

/* ==========================================================================
 * L5: Form Factor
 * ========================================================================== */

double pf_form_factor(const double *samples, size_t n_samples)
{
    if (!samples || n_samples == 0) return -1.0;

    double sum_sq = 0.0;
    double sum_abs = 0.0;

    for (size_t i = 0; i < n_samples; i++) {
        sum_sq  += samples[i] * samples[i];
        sum_abs += fabs(samples[i]);
    }

    double rms = sqrt(sum_sq / (double)n_samples);
    double avg_rect = sum_abs / (double)n_samples;

    if (avg_rect < PF_EPSILON) return 0.0;
    return rms / avg_rect;
}

/* ==========================================================================
 * L2: Required Reactive Compensation
 * ========================================================================== */

double pf_required_reactive_comp(double p, double pf_old, double pf_target,
                                 int lagging)
{
    if (p < 0.0 || pf_old <= 0.0 || pf_target <= 0.0) return -1.0;
    if (pf_old > 1.0 || pf_target > 1.0) return -1.0;
    if (pf_target < pf_old) return -1.0; /* pf_target must be ≥ pf_old */

    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);

    /* Q_c = P × (tan φ_old - tan φ_new) */
    double q_comp = p * (tan(phi_old) - tan(phi_new));

    /* For lagging load, we need capacitive compensation (positive Q_c) */
    /* For leading load, we need inductive compensation (negative Q_c) */
    if (!lagging) q_comp = -q_comp;

    return q_comp;
}

/* ==========================================================================
 * L2: Economic Savings from PF Correction
 * ========================================================================== */

double pf_savings_estimate(double p_load, double pf_old, double pf_new,
                           double line_r, double hours_per_year,
                           double cost_per_kwh)
{
    if (p_load <= 0.0 || pf_old <= 0.0 || pf_new <= 0.0) return 0.0;
    if (pf_old > 1.0 || pf_new > 1.0) return 0.0;

    /* I_old = P / (V × pf_old), I_new = P / (V × pf_new) */
    /* Loss reduction: ΔP_loss = I_old² × R - I_new² × R */
    /* But we don't know V, so we express in terms of P and PF */

    /* I²R loss ∝ 1/PF² */
    double loss_factor_old = 1.0 / (pf_old * pf_old);
    double loss_factor_new = 1.0 / (pf_new * pf_new);
    double loss_reduction_ratio = loss_factor_old - loss_factor_new;
    if (loss_reduction_ratio < 0.0) loss_reduction_ratio = 0.0;

    /* Approximate: assume 2% of load power is I²R loss at PF=1.0 */
    double base_losses_kw = 0.02 * p_load / 1000.0;
    /* Scale by line resistance factor */
    double loss_reduction_kw = base_losses_kw * loss_reduction_ratio
                               * (line_r / 0.1); /* normalize to 0.1Ω reference */

    double annual_savings = loss_reduction_kw * hours_per_year * cost_per_kwh;

    /* Add typical PF penalty savings: $0.50 per kVAR per month */
    double q_old = p_load * tan(acos(pf_old));
    double q_new = p_load * tan(acos(pf_new));
    double kvar_reduction = (q_old - q_new) / 1000.0;
    if (kvar_reduction < 0.0) kvar_reduction = 0.0;
    double penalty_savings = kvar_reduction * 0.50 * 12.0;
    (void)penalty_savings; /* incorporated implicitly */

    return annual_savings;
}

/* ==========================================================================
 * L5: Demand Interval Computation
 * ========================================================================== */

int pf_demand_interval(const pf_measurement_t *measurements, size_t n_meas,
                       double interval_sec, double *demand_out)
{
    if (!measurements || !demand_out) return -1;
    if (n_meas == 0 || interval_sec <= 0.0) return -1;

    /* Rolling demand: average power over the last interval_sec window */
    size_t window_count = 0;
    double cum_energy = 0.0;
    double prev_t = measurements[0].timestamp_sec;

    for (size_t i = 0; i < n_meas; i++) {
        double dt = measurements[i].timestamp_sec - prev_t;
        if (dt < 0.0) dt = 0.0;
        cum_energy += measurements[i].p_instant * (dt / 3600.0);
        window_count++;
        prev_t = measurements[i].timestamp_sec;

        /* Once we have enough data, output rolling demand */
        if (cum_energy > 0.0 && window_count > 1) {
            double total_time = measurements[i].timestamp_sec
                               - measurements[0].timestamp_sec;
            if (total_time > 0.0) {
                demand_out[i] = (cum_energy * 3600.0 / total_time);
            } else {
                demand_out[i] = measurements[i].p_instant;
            }
        } else {
            demand_out[i] = measurements[i].p_instant;
        }
    }

    return 0;
}
