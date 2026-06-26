/**
 * @file resonance.c
 * @brief Resonance analysis implementation
 *
 * Implements complete resonance analysis for RLC circuits:
 * - Series and parallel resonance computation
 * - Transfer function derivation
 * - Universal resonance curve
 * - Quality factor from bandwidth and energy
 * - Step response of resonant circuits
 * - Coupled resonators (magnetic coupling)
 * - Quartz crystal resonator model
 *
 * Reference: Hayt et al. (2019), Terman (1943), Sedra & Smith (2020)
 * Course: Berkeley EE16B, MIT 6.003, Stanford EE102A
 */

#include "resonance.h"
#include "frequency_response.h"
#include "transfer_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Series RLC Resonance
 * ============================================================================ */

resonance_result_t resonance_series(const rlc_params_t *params)
{
    resonance_result_t result;
    memset(&result, 0, sizeof(resonance_result_t));

    if (!params || params->L <= 0.0 || params->C <= 0.0) return result;

    double R = params->R;
    double L = params->L;
    double C = params->C;

    /* Natural (undamped) resonant frequency */
    double w0 = 1.0 / sqrt(L * C);
    double f0 = w0 / (2.0 * M_PI);

    result.resonant_freq_rad = w0;
    result.resonant_freq_hz = f0;

    /* Quality factor for series RLC */
    double Q = w0 * L / R;
    if (R <= 0.0) Q = INFINITY;  /* Lossless */

    result.quality_factor = Q;
    result.damping_factor = 1.0 / (2.0 * Q);

    /* Bandwidth */
    double BW = w0 / Q;
    result.bandwidth_hz = BW / (2.0 * M_PI);

    /* Half-power frequencies */
    if (Q > 0.5) {
        double sqrt_term = sqrt(1.0 + 1.0 / (4.0 * Q * Q));
        double fL = f0 * (sqrt_term - 1.0 / (2.0 * Q));
        double fH = f0 * (sqrt_term + 1.0 / (2.0 * Q));
        result.half_power_low_hz = fL;
        result.half_power_high_hz = fH;
    } else {
        /* Overdamped: no distinct half-power points */
        result.half_power_low_hz = 0.0;
        result.half_power_high_hz = 0.0;
    }

    /* Impedance at resonance: Z(jω₀) = R (purely resistive, minimum) */
    result.impedance_at_res = R;

    /* Peak magnitude depends on what's being measured.
     * For series RLC, the current is maximum at resonance:
     * |I_max| = V_s/R. The voltage across R is V_s.
     * Voltage across L or C: Q·V_s. */

    /* For the transfer function V_R/V_in (bandpass):
     * |H(jω₀)| = 1 (at resonance, V_R = V_in) */
    if (R > 0.0) {
        result.peak_magnitude = 1.0;
        result.peak_magnitude_db = 0.0;
    } else {
        result.peak_magnitude = INFINITY;
        result.peak_magnitude_db = INFINITY;
    }

    return result;
}

/* ============================================================================
 * Parallel RLC Resonance
 * ============================================================================ */

resonance_result_t resonance_parallel(const rlc_params_t *params)
{
    resonance_result_t result;
    memset(&result, 0, sizeof(resonance_result_t));

    if (!params || params->L <= 0.0 || params->C <= 0.0) return result;

    double R = params->R;
    double L = params->L;
    double C = params->C;

    /* Natural resonant frequency (same as series) */
    double w0 = 1.0 / sqrt(L * C);
    double f0 = w0 / (2.0 * M_PI);

    result.resonant_freq_rad = w0;
    result.resonant_freq_hz = f0;

    /* Quality factor for parallel RLC: Q = R/(ω₀L) = ω₀RC */
    double Q;
    if (R > 0.0) {
        Q = R / (w0 * L);
    } else {
        Q = INFINITY;  /* Ideal parallel LC (no R) */
    }

    result.quality_factor = Q;
    result.damping_factor = (Q > 0.0 && isfinite(Q)) ? 1.0 / (2.0 * Q) : 0.0;

    /* Bandwidth */
    double BW = (Q > 0.0 && isfinite(Q)) ? (w0 / Q) : 0.0;
    result.bandwidth_hz = BW / (2.0 * M_PI);

    /* Half-power frequencies */
    if (Q > 0.5) {
        double sqrt_term = sqrt(1.0 + 1.0 / (4.0 * Q * Q));
        double fL = f0 * (sqrt_term - 1.0 / (2.0 * Q));
        double fH = f0 * (sqrt_term + 1.0 / (2.0 * Q));
        result.half_power_low_hz = fL;
        result.half_power_high_hz = fH;
    }

    /* Impedance at resonance: Z(jω₀) = R (purely resistive, maximum) */
    result.impedance_at_res = R;
    result.peak_magnitude = 1.0;
    result.peak_magnitude_db = 0.0;

    return result;
}

/* ============================================================================
 * Transfer Functions for RLC Circuits
 * ============================================================================ */

tf_polynomial_t *resonance_series_tf(const rlc_params_t *params, int output)
{
    if (!params) return NULL;

    double R = params->R;
    double L = params->L;
    double C = params->C;

    /* Characteristic polynomial: s² + (R/L)s + 1/(LC) */
    double den[3] = { 1.0 / (L * C), R / L, 1.0 };

    /* LC = L·C, w0² = 1/LC */

    switch (output) {
        case 0: {  /* V_R/V_in: Bandpass */
            double num[2] = { 0.0, R / L };
            return tf_polynomial_create(num, 1, den, 2);
        }
        case 1: {  /* V_L/V_in: Highpass with resonance */
            double num[3] = { 0.0, 0.0, 1.0 };
            return tf_polynomial_create(num, 2, den, 2);
        }
        case 2: {  /* V_C/V_in: Lowpass */
            double num[1] = { 1.0 / (L * C) };
            return tf_polynomial_create(num, 0, den, 2);
        }
        case 3: {  /* V_L+V_C/V_in: Bandstop (notch) */
            double num[3] = { 1.0 / (L * C), 0.0, 1.0 };
            return tf_polynomial_create(num, 2, den, 2);
        }
        case 4: {  /* I/V_in = admittance: Bandpass */
            double num[2] = { 0.0, 1.0 / L };
            return tf_polynomial_create(num, 1, den, 2);
        }
        default:
            return NULL;
    }
}

tf_polynomial_t *resonance_parallel_tf(const rlc_params_t *params, int output)
{
    if (!params) return NULL;

    double R = params->R;
    double L = params->L;
    double C = params->C;

    /* Characteristic polynomial for admittance:
     * Y(s) = 1/R + 1/(sL) + sC = (s² + s/(RC) + 1/(LC)) / (s/C)
     * Impedance Z(s) = (s/C) / (s² + s/(RC) + 1/(LC))
     *              = s / (C·s² + s/R + C/(LC))
     *
     * Normalized: s² + s/(RC) + 1/(LC)
     */

    double den[3] = { 1.0 / (L * C), 1.0 / (R * C), 1.0 };

    switch (output) {
        case 0: {  /* V_out (impedance): Bandpass */
            double num[2] = { 0.0, 1.0 / C };
            return tf_polynomial_create(num, 1, den, 2);
        }
        case 1: {  /* I_R: same shape as V_out */
            double num[1] = { 1.0 / (R * C) };
            return tf_polynomial_create(num, 0, den, 2);
        }
        case 2: {  /* I_L: Lowpass */
            double num[1] = { 1.0 / (L * C) };
            return tf_polynomial_create(num, 0, den, 2);
        }
        case 3: {  /* I_C: Highpass */
            double num[3] = { 0.0, 0.0, 1.0 };
            return tf_polynomial_create(num, 2, den, 2);
        }
        default:
            return NULL;
    }
}

/* ============================================================================
 * Universal Resonance Curve
 * ============================================================================ */

/**
 * The universal resonance curve for a second-order bandpass:
 *
 * |H(Ω)| = 1/√(1 + Q²·(Ω - 1/Ω)²)
 *
 * where Ω = f/f₀ is the normalized frequency.
 *
 * At Ω = 1: |H| = 1 (peak)
 * At Ω = 0: |H| = 0
 * At Ω → ∞: |H| = 0
 *
 * Half-power: 1/√(1 + Q²·(Ω - 1/Ω)²) = 1/√2
 *            → Q·|Ω - 1/Ω| = 1
 *            → Ω = √(1 + 1/(4Q²)) ± 1/(2Q)
 *
 * For Q >> 1: Ω_L ≈ 1 - 1/(2Q), Ω_H ≈ 1 + 1/(2Q)
 */
double *resonance_universal_curve(const double *norm_freq,
                                    double Q, size_t n)
{
    if (!norm_freq || Q <= 0.0 || n == 0) return NULL;

    double *mag = (double *)malloc(n * sizeof(double));
    if (!mag) return NULL;

    for (size_t i = 0; i < n; i++) {
        double Omega = norm_freq[i];
        if (Omega <= 0.0) {
            mag[i] = 0.0;
            continue;
        }
        double detuning = Omega - 1.0 / Omega;
        double denom = 1.0 + Q * Q * detuning * detuning;
        mag[i] = 1.0 / sqrt(denom);
    }

    return mag;
}

/* ============================================================================
 * Quality Factor Analysis
 * ============================================================================ */

double resonance_q_from_bandwidth(double f0, double f_low, double f_high)
{
    if (f0 <= 0.0 || f_low <= 0.0 || f_high <= f_low) return 0.0;
    return f0 / (f_high - f_low);
}

double resonance_q_from_energy(const rlc_params_t *params,
                                 resonance_topology_t topology)
{
    if (!params) return 0.0;

    double R = params->R;
    double L = params->L;
    double C = params->C;
    double w0 = 1.0 / sqrt(L * C);

    if (R <= 0.0) return INFINITY;

    switch (topology) {
        case RESONANCE_SERIES_RLC:
            /* Q = ω₀L/R = 1/(ω₀RC) */
            return w0 * L / R;

        case RESONANCE_PARALLEL_RLC:
            /* Q = R/(ω₀L) = ω₀RC */
            return R / (w0 * L);

        default:
            return 0.0;
    }
}

double resonance_damping_from_q(double Q)
{
    if (Q <= 0.0) return INFINITY;
    return 1.0 / (2.0 * Q);
}

/* ============================================================================
 * Step Response
 * ============================================================================ */

/**
 * Step response of a resonant second-order system.
 *
 * Normalized transfer function for series RLC (V_C step response):
 *   H(s) = ω₀²/(s² + 2ζω₀s + ω₀²)
 *
 * For a unit step input:
 *
 * Underdamped (ζ < 1):
 *   v(t) = 1 - e^{-ζω₀t}·[cos(ω_d·t) + (ζ/√(1-ζ²))·sin(ω_d·t)]
 *   where ω_d = ω₀·√(1-ζ²)
 *
 * Critically damped (ζ = 1):
 *   v(t) = 1 - (1 + ω₀t)·e^{-ω₀t}
 *
 * Overdamped (ζ > 1):
 *   v(t) = 1 - (1/(τ₂-τ₁))·(τ₂·e^{-t/τ₁} - τ₁·e^{-t/τ₂})
 *   where τ₁,₂ = 1/(ω₀(ζ ± √(ζ²-1)))
 */
void resonance_step_response(const rlc_params_t *params,
                              const double *t, double *v_out,
                              size_t n, resonance_topology_t topology)
{
    if (!params || !t || !v_out || n == 0) return;

    double R = params->R;
    double L = params->L;
    double C = params->C;
    double w0 = 1.0 / sqrt(L * C);

    double zeta;
    if (topology == RESONANCE_SERIES_RLC) {
        zeta = R / (2.0 * w0 * L);
    } else {
        zeta = 1.0 / (2.0 * R * w0 * C);
    }

    if (R <= 0.0) {
        /* Lossless: sustained oscillation */
        for (size_t i = 0; i < n; i++) {
            v_out[i] = 1.0 - cos(w0 * t[i]);
        }
        return;
    }

    if (zeta < 1.0) {
        /* Underdamped */
        double wd = w0 * sqrt(1.0 - zeta * zeta);
        double alpha = zeta / sqrt(1.0 - zeta * zeta);
        for (size_t i = 0; i < n; i++) {
            double exp_term = exp(-zeta * w0 * t[i]);
            double cos_term = cos(wd * t[i]);
            double sin_term = sin(wd * t[i]);
            v_out[i] = 1.0 - exp_term * (cos_term + alpha * sin_term);
        }
    } else if (fabs(zeta - 1.0) < 1e-10) {
        /* Critically damped */
        for (size_t i = 0; i < n; i++) {
            double exp_term = exp(-w0 * t[i]);
            v_out[i] = 1.0 - (1.0 + w0 * t[i]) * exp_term;
        }
    } else {
        /* Overdamped */
        double sqrt_term = sqrt(zeta * zeta - 1.0);
        double tau1 = 1.0 / (w0 * (zeta - sqrt_term));
        double tau2 = 1.0 / (w0 * (zeta + sqrt_term));
        for (size_t i = 0; i < n; i++) {
            v_out[i] = 1.0 - (tau2 * exp(-t[i] / tau1)
                              - tau1 * exp(-t[i] / tau2))
                             / (tau2 - tau1);
        }
    }
}

/* ============================================================================
 * Coupled Resonators (L8 Advanced)
 * ============================================================================ */

/**
 * Two magnetically-coupled RLC resonators.
 *
 * Primary: L₁, C₁, R₁
 * Secondary: L₂, C₂, R₂
 * Mutual inductance: M
 * Coupling coefficient: k = M/√(L₁L₂)
 *
 * The frequency response is computed from the coupled differential
 * equations in the frequency domain:
 *
 * V₁ = (R₁ + jωL₁ + 1/(jωC₁))·I₁ + jωM·I₂
 * 0  = jωM·I₁ + (R₂ + jωL₂ + 1/(jωC₂))·I₂
 *
 * Solving for V₂/V₁ (voltage transfer to secondary):
 *   H(jω) = jωM·R₂ / (Z₁Z₂ + ω²M²)
 * where Z₁ = R₁ + jωL₁ + 1/(jωC₁), Z₂ = R₂ + jωL₂ + 1/(jωC₂)
 */
freq_response_t *resonance_coupled(const rlc_params_t *primary,
                                     const rlc_params_t *secondary,
                                     double M,
                                     const double *freq, size_t n_freq)
{
    if (!primary || !secondary || !freq || n_freq == 0) return NULL;

    freq_response_t *resp = freq_response_alloc(n_freq, freq[0],
                                                  freq[n_freq - 1],
                                                  FREQ_SCALE_LINEAR);
    if (!resp) return NULL;

    double R1 = primary->R, L1 = primary->L, C1 = primary->C;
    double R2 = secondary->R, L2 = secondary->L, C2 = secondary->C;

    for (size_t i = 0; i < n_freq; i++) {
        double f = freq[i];
        double w = 2.0 * M_PI * f;

        /* Impedances */
        double _Complex Z1 = R1 + I * (w * L1 - 1.0 / (w * C1));
        double _Complex Z2 = R2 + I * (w * L2 - 1.0 / (w * C2));
        double _Complex ZM = I * w * M;

        /* V₂/V₁ = ZM·R₂/(Z₁Z₂ + ZM²) */
        double _Complex denom = Z1 * Z2 + ZM * ZM;
        double _Complex H;
        if (cabs(denom) > 1e-15) {
            H = ZM * R2 / denom;
        } else {
            H = 0.0;
        }

        resp->points[i].frequency = f;
        resp->points[i].angular_freq = w;
        resp->points[i].real = creal(H);
        resp->points[i].imag = cimag(H);
        resp->points[i].magnitude = cabs(H);
        resp->points[i].phase_rad = carg(H);
        resp->points[i].magnitude_db = magnitude_to_db(cabs(H));
    }

    return resp;
}

/* ============================================================================
 * Crystal Resonator Model (L7 Application)
 * ============================================================================ */

/**
 * Quartz crystal equivalent circuit model.
 *
 * The crystal is modeled as a series RLC branch (motional arm:
 * L₁, C₁, R₁) in parallel with a capacitor C₀ (static/holder
 * capacitance).
 *
 * The impedance is:
 *   Z(s) = (s²L₁C₁ + sR₁C₁ + 1)
 *        / (s³L₁C₁C₀ + s²R₁C₁C₀ + s(C₁+C₀))
 *
 * Series resonance: f_s = 1/(2π√(L₁C₁))
 * Parallel resonance: f_p = f_s·√(1 + C₁/C₀)
 *
 * Quality factor: Q = ω_s·L₁/R₁ = 1/(ω_s·R₁·C₁)
 *
 * Typical values for a 10 MHz AT-cut crystal:
 *   L₁ ≈ 10 mH, C₁ ≈ 25 fF, R₁ ≈ 10 Ω, C₀ ≈ 5 pF
 *   Q ≈ 62,800, f_p - f_s ≈ 2500 Hz (250 ppm)
 */
tf_polynomial_t *resonance_crystal_model(double C0, double L1, double C1,
                                           double R1,
                                           double *fs, double *fp, double *Q)
{
    if (C0 <= 0.0 || L1 <= 0.0 || C1 <= 0.0) return NULL;

    double ws = 1.0 / sqrt(L1 * C1);
    double fs_val = ws / (2.0 * M_PI);
    double fp_val = fs_val * sqrt(1.0 + C1 / C0);
    double Q_val = ws * L1 / R1;

    if (fs) *fs = fs_val;
    if (fp) *fp = fp_val;
    if (Q) *Q = Q_val;

    /* Z(s) = N(s)/D(s) where:
     * N(s) = s²L₁C₁ + sR₁C₁ + 1
     * D(s) = s³L₁C₁C₀ + s²R₁C₁C₀ + s(C₁+C₀)
     */

    /* Numerator coefficients: b₀ + b₁s + b₂s² */
    double num[3] = { 1.0, R1 * C1, L1 * C1 };

    /* Denominator coefficients: a₀ + a₁s + a₂s² + a₃s³
     * Actually a₀ = 0 (no constant term in D(s)):
     * D(s) = s·[(C₁+C₀) + sR₁C₁C₀ + s²L₁C₁C₀] */
    double den[4] = { 0.0, C1 + C0, R1 * C1 * C0, L1 * C1 * C0 };

    return tf_polynomial_create(num, 2, den, 3);
}
