/**
 * @file example_rlc_resonance.c
 * @brief Example: Series and Parallel RLC Resonance Analysis
 *
 * This example demonstrates complete resonance analysis:
 * 1. Series RLC: compute resonance parameters, Q, bandwidth, transfer functions
 * 2. Parallel RLC: same analysis
 * 3. Universal resonance curve comparison for different Q values
 * 4. Step response characterization (overshoot, settling time)
 * 5. Bode plot of the bandpass transfer function
 *
 * L6 Canonical Problem: RLC resonance is the quintessential frequency
 * response problem in circuit analysis. It illustrates:
 * - Reactance cancellation at resonance
 * - Voltage/current magnification
 * - Half-power bandwidth and quality factor
 * - Damping factor and its effect on frequency and time response
 *
 * Reference: Hayt, Kemmerly & Durbin (2019), Ch. 15
 * Course: Berkeley EE16B, MIT 6.003
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "frequency_response.h"
#include "transfer_function.h"
#include "bode_plot.h"
#include "resonance.h"

/* Print resonance analysis results */
static void print_resonance(const char *label, const resonance_result_t *r)
{
    printf("\n--- %s ---\n", label);
    printf("  Resonant frequency:  f₀ = %.2f Hz  (ω₀ = %.2f rad/s)\n",
           r->resonant_freq_hz, r->resonant_freq_rad);
    printf("  Quality factor:      Q  = %.3f\n", r->quality_factor);
    printf("  Damping factor:      ζ  = %.4f\n", r->damping_factor);
    printf("  -3 dB Bandwidth:     BW = %.2f Hz\n", r->bandwidth_hz);

    if (r->quality_factor > 0.5) {
        printf("  Half-power freq:     f_L = %.2f Hz,  f_H = %.2f Hz\n",
               r->half_power_low_hz, r->half_power_high_hz);
    } else {
        printf("  Half-power freq:     N/A (overdamped, ζ > %.2f)\n",
               1.0 / sqrt(2.0));
    }

    printf("  Impedance at f₀:     |Z| = %.2f Ω\n", r->impedance_at_res);

    /* Interpret the response type */
    if (r->damping_factor < 1.0) {
        double overshoot_pct = 100.0 * exp(-M_PI * r->damping_factor
                                   / sqrt(1.0 - r->damping_factor
                                           * r->damping_factor));
        printf("  Response type:       UNDERDAMPED\n");
        printf("  Step overshoot:      %.1f%%\n", overshoot_pct);
    } else if (r->damping_factor < 1.0001) {
        printf("  Response type:       CRITICALLY DAMPED\n");
    } else {
        printf("  Response type:       OVERDAMPED\n");
    }
}

int main(void)
{
    printf("╔═══════════════════════════════════╗\n");
    printf("║  RLC Resonance Analysis Example  ║\n");
    printf("╚═══════════════════════════════════╝\n");

    /* ================================================================
     * Part 1: Series RLC Resonance
     * ================================================================ */
    printf("\n══════════ Series RLC Resonance ══════════\n");

    /* Series RLC with moderate Q:
     * R = 10 Ω, L = 10 mH, C = 100 nF
     * f₀ = 1/(2π√(10e-3 × 100e-9))
     *    = 1/(2π√(1e-9)) = 1/(2π × 3.162e-5) ≈ 5033 Hz
     * Q = ω₀L/R = 2π×5033×10e-3/10 ≈ 31.6
     * BW = f₀/Q ≈ 159 Hz */

    rlc_params_t series_params = { 10.0, 10e-3, 100e-9, 1.0 };
    resonance_result_t series_res = resonance_series(&series_params);
    print_resonance("Series RLC (R=10Ω, L=10mH, C=100nF)", &series_res);

    /* Verify Q from bandwidth formula */
    double q_bw = resonance_q_from_bandwidth(series_res.resonant_freq_hz,
                                               series_res.half_power_low_hz,
                                               series_res.half_power_high_hz);
    printf("  Q from bandwidth:    %.3f (should match above)\n", q_bw);

    /* Q from energy definition */
    double q_energy = resonance_q_from_energy(&series_params,
                                                RESONANCE_SERIES_RLC);
    printf("  Q from energy:       %.3f\n", q_energy);

    /* ================================================================
     * Part 2: Effect of damping on resonance
     * ================================================================ */
    printf("\n═══════ Effect of Damping (Series RLC) ═══════\n");

    double R_values[] = { 1.0, 5.0, 10.0, 31.62, 100.0, 500.0 };
    int n_R = 6;

    printf("\n  %8s  %8s  %8s  %8s  %8s  %12s\n",
           "R (Ω)", "f₀ (Hz)", "Q", "ζ", "BW (Hz)", "Peak |Z| (Ω)");
    printf("  %8s  %8s  %8s  %8s  %8s  %12s\n",
           "------", "------", "------", "------", "------", "----------");

    for (int i = 0; i < n_R; i++) {
        rlc_params_t p = { R_values[i], 10e-3, 100e-9, 1.0 };
        resonance_result_t r = resonance_series(&p);
        printf("  %8.2f  %8.1f  %8.2f  %8.4f  %8.1f  %12.2f\n",
               R_values[i], r.resonant_freq_hz, r.quality_factor,
               r.damping_factor, r.bandwidth_hz, r.impedance_at_res);
    }

    /* ================================================================
     * Part 3: Series RLC Transfer Functions
     * ================================================================ */
    printf("\n═══════ Series RLC Transfer Functions ═══════\n");

    /* V_R/V_in: Bandpass */
    tf_polynomial_t *tf_vr = resonance_series_tf(&series_params, 0);
    /* V_C/V_in: Lowpass */
    tf_polynomial_t *tf_vc = resonance_series_tf(&series_params, 2);
    /* V_L/V_in: Highpass */
    tf_polynomial_t *tf_vl = resonance_series_tf(&series_params, 1);

    if (tf_vr && tf_vc && tf_vl) {
        double f0 = series_res.resonant_freq_hz;

        printf("\n  Magnitude at key frequencies:\n");
        printf("  %12s  %12s  %12s  %12s\n",
               "Frequency", "|V_R/V_in|", "|V_C/V_in|", "|V_L/V_in|");
        double test_f[] = { f0 * 0.1, f0 * 0.5, f0, f0 * 2.0, f0 * 10.0 };
        for (int i = 0; i < 5; i++) {
            double mag_vr = tf_magnitude_at(tf_vr, test_f[i]);
            double mag_vc = tf_magnitude_at(tf_vc, test_f[i]);
            double mag_vl = tf_magnitude_at(tf_vl, test_f[i]);
            printf("  %12.1f  %12.4f  %12.4f  %12.4f\n",
                   test_f[i], mag_vr, mag_vc, mag_vl);
        }

        /* Generate Bode plot of bandpass */
        bode_plot_t *bode = bode_compute(tf_vr, f0 * 0.01, f0 * 100.0, 50);
        if (bode) {
            double dom = bode_dominant_pole_freq(bode);
            printf("\n  Bandpass Bode analysis:\n");
            printf("  Center frequency: %.1f Hz\n", f0);
            printf("  Dominant pole:    %.1f Hz\n", dom);

            bode_free(bode);
        }
    }

    /* ================================================================
     * Part 4: Parallel RLC Resonance
     * ================================================================ */
    printf("\n═════════ Parallel RLC Resonance ═════════\n");

    /* Parallel RLC: R = 1 kΩ, L = 10 mH, C = 100 nF
     * f₀ same ≈ 5033 Hz
     * Q = R/(ω₀L) = 1000/(2π×5033×10e-3) ≈ 3.16
     * BW = f₀/Q ≈ 1592 Hz */

    rlc_params_t parallel_params = { 1000.0, 10e-3, 100e-9, 1.0 };
    resonance_result_t parallel_res = resonance_parallel(&parallel_params);
    print_resonance("Parallel RLC (R=1kΩ, L=10mH, C=100nF)", &parallel_res);

    /* ================================================================
     * Part 5: Universal Resonance Curve
     * ================================================================ */
    printf("\n════════ Universal Resonance Curve ════════\n");

    double Q_values[] = { 1.0, 5.0, 10.0, 50.0 };
    int n_Q = 4;

    printf("\n  Normalized magnitude |H(Ω)|/|H_max| at Ω = f/f₀:\n");
    printf("  %8s", "Ω");
    for (int q = 0; q < n_Q; q++)
        printf("  Q=%-8.1f", Q_values[q]);
    printf("\n");

    /* Sample a few key frequencies */
    double sample_Omega[] = { 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2 };
    for (int s = 0; s < 7; s++) {
        double Omega = sample_Omega[s];
        printf("  %8.3f", Omega);
        for (int q = 0; q < n_Q; q++) {
            double *curve = resonance_universal_curve(&Omega, Q_values[q], 1);
            printf("  %8.4f", curve[0]);
            free(curve);
        }
        printf("\n");
    }

    /* ================================================================
     * Part 6: Step Response
     * ================================================================ */
    printf("\n═══════════ Step Response Analysis ═══════════\n");

    rlc_params_t step_params = { 10.0, 10e-3, 100e-9, 1.0 };
    resonance_result_t step_res = resonance_series(&step_params);

    int n_t = 200;
    double t_max = 10.0 / (step_res.damping_factor * step_res.resonant_freq_rad);
    double *t = (double *)malloc(n_t * sizeof(double));
    double *v = (double *)malloc(n_t * sizeof(double));

    for (int i = 0; i < n_t; i++) {
        t[i] = (double)i * t_max / (double)(n_t - 1);
    }

    resonance_step_response(&step_params, t, v, n_t,
                              RESONANCE_SERIES_RLC);

    /* Find overshoot */
    double max_v = 0.0;
    double t_peak = 0.0;
    for (int i = 0; i < n_t; i++) {
        if (v[i] > max_v) {
            max_v = v[i];
            t_peak = t[i];
        }
    }

    double overshoot = (max_v - 1.0) * 100.0;
    printf("  Damping: ζ = %.3f\n", step_res.damping_factor);
    printf("  Peak overshoot: %.1f%% at t = %.3f ms\n",
           overshoot, t_peak * 1000.0);
    printf("  Theoretical overshoot: %.1f%%\n",
           100.0 * exp(-M_PI * step_res.damping_factor
                       / sqrt(1.0 - step_res.damping_factor
                               * step_res.damping_factor)));

    /* Find settling time (within 2% of final value) */
    for (int i = n_t - 1; i >= 0; i--) {
        if (fabs(v[i] - 1.0) > 0.02) {
            printf("  2%% settling time: %.3f ms\n", t[i + 1] * 1000.0);
            break;
        }
    }

    free(t);
    free(v);

    /* Cleanup */
    tf_polynomial_free(tf_vr);
    tf_polynomial_free(tf_vc);
    tf_polynomial_free(tf_vl);

    printf("\nResonance analysis complete.\n");
    return 0;
}
