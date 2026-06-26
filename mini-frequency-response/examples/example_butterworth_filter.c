/**
 * @file example_butterworth_filter.c
 * @brief Example: Complete Butterworth lowpass filter design and analysis
 *
 * This example demonstrates:
 * 1. Designing a 5th-order Butterworth lowpass filter (f_c = 1 kHz)
 * 2. Computing the transfer function
 * 3. Generating Bode plots (exact and asymptotic)
 * 4. Finding poles and zeros
 * 5. Computing g-values for LC ladder realization
 * 6. Denormalizing to actual component values
 *
 * L6 Canonical Problem: End-to-end filter design from specification
 * to component values, with frequency response verification.
 *
 * Reference: Zverev (1967), Sedra & Smith (2020) Ch. 16
 * Course: Berkeley EE105, Stanford EE247
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "frequency_response.h"
#include "transfer_function.h"
#include "bode_plot.h"
#include "filter_design.h"
#include "resonance.h"
#include "stability.h"

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Butterworth Lowpass Filter Design Example  ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* Step 1: Design specification */
    filter_spec_t spec;
    spec.approx = APPROX_BUTTERWORTH;
    spec.f_pass = 1000.0;           /* Passband edge: 1 kHz */
    spec.f_stop = 5000.0;           /* Stopband edge: 5 kHz */
    spec.a_pass = 3.0;              /* Max passband ripple: 3 dB */
    spec.a_stop = 40.0;             /* Min stopband attenuation: 40 dB */
    spec.dc_gain = 1.0;
    spec.source_impedance = 50.0;
    spec.load_impedance = 50.0;

    printf("Filter Specification:\n");
    printf("  Type:      Butterworth (maximally flat)\n");
    printf("  Passband:  0 - %.0f Hz (≤ %.1f dB ripple)\n",
           spec.f_pass, spec.a_pass);
    printf("  Stopband:  ≥ %.0f Hz (≥ %.1f dB attenuation)\n",
           spec.f_stop, spec.a_stop);
    printf("  Impedance: %.0f Ω\n\n", spec.source_impedance);

    /* Step 2: Compute required order */
    int order = filter_order_butterworth(&spec);
    printf("Computed filter order: N = %d\n", order);
    printf("  (Required for %.0f dB attenuation at %.0f Hz)\n\n",
           spec.a_stop, spec.f_stop);

    /* Step 3: Generate prototype (normalized ω_c = 1 rad/s) */
    tf_polynomial_t *proto = filter_butterworth_prototype(order);
    if (!proto) {
        printf("Error: Failed to generate prototype!\n");
        return 1;
    }

    printf("Butterworth Polynomial B_%d(s) coefficients:\n", order);
    printf("  Denominator (ascending powers):\n");
    for (size_t i = 0; i <= proto->den_order; i++) {
        printf("    a_%zu = %.6f\n", i, proto->den[i]);
    }

    /* Step 4: Frequency-transform to 1 kHz cutoff */
    double wc = 2.0 * M_PI * spec.f_pass;
    tf_polynomial_t *lp = filter_lp_to_lp(proto, wc);
    if (!lp) {
        printf("Error: Frequency transformation failed!\n");
        tf_polynomial_free(proto);
        return 1;
    }

    /* Step 5: Pole-zero analysis */
    tf_pole_zero_t *pz = tf_to_pole_zero(lp);
    if (pz) {
        printf("\nPoles of the %.0f Hz Butterworth filter:\n", spec.f_pass);
        for (size_t i = 0; i < pz->num_poles; i++) {
            double re = creal(pz->poles[i]);
            double im = cimag(pz->poles[i]);
            double freq = cabs(pz->poles[i]) / (2.0 * M_PI);
            printf("  p_%zu = %.1f + j(%.1f)  →  |p| = %.1f Hz\n",
                   i + 1, re, im, freq);
        }

        int is_mp = tf_is_minimum_phase(pz);
        printf("  Minimum-phase: %s\n", is_mp ? "YES" : "NO");
    }

    /* Step 6: Stability check */
    int stable = stability_is_hurwitz(lp);
    printf("\nStability (Routh-Hurwitz): %s\n",
           stable == 1 ? "STABLE (all poles in LHP)" :
           stable == 0 ? "UNSTABLE" : "ERROR");

    /* Step 7: Bode plot computation */
    printf("\nComputing Bode plot...\n");
    bode_plot_t *bode = bode_compute(lp, 10.0, 100000.0, 50);
    if (bode) {
        printf("  Frequency points: %zu\n",
               bode->magnitude_response.num_points);
        printf("  DC gain: %.2f dB\n", bode->dc_gain_db);
        printf("  HF slope: %.1f dB/dec\n", bode->hf_slope_db_dec);

        /* Dominant pole */
        double dom_pole = bode_dominant_pole_freq(bode);
        printf("  -3 dB cutoff: %.2f Hz\n", dom_pole);

        /* Gain-bandwidth product */
        double gbwp = bode_gain_bandwidth_product(bode);
        printf("  GBWP: %.2f Hz\n", gbwp);

        /* Check for peaking (Butterworth has none) */
        double peak_f, peak_db;
        int has_peak = bode_detect_peaking(bode, &peak_f, &peak_db);
        printf("  Resonance peaking: %s\n",
               has_peak ? "DETECTED" : "NONE (maximally flat)");
    }

    /* Step 8: Asymptotic Bode */
    bode_plot_t *bode_asym = bode_asymptotic(pz, 10.0, 100000.0, 50);
    if (bode_asym) {
        printf("\nAsymptotic Bode approximation:\n");
        printf("  DC gain: %.2f dB\n", bode_asym->dc_gain_db);
        printf("  HF slope: %.1f dB/dec (expected: %.0f dB/dec for N=%d)\n",
               bode_asym->hf_slope_db_dec, -20.0 * order, order);
    }

    /* Step 9: LC ladder g-values */
    size_t n_elem = 0;
    double *g = filter_g_values_butterworth(order, &n_elem);
    if (g) {
        printf("\nNormalized g-values (LC ladder prototype, R₀=1Ω, ω_c=1):\n");
        printf("  g₀ = %.6f (source resistance)\n", g[0]);
        for (size_t i = 1; i < n_elem - 1; i++) {
            printf("  g_%zu = %.6f (%s)\n", i, g[i],
                   i % 2 == 1 ? "series L / shunt C" : "shunt C / series L");
        }
        printf("  g_%zu = %.6f (load resistance)\n", n_elem - 1, g[n_elem - 1]);

        /* Denormalize to actual values */
        double *L = NULL, *C = NULL;
        double R_load = 0.0;
        filter_denormalize(g, n_elem, spec.source_impedance, wc,
                            &L, &C, &R_load);

        printf("\nDenormalized component values (R₀=%.0fΩ, f_c=%.0fHz):\n",
               spec.source_impedance, spec.f_pass);
        for (size_t i = 0; i < (n_elem - 2) / 2; i++) {
            if (L && L[i] > 0.0)
                printf("  L_%zu = %.6e H = %.3f µH\n",
                       i + 1, L[i], L[i] * 1e6);
        }
        for (size_t i = 0; i < (n_elem - 2) / 2; i++) {
            if (C && C[i] > 0.0)
                printf("  C_%zu = %.6e F = %.3f nF\n",
                       i + 1, C[i], C[i] * 1e9);
        }
        printf("  R_load = %.2f Ω\n", R_load);

        free(L);
        free(C);
        free(g);
    }

    /* Step 10: Magnitude response at key frequencies */
    printf("\nMagnitude response at key frequencies:\n");
    double test_freqs[] = { 100.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0 };
    for (int i = 0; i < 6; i++) {
        double mag = tf_magnitude_at(lp, test_freqs[i]);
        double mag_db = magnitude_to_db(mag);
        printf("  @ %8.1f Hz: |H| = %.6f  = %+.2f dB\n",
               test_freqs[i], mag, mag_db);
    }

    /* Cleanup */
    bode_free(bode);
    bode_free(bode_asym);
    tf_pole_zero_free(pz);
    tf_polynomial_free(lp);
    tf_polynomial_free(proto);

    printf("\nFilter design complete.\n");
    return 0;
}
