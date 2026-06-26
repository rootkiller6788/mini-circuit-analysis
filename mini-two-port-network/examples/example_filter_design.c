/**
 * @file example_filter_design.c
 * @brief Example: Butterworth & Chebyshev LC Ladder Filter Design
 *
 * Demonstrates complete LC ladder filter synthesis using two-port
 * ABCD cascade analysis. This is the standard method for designing
 * passive RF/microwave filters.
 *
 * Workflow:
 *   1. Choose filter type (Butterworth/Chebyshev) and order
 *   2. Look up normalized g-values (low-pass prototype)
 *   3. Build the ABCD cascade from alternating series/shunt elements
 *   4. Compute the frequency response (S21 in dB vs frequency)
 *   5. Apply frequency transformation to scale to desired cutoff
 *
 * This example covers:
 *   - Butterworth 3rd-order LP at 100 MHz
 *   - Chebyshev 3rd-order LP at 100 MHz with 0.5 dB ripple
 *   - Frequency response comparison
 *
 * Course: Georgia Tech ECE 6350 — Filter synthesis
 * Ref: Matthaei, Young, Jones, "Microwave Filters..."
 */

#include <stdio.h>
#include <math.h>
#include "../include/two_port.h"
#include "../include/abcd_params.h"
#include "../include/conversion.h"
#include "../include/s_params.h"
#include "../include/filter_design.h"

int main(void) {
    printf("========================================\n");
    printf("LC Ladder Filter Design & Analysis\n");
    printf("========================================\n\n");

    double fc = 100e6;              /* Cutoff frequency: 100 MHz */
    double omega_c = 2.0 * M_PI * fc;
    double z0 = 50.0;               /* System impedance */
    int order = 3;                  /* 3rd-order filter */

    /* ============================================================
     * Part 1: Butterworth Low-Pass Filter
     * ============================================================ */
    printf("--- Butterworth Low-Pass (N=%d) ---\n", order);

    double bw_g[3];
    filter_butterworth_g_values(order, bw_g);
    printf("Normalized g-values:\n");
    for (int i = 0; i < order; i++) {
        printf("  g%d = %.4f %s\n", i+1, bw_g[i],
               (i % 2 == 0) ? "H (series L)" : "F (shunt C)");
    }

    /* Compute actual component values at fc=100MHz, Z0=50Ω */
    printf("\nDenormalized component values at %.0f MHz:\n", fc / 1e6);
    for (int i = 0; i < order; i++) {
        if (i % 2 == 0) {
            /* Series L: L = g_i * Z0 / ωc */
            double l_val = bw_g[i] * z0 / omega_c;
            printf("  L%d = %.2f nH\n", i+1, l_val * 1e9);
        } else {
            /* Shunt C: C = g_i / (Z0 * ωc) */
            double c_val = bw_g[i] / (z0 * omega_c);
            printf("  C%d = %.2f pF\n", i+1, c_val * 1e12);
        }
    }

    /* Build and analyze the Butterworth filter */
    printf("\nFrequency Response (Butterworth):\n");
    printf("  Freq (MHz)    |S21| (dB)    Phase (°)\n");
    printf("  ----------    ----------    ---------\n");

    double test_freqs[] = {0.1, 10.0, 50.0, 100.0, 150.0, 200.0, 300.0, 500.0};
    for (int i = 0; i < 8; i++) {
        double f = test_freqs[i] * 1e6;
        double omega = 2.0 * M_PI * f;

        /* Scale: normalize to prototype (ω/ωc) */
        double omega_norm = omega / omega_c;

        /* Build ladder */
        matrix2x2_t abcd = filter_build_lc_ladder(bw_g, order, omega_norm, 0);
        matrix2x2_t s = convert_abcd_to_s(abcd, z0);
        double s21_db = sparams_s21_db(s);
        double phase = sparams_transmission_phase_deg(s);

        printf("  %8.1f      %7.2f      %8.1f\n", f / 1e6, s21_db, phase);
    }

    /* ============================================================
     * Part 2: Chebyshev Low-Pass Filter (0.5 dB ripple)
     * ============================================================ */
    double ripple_db = 0.5;
    printf("\n--- Chebyshev-I Low-Pass (N=%d, %.1f dB ripple) ---\n",
           order, ripple_db);

    double ch_g[3];
    filter_chebyshev_g_values(order, ripple_db, ch_g);
    printf("Normalized g-values:\n");
    for (int i = 0; i < order; i++) {
        printf("  g%d = %.4f %s\n", i+1, ch_g[i],
               (i % 2 == 0) ? "H" : "F");
    }

    printf("\nFrequency Response (Chebyshev-I):\n");
    printf("  Freq (MHz)    |S21| (dB)    Phase (°)\n");
    printf("  ----------    ----------    ---------\n");

    for (int i = 0; i < 8; i++) {
        double f = test_freqs[i] * 1e6;
        double omega = 2.0 * M_PI * f;
        double omega_norm = omega / omega_c;

        matrix2x2_t abcd = filter_build_lc_ladder(ch_g, order, omega_norm, 0);
        matrix2x2_t s = convert_abcd_to_s(abcd, z0);
        double s21_db = sparams_s21_db(s);
        double phase = sparams_transmission_phase_deg(s);

        printf("  %8.1f      %7.2f      %8.1f\n", f / 1e6, s21_db, phase);
    }

    /* ============================================================
     * Part 3: High-Pass Transformation Example
     * ============================================================ */
    printf("\n--- HP Transformation (100 MHz cutoff) ---\n");
    /* HP: s → ωc/s, so at prototype ω=1, the HP cutoff is at ωc */
    printf("Normalized prototype freq | HP response at %.0f MHz\n", fc/1e6);

    for (int i = 0; i < 5; i++) {
        double omega_norm = 0.2 * (i + 1);  /* 0.2, 0.4, 0.6, 0.8, 1.0 */
        double omega_hp_equiv = filter_lp_to_hp(omega_norm, omega_c);

        /* Compute LP prototype response at the equivalent frequency */
        complex_t h = filter_butterworth_lp(order, omega_norm);
        double mag_db = 20.0 * log10(complex_mag(h));

        printf("  ω/ωc = %.1f → f = %.1f MHz → |H| = %.1f dB\n",
               omega_norm, omega_hp_equiv / (2.0 * M_PI) / 1e6, mag_db);
    }

    /* ============================================================
     * Part 4: Bessel Filter Comparison (linear phase)
     * ============================================================ */
    printf("\n--- Bessel Filter Group Delay Flatness ---\n");
    printf("  Freq         |H_Bessel|    |H_Butter|    τg_Bessel\n");
    printf("  ----         ----------    ----------    ----------\n");

    for (int i = 0; i < 6; i++) {
        double omega = 0.2 * (i + 1);  /* normalized */
        complex_t h_bessel = filter_bessel_lp(3, omega);
        complex_t h_butter = filter_butterworth_lp(3, omega);
        /* Group delay approximated by phase/omega */
        double tau_bessel = -complex_arg(h_bessel) / (omega * omega_c);
        if (omega < 0.01) tau_bessel = 0;

        printf("  ω=%.1f        %.3f        %.3f        %.3f ns\n",
               omega, complex_mag(h_bessel), complex_mag(h_butter),
               tau_bessel * 1e9);
    }

    printf("\n========================================\n");
    printf("Filter Design Example Complete\n");
    printf("========================================\n");

    return 0;
}
