/**
 * @file example_harmonic_analysis.c
 * @brief Harmonic power analysis of a non-linear load
 *
 * Demonstrates complete harmonic power analysis per IEEE 1459-2010
 * for a single-phase rectifier load with significant harmonic content.
 *
 * This example (L6) covers: DFT-based harmonic decomposition,
 * THD computation, displacement vs. distortion PF, K-factor,
 * and IEEE 519 compliance assessment.
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/power_factor.h"
#include "../include/harmonic_power.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Generate a rectified current waveform (typical of SMPS with PFC) */
static double rectifier_current(double t, double f, double i_fund_peak)
{
    /* Simplified full-wave rectifier current:
     * Fundamental + 3rd (30%) + 5th (15%) + 7th (8%) + 9th (5%) */
    double omega = 2.0 * M_PI * f;
    double i = i_fund_peak * (
        sin(omega * t)
        + 0.30 * sin(3.0 * omega * t + M_PI)
        + 0.15 * sin(5.0 * omega * t)
        + 0.08 * sin(7.0 * omega * t + M_PI)
        + 0.05 * sin(9.0 * omega * t)
    );
    return i;
}

int main(void)
{
    printf("========================================\n");
    printf("  Harmonic Power Analysis — Rectifier Load\n");
    printf("  IEEE 1459-2010 / IEEE 519-2014\n");
    printf("========================================\n\n");

    /* Configuration */
    double f_fund   = 60.0;
    double f_sample = 6000.0;  /* 100 samples per cycle */
    size_t n_cycles = 16;
    size_t N = (size_t)(n_cycles * f_sample / f_fund); /* 1600 samples */
    double dt = 1.0 / f_sample;

    /* Generate waveforms */
    double *v_waveform = (double*)malloc(N * sizeof(double));
    double *i_waveform = (double*)malloc(N * sizeof(double));
    if (!v_waveform || !i_waveform) {
        printf("Memory allocation failed\n");
        return 1;
    }

    double v_peak = 170.0;  /* 120V RMS */
    double i_fund_peak = 10.0;

    for (size_t n = 0; n < N; n++) {
        double t = n * dt;
        v_waveform[n] = v_peak * sin(2.0 * M_PI * f_fund * t);
        i_waveform[n] = rectifier_current(t, f_fund, i_fund_peak);
    }

    /* ==================================================================
     * 1. Time-domain PF computation
     * ================================================================== */
    pf_single_phase_t time_power;
    pf_compute_from_samples(v_waveform, i_waveform, N, dt, &time_power);

    printf("Time-Domain Analysis:\n");
    printf("  V_rms = %.1f V\n", time_power.v_rms);
    printf("  I_rms = %.2f A\n", time_power.i_rms);
    printf("  P     = %.1f W\n", time_power.p_real);
    printf("  S     = %.1f VA\n", time_power.s_apparent);
    printf("  PF    = %.4f\n\n", time_power.pf);

    /* ==================================================================
     * 2. Frequency-domain harmonic spectrum analysis
     * ================================================================== */
    hp_spectrum_t spectrum;
    hp_analyze_spectrum(v_waveform, i_waveform, N,
                         f_sample, f_fund, &spectrum);

    printf("Harmonic Spectrum (first 10 harmonics):\n");
    printf("  Order | V_mag [V] | V_phase  | I_mag [A] | I_phase  | P_h [W]\n");
    printf("  ------|-----------|----------|-----------|----------|--------\n");

    uint32_t n_show = spectrum.num_harmonics;
    if (n_show > 10) n_show = 10;

    for (uint32_t k = 0; k < n_show; k++) {
        printf("  %5u | %9.3f | %8.1f | %9.3f | %8.1f | %7.1f\n",
               spectrum.harmonics[k].order,
               spectrum.harmonics[k].v_mag,
               spectrum.harmonics[k].v_phase_deg,
               spectrum.harmonics[k].i_mag,
               spectrum.harmonics[k].i_phase_deg,
               spectrum.harmonics[k].p_h);
    }

    /* ==================================================================
     * 3. THD and Distortion Analysis
     * ================================================================== */
    printf("\nDistortion Analysis (IEEE 1459-2010):\n");
    printf("  THDv = %.2f %%\n", spectrum.thd_v_percent);
    printf("  THDi = %.2f %%\n", spectrum.thd_i_percent);
    printf("  Displacement PF (cos φ₁) = %.4f\n", spectrum.pf_displacement);
    printf("  Distortion PF           = %.4f\n", spectrum.pf_distortion);
    printf("  True PF  = displacement × distortion = %.4f\n",
           spectrum.pf_true);
    printf("  PF reduction due to harmonics = %.2f %%\n\n",
           (1.0 - spectrum.pf_distortion) * 100.0);

    /* ==================================================================
     * 4. K-Factor Analysis (Transformer Derating)
     * ================================================================== */
    double *i_mags = (double*)malloc(spectrum.num_harmonics * sizeof(double));
    for (uint32_t k = 0; k < spectrum.num_harmonics; k++) {
        i_mags[k] = spectrum.harmonics[k].i_mag;
    }

    double k = hp_k_factor(i_mags, spectrum.num_harmonics);
    double derating = hp_transformer_derating(k, 0.08);

    printf("Transformer K-Factor Analysis (IEEE C57.110):\n");
    printf("  K-factor    = %.2f\n", k);
    printf("  Derating    = %.2f (can use %.0f%% of rated kVA)\n",
           derating, derating * 100.0);
    printf("  Recommended = K-%.0f transformer\n\n", ceil(k));

    /* ==================================================================
     * 5. IEEE 519 Compliance Check
     * ================================================================== */
    hp_limits_t limits;
    hp_init_ieee519_limits(&limits, HP_STD_IEEE519, 5000.0, 100.0, 480.0);

    int violations;
    int compliant = hp_check_ieee519(&spectrum, &limits, &violations);

    printf("IEEE 519-2014 Compliance:\n");
    printf("  TDD limit = %.1f %%\n", limits.tdd_limit_percent);
    printf("  TDD actual = %.1f %%\n", spectrum.tdd_i_percent);
    if (compliant == 0) {
        printf("  Status: COMPLIANT ✅\n");
    } else {
        printf("  Status: NON-COMPLIANT ❌ (%d violations)\n", violations);
    }

    /* ==================================================================
     * 6. Parseval Verification
     * ================================================================== */
    double *p_harmonics = (double*)malloc(spectrum.num_harmonics * sizeof(double));
    for (uint32_t k = 0; k < spectrum.num_harmonics; k++) {
        p_harmonics[k] = spectrum.harmonics[k].p_h;
    }

    int parseval_ok = hp_verify_parseval_power(time_power.p_real,
                                                p_harmonics,
                                                spectrum.num_harmonics, 1.0);
    printf("\nParseval's Theorem Verification:\n");
    printf("  Time-domain P = %.1f W\n", time_power.p_real);
    printf("  Freq-domain ΣP_h = %.1f W\n",
           spectrum.harmonics[0].p_h + spectrum.harmonics[2].p_h
           + spectrum.harmonics[4].p_h);
    printf("  Parseval: %s\n\n", (parseval_ok == 0) ? "HOLDS ✓" : "VIOLATED ✗");

    /* ==================================================================
     * 7. Mitigation Recommendations
     * ================================================================== */
    char recommendation[1024];
    hp_recommend_mitigation(&spectrum, recommendation, sizeof(recommendation));
    printf("%s\n", recommendation);

    printf("========================================\n");

    free(v_waveform);
    free(i_waveform);
    free(i_mags);
    free(p_harmonics);
    return 0;
}
