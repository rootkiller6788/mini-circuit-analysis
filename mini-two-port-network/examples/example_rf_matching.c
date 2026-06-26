/**
 * @file example_rf_matching.c
 * @brief Example: RF Impedance Matching Network Design
 *
 * Demonstrates the complete RF amplifier matching workflow using
 * S-parameters and Smith chart concepts (implemented algebraically).
 *
 * This example models a typical LNA (Low Noise Amplifier) design
 * at 2.4 GHz (WiFi/Bluetooth ISM band), covering:
 *   1. Transistor S-parameters at operating frequency
 *   2. Stability analysis (Rollett K, μ-factor, stability circles)
 *   3. Stabilization if needed
 *   4. Conjugate impedance matching for maximum gain
 *   5. L-network and π-network matching circuit design
 *   6. Noise figure trade-off analysis
 *
 * This is the standard workflow taught in Stanford EE359 and
 * Georgia Tech ECE 6350 for microwave amplifier design.
 *
 * Reference: Gonzalez, "Microwave Transistor Amplifiers", Ch. 3
 * Course: Stanford EE359 — RF Amplifier Design
 */

#include <stdio.h>
#include <math.h>
#include "../include/two_port.h"
#include "../include/s_params.h"
#include "../include/conversion.h"
#include "../include/stability.h"
#include "../include/abcd_params.h"
#include "../include/network_synthesis.h"

int main(void) {
    printf("========================================\n");
    printf("2.4 GHz LNA Matching Network Design\n");
    printf("========================================\n\n");

    double freq = 2.4e9;           /* Operating frequency: 2.4 GHz */
    double omega = 2.0 * M_PI * freq;
    double z0 = 50.0;              /* System impedance */

    /* ============================================================
     * Step 1: Define the transistor S-parameters (typical GaAs pHEMT)
     * ============================================================ */
    /* ATF-54143-like pHEMT at 2.4 GHz, Vds=3V, Ids=60mA */
    matrix2x2_t s_transistor = sparams_create(
        complex_make(0.65, -0.35),   /* s11: mag=0.74, ang=-28° */
        complex_make(0.05, 0.02),    /* s12: mag=0.054, ang=22° (isolation) */
        complex_make(3.2, -1.1),     /* s21: mag=3.38, ang=-19° (gain) */
        complex_make(0.45, -0.25)    /* s22: mag=0.51, ang=-29° */
    );

    printf("Transistor S-Parameters at 2.4 GHz:\n");
    printf("  s11 = %.3f ∠ %.1f°  (input reflection)\n",
           complex_mag(s_transistor.m11),
           complex_arg(s_transistor.m11) * 180.0 / M_PI);
    printf("  s12 = %.3f ∠ %.1f°  (reverse isolation)\n",
           complex_mag(s_transistor.m12),
           complex_arg(s_transistor.m12) * 180.0 / M_PI);
    printf("  s21 = %.3f ∠ %.1f°  (forward gain)\n",
           complex_mag(s_transistor.m21),
           complex_arg(s_transistor.m21) * 180.0 / M_PI);
    printf("  s22 = %.3f ∠ %.1f°  (output reflection)\n\n",
           complex_mag(s_transistor.m22),
           complex_arg(s_transistor.m22) * 180.0 / M_PI);

    /* Basic S-parameter metrics */
    printf("Performance Metrics:\n");
    printf("  Forward Gain = %.1f dB\n", sparams_s21_db(s_transistor));
    printf("  Input Return Loss = %.1f dB\n",
           sparams_input_return_loss_db(s_transistor));
    printf("  Output Return Loss = %.1f dB\n",
           sparams_output_return_loss_db(s_transistor));
    printf("  Reverse Isolation = %.1f dB\n",
           sparams_isolation_db(s_transistor));
    printf("  MSG = %.1f dB\n\n",
           10.0 * log10(sparams_msg(s_transistor)));

    /* ============================================================
     * Step 2: Stability Analysis
     * ============================================================ */
    double k = sparams_rollett_k(s_transistor);
    double delta_mag = sparams_delta_mag(s_transistor);
    double mu = sparams_mu_factor(s_transistor);

    printf("Stability Analysis:\n");
    printf("  Rollett K = %.3f\n", k);
    printf("  |Δ| = %.4f\n", delta_mag);
    printf("  μ-factor = %.3f\n\n", mu);

    stability_t stab = stability_classify(s_transistor, PARAM_S, z0);
    switch (stab) {
        case STABLE_UNCONDITIONAL:
            printf("  Status: UNCONDITIONALLY STABLE (K>1, |Δ|<1)\n\n");
            break;
        case STABLE_CONDITIONAL:
            printf("  Status: CONDITIONALLY STABLE\n");
            /* Show stability circles */
            {
                complex_t cs, cl;
                double rs, rl;
                stability_circles(s_transistor, &cs, &rs, &cl, &rl);
                printf("  Source stability circle: center=%.2f∠%.1f°, radius=%.3f\n",
                       complex_mag(cs), complex_arg(cs)*180/M_PI, rs);
                printf("  Load stability circle:  center=%.2f∠%.1f°, radius=%.3f\n",
                       complex_mag(cl), complex_arg(cl)*180/M_PI, rl);
            }
            /* Find stabilizing resistance */
            {
                double r_stab = stability_stabilizing_series_r(s_transistor, z0);
                printf("  Recommended series R for stabilization: %.1f Ω\n\n", r_stab);
            }
            break;
        default:
            printf("  Status: UNSTABLE — requires redesign\n\n");
            break;
    }

    /* ============================================================
     * Step 3: Conjugate Matching for Maximum Gain
     * ============================================================ */
    complex_t gms, gml;
    int match_ok = sparams_conjugate_match(s_transistor, &gms, &gml);

    if (match_ok == 0) {
        printf("Simultaneous Conjugate Match:\n");
        printf("  ΓMS = %.3f ∠ %.1f°\n",
               complex_mag(gms), complex_arg(gms) * 180.0 / M_PI);
        printf("       → ZS_opt = %.1f + j%.1f Ω\n",
               impedance_from_gamma(gms, z0).real,
               impedance_from_gamma(gms, z0).imag);
        printf("  ΓML = %.3f ∠ %.1f°\n",
               complex_mag(gml), complex_arg(gml) * 180.0 / M_PI);
        printf("       → ZL_opt = %.1f + j%.1f Ω\n\n",
               impedance_from_gamma(gml, z0).real,
               impedance_from_gamma(gml, z0).imag);

        /* Maximum Available Gain */
        double mag = sparams_mag(s_transistor);
        printf("  MAG = %.1f dB (maximum available gain with perfect matching)\n\n",
               10.0 * log10(mag));
    } else {
        printf("  Device is not unconditionally stable — conjugate match undefined.\n");
        printf("  Use MSG = %.1f dB as gain estimate.\n\n",
               10.0 * log10(sparams_msg(s_transistor)));
    }

    /* ============================================================
     * Step 4: Noise Figure Analysis
     * ============================================================ */
    double nfmin = 1.6;          /* NFmin (linear) ≈ 2.0 dB */
    double rn = 0.2;             /* Normalized noise resistance Rn/Z0 */
    complex_t gopt = complex_make(0.45, 0.15);  /* Γopt for minimum noise */

    /* Compute noise figure for the conjugate gain match source */
    if (match_ok == 0) {
        double nf_at_gain_match = sparams_noise_figure(nfmin, rn, gopt, gms);
        printf("Noise Analysis:\n");
        printf("  NFmin = %.2f dB\n", 10.0 * log10(nfmin));
        printf("  Γopt  = %.3f ∠ %.1f°\n",
               complex_mag(gopt), complex_arg(gopt) * 180.0 / M_PI);
        printf("  NF at gain match = %.2f dB\n", 10.0 * log10(nf_at_gain_match));
        printf("  (Trade-off: gain match vs noise match)\n\n");
    }

    /* ============================================================
     * Step 5: L-Network Matching Circuit Design
     * ============================================================ */
    printf("Input Matching Network (L-section):\n");

    complex_t z_antenna = complex_make(50.0, 0.0);  /* Antenna = 50 Ω */
    complex_t z_amp_input;
    if (match_ok == 0) {
        z_amp_input = impedance_from_gamma(gms, z0);
    } else {
        /* Use 50Ω as approximation */
        z_amp_input = complex_make(50.0, 0.0);
    }

    double l_match, c_match;
    int l_ret = synthesize_l_match(z_antenna, z_amp_input, omega, &l_match, &c_match);

    if (l_ret == 0) {
        if (l_match > 0) printf("  Series L = %.2f nH\n", l_match * 1e9);
        if (c_match > 0) printf("  Shunt C  = %.2f pF\n", c_match * 1e12);
    }

    /* ============================================================
     * Step 6: π-Network Matching (broader bandwidth)
     * ============================================================ */
    printf("\nOutput Matching Network (π-section, Q=3):\n");

    complex_t z_load = complex_make(50.0, 0.0);
    complex_t z_amp_output;
    if (match_ok == 0) {
        z_amp_output = impedance_from_gamma(gml, z0);
    } else {
        z_amp_output = complex_make(50.0, 0.0);
    }

    complex_t ya_pi, yb_pi, yc_pi;
    int pi_ret = synthesize_pi_match(z_amp_output, z_load, omega,
                                      3.0, &ya_pi, &yb_pi, &yc_pi);
    if (pi_ret == 0) {
        printf("  Input shunt:  C = %.2f pF\n", ya_pi.imag / omega * 1e12);
        printf("  Series:       L = %.2f nH\n", -1.0 / (yb_pi.imag * omega) * 1e9);
        printf("  Output shunt: C = %.2f pF\n", yc_pi.imag / omega * 1e12);
    }

    /* ============================================================
     * Step 7: Complete Chain Simulation
     * ============================================================ */
    printf("\nComplete RF Chain Simulation (Input Match → LNA → Output Match):\n");

    /* Build the complete chain in ABCD */
    matrix2x2_t abcd_input_match = abcd_l_network(
        complex_make(0.0, omega * l_match),  /* Series L */
        complex_make(0.0, omega * c_match)    /* Shunt C */
    );

    /* LNA ABCD from S-parameters */
    matrix2x2_t abcd_lna = convert_s_to_abcd(s_transistor, z0);

    /* Output π-network ABCD */
    matrix2x2_t abcd_output_match = abcd_pi_network(ya_pi, yb_pi, yc_pi);

    /* Cascade: Input Match → LNA → Output Match */
    matrix2x2_t abcd_chain = matrix2x2_mul(
        abcd_input_match,
        matrix2x2_mul(abcd_lna, abcd_output_match)
    );

    matrix2x2_t s_chain = convert_abcd_to_s(abcd_chain, z0);
    printf("  System Gain = %.1f dB\n", sparams_s21_db(s_chain));
    printf("  Input RL    = %.1f dB\n", sparams_input_return_loss_db(s_chain));
    printf("  Output RL   = %.1f dB\n", sparams_output_return_loss_db(s_chain));

    printf("\n========================================\n");
    printf("LNA Design Example Complete\n");
    printf("========================================\n");

    return 0;
}
