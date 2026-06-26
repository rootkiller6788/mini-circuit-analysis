/**
 * @file example_cascade_analysis.c
 * @brief Example: Cascade Analysis of a Multi-Stage RF System
 *
 * Demonstrates ABCD cascade analysis of a complete RF receive chain:
 *   Antenna → BPF → LNA → Attenuator → Mixer → IF Amplifier
 *
 * Each component is modeled with its ABCD matrix, and the cascade
 * is computed by matrix multiplication. The overall gain, noise figure,
 * and input/output match are then derived from the combined S-parameters.
 *
 * This example also demonstrates de-embedding: removing the effect
 * of the input matching network to find the raw LNA S-parameters.
 *
 * Course: TU Munich HF Engineering — System cascade analysis
 *         ETH 227-0455 — RF system design
 */

#include <stdio.h>
#include <math.h>
#include "../include/two_port.h"
#include "../include/abcd_params.h"
#include "../include/conversion.h"
#include "../include/s_params.h"

int main(void) {
    printf("========================================\n");
    printf("Multi-Stage RF Receiver Cascade Analysis\n");
    printf("========================================\n\n");

    double freq = 1.57542e9;  /* GPS L1 frequency */
    double omega = 2.0 * M_PI * freq;
    double z0 = 50.0;

    /* ============================================================
     * Component 1: Band-Pass Filter (SAW filter)
     * ============================================================ */
    /* Ideal SAW filter: Insertion Loss = 2.5 dB, Return Loss = 15 dB
       s11 = 10^(-15/20) ≈ 0.178, s21 = 10^(-2.5/20) ≈ 0.75 */
    matrix2x2_t s_bpf = sparams_create(
        complex_make(0.178, 0.0),
        complex_make(0.75, 0.0),    /* s12 = s21 for passive reciprocal */
        complex_make(0.75, 0.0),
        complex_make(0.178, 0.0)
    );
    matrix2x2_t abcd_bpf = convert_s_to_abcd(s_bpf, z0);

    printf("Stage 1: GPS SAW Filter\n");
    printf("  IL = %.1f dB, RL = %.1f dB\n",
           2.5, 15.0);

    /* ============================================================
     * Component 2: LNA (Low Noise Amplifier)
     * ============================================================ */
    /* NF = 1.2 dB, Gain = 18 dB, decent match
       s21 = 10^(18/20) ≈ 7.94 */
    matrix2x2_t s_lna = sparams_create(
        complex_make(0.15, -0.10),   /* Input RL ≈ 15 dB */
        complex_make(0.01, 0.0),     /* Isolation ≈ 40 dB */
        complex_make(7.94, 0.0),     /* Gain = 18 dB */
        complex_make(0.12, -0.08)    /* Output RL ≈ 17 dB */
    );
    matrix2x2_t abcd_lna = convert_s_to_abcd(s_lna, z0);

    printf("Stage 2: LNA\n");
    printf("  Gain = 18.0 dB, NF = 1.2 dB\n");

    /* ============================================================
     * Component 3: Digital Step Attenuator
     * ============================================================ */
    /* ATT = 6 dB, well matched */
    matrix2x2_t s_att = sparams_attenuator(6.0);
    matrix2x2_t abcd_att = convert_s_to_abcd(s_att, z0);

    printf("Stage 3: Digital Attenuator\n");
    printf("  ATT = 6.0 dB\n");

    /* ============================================================
     * Component 4: Mixer (passive, double-balanced)
     * ============================================================ */
    /* Conversion loss = 7 dB, reasonable match */
    matrix2x2_t s_mixer = sparams_create(
        complex_make(0.10, 0.05),
        complex_make(0.446, 0.0),    /* 10^(-7/20) = 0.446 */
        complex_make(0.446, 0.0),
        complex_make(0.12, 0.06)
    );
    matrix2x2_t abcd_mixer = convert_s_to_abcd(s_mixer, z0);

    printf("Stage 4: Passive Mixer\n");
    printf("  Conversion Loss = 7.0 dB\n");

    /* ============================================================
     * Component 5: IF Amplifier
     * ============================================================ */
    /* Gain = 20 dB, NF = 3 dB */
    matrix2x2_t s_ifamp = sparams_create(
        complex_make(0.08, -0.04),
        complex_make(0.005, 0.0),
        complex_make(10.0, 0.0),    /* 20 dB gain */
        complex_make(0.10, -0.06)
    );
    matrix2x2_t abcd_ifamp = convert_s_to_abcd(s_ifamp, z0);

    printf("Stage 5: IF Amplifier\n");
    printf("  Gain = 20.0 dB, NF = 3.0 dB\n\n");

    /* ============================================================
     * Cascade All Stages
     * ============================================================ */
    matrix2x2_t stages[] = {
        abcd_bpf, abcd_lna, abcd_att, abcd_mixer, abcd_ifamp
    };
    matrix2x2_t abcd_total = abcd_cascade_n(stages, 5);
    matrix2x2_t s_total = convert_abcd_to_s(abcd_total, z0);

    printf("=== Cascade Results ===\n\n");
    printf("Parameter          Per-Stage          Cumulative\n");
    printf("---------          ---------          ----------\n");

    /* Gain breakdown */
    printf("Gain (dB):\n");
    printf("  BPF              -2.5 dB\n");
    printf("  LNA              +18.0 dB           %.1f dB\n",
           sparams_s21_db(convert_abcd_to_s(
               abcd_cascade_n(stages, 2), z0)));
    printf("  Attenuator       -6.0 dB            %.1f dB\n",
           sparams_s21_db(convert_abcd_to_s(
               abcd_cascade_n(stages, 3), z0)));
    printf("  Mixer            -7.0 dB            %.1f dB\n",
           sparams_s21_db(convert_abcd_to_s(
               abcd_cascade_n(stages, 4), z0)));
    printf("  IF Amp           +20.0 dB           %.1f dB\n",
           sparams_s21_db(s_total));

    printf("\nOverall System Performance:\n");
    printf("  Total Gain       = %.1f dB\n", sparams_s21_db(s_total));
    printf("  Input Return Loss  = %.1f dB\n",
           sparams_input_return_loss_db(s_total));
    printf("  Output Return Loss = %.1f dB\n",
           sparams_output_return_loss_db(s_total));

    /* ============================================================
     * De-embedding Demo
     * ============================================================ */
    printf("\n=== De-Embedding Demo ===\n\n");

    /* Simulate: Measure total chain → extract LNA alone */
    matrix2x2_t abcd_left_fixture = abcd_bpf;
    matrix2x2_t abcd_right_fixture = abcd_cascade_n(stages + 2, 3);
    matrix2x2_t abcd_dut = abcd_deembed(abcd_total,
                                         abcd_left_fixture,
                                         abcd_right_fixture);

    matrix2x2_t s_extracted = convert_abcd_to_s(abcd_dut, z0);
    printf("Extracted LNA S-parameters (from de-embedding):\n");
    printf("  s11 = %.3f ∠ %.1f°  (original: %.3f ∠ %.1f°)\n",
           complex_mag(s_extracted.m11),
           complex_arg(s_extracted.m11) * 180.0 / M_PI,
           complex_mag(s_lna.m11),
           complex_arg(s_lna.m11) * 180.0 / M_PI);
    printf("  s21 = %.2f ∠ %.1f°  (original: %.2f ∠ %.1f°)\n",
           complex_mag(s_extracted.m21),
           complex_arg(s_extracted.m21) * 180.0 / M_PI,
           complex_mag(s_lna.m21),
           complex_arg(s_lna.m21) * 180.0 / M_PI);

    printf("\n========================================\n");
    printf("Cascade Analysis Example Complete\n");
    printf("========================================\n");

    return 0;
}
