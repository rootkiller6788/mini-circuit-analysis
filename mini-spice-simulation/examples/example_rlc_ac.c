/**
 * @file example_rlc_ac.c
 * @brief AC frequency analysis of an RLC bandpass filter
 *
 * L6 Canonical Problem: RLC filter frequency response.
 *
 * Circuit:
 *   Vin (AC 1V) → R1 (100Ω) → L1 (10mH) → node "out"
 *   C1 (1µF) from "out" to GND
 *
 * This forms a series RLC with output taken across C1.
 * Transfer function H(s) = 1/(LC*s^2 + RC*s + 1)
 * Resonant frequency: f0 = 1/(2π√(LC)) ≈ 1592 Hz
 * Quality factor: Q = (1/R)*√(L/C) = (1/100)*√(0.01/1e-6) ≈ 1.0
 *
 * Reference: Hayt, Kemmerly & Durbin, "Engineering Circuit Analysis"
 *            9th ed. (2019) §16 — Frequency Response
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spice_core.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Example: RLC Bandpass Filter (AC Analysis)  ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    const char *cir_file = "__example_rlc_ac.cir";
    FILE *fp = fopen(cir_file, "w");
    if (!fp) { fprintf(stderr, "Cannot create netlist file\n"); return 1; }

    /* Series RLC bandpass filter */
    fprintf(fp, "RLC Bandpass Filter — AC Analysis\n");
    fprintf(fp, "Vin in 0 DC 0 AC 1 0\n");
    fprintf(fp, "R1 in mid 100\n");
    fprintf(fp, "L1 mid out 10mH\n");
    fprintf(fp, "C1 out 0 1u\n");
    fprintf(fp, ".AC DEC 10 10 100k\n");
    fprintf(fp, ".END\n");
    fclose(fp);

    spice_simulator_t sim;
    spice_simulator_init(&sim);

    printf("Loading netlist...\n");
    if (spice_simulator_load(&sim, cir_file) != 0) {
        fprintf(stderr, "Failed to load netlist\n");
        remove(cir_file);
        return 1;
    }

    spice_netlist_print_summary(&sim.netlist);

    /* RLC theory */
    double R = 100.0, L = 0.01, C = 1e-6;
    double f0 = 1.0 / (2.0 * M_PI * sqrt(L * C));
    double Q = (1.0 / R) * sqrt(L / C);
    printf("\n  Theoretical parameters:\n");
    printf("    R = %.0f Ω, L = %.0f mH, C = %.1f µF\n", R, L*1000, C*1e6);
    printf("    Resonant frequency f0 = %.1f Hz\n", f0);
    printf("    Quality factor Q = %.2f\n", Q);
    printf("    Bandwidth BW = f0/Q = %.1f Hz\n\n", f0 / Q);

    printf("Running AC analysis...\n");

    /* Run simulation */
    spice_simulator_run(&sim);

    if (sim.ac_result && sim.ac_result->converged) {
        printf("AC analysis converged: %d frequency points\n\n",
               sim.ac_result->num_freqs);

        spice_ac_export_csv(sim.ac_result, "__example_rlc_ac_output.csv");
        printf("Data exported to __example_rlc_ac_output.csv\n");

        /* Find the peak (resonance) */
        double max_mag = 0.0;
        double peak_freq = 0.0;

        for (int32_t f = 0; f < sim.ac_result->num_freqs; f++) {
            /* Find V(out) magnitude */
            /* "out" node — need to find its index */
            int out_idx = 0;
            for (int32_t n = 1; n <= sim.netlist.num_nodes; n++) {
                if (strcmp(sim.netlist.nodes[n].name, "out") == 0) {
                    out_idx = n;
                    break;
                }
            }
            if (out_idx > 0) {
                double mag = cabs(sim.ac_result->node_voltages[f][out_idx]);
                if (mag > max_mag) {
                    max_mag = mag;
                    peak_freq = sim.ac_result->frequencies[f];
                }
            }
        }

        printf("  Peak response: |V(out)| = %.4f at f = %.1f Hz\n",
               max_mag, peak_freq);
        printf("  Theoretical f0 = %.1f Hz (error = %.1f%%)\n",
               f0, fabs(peak_freq - f0) / f0 * 100.0);

        /* Print a few data points near resonance */
        printf("\n  %-12s %-14s %-14s\n", "Freq (Hz)", "|Vout|", "Phase (°)");
        printf("  %-12s %-14s %-14s\n", "────────", "──────", "────────");

        int out_idx = 0;
        for (int32_t n = 1; n <= sim.netlist.num_nodes; n++) {
            if (strcmp(sim.netlist.nodes[n].name, "out") == 0) {
                out_idx = n; break;
            }
        }

        /* Show points near resonance */
        for (int32_t f = 0; f < sim.ac_result->num_freqs; f++) {
            double freq = sim.ac_result->frequencies[f];
            if (freq >= f0 * 0.3 && freq <= f0 * 3.0) {
                double mag = cabs(sim.ac_result->node_voltages[f][out_idx]);
                double phase = atan2(cimag(sim.ac_result->node_voltages[f][out_idx]),
                                     creal(sim.ac_result->node_voltages[f][out_idx]))
                               * 180.0 / M_PI;
                printf("  %-12.1f %-14.6f %-14.2f\n", freq, mag, phase);
            }
        }
    }

    spice_simulator_print_results(&sim);
    spice_simulator_destroy(&sim);
    remove(cir_file);

    printf("\nExample complete.\n\n");
    return 0;
}
