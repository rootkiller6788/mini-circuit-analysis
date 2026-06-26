/**
 * @file example_rc_transient.c
 * @brief Transient analysis of an RC charging circuit
 *
 * L6 Canonical Problem: RC circuit step response.
 * V = V0 * (1 - exp(-t/τ)) where τ = R*C.
 *
 * Circuit:
 *   V1 = 5V DC step from t=0
 *   R1 = 1kΩ in series with C1
 *   C1 = 1µF to GND
 *
 * τ = 1e3 * 1e-6 = 1ms
 *
 * At t = τ:  V = 5 * (1 - 1/e) = 3.161 V
 * At t = 5τ: V = 5 * (1 - e^{-5}) ≈ 4.966 V (≈99.3% of final value)
 *
 * Reference: Hayt, Kemmerly & Durbin, "Engineering Circuit Analysis"
 *            9th ed. (2019) §8.2 — The Source-Free RC Circuit
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spice_core.h"

int main(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Example: RC Circuit Transient Response      ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    const char *cir_file = "__example_rc_tran.cir";
    FILE *fp = fopen(cir_file, "w");
    if (!fp) { fprintf(stderr, "Cannot create netlist file\n"); return 1; }

    fprintf(fp, "RC Circuit — Transient Analysis\n");
    fprintf(fp, "V1 1 0 DC 5\n");
    fprintf(fp, "R1 1 2 1000\n");
    fprintf(fp, "C1 2 0 1u IC=0\n");
    fprintf(fp, ".TRAN 1e-4 5e-3\n");
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

    printf("\nRunning transient analysis...\n");
    printf("  R = 1 kΩ, C = 1 µF, τ = RC = 1 ms\n");
    printf("  Time step = 100 µs, Stop time = 5 ms\n\n");

    /* Run transient */
    if (spice_simulator_run_tran(&sim) != 0) {
        fprintf(stderr, "Warning: transient analysis had convergence issues\n");
    }

    if (sim.tran_result && sim.tran_result->num_steps > 0) {
        printf("Simulation completed: %d time steps\n\n", sim.tran_result->num_steps);

        /* Export data */
        spice_tran_export_csv(sim.tran_result, "__example_rc_output.csv");
        printf("Data exported to __example_rc_output.csv\n");

        /* Sample some time points and verify against theory */
        double R = 1000.0, C = 1e-6, tau = R * C;
        double V0 = 5.0;

        printf("\n  %-12s %-12s %-12s %s\n", "Time(ms)", "V_cap(V)", "Theory(V)", "Match?");
        printf("  %-12s %-12s %-12s %s\n", "────────", "────────", "────────", "─────");

        /* Check at key time points */
        double check_times[] = {0.001, 0.002, 0.003, 0.005}; /* τ, 2τ, 3τ, 5τ */
        for (int ci = 0; ci < 4; ci++) {
            double t_target = check_times[ci];
            double v_measured = 0.0;
            int found = 0;

            for (int32_t s = 0; s < sim.tran_result->num_steps; s++) {
                if (sim.tran_result->time_points[s] >= t_target - 1e-6) {
                    /* node 2 (index 2) should be the capacitor node */
                    if (sim.netlist.num_nodes >= 2)
                        v_measured = sim.tran_result->node_voltages[s][2];
                    found = 1;
                    break;
                }
            }

            if (found) {
                double v_theory = V0 * (1.0 - exp(-t_target / tau));
                int match = (fabs(v_measured - v_theory) < 0.1);
                printf("  %-12.3f %-12.6f %-12.6f %s\n",
                       t_target * 1000, v_measured, v_theory,
                       match ? "✓" : "✗");
            }
        }
    }

    /* Print full results */
    spice_simulator_print_results(&sim);

    spice_simulator_destroy(&sim);
    remove(cir_file);

    printf("\nExample complete.\n\n");
    return 0;
}
