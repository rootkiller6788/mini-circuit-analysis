/**
 * @file example_divider_dc.c
 * @brief DC analysis of a resistive voltage divider
 *
 * L6 Canonical Problem: Find the DC operating point of a two-resistor
 * voltage divider circuit. Demonstrates netlist parsing, MNA assembly,
 * Newton-Raphson DC solve, and result querying.
 *
 * Circuit:
 *   V1 (10V DC) between node "src" and GND
 *   R1 (1kΩ) between "src" and "out"
 *   R2 (2kΩ) between "out" and GND
 *
 * Expected output: V(out) = 10 * 2000 / (1000 + 2000) = 6.667 V
 *
 * Reference: Sedra & Smith (2020) §1.6, Voltage Dividers
 */

#include <stdio.h>
#include <stdlib.h>
#include "spice_core.h"

int main(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  Example: Resistive Voltage Divider (DC)     ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* Create a temporary SPICE netlist */
    const char *cir_file = "__example_divider.cir";
    FILE *fp = fopen(cir_file, "w");
    if (!fp) { fprintf(stderr, "Cannot create netlist file\n"); return 1; }

    fprintf(fp, "Resistive Voltage Divider — DC Analysis\n");
    fprintf(fp, "V1 src 0 DC 10\n");
    fprintf(fp, "R1 src out 1000\n");
    fprintf(fp, "R2 out 0 2000\n");
    fprintf(fp, ".OP\n");
    fprintf(fp, ".END\n");
    fclose(fp);

    /* Initialize and load */
    spice_simulator_t sim;
    spice_simulator_init(&sim);

    printf("Loading netlist: %s\n", cir_file);
    if (spice_simulator_load(&sim, cir_file) != 0) {
        fprintf(stderr, "Failed to load netlist\n");
        remove(cir_file);
        return 1;
    }

    spice_netlist_print_summary(&sim.netlist);

    /* Run DC analysis */
    printf("\nRunning DC operating point analysis...\n");
    if (spice_simulator_run_dc(&sim) != 0) {
        fprintf(stderr, "DC analysis failed\n");
        spice_simulator_destroy(&sim);
        remove(cir_file);
        return 1;
    }

    /* Query results */
    double vsrc, vout;
    if (spice_get_node_voltage(&sim, "src", &vsrc) == 0) {
        printf("  V(src) = %+.6f V\n", vsrc);
    }
    if (spice_get_node_voltage(&sim, "out", &vout) == 0) {
        printf("  V(out) = %+.6f V\n", vout);

        /* Verify using voltage divider formula */
        double vdiv = vsrc * 2000.0 / (1000.0 + 2000.0);
        printf("\n  Expected (voltage divider): %+.6f V\n", vdiv);
        printf("  Error: %.2e V\n", vout - vdiv);
    }

    /* Print all results */
    spice_simulator_print_results(&sim);

    /* Cleanup */
    spice_simulator_destroy(&sim);
    remove(cir_file);

    printf("\nExample complete.\n\n");
    return 0;
}
