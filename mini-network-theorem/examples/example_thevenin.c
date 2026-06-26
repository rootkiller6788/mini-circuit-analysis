/**
 * @file example_thevenin.c
 * @brief End-to-End Example: Thevenin Equivalent of a Multi-Source Circuit
 *
 * This example demonstrates:
 *   1. Building a circuit topology programmatically
 *   2. Computing the Thevenin equivalent at specified terminals
 *   3. Computing the Norton equivalent via source transformation
 *   4. Verifying the Maximum Power Transfer Theorem
 *   5. Computing power delivered to an arbitrary load
 *
 * Circuit analyzed:
 *       R1=4ohm    Node1    R2=6ohm
 *   +──/\/\/\/───┬──────/\/\/\/───+
 *   |            │                |
 *  V1=12V       R3=12ohm        V2=6V
 *   |            │                |
 *   +────────────┴────────────────+
 *              GND (Node0)
 *
 * We want the Thevenin equivalent at terminals (Node1, GND).
 */

#include "network_theorem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("========================================\n");
    printf(" Example: Thevenin & Norton Equivalents\n");
    printf("========================================\n\n");

    /* Build circuit with 2 nodes (0=GND, 1=Node1), 3 branches, 4 elements */
    CircuitTopology *ckt = topology_create(2, 3, 5, 0);
    if (!ckt) {
        printf("ERROR: Failed to allocate topology\n");
        return 1;
    }

    /* Element 0: Voltage source V1=12V from GND to Node1 */
    ckt->elements[0].type = ELEM_VOLTAGE_SOURCE;
    ckt->elements[0].node_from = 0;
    ckt->elements[0].node_to = 1;
    ckt->elements[0].value = 12.0;

    /* Element 1: Resistor R1=4ohm from Node1 to GND */
    ckt->elements[1].type = ELEM_RESISTOR;
    ckt->elements[1].node_from = 1;
    ckt->elements[1].node_to = 0;
    ckt->elements[1].value = 4.0;

    /* Element 2: Voltage source V2=6V from GND to Node1 */
    ckt->elements[2].type = ELEM_VOLTAGE_SOURCE;
    ckt->elements[2].node_from = 0;
    ckt->elements[2].node_to = 1;
    ckt->elements[2].value = 6.0;

    /* Element 3: Resistor R2=6ohm from Node1 to GND */
    ckt->elements[3].type = ELEM_RESISTOR;
    ckt->elements[3].node_from = 1;
    ckt->elements[3].node_to = 0;
    ckt->elements[3].value = 6.0;

    /* Element 4: Resistor R3=12ohm from Node1 to GND */
    ckt->elements[4].type = ELEM_RESISTOR;
    ckt->elements[4].node_from = 1;
    ckt->elements[4].node_to = 0;
    ckt->elements[4].value = 12.0;

    printf("Circuit built: 2 voltage sources, 3 resistors in parallel\n");
    printf("  V1 = 12V, R1 = 4 ohm\n");
    printf("  V2 = 6V,  R2 = 6 ohm\n");
    printf("  R3 = 12 ohm\n\n");

    /* --- Step 1: Manual analysis using Millman's Theorem ---
     * Three parallel branches from Node1 to GND:
     *   Branch 1: V1=12V in series with R1=4ohm → I_eq1 = 12/4 = 3A
     *   Branch 2: V2=6V in series with R2=6ohm  → I_eq2 = 6/6 = 1A
     *   Branch 3: R3=12ohm (no source)           → I_eq3 = 0
     *
     * V_common = (12/4 + 6/6 + 0) / (1/4 + 1/6 + 1/12)
     *          = (3 + 1 + 0) / (0.25 + 0.1667 + 0.0833)
     *          = 4 / 0.5 = 8V
     *
     * Z_th = 4 || 6 || 12 = 1 / (1/4 + 1/6 + 1/12) = 1 / 0.5 = 2 ohm
     */
    double v_branch[] = {12.0, 6.0, 0.0};
    double r_branch[] = {4.0, 6.0, 12.0};
    MillmanResult mr;
    if (millman_solve(v_branch, r_branch, 3, &mr) == 0) {
        printf("Millman's Theorem Analysis:\n");
        printf("  V_common = %.4f V (expected 8.0 V)\n", mr.common_node_voltage);
        free(mr.branch_voltages);
        free(mr.branch_impedances);
    }

    /* Equivalent resistance: 4 || 6 || 12 */
    ComplexImpedance z1 = impedance_resistor(4.0);
    ComplexImpedance z2 = impedance_resistor(6.0);
    ComplexImpedance z3 = impedance_resistor(12.0);
    ComplexImpedance zp12 = impedance_parallel(z1, z2);
    ComplexImpedance zp123 = impedance_parallel(zp12, z3);
    printf("  Z_th = %.4f ohm (expected 2.0 ohm)\n\n", zp123.real);

    /* --- Step 2: Thevenin & Norton Computation --- */
    TheveninEquivalent thev;
    if (compute_thevenin(ckt, 1, 0, &thev) == 0) {
        printf("Thevenin Equivalent:\n");
        printf("  V_th = %.4f V\n", thev.V_th.dc_offset);
        printf("  Z_th = %.4f + j%.4f ohm\n", thev.Z_th.real, thev.Z_th.imag);
    }

    /* Compute Norton via conversion */
    NortonEquivalent nort;
    thevenin_to_norton(&thev, &nort);
    printf("\nNorton Equivalent (from Thevenin):\n");
    printf("  I_n = %.4f A\n", nort.I_n.dc_offset);
    printf("  Y_n = %.4f + j%.4f S\n", nort.Y_n.real, nort.Y_n.imag);

    /* --- Step 3: Maximum Power Transfer --- */
    MaxPowerTransferResult mpt;
    if (max_power_transfer(&thev, &mpt) == 0) {
        printf("\nMaximum Power Transfer Analysis:\n");
        printf("  Source: Z_s = %.2f + j%.2f ohm\n",
               mpt.source_impedance_real, mpt.source_impedance_imag);
        printf("  Optimal load: Z_L = %.2f + j%.2f ohm\n",
               mpt.optimal_load_real, mpt.optimal_load_imag);
        printf("  Maximum power: P_max = %.4f W\n", mpt.max_power);
        printf("  Matching type: %s\n", mpt.is_ac ? "AC conjugate match" : "DC resistive match");
    }

    /* --- Step 4: Power delivered to non-optimal loads --- */
    double loads[] = {1.0, 2.0, 4.0, 10.0};
    printf("\nPower delivered to various loads (Vth=%.2fV, Rth=%.2f ohm):\n",
           thev.V_th.dc_offset, thev.Z_th.real);
    for (int i = 0; i < 4; i++) {
        double RL = loads[i];
        double P = (thev.V_th.dc_offset * thev.V_th.dc_offset * RL)
                 / ((thev.Z_th.real + RL) * (thev.Z_th.real + RL));
        printf("  RL = %5.1f ohm → P = %8.4f W\n", RL, P);
    }

    /* --- Step 5: Power calculations --- */
    printf("\nPower Analysis:\n");
    double P = power_dissipated(thev.V_th.dc_offset, 0.0, thev.Z_th.real + 2.0);
    printf("  Power in 2 ohm load: %.4f W\n", P);
    double pf = power_factor(0.0, 0.0);
    printf("  Power factor (DC): %.4f\n", pf);
    double S = apparent_power(thev.V_th.dc_offset, thev.V_th.dc_offset / (thev.Z_th.real + 2.0));
    printf("  Apparent power (DC = real): %.4f VA\n", S);

    /* Cleanup */
    topology_free(ckt);

    printf("\n========================================\n");
    printf(" Example Complete\n");
    printf("========================================\n");
    return 0;
}
