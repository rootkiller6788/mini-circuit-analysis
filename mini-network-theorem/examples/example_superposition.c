/**
 * @file example_superposition.c
 * @brief End-to-End Example: Superposition Theorem in Multi-Source Circuit
 *
 * Demonstrates the Superposition Theorem by analyzing a circuit with
 * both a voltage source and a current source, showing that the total
 * response equals the sum of individual responses.
 *
 * Circuit:
 *       R1=3ohm     Node1     R3=2ohm
 *   +──/\/\/\/───┬──────────/\/\/\/───+
 *   |            │                    |
 *  V1=9V        R2=6ohm            I1=2A
 *   |            │                    |
 *   +────────────┴────────────────────+
 *              GND (Node0)
 *
 * L7 Application context:
 *   This circuit pattern appears in sensor signal conditioning
 *   (Boeing strain gauge bridge excitation with offset compensation)
 *   and smart grid distributed generation analysis.
 */

#include "network_theorem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("========================================\n");
    printf(" Example: Superposition Theorem\n");
    printf("========================================\n\n");

    /* Build circuit */
    CircuitTopology *ckt = topology_create(2, 3, 3, 0);
    if (!ckt) { printf("ERROR: allocation failed\n"); return 1; }

    /* V1: 9V from GND(0) to Node1(1) */
    ckt->elements[0].type = ELEM_VOLTAGE_SOURCE;
    ckt->elements[0].node_from = 0;
    ckt->elements[0].node_to = 1;
    ckt->elements[0].value = 9.0;

    /* R1: 3ohm from Node1 to GND */
    ckt->elements[1].type = ELEM_RESISTOR;
    ckt->elements[1].node_from = 1;
    ckt->elements[1].node_to = 0;
    ckt->elements[1].value = 3.0;

    /* I1: 2A from GND to Node1 */
    ckt->elements[2].type = ELEM_CURRENT_SOURCE;
    ckt->elements[2].node_from = 0;
    ckt->elements[2].node_to = 1;
    ckt->elements[2].value = 2.0;

    printf("Circuit: V1=9V source and I1=2A source with R1=3ohm\n\n");

    /* --- Manual Superposition Analysis --- */
    printf("Manual Superposition Analysis:\n");

    /* With V1 alone (I1 deactivated → open circuit):
     * V_node1 = 9V (voltage source directly across R1)
     */
    printf("  V1 alone (I1 open): V_node1 = 9.000 V\n");

    /* With I1 alone (V1 deactivated → short circuit):
     * R1 sees 2A flowing through it from node1 to GND
     * V_node1 = I1 * R1 = 2 * 3 = 6V
     */
    printf("  I1 alone (V1 short): V_node1 = 2A * 3ohm = 6.000 V\n");

    /* Superposition total */
    printf("  Total (superposition): V_node1 = 9.0 + 6.0 = 15.000 V\n\n");

    /* --- Automated Superposition Solver --- */
    SuperpositionResult spr;
    if (superposition_solve(ckt, 1, &spr) == 0) {
        printf("Automated Superposition Solver Result:\n");
        printf("  Total V_node1 = %.4f V\n", spr.total_response);
        printf("  Number of independent sources: %u\n", spr.num_sources);
        for (uint32_t i = 0; i < spr.num_sources; i++) {
            printf("    Source %u contribution: %+.4f V\n",
                   i, spr.individual_contributions[i]);
        }
        free(spr.individual_contributions);
    }

    /* --- Linearity Verification --- */
    double error = superposition_verify(ckt, 1);
    printf("\nLinearity verification (|V_all - V_superpos|): %.2e\n", error);
    if (error < 1e-6) {
        printf("  ✓ Superposition holds (circuit is linear)\n");
    } else {
        printf("  ⚠ Superposition error detected\n");
    }

    /* --- L7 Application: Sensor Bridge Offset Compensation ---
     * In Boeing strain gauge applications, multiple excitation sources
     * are used to compensate for thermal offset. The superposition
     * principle allows independent analysis of each excitation.
     */
    printf("\nL7 Application: Sensor Bridge Offset Compensation\n");
    printf("  Context: Boeing 787 strain gauge system\n");
    printf("  Multiple excitation sources for thermal drift compensation\n");
    printf("  Superposition enables independent calibration of each channel\n");

    /* --- Power Analysis --- */
    printf("\nPower Analysis at Node1:\n");
    double V_node = 15.0;
    double I_total = 0.0;
    if (ckt->elements[0].type == ELEM_VOLTAGE_SOURCE) {
        I_total = V_node / ckt->elements[1].value + ckt->elements[2].value;
    }
    double P_total = V_node * I_total;
    printf("  Voltage: %.2f V\n", V_node);
    printf("  Power delivered to R1: %.2f W (V^2/R = 225/3)\n",
           V_node * V_node / 3.0);
    printf("  Total power: %.2f W\n", P_total);

    topology_free(ckt);

    printf("\n========================================\n");
    printf(" Example Complete\n");
    printf("========================================\n");
    return 0;
}
