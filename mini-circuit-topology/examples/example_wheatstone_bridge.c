#include "../include/circuit_topology.h"
#include "../include/circuit_mna.h"
#include "../include/circuit_analysis.h"
#include <stdio.h>

int main(void) {
    ct_circuit_t bal;
    ct_circuit_init(&bal, "Balanced Bridge");
    ct_add_node(&bal, "Vexc", 0);
    ct_add_node(&bal, "n2", 0);
    ct_add_node(&bal, "n3", 0);
    ct_add_branch(&bal, CT_ELEM_VSOURCE, 1, 0, 5.0, 0.0, "Vexc");
    ct_add_branch(&bal, CT_ELEM_RESISTOR, 1, 2, 350.0, 0.0, "R1");
    ct_add_branch(&bal, CT_ELEM_RESISTOR, 1, 3, 350.0, 0.0, "R2");
    ct_add_branch(&bal, CT_ELEM_RESISTOR, 2, 0, 350.0, 0.0, "R3");
    ct_add_branch(&bal, CT_ELEM_RESISTOR, 3, 0, 350.0, 0.0, "R4");
    ct_dc_result_t r;
    ct_mna_dc_solve(&bal, &r);
    double vout = r.node_voltages[2] - r.node_voltages[3];
    printf("Balanced: Vout=%.6f V\n", vout);
    printf("V(n2)=%.4f V, V(n3)=%.4f V\n", r.node_voltages[2], r.node_voltages[3]);
    ct_circuit_t unb;
    ct_circuit_init(&unb, "Unbalanced");
    ct_add_node(&unb, "Vexc", 0);
    ct_add_node(&unb, "n2", 0);
    ct_add_node(&unb, "n3", 0);
    ct_add_branch(&unb, CT_ELEM_VSOURCE, 1, 0, 5.0, 0.0, "Vexc");
    ct_add_branch(&unb, CT_ELEM_RESISTOR, 1, 2, 350.0, 0.0, "R1");
    ct_add_branch(&unb, CT_ELEM_RESISTOR, 1, 3, 350.0, 0.0, "R2");
    ct_add_branch(&unb, CT_ELEM_RESISTOR, 2, 0, 353.5, 0.0, "R3");
    ct_add_branch(&unb, CT_ELEM_RESISTOR, 3, 0, 350.0, 0.0, "R4");
    ct_dc_result_t r2;
    ct_mna_dc_solve(&unb, &r2);
    double vout2 = r2.node_voltages[2] - r2.node_voltages[3];
    printf("Unbalanced (+1%% R3): Vout=%.6f V\n", vout2);
    printf("Sensitivity: %.6f V/Ohm\n", (vout2 - vout) / 3.5);
    printf("Done.\n");
    return 0;
}
