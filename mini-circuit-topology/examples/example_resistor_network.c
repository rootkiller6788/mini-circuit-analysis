#include "../include/circuit_topology.h"
#include "../include/circuit_mna.h"
#include "../include/circuit_analysis.h"
#include <stdio.h>

int main(void) {
    ct_circuit_t c;
    ct_circuit_init(&c, "Resistor Network");
    ct_add_node(&c, "n1", 0);
    ct_add_node(&c, "n2", 0);
    ct_add_branch(&c, CT_ELEM_VSOURCE, 1, 0, 12.0, 0.0, "V1");
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, 1000.0, 0.0, "R1");
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 0, 2000.0, 0.0, "R2");
    ct_dc_result_t r;
    ct_mna_dc_solve(&c, &r);
    printf("V(n1)=%.4f V, V(n2)=%.4f\n", r.node_voltages[1], r.node_voltages[2]);
    printf("P_diss=%.4f W\n", r.total_power_dissipated);
    double V_th, R_th;
    ct_thevenin_equivalent(&c, 2, 0, &V_th, &R_th);
    printf("Thevenin: Vth=%.4f V, Rth=%.2f Ohm\n", V_th, R_th);
    printf("Done.\n");
    return 0;
}
