#include "../include/circuit_topology.h"
#include "../include/circuit_mna.h"
#include "../include/circuit_analysis.h"
#include <stdio.h>
#include <math.h>

#define M_PI 3.14159265358979323846

int main(void) {
    double R = 100.0;
    double L = 10e-3;
    double Cv = 100e-9;
    double f_res = 1.0 / (2.0 * M_PI * sqrt(L * Cv));
    double Q = (1.0 / R) * sqrt(L / Cv);
    printf("=== RLC Circuit AC Analysis ===\n");
    printf("R=%.1f Ohm, L=%.2f mH, C=%.1f nF\n", R, L*1e3, Cv*1e9);
    printf("f0=%.2f kHz, Q=%.2f\n\n", f_res/1e3, Q);
    ct_circuit_t c;
    ct_circuit_init(&c, "Series RLC");
    ct_add_node(&c, "n1", 0);
    ct_add_node(&c, "n2", 0);
    ct_add_branch(&c, CT_ELEM_VSOURCE, 1, 0, 1.0, 0.0, "Vin");
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, R, 0.0, "R1");
    ct_add_branch(&c, CT_ELEM_INDUCTOR, 2, 0, L, 0.0, "L1");
    ct_add_branch(&c, CT_ELEM_CAPACITOR, 2, 0, Cv, 0.0, "C1");
    double freqs[] = {f_res*0.1, f_res*0.5, f_res, f_res*2, f_res*10};
    printf("Frequency Sweep:\n");
    for (int i = 0; i < 5; i++) {
        ct_ac_result_t ac;
        ct_mna_ac_solve(&c, freqs[i], &ac);
        printf("  %.1f Hz: |Vout|=%.4f V, phase=%.1f deg\n",
               freqs[i], ac.magnitude_v[2], ac.phase_v[2]);
    }
    printf("Done.\n");
    return 0;
}
