/**
 * @file example_bridge.c
 * @brief End-to-End Example: Wheatstone Bridge and Delta-Wye Transformation
 *
 * Demonstrates:
 *   1. Wheatstone bridge analysis (balanced and unbalanced)
 *   2. Delta-Wye transformation for bridge simplification
 *   3. Two-port parameter analysis of the bridge
 *   4. Reciprocity verification
 *
 * L7 Application: Strain gauge measurement system (Boeing 787)
 *   - Wheatstone bridge is the standard circuit for precision
 *     resistance measurement in strain gauges
 *   - Temperature compensation uses dummy gauge in adjacent arm
 *   - Smart grid application: power line sag monitoring
 *
 * L8 Advanced: Sensitivity analysis of the bridge circuit
 */

#include "network_theorem.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("========================================\n");
    printf(" Example: Wheatstone Bridge Analysis\n");
    printf("========================================\n\n");

    /* --- Part 1: Balanced Bridge --- */
    printf("--- Balanced Wheatstone Bridge ---\n");
    printf("R1=120ohm, R2=120ohm, R3=350ohm, R4=350ohm, Vex=5V\n\n");

    WheatstoneBridge bridge;
    wheatstone_analyze(120.0, 120.0, 350.0, 350.0, 1e6, 5.0, &bridge);

    printf("Bridge status: %s\n", bridge.is_balanced ? "BALANCED" : "UNBALANCED");
    printf("  V_left  = 5 * 350/(120+350) = %.4f V\n",
           5.0 * 350.0 / 470.0);
    printf("  V_right = 5 * 350/(120+350) = %.4f V\n",
           5.0 * 350.0 / 470.0);
    printf("  V_output = V_left - V_right = %.6f V\n", bridge.V_output);
    printf("  Galvanometer current: %.6f uA\n\n", bridge.I_g * 1e6);

    /* --- Part 2: Unbalanced Bridge (simulating strain) --- */
    printf("--- Unbalanced Bridge (Delta R = +2 ohm on R1) ---\n");
    printf("R1=122ohm, R2=120ohm, R3=350ohm, R4=350ohm, Vex=5V\n\n");

    wheatstone_analyze(122.0, 120.0, 350.0, 350.0, 1e6, 5.0, &bridge);

    printf("Bridge status: %s\n", bridge.is_balanced ? "BALANCED" : "UNBALANCED");
    printf("  V_left  = 5 * 350/(122+350) = %.4f V\n",
           5.0 * 350.0 / 472.0);
    printf("  V_right = 5 * 350/(120+350) = %.4f V\n",
           5.0 * 350.0 / 470.0);
    printf("  V_output = %.6f V\n", bridge.V_output);
    printf("  Sensitivity: %.2f uV per ohm change\n",
           bridge.V_output * 1e6 / 2.0);

    /* --- Part 3: Bridge Thevenin Equivalent --- */
    printf("\n--- Bridge Thevenin Equivalent at Output ---\n");

    /* R_th of bridge output:
     * With V_source shorted: looking into the output terminals:
     *   Left side: R1||R3 = 122||350 = 122*350/(122+350) = 90.47 ohm
     *   Right side: R2||R4 = 120||350 = 120*350/(120+350) = 89.36 ohm
     *   R_th = 90.47 + 89.36 = 179.83 ohm
     */
    double R_left = combine_parallel_resistors(122.0, 350.0);
    double R_right = combine_parallel_resistors(120.0, 350.0);
    double R_th_bridge = combine_series_resistors(R_left, R_right);

    printf("  R_th = R1||R3 + R2||R4\n");
    printf("       = %.2f||%.0f + %.0f||%.0f\n", 122.0, 350.0, 120.0, 350.0);
    printf("       = %.4f + %.4f = %.4f ohm\n", R_left, R_right, R_th_bridge);
    printf("  V_th = V_output = %.6f V\n", bridge.V_output);

    /* Maximum power to a load at the bridge output */
    double P_max_bridge = bridge.V_output * bridge.V_output
                         / (4.0 * R_th_bridge);
    printf("  P_max (matched load) = %.6f uW\n", P_max_bridge * 1e6);

    /* --- Part 4: Delta-Wye Transformation --- */
    printf("\n--- Delta-Wye Transformation for Bridge Simplification ---\n");

    /* The bridge contains a delta: R1(122), R3(350), and the upper
     * diagonal between the two midpoints. We can transform this delta
     * to a wye for easier analysis. */

    /* Consider delta formed by R1, R2, and the top wire (R_top ≈ 0):
     * This is actually a simpler case. Let's demonstrate with a
     * numeric example of the Y-Delta transformation. */
    DeltaNetwork delta_net = {10.0, 20.0, 30.0};
    WyeNetwork wye_net;
    delta_to_wye(&delta_net, &wye_net);

    printf("  Delta network: R12=%.1f, R23=%.1f, R31=%.1f ohm\n",
           delta_net.R12, delta_net.R23, delta_net.R31);
    printf("  Equivalent Wye: R1=%.4f, R2=%.4f, R3=%.4f ohm\n",
           wye_net.R1, wye_net.R2, wye_net.R3);

    /* Verify equivalence: R_12(wye) = R1 + R2 =? R12||(R23+R31) */
    double R12_wye = wye_net.R1 + wye_net.R2;
    double R12_delta_check = combine_parallel_resistors(
        delta_net.R12,
        combine_series_resistors(delta_net.R23, delta_net.R31));
    printf("  Verification: R12 via Wye = %.4f, via Delta = %.4f ohm\n",
           R12_wye, R12_delta_check);

    /* Round-trip: Wye -> Delta -> Wye */
    DeltaNetwork delta2;
    wye_to_delta(&wye_net, &delta2);
    WyeNetwork wye2;
    delta_to_wye(&delta2, &wye2);
    printf("  Round-trip: original R1=%.4f, recovered R1=%.4f ohm\n",
           wye_net.R1, wye2.R1);

    /* --- Part 5: Two-Port Analysis of Bridge --- */
    printf("\n--- Two-Port Parameter Analysis ---\n");

    /* Model the bridge as a two-port with:
     * Port 1: excitation input
     * Port 2: bridge output */

    /* Approximate Z-parameters for the bridge network */
    ZParameters zp;
    zp.z11.real = 122.0 + 350.0;  /* Total input impedance (approx) */
    zp.z11.imag = 0.0;
    zp.z12.real = 0.0;            /* Transfer impedance */
    zp.z12.imag = 0.0;
    zp.z21.real = 0.0;
    zp.z21.imag = 0.0;
    zp.z22.real = R_th_bridge;
    zp.z22.imag = 0.0;

    printf("  Z-parameters (approximate DC):\n");
    printf("    z11 = %.2f ohm (input impedance)\n", zp.z11.real);
    printf("    z22 = %.2f ohm (output impedance)\n", zp.z22.real);

    /* Convert to ABCD for cascading */
    ABCDParameters abcd;
    if (z_to_abcd(&zp, &abcd) == 0) {
        printf("  ABCD parameters:\n");
        printf("    A = %.4f, B = %.4f ohm\n", abcd.A_real, abcd.B_real);
        printf("    C = %.4f S, D = %.4f\n", abcd.C_real, abcd.D_real);
    }

    /* Reciprocity check */
    ReciprocityParams rp;
    verify_reciprocity(&zp, &rp);
    printf("  Reciprocity: %s (error = %.2e)\n",
           rp.is_reciprocal ? "YES" : "NO", rp.reciprocity_error);

    /* --- L7 Application Context --- */
    printf("\n--- L7 Application: Boeing 787 Strain Gauge System ---\n");
    printf("  The Wheatstone bridge is deployed in hundreds of locations\n");
    printf("  on the Boeing 787 airframe for structural health monitoring.\n");
    printf("  Each bridge uses 350-ohm foil strain gauges with 5V excitation.\n");
    printf("  Typical strain resolution: 1 microstrain (1e-6)\n");
    printf("  Corresponding resistance change: 0.00024 ohm (GF=2.0)\n");
    printf("  Bridge output for 1 microstrain: %.3f uV\n\n",
           (5.0 * 0.00024 / 350.0) * 1e6);

    /* --- L8 Advanced: Sensitivity Analysis --- */
    printf("--- L8 Advanced: Bridge Sensitivity Analysis ---\n");
    printf("  Sensitivity S = d(Vout)/d(R1) at balance point:\n");
    printf("  S = Vex * R3 / (R1+R3)^2\n");
    double S_sensitivity = 5.0 * 350.0 / (470.0 * 470.0);
    printf("  S = 5 * 350 / %.0f^2 = %.6f V/ohm\n", 470.0, S_sensitivity);
    printf("  For 1ppm resistance change: dV = %.3f uV\n",
           S_sensitivity * 120.0 * 1e-6 * 1e6);

    printf("\n========================================\n");
    printf(" Example Complete\n");
    printf("========================================\n");
    return 0;
}
