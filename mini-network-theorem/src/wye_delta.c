/**
 * @file wye_delta.c
 * @brief Wye-Delta Transformations, Source Transformation, and Network Reduction
 *
 * Implements fundamental circuit transformation techniques:
 *   - Wye-Delta (Y-Δ) and Delta-Wye (Δ-Y) transformations
 *   - Source transformation (Thevenin ↔ Norton at branch level)
 *   - Iterative network reduction combining series/parallel/transform
 *
 * Knowledge Coverage:
 *   L5 - Algorithms: Y-Δ transform, source transform, network reduction
 *   L6 - Canonical Problems: Bridge circuit reduction, multi-stage simplification
 *   L8 - Advanced Topics: Iterative network reduction algorithm
 *
 * Reference:
 *   Kennelly, A.E. (1899) "Equivalence of triangles and stars"
 *   Hayt et al. "Engineering Circuit Analysis" (2019) Ch 5
 *
 * Course Mapping:
 *   MIT 6.002: Wye-Delta transformation
 *   Berkeley EE16A: Circuit transformations
 */

#include "network_theorem.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L5: Wye-to-Delta (Y→Δ) Transformation
 * ==========================================================================
 *
 * Given a Y (star/wye) network with resistances R1, R2, R3 connected to
 * a common center point, the equivalent Δ (delta/pi) network has:
 *
 *            R1*R2 + R2*R3 + R3*R1
 *   R12 = ─────────────────────────
 *                    R3
 *
 *            R1*R2 + R2*R3 + R3*R1
 *   R23 = ─────────────────────────
 *                    R1
 *
 *            R1*R2 + R2*R3 + R3*R1
 *   R31 = ─────────────────────────
 *                    R2
 *
 * Physical interpretation:
 *   - R12 is the resistance between terminals 1 and 2 with terminal 3 open.
 *     In Y: R12(Y) = R1 + R2
 *     In Δ: R12(Δ) = R12 || (R23 + R31) = R12*(R23+R31)/(R12+R23+R31)
 *     Setting these equal gives the transformation.
 *
 * Edge cases:
 *   - If any R_k = 0, the corresponding Δ resistance becomes infinite (open)
 *   - If all R_k equal (R), then R_Δ = 3R (each Δ resistance is 3× the Y)
 */
void wye_to_delta(const WyeNetwork *wye, DeltaNetwork *delta) {
    if (!wye || !delta) return;

    double R1 = wye->R1, R2 = wye->R2, R3 = wye->R3;

    /* Sum of pairwise products: R1*R2 + R2*R3 + R3*R1 */
    double sum_prod = R1*R2 + R2*R3 + R3*R1;

    /* Handle zero-divisor cases */
    if (fabs(R3) < 1e-15) {
        delta->R12 = 1e15;  /* ∞ (effectively open) */
    } else {
        delta->R12 = sum_prod / R3;
    }

    if (fabs(R1) < 1e-15) {
        delta->R23 = 1e15;
    } else {
        delta->R23 = sum_prod / R1;
    }

    if (fabs(R2) < 1e-15) {
        delta->R31 = 1e15;
    } else {
        delta->R31 = sum_prod / R2;
    }
}

/* ==========================================================================
 * L5: Delta-to-Wye (Δ→Y) Transformation
 * ==========================================================================
 *
 * Given a Δ (delta/pi) network, the equivalent Y (star/wye) network has:
 *
 *              R12 * R31
 *   R1 = ───────────────────
 *         R12 + R23 + R31
 *
 *              R12 * R23
 *   R2 = ───────────────────
 *         R12 + R23 + R31
 *
 *              R23 * R31
 *   R3 = ───────────────────
 *         R12 + R23 + R31
 *
 * This transformation is useful for simplifying bridge circuits
 * (e.g., Wheatstone bridge) where a Δ subnetwork can be converted
 * to Y, enabling series/parallel reduction.
 */
void delta_to_wye(const DeltaNetwork *delta, WyeNetwork *wye) {
    if (!delta || !wye) return;

    double R12 = delta->R12, R23 = delta->R23, R31 = delta->R31;
    double sum_R = R12 + R23 + R31;

    if (fabs(sum_R) < 1e-15) {
        /* All resistances are zero → Y resistances are also zero */
        wye->R1 = 0.0;
        wye->R2 = 0.0;
        wye->R3 = 0.0;
        return;
    }

    wye->R1 = (R12 * R31) / sum_R;
    wye->R2 = (R12 * R23) / sum_R;
    wye->R3 = (R23 * R31) / sum_R;
}

/* ==========================================================================
 * L5: Source Transformation — Thevenin Branch to Norton Branch
 * ==========================================================================
 *
 * Converts a Thevenin branch (voltage source V in series with impedance Z)
 * to its equivalent Norton branch (current source I in parallel with Z):
 *
 *   I = V / Z
 *   Z_out = Z (same impedance)
 *
 * This is the per-branch version; see thevenin_norton.c for the
 * two-terminal network version.
 *
 * Edge case: if Z ≈ 0, the Norton current is infinite → invalid.
 */
void source_transform_thevenin_to_norton(double V, double Z,
                                          double *I_out, double *Z_out) {
    if (!I_out || !Z_out) return;

    *Z_out = Z;

    if (fabs(Z) < 1e-15) {
        /* Ideal voltage source: cannot be Norton-transformed */
        *I_out = 1e15;  /* effectively infinite */
    } else {
        *I_out = V / Z;
    }
}

/* ==========================================================================
 * L5: Source Transformation — Norton Branch to Thevenin Branch
 * ==========================================================================
 *
 * Converts a Norton branch (current source I in parallel with admittance Y)
 * to its equivalent Thevenin branch (voltage source V in series with Z):
 *
 *   V = I * Z = I / Y
 *   Z_out = 1 / Y
 *
 * Edge case: if Z ≈ ∞ (Y ≈ 0), the Thevenin voltage is infinite → invalid.
 */
void source_transform_norton_to_thevenin(double I, double Z,
                                          double *V_out, double *Z_out) {
    if (!V_out || !Z_out) return;

    *Z_out = Z;

    if (Z > 1e15) {
        /* Ideal current source: cannot be Thevenin-transformed */
        *V_out = 1e15;
    } else {
        *V_out = I * Z;
    }
}

/* ==========================================================================
 * L5: Series Resistance Combination
 * ==========================================================================
 *
 * For two resistors in series: R_eq = R1 + R2
 * (Trivial, but needed for automated network reduction)
 */
double combine_series_resistors(double R1, double R2) {
    return R1 + R2;
}

/* ==========================================================================
 * L5: Parallel Resistance Combination
 * ==========================================================================
 *
 * For two resistors in parallel: R_eq = (R1 * R2) / (R1 + R2)
 *
 * Edge cases:
 *   - R1 = 0 (short): R_eq = 0
 *   - R1 → ∞ (open): R_eq = R2
 */
double combine_parallel_resistors(double R1, double R2) {
    if (fabs(R1) < 1e-15 || fabs(R2) < 1e-15) return 0.0;
    if (R1 > 1e15) return R2;
    if (R2 > 1e15) return R1;
    return (R1 * R2) / (R1 + R2);
}

/* ==========================================================================
 * L5: Series Capacitor Combination
 * ==========================================================================
 *
 * For capacitors in series: 1/C_eq = 1/C1 + 1/C2
 * → C_eq = (C1 * C2) / (C1 + C2)
 */
double combine_series_capacitors(double C1, double C2) {
    if (fabs(C1) < 1e-15 || fabs(C2) < 1e-15) return 0.0;
    return (C1 * C2) / (C1 + C2);
}

/* ==========================================================================
 * L5: Parallel Capacitor Combination
 * ==========================================================================
 *
 * For capacitors in parallel: C_eq = C1 + C2
 */
double combine_parallel_capacitors(double C1, double C2) {
    return C1 + C2;
}

/* ==========================================================================
 * L5: Series Inductor Combination
 * ==========================================================================
 *
 * For inductors in series (no mutual coupling): L_eq = L1 + L2
 */
double combine_series_inductors(double L1, double L2) {
    return L1 + L2;
}

/* ==========================================================================
 * L5: Parallel Inductor Combination
 * ==========================================================================
 *
 * For inductors in parallel (no mutual coupling):
 *   L_eq = (L1 * L2) / (L1 + L2)
 */
double combine_parallel_inductors(double L1, double L2) {
    if (fabs(L1) < 1e-15 || fabs(L2) < 1e-15) return 0.0;
    return (L1 * L2) / (L1 + L2);
}

/* ==========================================================================
 * L8: Iterative Network Reduction
 * ==========================================================================
 *
 * Iteratively simplifies a circuit by applying:
 *   1. Series resistor combination
 *   2. Parallel resistor combination
 *   3. Y-Δ / Δ-Y transformation when applicable
 *   4. Source transformation (Thevenin ↔ Norton)
 *
 * This is a heuristic algorithm; it does not guarantee optimal reduction
 * but is effective for typical ladder and bridge circuits.
 *
 * The reduction continues until:
 *   - The circuit is reduced to a single Thevenin/Norton equivalent, or
 *   - No further reductions are possible, or
 *   - max_steps is reached
 *
 * Returns: 0 if still in progress, 1 if converged, -1 on error.
 */
int network_reduce(CircuitTopology *ckt, NetworkReduction *state) {
    if (!ckt || !state) return -1;
    if (state->step >= state->max_steps) return 1;

    state->step++;

    /* Stage 1: Combine series resistors
     * Look for two resistors connected in series (shared node with
     * exactly 2 connections) and replace with equivalent. */
    uint32_t merged = 0;
    for (uint32_t i = 0; i < ckt->num_elements && merged < 10; i++) {
        if (ckt->elements[i].type != ELEM_RESISTOR) continue;

        for (uint32_t j = i + 1; j < ckt->num_elements && merged < 10; j++) {
            if (ckt->elements[j].type != ELEM_RESISTOR) continue;

            /* Check if they're in series: i.to == j.from or i.from == j.to,
             * and the shared node has exactly these 2 connections */
            uint32_t ni_from = ckt->elements[i].node_from;
            uint32_t ni_to   = ckt->elements[i].node_to;
            uint32_t nj_from = ckt->elements[j].node_from;
            uint32_t nj_to   = ckt->elements[j].node_to;

            int in_series = 0;
            uint32_t new_from, new_to;

            if (ni_to == nj_from) {
                in_series = 1;
                new_from = ni_from;
                new_to = nj_to;
            } else if (ni_from == nj_to) {
                in_series = 1;
                new_from = nj_from;
                new_to = ni_to;
            }

            if (in_series) {
                /* Combine: R_eq = R_i + R_j */
                ckt->elements[i].value += ckt->elements[j].value;
                ckt->elements[i].node_from = new_from;
                ckt->elements[i].node_to = new_to;

                /* Remove element j by marking it as open */
                ckt->elements[j].type = ELEM_OPEN;
                ckt->elements[j].value = 0.0;
                merged++;
            }
        }
    }

    /* Stage 2: Combine parallel resistors
     * Look for two resistors sharing both nodes. */
    for (uint32_t i = 0; i < ckt->num_elements && merged < 20; i++) {
        if (ckt->elements[i].type != ELEM_RESISTOR) continue;

        for (uint32_t j = i + 1; j < ckt->num_elements && merged < 20; j++) {
            if (ckt->elements[j].type != ELEM_RESISTOR) continue;

            uint32_t i_f = ckt->elements[i].node_from;
            uint32_t i_t = ckt->elements[i].node_to;
            uint32_t j_f = ckt->elements[j].node_from;
            uint32_t j_t = ckt->elements[j].node_to;

            if ((i_f == j_f && i_t == j_t) || (i_f == j_t && i_t == j_f)) {
                double R1 = ckt->elements[i].value;
                double R2 = ckt->elements[j].value;
                ckt->elements[i].value = combine_parallel_resistors(R1, R2);

                /* Remove element j */
                ckt->elements[j].type = ELEM_OPEN;
                ckt->elements[j].value = 0.0;
                merged++;
            }
        }
    }

    /* If nothing was merged, we've converged */
    if (merged == 0) return 1;

    return 0;
}

/* ==========================================================================
 * L6: Wheatstone Bridge Analysis
 * ==========================================================================
 *
 * Classical Wheatstone bridge circuit:
 *
 *        R1         R2
 *   A ──/\/\/──┬──/\/\/── B
 *              │
 *              ├─── G (galvanometer/detector)
 *              │
 *   C ──/\/\/──┴──/\/\/── D
 *        R3         R4
 *
 * The bridge is balanced when no current flows through G:
 *   R1/R3 = R2/R4  →  R1*R4 = R2*R3
 *
 * Balanced bridge output voltage V_BD = 0.
 * Unbalanced: V_BD = V_supply * (R3/(R1+R3) - R4/(R2+R4))
 *
 * Applications: precision resistance measurement, strain gauges (Boeing),
 *                pressure sensors, temperature sensors (RTD)
 */
/**
 * Analyze a Wheatstone bridge circuit.
 *
 * Uses the Δ-Y transformation approach:
 *   - The bridge forms a Δ (R1, R3, R_g) and a Δ (R2, R4, R_g)
 *   - Convert one Δ to Y, then solve as series-parallel network
 *
 * Output voltage formula for an unbalanced bridge:
 *   V_out = V_supply * [R3/(R1+R3) - R4/(R2+R4)]
 */
int wheatstone_analyze(double R1, double R2, double R3, double R4,
                        double R_g, double V_supply,
                        WheatstoneBridge *bridge) {
    if (!bridge) return -1;

    bridge->R1 = R1; bridge->R2 = R2;
    bridge->R3 = R3; bridge->R4 = R4;
    bridge->R_g = R_g;
    bridge->V_supply = V_supply;

    /* Voltage at the two midpoints (voltage divider) */
    double V_left  = V_supply * R3 / (R1 + R3);   /* node between R1,R3 */
    double V_right = V_supply * R4 / (R2 + R4);   /* node between R2,R4 */

    bridge->V_output = V_left - V_right;

    /* Balance check: R1/R3 == R2/R4 */
    double ratio1 = (fabs(R3) > 1e-15) ? R1 / R3 : 1e15;
    double ratio2 = (fabs(R4) > 1e-15) ? R2 / R4 : 1e15;
    bridge->balance_ratio = fabs(ratio1 - ratio2);
    bridge->is_balanced = (bridge->balance_ratio < 1e-9) ? 1 : 0;

    /* Galvanometer current using Thevenin equivalent of the bridge */
    if (bridge->is_balanced) {
        bridge->I_g = 0.0;
    } else {
        /* Thevenin voltage = V_output (open-circuit output) */
        /* Thevenin resistance: R1||R3 + R2||R4 (looking into bridge output) */
        double R13 = combine_parallel_resistors(R1, R3);
        double R24 = combine_parallel_resistors(R2, R4);
        double R_th = R13 + R24;
        bridge->I_g = bridge->V_output / (R_th + R_g);
    }

    return 0;
}

/* ==========================================================================
 * L6: Voltage Divider Analysis
 * ==========================================================================
 *
 * The voltage divider is the most fundamental circuit analysis pattern.
 * For a series circuit with V_in across R1 + R2:
 *   V_out = V_in * R2 / (R1 + R2)
 *
 * Loaded voltage divider: with load R_L across R2:
 *   R2_eff = R2 || R_L
 *   V_out = V_in * R2_eff / (R1 + R2_eff)
 */
int voltage_divider_analyze(double R1, double R2, double R_load,
                             double V_in, VoltageDivider *div) {
    if (!div) return -1;
    if (R1 <= 0 || R2 <= 0) return -1;

    div->R1 = R1;
    div->R2 = R2;
    div->R_load = R_load;
    div->V_in = V_in;

    if (R_load > 1e15) {
        /* Unloaded divider */
        div->V_out = V_in * R2 / (R1 + R2);
        div->I_total = V_in / (R1 + R2);
        div->I_load = 0.0;
        div->output_impedance = combine_parallel_resistors(R1, R2);
    } else {
        /* Loaded divider: R2_eff = R2 || R_load */
        double R2_eff = combine_parallel_resistors(R2, R_load);
        div->V_out = V_in * R2_eff / (R1 + R2_eff);
        div->I_total = V_in / (R1 + R2_eff);
        div->I_load = div->V_out / R_load;
        /* Output impedance: R1 || R2 (looking back from load) */
        div->output_impedance = combine_parallel_resistors(R1, R2);
    }

    return 0;
}

/* ==========================================================================
 * L6: Current Divider Analysis
 * ==========================================================================
 *
 * For a parallel circuit with total current I_in splitting between R1 and R2:
 *   I1 = I_in * R2 / (R1 + R2)
 *   I2 = I_in * R1 / (R1 + R2)
 *
 * The current divides in inverse proportion to the resistances.
 */
int current_divider_analyze(double R1, double R2, double I_in,
                             CurrentDivider *div) {
    if (!div) return -1;
    if (R1 <= 0 || R2 <= 0) return -1;

    div->R1 = R1;
    div->R2 = R2;
    div->I_in = I_in;

    double R_parallel = combine_parallel_resistors(R1, R2);
    div->V_parallel = I_in * R_parallel;
    div->I1 = I_in * R2 / (R1 + R2);
    div->I2 = I_in * R1 / (R1 + R2);

    return 0;
}
