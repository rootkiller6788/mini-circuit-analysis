/**
 * @file thevenin_norton.c
 * @brief Thevenin and Norton Equivalent Circuit Computation
 *
 * Implements the two fundamental equivalent circuit theorems:
 *   - Thevenin's Theorem (Leon Thevenin, 1883): Any linear two-terminal
 *     network can be replaced by V_th in series with Z_th.
 *   - Norton's Theorem (Edward Norton, 1926): Any linear two-terminal
 *     network can be replaced by I_n in parallel with Y_n.
 *
 * Also provides source transformation between Thevenin and Norton forms,
 * and computation of equivalent impedance by source deactivation.
 *
 * Knowledge Coverage:
 *   L4 - Fundamental Laws: Thevenin & Norton theorems with full computation
 *   L2 - Core Concepts: Source transformation, equivalent circuits
 *
 * Reference:
 *   Thevenin, L. (1883) "Sur un nouveau theoreme d'electricite dynamique"
 *   Norton, E.L. (1926) "Design of finite networks..."
 *   Johnson, D.E. et al. "Electric Circuit Analysis" (1997) Ch 5
 */

#include "network_theorem.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L4: Thevenin Equivalent Computation
 * ==========================================================================
 *
 * Algorithm:
 *   1. V_th = open-circuit voltage V_oc at terminal pair (a,b)
 *      - Insert an open circuit between terminals a and b
 *      - Solve the nodal equations for all node voltages
 *      - V_th = V[a] - V[b]
 *
 *   2. Z_th = equivalent impedance seen from terminals (a,b)
 *      - Deactivate all independent sources:
 *        * Voltage sources → short circuit (0 Ω)
 *        * Current sources → open circuit (∞ Ω)
 *      - Apply a test current I_test = 1A at terminals
 *      - Solve for V_test across terminals
 *      - Z_th = V_test / I_test = V_test
 *
 * Handles circuits with:
 *   - Resistors, capacitors, inductors (DC analysis: C→open, L→short)
 *   - Independent voltage and current sources
 *   - Dependent sources (simplified: treated as active elements)
 *
 * Returns: 0 on success, -1 if singular circuit matrix.
 */
int compute_thevenin(const CircuitTopology *ckt,
                     uint32_t terminal_a, uint32_t terminal_b,
                     TheveninEquivalent *result) {
    if (!ckt || !result) return -1;
    uint32_t n = ckt->num_nodes;
    if (n == 0 || terminal_a >= n || terminal_b >= n) return -1;

    /* --- Step 1: Compute Open-Circuit Voltage V_th --- */
    double *Y = build_nodal_admittance_matrix(ckt, &n);
    if (!Y) return -1;

    double *I_vec = (double *)calloc(n, sizeof(double));
    double *V = (double *)calloc(n, sizeof(double));
    if (!I_vec || !V) { free(Y); free(I_vec); free(V); return -1; }

    build_current_vector(ckt, I_vec, n);

    if (solve_nodal_voltages(Y, I_vec, n, V) != 0) {
        free(Y); free(I_vec); free(V); return -1;
    }

    result->V_th.dc_offset = V[terminal_a] - V[terminal_b];
    result->V_th.ac_amplitude = 0.0;
    result->V_th.angular_freq = 0.0;
    result->V_th.phase = 0.0;
    result->frequency = 0.0;

    /* --- Step 2: Compute Equivalent Impedance Z_th --- */
    /* Deactivate independent sources by zeroing them in the Y matrix.
     * Rebuild Y with sources deactivated (already done by build_nodal
     * which skips voltage sources). Current sources are open → no contribution.
     *
     * Apply test current of 1.0A from terminal_a to terminal_b.
     */
    memset(I_vec, 0, n * sizeof(double));
    I_vec[terminal_a] -= 1.0;
    I_vec[terminal_b] += 1.0;

    /* Rebuild Y for deactivated sources (same as before since we skip Vs) */
    double *Y_deact = build_nodal_admittance_matrix(ckt, &n);
    if (!Y_deact) { free(Y); free(I_vec); free(V); return -1; }

    double *V_test = (double *)calloc(n, sizeof(double));
    if (!V_test) {
        free(Y); free(Y_deact); free(I_vec); free(V); return -1;
    }

    int ret = solve_nodal_voltages(Y_deact, I_vec, n, V_test);
    if (ret != 0) {
        free(Y); free(Y_deact); free(I_vec); free(V); free(V_test);
        return -1;
    }

    double V_ab = V_test[terminal_a] - V_test[terminal_b];
    result->Z_th.real = V_ab;      /* Since I_test = 1A: Z = V/I = V */
    result->Z_th.imag = 0.0;

    free(Y);
    free(Y_deact);
    free(I_vec);
    free(V);
    free(V_test);
    return 0;
}

/* ==========================================================================
 * L4: Norton Equivalent Computation
 * ==========================================================================
 *
 * Algorithm:
 *   1. I_n = short-circuit current I_sc at terminal pair
 *      - Connect a short between terminals a and b
 *      - Solve for the current flowing through the short
 *
 *   2. Y_n = 1 / Z_th (or compute directly as equivalent admittance)
 *
 * For a short circuit between a and b, we merge nodes a and b
 * (set V_a = V_b) and solve the reduced system.
 *
 * Simpler approach: compute Thevenin first, then convert.
 */
int compute_norton(const CircuitTopology *ckt,
                   uint32_t terminal_a, uint32_t terminal_b,
                   NortonEquivalent *result) {
    if (!ckt || !result) return -1;

    /* Compute Thevenin equivalent first */
    TheveninEquivalent thev;
    int ret = compute_thevenin(ckt, terminal_a, terminal_b, &thev);
    if (ret != 0) return -1;

    /* Convert Thevenin → Norton */
    thevenin_to_norton(&thev, result);
    return 0;
}

/* ==========================================================================
 * L2: Source Transformation — Thevenin to Norton
 * ==========================================================================
 *
 * Given V_th in series with Z_th:
 *   I_n  = V_th / Z_th
 *   Y_n  = 1 / Z_th
 *
 * The internal impedance/admittance is the same (Z_th = 1/Y_n).
 *
 * Handles edge cases:
 *   - Z_th = 0: I_n is infinite (ideal voltage source cannot be Norton-transformed)
 *   - Z_th → ∞: I_n = 0 (open circuit)
 */
void thevenin_to_norton(const TheveninEquivalent *thev, NortonEquivalent *nort) {
    if (!thev || !nort) return;

    double R = thev->Z_th.real;
    double X = thev->Z_th.imag;

    /* Compute admittance: Y = 1 / Z = 1/(R + jX) = (R - jX) / (R² + X²) */
    double denom = R * R + X * X;
    if (denom < 1e-30) {
        /* Z_th ≈ 0: short circuit, infinite admittance */
        nort->Y_n.real = 1e15;
        nort->Y_n.imag = 0.0;
        nort->I_n.dc_offset = 0.0;
        nort->I_n.ac_amplitude = 0.0;
        nort->I_n.angular_freq = 0.0;
        nort->I_n.phase = 0.0;
        nort->frequency = thev->frequency;
        return;
    }

    nort->Y_n.real =  R / denom;
    nort->Y_n.imag = -X / denom;

    /* I_n = V_th / Z_th = V_th * Y_n
     * For DC: I_dc = V_dc / R
     * For complex: I = V * Y = V * (G + jB) */
    double Vdc = thev->V_th.dc_offset;
    nort->I_n.dc_offset = Vdc * nort->Y_n.real;
    nort->I_n.ac_amplitude = 0.0;
    nort->I_n.angular_freq = 0.0;
    nort->I_n.phase = 0.0;
    nort->frequency = thev->frequency;
}

/* ==========================================================================
 * L2: Source Transformation — Norton to Thevenin
 * ==========================================================================
 *
 * Given I_n in parallel with Y_n:
 *   V_th = I_n / Y_n = I_n * Z_th
 *   Z_th = 1 / Y_n
 */
void norton_to_thevenin(const NortonEquivalent *nort, TheveninEquivalent *thev) {
    if (!nort || !thev) return;

    double G = nort->Y_n.real;
    double B = nort->Y_n.imag;

    /* Z = 1 / Y = 1/(G + jB) = (G - jB) / (G² + B²) */
    double denom = G * G + B * B;
    if (denom < 1e-30) {
        /* Y_n ≈ 0: open circuit, infinite impedance */
        thev->Z_th.real = 1e15;
        thev->Z_th.imag = 0.0;
        thev->V_th.dc_offset = 0.0;
        thev->V_th.ac_amplitude = 0.0;
        thev->V_th.angular_freq = 0.0;
        thev->V_th.phase = 0.0;
        thev->frequency = nort->frequency;
        return;
    }

    thev->Z_th.real =  G / denom;
    thev->Z_th.imag = -B / denom;

    /* V_th = I_n * Z_th = I_n * (R + jX) */
    double Idc = nort->I_n.dc_offset;
    thev->V_th.dc_offset = Idc * thev->Z_th.real;
    thev->V_th.ac_amplitude = 0.0;
    thev->V_th.angular_freq = 0.0;
    thev->V_th.phase = 0.0;
    thev->frequency = nort->frequency;
}

/* ==========================================================================
 * L2: Core Impedance Computation Functions
 * ========================================================================== */

/** Impedance of an ideal resistor: Z = R + j0 */
ComplexImpedance impedance_resistor(double R) {
    ComplexImpedance z;
    z.real = R;
    z.imag = 0.0;
    return z;
}

/** Impedance of an ideal capacitor: Z = 1/(jωC) = 0 - j/(ωC) */
ComplexImpedance impedance_capacitor(double C, double frequency) {
    ComplexImpedance z;
    if (frequency <= 0.0 || C <= 0.0) {
        /* DC or invalid: capacitor is open circuit */
        z.real = 1e15;
        z.imag = 0.0;
        return z;
    }
    double omega = 2.0 * 3.14159265358979323846 * frequency;
    z.real = 0.0;
    z.imag = -1.0 / (omega * C);
    return z;
}

/** Impedance of an ideal inductor: Z = jωL */
ComplexImpedance impedance_inductor(double L, double frequency) {
    ComplexImpedance z;
    double omega = 2.0 * 3.14159265358979323846 * frequency;
    z.real = 0.0;
    z.imag = omega * L;
    return z;
}

/** Series combination: Z_eq = Z1 + Z2 */
ComplexImpedance impedance_series(ComplexImpedance z1, ComplexImpedance z2) {
    ComplexImpedance z;
    z.real = z1.real + z2.real;
    z.imag = z1.imag + z2.imag;
    return z;
}

/**
 * Parallel combination: Z_eq = (Z1 * Z2) / (Z1 + Z2)
 *
 * Implementation uses complex arithmetic:
 *   Z_eq = (Z1 * Z2) / (Z1 + Z2)
 *
 * Let Z1 = R1 + jX1, Z2 = R2 + jX2
 * Numerator:   (R1+jX1)*(R2+jX2) = (R1*R2-X1*X2) + j(R1*X2+R2*X1)
 * Denominator: (R1+R2) + j(X1+X2)
 *
 * Handles edge cases:
 *   - Z1 = 0 (short): Z_eq = 0
 *   - Z1 → ∞ (open): Z_eq = Z2
 *   - Z2 → ∞ (open): Z_eq = Z1
 */
ComplexImpedance impedance_parallel(ComplexImpedance z1, ComplexImpedance z2) {
    ComplexImpedance z;

    /* Check for short circuits */
    double mag1 = z1.real * z1.real + z1.imag * z1.imag;
    double mag2 = z2.real * z2.real + z2.imag * z2.imag;

    /* If either is effectively a short */
    if (mag1 < 1e-30) {
        z.real = 0.0;
        z.imag = 0.0;
        return z;
    }
    if (mag2 < 1e-30) {
        z.real = 0.0;
        z.imag = 0.0;
        return z;
    }

    /* If either is effectively an open */
    if (mag1 > 1e30) { z = z2; return z; }
    if (mag2 > 1e30) { z = z1; return z; }

    /* Z1 + Z2 */
    double sum_real = z1.real + z2.real;
    double sum_imag = z1.imag + z2.imag;

    /* Z1 * Z2 */
    double prod_real = z1.real * z2.real - z1.imag * z2.imag;
    double prod_imag = z1.real * z2.imag + z1.imag * z2.real;

    /* Z_eq = prod / sum = prod * conj(sum) / |sum|^2 */
    double denom = sum_real * sum_real + sum_imag * sum_imag;
    if (denom < 1e-30) {
        z.real = 1e15;
        z.imag = 0.0;
        return z;
    }

    z.real = (prod_real * sum_real + prod_imag * sum_imag) / denom;
    z.imag = (prod_imag * sum_real - prod_real * sum_imag) / denom;
    return z;
}

/** Convert impedance to admittance: Y = 1/Z = (R - jX) / (R² + X²) */
ComplexAdmittance impedance_to_admittance(ComplexImpedance z) {
    ComplexAdmittance y;
    double denom = z.real * z.real + z.imag * z.imag;

    if (denom < 1e-30) {
        /* Z ≈ 0: admittance is infinite (short circuit) */
        y.real = 1e15;
        y.imag = 0.0;
        return y;
    }

    y.real =  z.real / denom;   /* G = R / |Z|² */
    y.imag = -z.imag / denom;   /* B = -X / |Z|² */
    return y;
}

/** Convert admittance to impedance: Z = 1/Y = (G - jB) / (G² + B²) */
ComplexImpedance admittance_to_impedance(ComplexAdmittance y) {
    ComplexImpedance z;
    double denom = y.real * y.real + y.imag * y.imag;

    if (denom < 1e-30) {
        /* Y ≈ 0: impedance is infinite (open circuit) */
        z.real = 1e15;
        z.imag = 0.0;
        return z;
    }

    z.real =  y.real / denom;   /* R = G / |Y|² */
    z.imag = -y.imag / denom;   /* X = -B / |Y|² */
    return z;
}

/* ==========================================================================
 * L2: Impedance Utility Functions
 * ========================================================================== */

/** Magnitude: |Z| = sqrt(R² + X²) */
double impedance_magnitude(ComplexImpedance z) {
    return sqrt(z.real * z.real + z.imag * z.imag);
}

/** Phase angle: φ = atan2(X, R) in radians */
double impedance_phase(ComplexImpedance z) {
    return atan2(z.imag, z.real);
}

/** Power dissipated in resistance: P = V²/R = I²*R = V*I (W) */
double power_dissipated(double voltage, double current, double resistance) {
    /* Three equivalent formulations; use the most numerically stable */
    if (fabs(resistance) > 1e-15) {
        return voltage * voltage / resistance;
    }
    return voltage * current;
}

/** Reactive power: Q = I² * X (VAR) */
double reactive_power(double current, double reactance) {
    return current * current * reactance;
}

/** Apparent power: S = V_rms * I_rms (VA) */
double apparent_power(double v_rms, double i_rms) {
    return v_rms * i_rms;
}

/** Power factor: pf = cos(θ_v - θ_i), range [-1, 1] */
double power_factor(double v_phase, double i_phase) {
    return cos(v_phase - i_phase);
}

/* ==========================================================================
 * L1: Topology Memory Management
 * ========================================================================== */

CircuitTopology *topology_create(uint32_t n_nodes, uint32_t n_branches,
                                  uint32_t n_elements, uint32_t n_meshes) {
    CircuitTopology *ckt = (CircuitTopology *)malloc(sizeof(CircuitTopology));
    if (!ckt) return NULL;

    ckt->num_nodes    = n_nodes;
    ckt->num_branches = n_branches;
    ckt->num_elements = n_elements;
    ckt->num_meshes   = n_meshes;
    ckt->ground_node  = 0;

    ckt->nodes    = n_nodes    ? (CircuitNode *)calloc(n_nodes, sizeof(CircuitNode)) : NULL;
    ckt->branches = n_branches ? (CircuitBranch *)calloc(n_branches, sizeof(CircuitBranch)) : NULL;
    ckt->elements = n_elements ? (CircuitElement *)calloc(n_elements, sizeof(CircuitElement)) : NULL;
    ckt->meshes   = n_meshes   ? (CircuitMesh *)calloc(n_meshes, sizeof(CircuitMesh)) : NULL;

    /* Initialize node IDs */
    for (uint32_t i = 0; i < n_nodes; i++) {
        ckt->nodes[i].id = i;
        ckt->nodes[i].voltage = 0.0;
        ckt->nodes[i].is_ground = (i == 0) ? 1 : 0;
    }

    return ckt;
}

void topology_free(CircuitTopology *ckt) {
    if (!ckt) return;
    /* Free mesh branch arrays */
    if (ckt->meshes) {
        for (uint32_t i = 0; i < ckt->num_meshes; i++) {
            free(ckt->meshes[i].branch_ids);
        }
    }
    free(ckt->nodes);
    free(ckt->branches);
    free(ckt->elements);
    free(ckt->meshes);
    free(ckt);
}
