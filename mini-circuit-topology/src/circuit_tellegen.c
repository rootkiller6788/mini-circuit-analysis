/**
 * @file circuit_tellegen.c
 * @brief Tellegen's Theorem and Kirchhoff's Laws Verification
 *
 * Implements verification of KCL, KVL, and Tellegen's Theorem for
 * circuit topology validation. These are fundamental consistency checks
 * used in circuit simulation to detect numerical errors.
 *
 * Tellegen's Theorem (B.D.H. Tellegen, 1952):
 *   "For any lumped electrical network, for any sets of branch voltages
 *    {v_k} satisfying KVL and branch currents {i_k} satisfying KCL,
 *    the sum of products is identically zero: sum(v_k * i_k) = 0."
 *
 * This theorem is topological in nature — it depends only on the
 * interconnection pattern (incidence matrix), not on the branch
 * constitutive relations. It is a direct consequence of the
 * orthogonality of the cut-set space and the loop space of the
 * circuit graph.
 *
 * References:
 *   - B.D.H. Tellegen, "A General Network Theorem, with Applications",
 *     Philips Research Reports, vol. 7, pp. 259-269, 1952.
 *   - P. Penfield, R. Spence, S. Duinker, "Tellegen's Theorem and
 *     Electrical Networks", MIT Press, 1970.
 *   - ETH 227-0455 / TU Munich High-Frequency Engineering
 *
 * Knowledge coverage:
 *   L4 (Fundamental Laws): Tellegen's Theorem, KCL, KVL
 */

#include "circuit_topology.h"
#include "circuit_graph.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * L4: KCL Verification
 * ========================================================================== */

int ct_verify_kcl(const ct_circuit_t *circuit, const double *branch_currents,
                  double tolerance, int32_t *violated_node)
{
    if (!circuit || !branch_currents) return -1;
    if (violated_node) *violated_node = -1;

    /* Build incidence matrix to compute KCL */
    ct_incidence_t A;
    ct_build_incidence_matrix(circuit, &A);

    /* For each non-ground node, compute sum of currents */
    for (int32_t row = 0; row < A.rows; row++) {
        double sum = 0.0;
        for (int32_t col = 0; col < A.cols; col++) {
            /* A[row][col] = +1 means current leaves node
             * KCL: sum of currents LEAVING a node = 0
             * Current in branch = courant flowing FROM node_from TO node_to
             * If A=+1, branch leaves node → current flows away → subtract from sum
             * If A=-1, branch enters node → current flows toward → add to sum */
            sum += A.data[row][col] * branch_currents[col];
        }

        /* For KCL, sum of currents entering = sum leaving,
         * so net sum (with sign convention) should be zero */
        if (fabs(sum) > tolerance) {
            if (violated_node) *violated_node = row + 1;  /* Node ID = row + 1 */
            return -1;
        }
    }

    return 0;
}

/* ==========================================================================
 * L4: KVL Verification
 * ========================================================================== */

int ct_verify_kvl(const ct_circuit_t *circuit, const double *branch_voltages,
                  double tolerance, int32_t *violated_loop)
{
    if (!circuit || !branch_voltages) return -1;
    if (violated_loop) *violated_loop = -1;

    /* Select a spanning tree and build the fundamental loop matrix */
    ct_tree_t tree;
    if (ct_select_tree(circuit, &tree, 0) != 0) {
        /* Disconnected circuit — KVL is trivially satisfied */
        return 0;
    }

    ct_loop_matrix_t B;
    ct_build_loop_matrix(circuit, &tree, &B);

    /* For each fundamental loop, compute sum of voltages */
    for (int32_t row = 0; row < B.rows; row++) {
        double sum = 0.0;
        for (int32_t col = 0; col < B.cols; col++) {
            sum += B.data[row][col] * branch_voltages[col];
        }

        if (fabs(sum) > tolerance) {
            if (violated_loop) *violated_loop = row;
            return -1;
        }
    }

    return 0;
}

/* ==========================================================================
 * L4: Tellegen's Theorem Verification
 * ========================================================================== */

int ct_verify_tellegen(const ct_circuit_t *circuit,
                       const double *branch_voltages,
                       const double *branch_currents,
                       double tolerance)
{
    if (!circuit || !branch_voltages || !branch_currents) return -1;

    int32_t b = circuit->num_branches;
    double total_power = 0.0;

    /* Compute sum(v_k * i_k) over all branches */
    for (int32_t k = 0; k < b; k++) {
        total_power += branch_voltages[k] * branch_currents[k];
    }

    /* Tellegen's theorem: this sum must be zero (within tolerance) */
    if (fabs(total_power) > tolerance) {
        return -1;
    }

    return 0;
}

/* ==========================================================================
 * L4: Tellegen's Theorem — Cross-Circuit Verification
 * ========================================================================== */

/**
 * @brief Verify Tellegen's theorem in its strong form: cross-product
 * between voltages from one circuit and currents from another circuit
 * with the SAME topology.
 *
 * The strong form of Tellegen's theorem states:
 *   sum(v'_k * i''_k) = 0
 *   sum(v''_k * i'_k) = 0
 *
 * where (v', i') satisfy KVL/KCL for circuit A, and
 * (v'', i'') satisfy KVL/KCL for circuit B with identical topology.
 *
 * This is the basis of the adjoint network method for sensitivity analysis
 * and is used in SPICE for computing noise contributions.
 *
 * @param circuit         Circuit topology (common to both)
 * @param v1, i1          Branch variables from first operating condition
 * @param v2, i2          Branch variables from second operating condition
 * @param sum_v1i2        Output: sum(v1_k * i2_k)
 * @param sum_v2i1        Output: sum(v2_k * i1_k)
 * @param tolerance       Numerical tolerance
 * @return 0 if both cross-sums are zero, -1 otherwise
 */
int ct_verify_tellegen_cross(const ct_circuit_t *circuit,
                             const double *v1, const double *i1,
                             const double *v2, const double *i2,
                             double *sum_v1i2, double *sum_v2i1,
                             double tolerance)
{
    if (!circuit || !v1 || !i1 || !v2 || !i2 || !sum_v1i2 || !sum_v2i1)
        return -1;

    int32_t b = circuit->num_branches;
    double sv1i2 = 0.0;
    double sv2i1 = 0.0;

    for (int32_t k = 0; k < b; k++) {
        sv1i2 += v1[k] * i2[k];
        sv2i1 += v2[k] * i1[k];
    }

    *sum_v1i2 = sv1i2;
    *sum_v2i1 = sv2i1;

    if (fabs(sv1i2) > tolerance || fabs(sv2i1) > tolerance) {
        return -1;
    }

    return 0;
}

/* ==========================================================================
 * L4: Power Conservation Check (Energy Balance)
 * ========================================================================== */

/**
 * @brief Verify power conservation in the circuit.
 *
 * For any circuit: Total Power Supplied = Total Power Dissipated.
 * This is a consequence of Tellegen's theorem applied to a single set
 * of branch variables:
 *   sum(V_s * I_s) = sum(R_k * I_k^2) + sum(G_k * V_k^2)
 *
 * where the left sum is over sources and the right sum is over
 * dissipative elements.
 *
 * @param circuit     Circuit topology
 * @param v           Branch voltages
 * @param i           Branch currents
 * @param p_supplied  Output: total power supplied by sources [W]
 * @param p_dissipated Output: total power dissipated [W]
 * @param tolerance   Numerical tolerance for equality check
 * @return 0 if power balances, -1 if imbalance exceeds tolerance
 */
int ct_verify_power_balance(const ct_circuit_t *circuit,
                            const double *v, const double *i,
                            double *p_supplied, double *p_dissipated,
                            double tolerance)
{
    if (!circuit || !v || !i || !p_supplied || !p_dissipated)
        return -1;

    double ps = 0.0;  /* Power supplied */
    double pd = 0.0;  /* Power dissipated */

    for (int32_t k = 0; k < circuit->num_branches; k++) {
        double p = v[k] * i[k];  /* Power absorbed by branch k */

        if (circuit->branches[k].is_active) {
            /* Active element: negative absorbed power = supplied */
            ps += -p;  /* If p < 0, source is delivering power */
        } else {
            /* Passive element: always dissipates (p >= 0 for resistors) */
            if (p > 0) {
                pd += p;
            }
        }
    }

    /* For active elements that actually consume power (e.g., op-amp biasing),
     * adjust the balance. The fundamental check is that net power = 0. */
    double net_power = 0.0;
    for (int32_t k = 0; k < circuit->num_branches; k++) {
        net_power += v[k] * i[k];
    }

    *p_supplied = ps;
    *p_dissipated = pd;

    if (fabs(net_power) > tolerance) {
        return -1;
    }

    return 0;
}

/* ==========================================================================
 * L4: Sensitivity via Adjoint Network (Tellegen Application)
 * ========================================================================== */

/**
 * @brief Compute the sensitivity of a node voltage to a branch parameter
 * using the adjoint network method.
 *
 * The adjoint network method (Director & Rohrer, 1969) uses Tellegen's
 * theorem to compute ALL sensitivities with just TWO circuit analyses
 * (original + adjoint), regardless of the number of parameters.
 *
 * For a resistor R_k between nodes a and b:
 *   dV_out / dR_k = -I_k * I_adj_k
 *
 * where I_k is the current in the original circuit and I_adj_k is the
 * current in the adjoint circuit (which has the same topology but
 * different excitations).
 *
 * This is the foundation of SPICE's .SENS analysis and is critical
 * for circuit optimization and yield analysis.
 *
 * Reference:
 *   - S.W. Director, R.A. Rohrer, "Automated Network Design — The
 *     Frequency-Domain Case", IEEE Trans. CT, vol. 16, pp. 330-337, 1969.
 *
 * @param circuit        Circuit topology
 * @param original_v     Branch voltages from original analysis
 * @param original_i     Branch currents from original analysis
 * @param adjoint_v      Branch voltages from adjoint analysis
 * @param adjoint_i      Branch currents from adjoint analysis
 * @param sensitivities  Output: d(V_out)/d(parameter_k) for each branch
 * @param num_branches   Number of branches to compute sensitivity for
 * @return 0 on success, -1 on error
 */
int ct_adjoint_sensitivity(const ct_circuit_t *circuit,
                           const double *original_v, const double *original_i,
                           const double *adjoint_v, const double *adjoint_i,
                           double *sensitivities, int32_t num_branches)
{
    if (!circuit || !original_v || !original_i ||
        !adjoint_v || !adjoint_i || !sensitivities)
        return -1;

    int32_t b = circuit->num_branches;
    if (num_branches > b) num_branches = b;

    for (int32_t k = 0; k < num_branches; k++) {
        const ct_branch_t *br = &circuit->branches[k];

        switch (br->elem_type) {
        case CT_ELEM_RESISTOR:
            /* dV_out/dR_k = -i_k * i_adj_k (Director-Rohrer formula) */
            sensitivities[k] = -original_i[k] * adjoint_i[k];
            break;

        case CT_ELEM_CAPACITOR:
            /* dV_out/dC_k = j*omega * v_k * v_adj_k (for AC) */
            sensitivities[k] = original_v[k] * adjoint_v[k];
            break;

        case CT_ELEM_VSOURCE:
            /* dV_out/dV_k = -i_adj_k (current through source in adjoint) */
            sensitivities[k] = -adjoint_i[k];
            break;

        default:
            sensitivities[k] = 0.0;
            break;
        }
    }

    return 0;
}

/* ==========================================================================
 * L4: Reciprocity Check
 * ========================================================================== */

/**
 * @brief Verify reciprocity for a two-port network.
 *
 * A network is reciprocal if the ratio of response to excitation is
 * unchanged when the positions of excitation and response are interchanged.
 *
 * For impedance parameters: Z_12 = Z_21
 * For admittance parameters: Y_12 = Y_21
 *
 * Reciprocity holds for all networks composed exclusively of linear,
 * time-invariant, passive, bilateral elements (R, L, C, transformers).
 * It does NOT hold for networks containing gyrators, circulators, or
 * active devices (transistors, op-amps).
 *
 * Lorentz Reciprocity Theorem (EM field analog):
 *   integral(E1 · J2) = integral(E2 · J1)
 *
 * @param circuit   Circuit topology
 * @param z12       Transfer impedance port1→port2 [Ohm]
 * @param z21       Transfer impedance port2→port1 [Ohm]
 * @param tolerance Numerical tolerance
 * @return 1 if reciprocal, 0 if not, -1 on error
 */
int ct_check_reciprocity(const ct_circuit_t *circuit,
                         double z12, double z21, double tolerance)
{
    if (!circuit) return -1;

    /* Check if circuit contains non-reciprocal elements */
    for (int32_t k = 0; k < circuit->num_branches; k++) {
        if (!circuit->branches[k].is_reciprocal) {
            /* Non-reciprocal element found — reciprocity may not hold */
            /* But still check if Z12 ≈ Z21 numerically */
        }
    }

    if (fabs(z12 - z21) < tolerance) {
        return 1;  /* Reciprocal */
    }

    return 0;  /* Non-reciprocal */
}
