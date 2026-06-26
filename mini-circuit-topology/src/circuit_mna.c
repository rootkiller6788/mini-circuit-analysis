/**
 * @file circuit_mna.c
 * @brief Modified Nodal Analysis (MNA) — Implementation
 *
 * Implements the MNA formulation used in SPICE-class circuit simulators.
 * Core algorithms: element stamping, matrix assembly, LU decomposition
 * with partial pivoting, forward/back substitution, DC/AC solution.
 *
 * References:
 *   - C.W. Ho et al., "The Modified Nodal Approach to Network Analysis",
 *     IEEE Trans. CAS, vol. 22, no. 6, pp. 504-509, 1975.
 *   - L.W. Nagel, "SPICE2: A Computer Program to Simulate Semiconductor
 *     Circuits", UCB/ERL M520, 1975.
 *   - J. Vlach, K. Singhal, "Computer Methods for Circuit Analysis and
 *     Design", 2nd ed., 1994.
 */

#include "circuit_topology.h"
#include "circuit_mna.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>

/* ==========================================================================
 * L5: MNA Initialization
 * ========================================================================== */

int ct_mna_init(const ct_circuit_t *circuit, ct_mna_system_t *mna)
{
    if (!circuit || !mna) return -1;

    memset(mna, 0, sizeof(ct_mna_system_t));

    int32_t n = circuit->num_nodes - 1;  /* Exclude ground */
    int32_t nv = circuit->num_vsrc;

    mna->num_nodes = n;
    mna->num_vsrc = nv;
    mna->size = n + nv;

    /* Set up node-to-MNA-row mapping */
    /* Ground node (0) is NOT mapped to any MNA row */
    mna->node_map[0] = -1;
    for (int32_t i = 1; i < circuit->num_nodes; i++) {
        int32_t mna_row = i - 1;
        mna->node_map[i] = mna_row;
        mna->row_type[mna_row] = 0;      /* Node voltage unknown */
        mna->row_index[mna_row] = i;     /* Node ID */
    }

    /* Set up voltage-source-to-MNA-row mapping */
    int32_t vsrc_row = n;
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        if (circuit->branches[b].is_voltage_source) {
            mna->vsrc_map[b] = vsrc_row;
            mna->row_type[vsrc_row] = 1;  /* V-source current unknown */
            mna->row_index[vsrc_row] = b; /* Branch ID */
            vsrc_row++;
        } else {
            mna->vsrc_map[b] = -1;  /* Not a voltage source */
        }
    }

    return 0;
}

/* ==========================================================================
 * L5: MNA Element Stamping
 * ========================================================================== */

int ct_mna_stamp_resistor(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                          double resistance)
{
    if (!mna || resistance <= 0.0) return -1;

    double g = 1.0 / resistance;  /* Conductance [Siemens] */
    int32_t rp = mna->node_map[node_p];
    int32_t rn = mna->node_map[node_n];

    if (rp >= 0) {
        mna->G[rp][rp] += g;
        if (rn >= 0) mna->G[rp][rn] -= g;
    }
    if (rn >= 0) {
        mna->G[rn][rn] += g;
        if (rp >= 0) mna->G[rn][rp] -= g;
    }

    return 0;
}

int ct_mna_stamp_capacitor(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                           double capacitance, double omega)
{
    if (!mna || capacitance < 0.0) return -1;

    /* For DC: capacitor is open circuit, no stamp needed */
    if (omega == 0.0) return 0;

    /* For AC: admittance = j * omega * C */
    double complex y = omega * capacitance * I;

    int32_t rp = mna->node_map[node_p];
    int32_t rn = mna->node_map[node_n];

    if (rp >= 0) {
        mna->G[rp][rp] += y;
        if (rn >= 0) mna->G[rp][rn] -= y;
    }
    if (rn >= 0) {
        mna->G[rn][rn] += y;
        if (rp >= 0) mna->G[rn][rp] -= y;
    }

    return 0;
}

int ct_mna_stamp_inductor(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                          double inductance, double omega)
{
    if (!mna || inductance <= 0.0) return -1;

    /* For DC: ideal inductor is a short circuit. This requires special
     * handling (treat as voltage source with V=0). For now, skip DC stamp. */
    if (omega == 0.0) return 0;

    /* For AC: admittance = 1 / (j * omega * L) = -j / (omega * L) */
    double complex y = 1.0 / (omega * inductance * I);

    int32_t rp = mna->node_map[node_p];
    int32_t rn = mna->node_map[node_n];

    if (rp >= 0) {
        mna->G[rp][rp] += y;
        if (rn >= 0) mna->G[rp][rn] -= y;
    }
    if (rn >= 0) {
        mna->G[rn][rn] += y;
        if (rp >= 0) mna->G[rn][rp] -= y;
    }

    return 0;
}

int ct_mna_stamp_vsource(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                         int32_t vsrc_id, double voltage)
{
    if (!mna) return -1;

    int32_t rp = mna->node_map[node_p];
    int32_t rn = mna->node_map[node_n];
    int32_t rv = mna->vsrc_map[vsrc_id];

    if (rv < 0) return -1;  /* Not registered as voltage source */

    /* KVL constraint row: v(rp) - v(rn) = V_source */
    if (rp >= 0) mna->G[rv][rp] = 1.0;
    if (rn >= 0) mna->G[rv][rn] = -1.0;
    mna->rhs[rv] = voltage;

    /* Contribution of i_v to node KCL equations */
    if (rp >= 0) mna->G[rp][rv] = 1.0;
    if (rn >= 0) mna->G[rn][rv] = -1.0;

    return 0;
}

int ct_mna_stamp_isource(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                         double current)
{
    if (!mna) return -1;

    int32_t rp = mna->node_map[node_p];
    int32_t rn = mna->node_map[node_n];

    /* Current enters node_p, leaves node_n:
     * KCL at node_p: +I flows out → subtract from RHS (I enters the node)
     * KCL at node_n: -I flows out → add to RHS
     * Convention: positive current flows from n+ to n- internally.
     * For external source: current entering node_p means RHS[rp] += I */
    if (rp >= 0) mna->rhs[rp] += current;
    if (rn >= 0) mna->rhs[rn] -= current;

    return 0;
}

int ct_mna_stamp_vccs(ct_mna_system_t *mna, int32_t nc_p, int32_t nc_n,
                      int32_t nout_p, int32_t nout_n, double gm)
{
    if (!mna) return -1;

    int32_t rcp = mna->node_map[nc_p];
    int32_t rcn = mna->node_map[nc_n];
    int32_t rop = mna->node_map[nout_p];
    int32_t ron = mna->node_map[nout_n];

    /* i_out = gm * (v(nc_p) - v(nc_n))
     * Stamps into the output node pair rows and control node pair columns */

    if (rop >= 0 && rcp >= 0) mna->G[rop][rcp] += gm;
    if (rop >= 0 && rcn >= 0) mna->G[rop][rcn] -= gm;
    if (ron >= 0 && rcp >= 0) mna->G[ron][rcp] -= gm;
    if (ron >= 0 && rcn >= 0) mna->G[ron][rcn] += gm;

    return 0;
}

int ct_mna_stamp_cccs(ct_mna_system_t *mna, int32_t control_br,
                      int32_t nout_p, int32_t nout_n, double beta)
{
    if (!mna) return -1;

    int32_t rop = mna->node_map[nout_p];
    int32_t ron = mna->node_map[nout_n];
    int32_t rctrl = mna->vsrc_map[control_br];

    if (rctrl < 0) return -1;  /* Control branch must be a voltage source */

    /* i_out = beta * i_control
     * i_control is the current through the voltage source at row rctrl
     * Stamps beta into the output node rows at the i_control column */

    if (rop >= 0) mna->G[rop][rctrl] += beta;
    if (ron >= 0) mna->G[ron][rctrl] -= beta;

    return 0;
}

int ct_mna_stamp_opamp(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                       int32_t node_out)
{
    if (!mna) return -1;

    /* Ideal op-amp: v(node_p) = v(node_n), infinite input Z, zero output Z.
     * Implemented as a nullor:
     *   - Nullator: vp - vn = 0, ip = in = 0 (port 1-2)
     *   - Norator: arbitrary v and i (port 3-ground)
     *
     * MNA stamp adds one extra row/column for the constraint current.
     * Simplified: add constraint equation vp - vn = 0 as a new row,
     * and add the constraint current to the KCL at node_out.
     */

    int32_t rp = mna->node_map[node_p];
    int32_t rn = mna->node_map[node_n];
    int32_t ro = mna->node_map[node_out];

    /* Use the next available row for the constraint (extend MNA matrix) */
    int32_t rc = mna->size;
    if (rc >= CT_MNA_MAX_SIZE) return -1;

    mna->row_type[rc] = 2;         /* Op-amp constraint current */
    mna->row_index[rc] = node_out;

    /* Constraint: v(rp) - v(rn) = 0 */
    if (rp >= 0) mna->G[rc][rp] = 1.0;
    if (rn >= 0) mna->G[rc][rn] = -1.0;
    mna->rhs[rc] = 0.0;

    /* Constraint current flows into output node */
    if (ro >= 0) mna->G[ro][rc] = 1.0;
    if (rp >= 0) mna->G[rp][rc] = -1.0;
    if (rn >= 0) mna->G[rn][rc] = 1.0;

    mna->size = rc + 1;

    return 0;
}

/* ==========================================================================
 * L5: MNA Matrix Assembly
 * ========================================================================== */

int ct_mna_assemble(const ct_circuit_t *circuit, ct_mna_system_t *mna)
{
    if (!circuit || !mna) return -1;

    /* Initialize MNA system */
    if (ct_mna_init(circuit, mna) != 0) return -1;

    double omega = 2.0 * 3.14159265358979323846 * circuit->frequency;

    /* Stamp each branch */
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        const ct_branch_t *br = &circuit->branches[b];
        int32_t np = br->node_from;
        int32_t nn = br->node_to;

        switch (br->elem_type) {
        case CT_ELEM_RESISTOR:
            ct_mna_stamp_resistor(mna, np, nn, br->value);
            break;

        case CT_ELEM_CAPACITOR:
            ct_mna_stamp_capacitor(mna, np, nn, br->value, omega);
            break;

        case CT_ELEM_INDUCTOR:
            ct_mna_stamp_inductor(mna, np, nn, br->value, omega);
            break;

        case CT_ELEM_VSOURCE:
            ct_mna_stamp_vsource(mna, np, nn, b, br->value);
            break;

        case CT_ELEM_ISOURCE:
            ct_mna_stamp_isource(mna, np, nn, br->value);
            break;

        case CT_ELEM_VCCS:
            /* For VCCS, value = gm, value2 encoded as control node info */
            /* Simplified: use fixed control across nodes 0 and a user node */
            ct_mna_stamp_vccs(mna, 0, np, np, nn, br->value);
            break;

        case CT_ELEM_CCCS:
            ct_mna_stamp_cccs(mna, b, np, nn, br->value);
            break;

        case CT_ELEM_OPAMP:
            /* Op-amp: node_from = v+, node_to = v-, value2 encodes vout */
            ct_mna_stamp_opamp(mna, np, nn, (int32_t)br->value2);
            break;

        default:
            /* Nonlinear elements handled in DC solve via iteration */
            break;
        }
    }

    return 0;
}

/* ==========================================================================
 * L5: LU Decomposition with Partial Pivoting (Complex)
 * ========================================================================== */

/**
 * @brief Find pivot row with maximum magnitude in column k.
 */
static int32_t find_pivot(const double complex G[CT_MNA_MAX_SIZE][CT_MNA_MAX_SIZE],
                          int32_t size, int32_t k, const int32_t *pivot)
{
    int32_t max_row = k;
    double max_val = cabs(G[pivot[k]][k]);

    for (int32_t i = k + 1; i < size; i++) {
        double val = cabs(G[pivot[i]][k]);
        if (val > max_val) {
            max_val = val;
            max_row = i;
        }
    }
    return max_row;
}

int ct_mna_lu_factor(ct_mna_system_t *mna, int *singular)
{
    if (!mna) return -1;

    int32_t n = mna->size;
    double complex (*G)[CT_MNA_MAX_SIZE] = mna->G;

    /* Initialize pivot array */
    for (int32_t i = 0; i < n; i++) {
        mna->pivot[i] = i;
    }

    *singular = 0;

    for (int32_t k = 0; k < n; k++) {
        /* Partial pivoting: find row with max |G[p][k]| */
        int32_t p = find_pivot(G, n, k, mna->pivot);

        if (cabs(G[mna->pivot[p]][k]) < DBL_EPSILON * 1e3) {
            *singular = 1;
            return 0;  /* Matrix is singular but not an error per se */
        }

        /* Swap pivot rows */
        if (p != k) {
            int32_t tmp = mna->pivot[k];
            mna->pivot[k] = mna->pivot[p];
            mna->pivot[p] = tmp;
        }

        int32_t pk = mna->pivot[k];
        double complex pivot_val = G[pk][k];

        /* Compute multipliers and update submatrix */
        for (int32_t i = k + 1; i < n; i++) {
            int32_t pi = mna->pivot[i];
            double complex multiplier = G[pi][k] / pivot_val;
            G[pi][k] = multiplier;  /* Store L factor */

            for (int32_t j = k + 1; j < n; j++) {
                G[pi][j] -= multiplier * G[pk][j];
            }
        }
    }

    mna->is_factored = 1;
    return 0;
}

/* ==========================================================================
 * L5: Forward/Back Substitution (Complex)
 * ========================================================================== */

int ct_mna_solve(ct_mna_system_t *mna)
{
    if (!mna || !mna->is_factored) return -1;

    int32_t n = mna->size;
    double complex (*G)[CT_MNA_MAX_SIZE] = mna->G;
    double complex *x = mna->sol;
    double complex *b = mna->rhs;

    /* Step 1: Apply pivot permutation to RHS */
    double complex y[CT_MNA_MAX_SIZE];
    for (int32_t i = 0; i < n; i++) {
        y[i] = b[mna->pivot[i]];
    }

    /* Step 2: Forward substitution (L * z = y) */
    /* L has unit diagonal, stored in lower triangular part */
    for (int32_t i = 1; i < n; i++) {
        for (int32_t j = 0; j < i; j++) {
            y[i] -= G[mna->pivot[i]][j] * y[j];
        }
    }

    /* Step 3: Back substitution (U * x = z) */
    for (int32_t i = n - 1; i >= 0; i--) {
        double complex sum = y[i];
        for (int32_t j = i + 1; j < n; j++) {
            sum -= G[mna->pivot[i]][j] * x[j];
        }
        double complex pivot_val = G[mna->pivot[i]][i];
        if (cabs(pivot_val) < DBL_EPSILON * 1e3) {
            x[i] = 0.0;
        } else {
            x[i] = sum / pivot_val;
        }
    }

    return 0;
}

/* ==========================================================================
 * L6: DC Operating Point
 * ========================================================================== */

int ct_mna_dc_solve(const ct_circuit_t *circuit, ct_dc_result_t *result)
{
    if (!circuit || !result) return -1;

    memset(result, 0, sizeof(ct_dc_result_t));

    ct_mna_system_t mna;

    /* Assemble MNA system */
    if (ct_mna_assemble(circuit, &mna) != 0) return -1;

    /* Factor and solve */
    int singular;
    if (ct_mna_lu_factor(&mna, &singular) != 0) return -1;
    if (singular) {
        result->converged = 0;
        return -1;
    }

    if (ct_mna_solve(&mna) != 0) return -1;

    /* Extract node voltages */
    result->node_voltages[0] = 0.0;  /* Ground */
    for (int32_t i = 1; i < circuit->num_nodes; i++) {
        int32_t row = mna.node_map[i];
        if (row >= 0 && row < mna.num_nodes) {
            result->node_voltages[i] = creal(mna.sol[row]);
        }
    }

    /* Extract branch currents and voltages */
    double bv[CT_MAX_BRANCHES], bi[CT_MAX_BRANCHES];
    ct_mna_extract_branch_results(circuit, &mna, bv, bi);

    for (int32_t b = 0; b < circuit->num_branches; b++) {
        result->branch_voltages[b] = bv[b];
        result->branch_currents[b] = bi[b];
        result->branch_power[b] = bv[b] * bi[b];

        if (result->branch_power[b] > 0) {
            if (circuit->branches[b].is_active) {
                result->total_power_supplied += result->branch_power[b];
            } else {
                result->total_power_dissipated += result->branch_power[b];
            }
        }
    }

    result->converged = 1;
    result->iterations = 1;

    return 0;
}

/* ==========================================================================
 * L6: AC Analysis
 * ========================================================================== */

int ct_mna_ac_solve(const ct_circuit_t *circuit, double freq_hz,
                    ct_ac_result_t *result)
{
    if (!circuit || !result) return -1;

    memset(result, 0, sizeof(ct_ac_result_t));
    result->frequency = freq_hz;

    /* Create a temporary circuit copy with AC domain settings */
    ct_circuit_t ac_circuit;
    memcpy(&ac_circuit, circuit, sizeof(ct_circuit_t));
    ac_circuit.domain = CT_DOMAIN_AC;
    ac_circuit.frequency = freq_hz;

    ct_mna_system_t mna;
    if (ct_mna_assemble(&ac_circuit, &mna) != 0) return -1;

    int singular;
    if (ct_mna_lu_factor(&mna, &singular) != 0) return -1;
    if (singular) return -1;
    if (ct_mna_solve(&mna) != 0) return -1;

    /* Extract complex node voltages */
    result->node_voltages[0] = 0.0;
    for (int32_t i = 1; i < circuit->num_nodes; i++) {
        int32_t row = mna.node_map[i];
        if (row >= 0 && row < mna.num_nodes) {
            result->node_voltages[i] = mna.sol[row];
            result->magnitude_v[i] = cabs(mna.sol[row]);
            result->phase_v[i] = carg(mna.sol[row]) * 180.0 / 3.14159265358979323846;
        }
    }

    /* Compute voltage gain between first two non-ground nodes (simplified) */
    if (circuit->num_nodes >= 3) {
        double complex v_in = result->node_voltages[1];
        double complex v_out = result->node_voltages[2];
        if (cabs(v_in) > DBL_EPSILON) {
            result->voltage_gain = cabs(v_out) / cabs(v_in);
            result->phase_shift = (carg(v_out) - carg(v_in)) * 180.0 / 3.14159265358979323846;
            if (result->phase_shift > 180.0) result->phase_shift -= 360.0;
            if (result->phase_shift < -180.0) result->phase_shift += 360.0;
        }
    }

    return 0;
}

/* ==========================================================================
 * L5: Extract Branch Results from MNA Solution
 * ========================================================================== */

int ct_mna_extract_branch_results(const ct_circuit_t *circuit,
                                  const ct_mna_system_t *mna,
                                  double *branch_voltages,
                                  double *branch_currents)
{
    if (!circuit || !mna || !branch_voltages || !branch_currents)
        return -1;

    /* Reconstruct node voltages from MNA solution */
    double vn[CT_MAX_NODES];
    vn[0] = 0.0;
    for (int32_t i = 1; i < circuit->num_nodes; i++) {
        int32_t row = mna->node_map[i];
        if (row >= 0 && row < mna->num_nodes) {
            vn[i] = creal(mna->sol[row]);
        } else {
            vn[i] = 0.0;
        }
    }

    /* Extract voltage source currents from MNA solution */
    double ivsrc[CT_MAX_BRANCHES];
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        ivsrc[b] = 0.0;
        if (circuit->branches[b].is_voltage_source) {
            int32_t row = mna->vsrc_map[b];
            if (row >= 0 && row < mna->size) {
                ivsrc[b] = creal(mna->sol[row]);
            }
        }
    }

    /* Compute branch voltages and currents */
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        const ct_branch_t *br = &circuit->branches[b];
        double v_br = vn[br->node_from] - vn[br->node_to];
        branch_voltages[b] = v_br;

        /* Compute current based on element type */
        switch (br->elem_type) {
        case CT_ELEM_RESISTOR:
            branch_currents[b] = v_br / br->value;
            break;
        case CT_ELEM_VSOURCE:
            branch_currents[b] = ivsrc[b];
            break;
        case CT_ELEM_ISOURCE:
            branch_currents[b] = br->value;
            break;
        case CT_ELEM_CAPACITOR:
            branch_currents[b] = 0.0;  /* DC: open circuit */
            break;
        case CT_ELEM_INDUCTOR:
            branch_currents[b] = 0.0;  /* DC: short circuit */
            break;
        default:
            branch_currents[b] = 0.0;
            break;
        }
    }

    return 0;
}
