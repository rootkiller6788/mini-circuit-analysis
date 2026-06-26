/**
 * @file circuit_analysis.c
 * @brief Classical circuit analysis methods — Implementation
 *
 * Implements node-voltage method, mesh-current method, superposition
 * principle, Thevenin/Norton equivalents, maximum power transfer,
 * and Y-Delta transformations.
 *
 * References:
 *   - W.H. Hayt, J.E. Kemmerly, "Engineering Circuit Analysis", 9th ed.
 *   - C.A. Desoer, E.S. Kuh, "Basic Circuit Theory" (1969)
 *   - MIT 6.002 / Berkeley EE16A / Michigan EECS 215
 */

#include "circuit_topology.h"
#include "circuit_analysis.h"
#include "circuit_mna.h"
#include "circuit_graph.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* ==========================================================================
 * Internal: Dense Gaussian Elimination Solver (for small matrices)
 * ========================================================================== */

/**
 * @brief Solve Ax = b using Gaussian elimination with partial pivoting.
 *
 * Complexity: O(n^3). Suitable for matrices up to ~100x100.
 * For larger matrices, use LU decomposition from circuit_mna.c.
 *
 * @param A  Matrix (n x n, row-major, modified in-place)
 * @param b  RHS vector (n, modified in-place → becomes solution x)
 * @param n  Matrix dimension
 * @return 0 on success, -1 if singular
 */
static int gauss_solve(double *A, double *b, int n)
{
    if (n <= 0) return -1;
    if (n == 1) {
        if (fabs(A[0]) < DBL_EPSILON) return -1;
        b[0] /= A[0];
        return 0;
    }

    /* Forward elimination with partial pivoting */
    for (int k = 0; k < n; k++) {
        /* Find pivot */
        int max_row = k;
        double max_val = fabs(A[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double val = fabs(A[i * n + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        if (max_val < DBL_EPSILON * 1e3) return -1;  /* Singular */

        /* Swap rows */
        if (max_row != k) {
            for (int j = k; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[max_row * n + j];
                A[max_row * n + j] = tmp;
            }
            double tmp = b[k];
            b[k] = b[max_row];
            b[max_row] = tmp;
        }

        /* Eliminate */
        double pivot = A[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double factor = A[i * n + k] / pivot;
            for (int j = k; j < n; j++) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    /* Back substitution */
    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum / A[i * n + i];
    }

    return 0;
}

/* ==========================================================================
 * L2: Node-Voltage Method
 * ========================================================================== */

int ct_build_nodal_admittance(const ct_circuit_t *circuit,
                              double *Y_n, int32_t n)
{
    if (!circuit || !Y_n || n <= 0) return -1;

    /* Initialize Y_n to zeros */
    for (int i = 0; i < n * n; i++) Y_n[i] = 0.0;

    /* For each branch, stamp its admittance into Y_n */
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        const ct_branch_t *br = &circuit->branches[b];
        double y = 0.0;  /* Admittance */

        switch (br->elem_type) {
        case CT_ELEM_RESISTOR:
            if (br->value > 0.0) y = 1.0 / br->value;
            break;
        case CT_ELEM_VSOURCE:
            /* Voltage source: not directly stampable in nodal Y.
             * This method is limited to circuits without floating V-sources. */
            continue;
        case CT_ELEM_ISOURCE:
            /* Current source: only affects RHS, not Y_n */
            continue;
        default:
            continue;
        }

        int32_t rp = br->node_from - 1;  /* Node index (0-based, exclude ground) */
        int32_t rn = br->node_to - 1;

        if (rp >= 0) {
            Y_n[rp * n + rp] += y;
            if (rn >= 0) Y_n[rp * n + rn] -= y;
        }
        if (rn >= 0) {
            Y_n[rn * n + rn] += y;
            if (rp >= 0) Y_n[rn * n + rp] -= y;
        }
    }

    return 0;
}

int ct_node_voltage_solve(const ct_circuit_t *circuit,
                          double *node_voltages,
                          double *branch_currents,
                          double *branch_voltages)
{
    if (!circuit || !node_voltages || !branch_currents || !branch_voltages)
        return -1;

    int n = circuit->num_nodes - 1;  /* Number of independent nodes */
    if (n <= 0) return -1;

    /* Allocate Y_n and i_n */
    double *Y_n = (double *)calloc(n * n, sizeof(double));
    double *i_n = (double *)calloc(n, sizeof(double));
    if (!Y_n || !i_n) {
        free(Y_n);
        free(i_n);
        return -1;
    }

    /* Build nodal admittance matrix */
    ct_build_nodal_admittance(circuit, Y_n, n);

    /* Build current source vector */
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        const ct_branch_t *br = &circuit->branches[b];
        if (br->elem_type == CT_ELEM_ISOURCE) {
            int32_t rp = br->node_from - 1;
            int32_t rn = br->node_to - 1;
            if (rp >= 0) i_n[rp] += br->value;
            if (rn >= 0) i_n[rn] -= br->value;
        }
    }

    /* Solve Y_n * v_n = i_n */
    if (gauss_solve(Y_n, i_n, n) != 0) {
        free(Y_n);
        free(i_n);
        return -1;
    }

    /* Extract node voltages */
    node_voltages[0] = 0.0;
    for (int i = 0; i < n; i++) {
        node_voltages[i + 1] = i_n[i];  /* Solution is in i_n after solve */
    }

    /* Compute branch voltages and currents */
    for (int32_t b = 0; b < circuit->num_branches; b++) {
        const ct_branch_t *br = &circuit->branches[b];
        double v_br = node_voltages[br->node_from] - node_voltages[br->node_to];
        branch_voltages[b] = v_br;

        switch (br->elem_type) {
        case CT_ELEM_RESISTOR:
            branch_currents[b] = v_br / br->value;
            break;
        case CT_ELEM_ISOURCE:
            branch_currents[b] = br->value;
            break;
        default:
            branch_currents[b] = 0.0;
            break;
        }
    }

    free(Y_n);
    free(i_n);
    return 0;
}

/* ==========================================================================
 * L2: Mesh-Current Method
 * ========================================================================== */

int ct_identify_meshes(const ct_circuit_t *circuit, int32_t max_meshes,
                       int8_t *mesh_branches, int32_t *num_meshes)
{
    if (!circuit || !mesh_branches || !num_meshes) return -1;

    int32_t b = circuit->num_branches;
    int32_t n = circuit->num_nodes;

    /* Number of independent meshes for a connected planar graph */
    int32_t expected_meshes = b - n + 1;
    if (expected_meshes <= 0 || expected_meshes > max_meshes) {
        *num_meshes = 0;
        return -1;
    }

    memset(mesh_branches, 0, max_meshes * CT_MAX_BRANCHES * sizeof(int8_t));

    /* Simplified mesh identification using the cycle basis approach.
     * Each fundamental cycle from a DFS tree corresponds to a mesh
     * in a planar embedding. */
    ct_tree_t tree;
    if (ct_select_tree(circuit, &tree, 0) != 0) {
        *num_meshes = 0;
        return -1;
    }

    *num_meshes = tree.num_links;
    if (*num_meshes > max_meshes) *num_meshes = max_meshes;

    for (int32_t l = 0; l < *num_meshes; l++) {
        int32_t link_br = tree.link_branches[l];
        int32_t cycle[CT_MAX_BRANCHES];
        int32_t cycle_len;

        ct_fundamental_cycle(circuit, &tree, link_br, cycle, &cycle_len);

        for (int32_t c = 0; c < cycle_len; c++) {
            int32_t br = cycle[c];
            /* Determine orientation: +1 if branch direction aligns with mesh
             * (simplified: all +1 for initial implementation) */
            mesh_branches[l * CT_MAX_BRANCHES + br] = 1;
        }
    }

    return 0;
}

int ct_build_mesh_impedance(const ct_circuit_t *circuit,
                            double *Z_m, int32_t m)
{
    if (!circuit || !Z_m || m <= 0) return -1;

    for (int i = 0; i < m * m; i++) Z_m[i] = 0.0;

    int8_t mesh_branches_c[CT_MAX_BRANCHES * CT_MAX_BRANCHES];
    int32_t num_meshes;

    if (ct_identify_meshes(circuit, CT_MAX_BRANCHES,
                           mesh_branches_c, &num_meshes) != 0) {
        return -1;
    }

    if (num_meshes != m) return -1;

    /* Build Z_m: self-impedance on diagonal, mutual off-diagonal */
    for (int32_t mi = 0; mi < m; mi++) {
        for (int32_t mj = 0; mj < m; mj++) {
            double z_sum = 0.0;
            for (int32_t b = 0; b < circuit->num_branches; b++) {
                int8_t oi = mesh_branches_c[mi * CT_MAX_BRANCHES + b];
                int8_t oj = mesh_branches_c[mj * CT_MAX_BRANCHES + b];
                if (oi != 0 && oj != 0) {
                    const ct_branch_t *br = &circuit->branches[b];
                    double sign = (oi == oj) ? 1.0 : -1.0;
                    if (br->elem_type == CT_ELEM_RESISTOR) {
                        z_sum += sign * br->value;
                    }
                }
            }
            Z_m[mi * m + mj] = z_sum;
        }
    }

    return 0;
}

int ct_mesh_current_solve(const ct_circuit_t *circuit,
                          double *mesh_currents,
                          double *branch_currents,
                          double *branch_voltages,
                          int32_t *num_meshes)
{
    if (!circuit || !mesh_currents || !branch_currents ||
        !branch_voltages || !num_meshes) return -1;

    int32_t b = circuit->num_branches;
    int32_t n = circuit->num_nodes;
    int32_t m = b - n + 1;  /* Number of independent meshes */

    if (m <= 0 || m > CT_MAX_BRANCHES) return -1;

    double *Z_m = (double *)calloc(m * m, sizeof(double));
    double *v_m = (double *)calloc(m, sizeof(double));
    if (!Z_m || !v_m) {
        free(Z_m);
        free(v_m);
        return -1;
    }

    /* Build mesh impedance matrix */
    if (ct_build_mesh_impedance(circuit, Z_m, m) != 0) {
        free(Z_m);
        free(v_m);
        return -1;
    }

    /* Build mesh voltage source vector */
    int8_t mesh_branches_c[CT_MAX_BRANCHES * CT_MAX_BRANCHES];
    ct_identify_meshes(circuit, CT_MAX_BRANCHES, mesh_branches_c, num_meshes);

    for (int32_t mi = 0; mi < m; mi++) {
        v_m[mi] = 0.0;
        for (int32_t br = 0; br < b; br++) {
            const ct_branch_t *branch = &circuit->branches[br];
            if (branch->elem_type == CT_ELEM_VSOURCE) {
                int8_t orient = mesh_branches_c[mi * CT_MAX_BRANCHES + br];
                if (orient != 0) {
                    /* Voltage rise in mesh direction = +source value */
                    v_m[mi] += orient * branch->value;
                }
            }
        }
    }

    /* Solve Z_m * i_m = v_m */
    if (gauss_solve(Z_m, v_m, m) != 0) {
        free(Z_m);
        free(v_m);
        return -1;
    }

    /* Extract mesh currents */
    for (int32_t i = 0; i < m; i++) {
        mesh_currents[i] = v_m[i];  /* Solution in v_m after solve */
    }

    /* Compute branch currents from mesh currents */
    for (int32_t br = 0; br < b; br++) {
        branch_currents[br] = 0.0;
        for (int32_t mi = 0; mi < m; mi++) {
            int8_t orient = mesh_branches_c[mi * CT_MAX_BRANCHES + br];
            if (orient != 0) {
                branch_currents[br] += orient * mesh_currents[mi];
            }
        }
    }

    /* Compute branch voltages from currents */
    for (int32_t br = 0; br < b; br++) {
        const ct_branch_t *branch = &circuit->branches[br];
        switch (branch->elem_type) {
        case CT_ELEM_RESISTOR:
            branch_voltages[br] = branch_currents[br] * branch->value;
            break;
        case CT_ELEM_VSOURCE:
            branch_voltages[br] = branch->value;
            break;
        default:
            branch_voltages[br] = 0.0;
            break;
        }
    }

    *num_meshes = m;
    free(Z_m);
    free(v_m);
    return 0;
}

/* ==========================================================================
 * L2: Superposition Principle
 * ========================================================================== */

int ct_superposition_solve(const ct_circuit_t *circuit,
                           double *node_voltages,
                           double *branch_currents,
                           int32_t *num_sources_used)
{
    if (!circuit || !node_voltages || !branch_currents || !num_sources_used)
        return -1;

    int32_t n = circuit->num_nodes;
    int32_t b = circuit->num_branches;

    /* Initialize outputs to zero */
    for (int32_t i = 0; i < n; i++) node_voltages[i] = 0.0;
    for (int32_t i = 0; i < b; i++) branch_currents[i] = 0.0;

    int32_t src_count = 0;

    /* Process each independent source one at a time */
    for (int32_t src = 0; src < b; src++) {
        const ct_branch_t *br = &circuit->branches[src];

        if (br->elem_type != CT_ELEM_VSOURCE &&
            br->elem_type != CT_ELEM_ISOURCE) {
            continue;
        }

        /* Create reduced circuit with all other sources deactivated */
        ct_circuit_t reduced;
        ct_circuit_init(&reduced, "Superposition sub-circuit");

        /* Copy all nodes */
        for (int32_t i = 1; i < n; i++) {
            uint8_t flags = 0;
            if (circuit->nodes[i].is_terminal) flags |= 0x01;
            if (circuit->nodes[i].is_internal) flags |= 0x02;
            ct_add_node(&reduced, circuit->nodes[i].name, flags);
        }

        /* Copy branches with source deactivation */
        for (int32_t i = 0; i < b; i++) {
            const ct_branch_t *obr = &circuit->branches[i];
            if (i == src) {
                /* Keep this source active */
                ct_add_branch(&reduced, obr->elem_type,
                              obr->node_from, obr->node_to,
                              obr->value, obr->value2, obr->name);
            } else if (obr->elem_type == CT_ELEM_VSOURCE) {
                /* Deactivate voltage source → short circuit (0V source) */
                ct_add_branch(&reduced, CT_ELEM_VSOURCE,
                              obr->node_from, obr->node_to,
                              0.0, obr->value2, obr->name);
            } else if (obr->elem_type == CT_ELEM_ISOURCE) {
                /* Deactivate current source → open circuit (skip it) */
                continue;
            } else {
                /* Passive elements remain unchanged */
                ct_add_branch(&reduced, obr->elem_type,
                              obr->node_from, obr->node_to,
                              obr->value, obr->value2, obr->name);
            }
        }

        /* Solve the reduced circuit */
        double partial_v[CT_MAX_NODES];
        double partial_i[CT_MAX_BRANCHES];
        double partial_bv[CT_MAX_BRANCHES];

        if (ct_node_voltage_solve(&reduced, partial_v, partial_i, partial_bv) == 0) {
            /* Accumulate results */
            for (int32_t i = 0; i < n; i++) {
                node_voltages[i] += partial_v[i];
            }
            for (int32_t i = 0; i < b; i++) {
                branch_currents[i] += partial_i[i];
            }
            src_count++;
        }
    }

    *num_sources_used = src_count;
    return 0;
}

/* ==========================================================================
 * L6: Thevenin and Norton Equivalents
 * ========================================================================== */

int ct_thevenin_equivalent(const ct_circuit_t *circuit,
                           int32_t node_a, int32_t node_b,
                           double *V_th, double *R_th)
{
    if (!circuit || !V_th || !R_th) return -1;

    /* Step 1: Compute V_th = open-circuit voltage between node_a and node_b */
    ct_dc_result_t dc_result;
    if (ct_mna_dc_solve(circuit, &dc_result) != 0) return -1;

    *V_th = dc_result.node_voltages[node_a] - dc_result.node_voltages[node_b];

    /* Step 2: Compute R_th = V_th / I_sc
     * Create a modified circuit with a short between node_a and node_b (0V source)
     * and measure the current. */
    ct_circuit_t sc_circuit;
    memcpy(&sc_circuit, circuit, sizeof(ct_circuit_t));

    /* Add a 0V test source between node_a and node_b to measure I_sc */
    int32_t test_br = ct_add_branch(&sc_circuit, CT_ELEM_VSOURCE,
                                     node_a, node_b, 0.0, 0.0, "Isc_test");

    ct_dc_result_t sc_result;
    if (ct_mna_dc_solve(&sc_circuit, &sc_result) != 0) return -1;

    double I_sc = sc_result.branch_currents[test_br];

    if (fabs(I_sc) > DBL_EPSILON) {
        *R_th = *V_th / I_sc;
        if (*R_th < 0) *R_th = -(*R_th);  /* R_th is always positive for passive nets */
    } else {
        *R_th = 1e12;  /* Effectively infinite (open circuit) */
    }

    return 0;
}

int ct_norton_equivalent(const ct_circuit_t *circuit,
                         int32_t node_a, int32_t node_b,
                         double *I_n, double *R_n)
{
    if (!circuit || !I_n || !R_n) return -1;

    double V_th, R_th;
    if (ct_thevenin_equivalent(circuit, node_a, node_b, &V_th, &R_th) != 0)
        return -1;

    if (R_th > DBL_EPSILON) {
        *I_n = V_th / R_th;
    } else {
        *I_n = 0.0;
    }
    *R_n = R_th;

    return 0;
}

int ct_max_power_transfer(double V_th, double R_th,
                          double *P_max, double *R_opt)
{
    if (!P_max || !R_opt || R_th <= 0.0) return -1;

    *R_opt = R_th;
    *P_max = (V_th * V_th) / (4.0 * R_th);

    return 0;
}

/* ==========================================================================
 * L3: Y-Delta Transformations
 * ========================================================================== */

int ct_y_to_delta(double R_a, double R_b, double R_c,
                  double *R_ab, double *R_bc, double *R_ca)
{
    if (!R_ab || !R_bc || !R_ca) return -1;
    if (R_a <= 0.0 || R_b <= 0.0 || R_c <= 0.0) return -1;

    double sum_prod = R_a * R_b + R_b * R_c + R_c * R_a;

    *R_ab = sum_prod / R_c;
    *R_bc = sum_prod / R_a;
    *R_ca = sum_prod / R_b;

    return 0;
}

int ct_delta_to_y(double R_ab, double R_bc, double R_ca,
                  double *R_a, double *R_b, double *R_c)
{
    if (!R_a || !R_b || !R_c) return -1;
    if (R_ab <= 0.0 || R_bc <= 0.0 || R_ca <= 0.0) return -1;

    double sum = R_ab + R_bc + R_ca;

    *R_a = R_ab * R_ca / sum;
    *R_b = R_ab * R_bc / sum;
    *R_c = R_bc * R_ca / sum;

    return 0;
}
