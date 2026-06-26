/**
 * @file superposition.c
 * @brief Superposition Theorem Implementation
 *
 * The Superposition Theorem states that in any linear circuit containing
 * multiple independent sources, the response (voltage or current) at any
 * point equals the algebraic sum of responses caused by each independent
 * source acting alone, with all other independent sources deactivated.
 *
 * Source deactivation rules:
 *   - Independent voltage sources → short circuit (0 V)
 *   - Independent current sources → open circuit (0 A)
 *   - Dependent sources remain active (they depend on circuit variables)
 *
 * Knowledge Coverage:
 *   L4 - Fundamental Laws: Superposition Theorem
 *   L5 - Algorithms: Iterative source deactivation solver
 *   L6 - Canonical Problems: Multi-source circuit analysis
 *
 * Reference:
 *   Desoer & Kuh "Basic Circuit Theory" (1969) Ch 10
 *   Alexander & Sadiku "Fundamentals of Electric Circuits" (2017) Ch 4
 *
 * Course Mapping:
 *   MIT 6.002: Linearity and superposition
 *   Berkeley EE16A: Superposition principle
 */

#include "network_theorem.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L5: Count Independent Sources in Circuit
 * ==========================================================================
 *
 * Scans the element list and counts voltage and current sources.
 * Returns the total count of independent sources.
 */
static uint32_t count_independent_sources(const CircuitTopology *ckt) {
    uint32_t count = 0;
    for (uint32_t i = 0; i < ckt->num_elements; i++) {
        ElementType t = ckt->elements[i].type;
        if (t == ELEM_VOLTAGE_SOURCE || t == ELEM_CURRENT_SOURCE) {
            count++;
        }
    }
    return count;
}

/* ==========================================================================
 * L5: Clone Circuit with One Source Active
 * ==========================================================================
 *
 * Creates a copy of the circuit topology where only source at index
 * 'active_idx' remains active; all other independent sources are
 * deactivated (V→0, I→0).
 *
 * Returns a newly allocated CircuitTopology that must be freed by caller.
 */
static CircuitTopology *clone_with_one_active_source(const CircuitTopology *ckt,
                                                      uint32_t active_idx) {
    if (!ckt) return NULL;

    CircuitTopology *clone = topology_create(ckt->num_nodes, ckt->num_branches,
                                              ckt->num_elements, ckt->num_meshes);
    if (!clone) return NULL;

    /* Copy nodes */
    memcpy(clone->nodes, ckt->nodes, ckt->num_nodes * sizeof(CircuitNode));

    /* Copy branches */
    memcpy(clone->branches, ckt->branches, ckt->num_branches * sizeof(CircuitBranch));

    /* Copy meshes (shallow copy of branch_ids — careful with double-free) */
    for (uint32_t m = 0; m < ckt->num_meshes; m++) {
        clone->meshes[m].id = ckt->meshes[m].id;
        clone->meshes[m].num_branches = ckt->meshes[m].num_branches;
        clone->meshes[m].mesh_current = 0.0;
        if (ckt->meshes[m].num_branches > 0) {
            clone->meshes[m].branch_ids = (uint32_t *)malloc(
                ckt->meshes[m].num_branches * sizeof(uint32_t));
            if (clone->meshes[m].branch_ids) {
                memcpy(clone->meshes[m].branch_ids, ckt->meshes[m].branch_ids,
                       ckt->meshes[m].num_branches * sizeof(uint32_t));
            }
        } else {
            clone->meshes[m].branch_ids = NULL;
        }
    }

    /* Copy elements but deactivate non-active independent sources */
    memcpy(clone->elements, ckt->elements, ckt->num_elements * sizeof(CircuitElement));

    uint32_t src_count = 0;
    for (uint32_t i = 0; i < ckt->num_elements; i++) {
        ElementType t = ckt->elements[i].type;

        if (t == ELEM_VOLTAGE_SOURCE || t == ELEM_CURRENT_SOURCE) {
            if (src_count != active_idx) {
                /* Deactivate this source */
                if (t == ELEM_VOLTAGE_SOURCE) {
                    /* Voltage source → short circuit (set value=0, keep as voltage source
                     * with 0V, which is equivalent to a short) */
                    clone->elements[i].value = 0.0;
                } else {
                    /* Current source → open circuit */
                    clone->elements[i].type = ELEM_OPEN;
                }
            }
            src_count++;
        }
    }

    return clone;
}

/* ==========================================================================
 * L4: Superposition Solver
 * ==========================================================================
 *
 * For a linear circuit with N independent sources, computes the voltage
 * at a specified target node as:
 *   V_target = Σ_{k=1}^{N} V_target(source_k acting alone)
 *
 * Algorithm:
 *   1. Count independent sources in the circuit
 *   2. For each source k:
 *      a. Clone the circuit with only source k active
 *      b. Solve nodal equations for this reduced circuit
 *      c. Record the voltage at the target node
 *   3. Sum all individual contributions
 *
 * Complexity: O(N * n³) where N = #sources, n = #nodes
 *
 * The result validates that superposition holds for linear circuits:
 * the sum of individual responses equals the response with all sources active.
 *
 * Returns: 0 on success, -1 on failure.
 */
int superposition_solve(const CircuitTopology *ckt,
                        uint32_t target_node,
                        SuperpositionResult *result) {
    if (!ckt || !result) return -1;
    if (target_node >= ckt->num_nodes) return -1;

    uint32_t num_sources = count_independent_sources(ckt);
    if (num_sources == 0) {
        result->total_response = 0.0;
        result->individual_contributions = NULL;
        result->num_sources = 0;
        result->is_voltage = 1;
        return 0;
    }

    result->num_sources = num_sources;
    result->is_voltage = 1;
    result->individual_contributions = (double *)calloc(num_sources, sizeof(double));
    if (!result->individual_contributions) return -1;

    double total = 0.0;

    /* Iterate over each independent source */
    for (uint32_t active = 0; active < num_sources; active++) {
        CircuitTopology *sub_ckt = clone_with_one_active_source(ckt, active);
        if (!sub_ckt) {
            free(result->individual_contributions);
            result->individual_contributions = NULL;
            return -1;
        }

        /* Solve nodal equations for this sub-circuit */
        uint32_t n = sub_ckt->num_nodes;
        double *Y = build_nodal_admittance_matrix(sub_ckt, &n);
        if (!Y) {
            topology_free(sub_ckt);
            free(result->individual_contributions);
            result->individual_contributions = NULL;
            return -1;
        }

        double *I_vec = (double *)calloc(n, sizeof(double));
        double *V = (double *)calloc(n, sizeof(double));
        if (!I_vec || !V) {
            free(Y); free(I_vec); free(V);
            topology_free(sub_ckt);
            free(result->individual_contributions);
            result->individual_contributions = NULL;
            return -1;
        }

        build_current_vector(sub_ckt, I_vec, n);

        int ret = solve_nodal_voltages(Y, I_vec, n, V);
        if (ret == 0) {
            double contribution = V[target_node];
            result->individual_contributions[active] = contribution;
            total += contribution;
        } else {
            result->individual_contributions[active] = 0.0;
        }

        free(Y);
        free(I_vec);
        free(V);
        topology_free(sub_ckt);
    }

    result->total_response = total;
    return 0;
}

/* ==========================================================================
 * L6: Superposition Verification
 * ==========================================================================
 *
 * Verifies the superposition principle by comparing the sum of individual
 * source responses against the response computed with all sources active
 * simultaneously.
 *
 * This is a canonical problem: does superposition hold for a given circuit?
 * For linear circuits, the difference should be zero (within numerical tolerance).
 *
 * Returns the absolute difference between the two methods.
 * A value near 0 confirms linearity; a large value suggests nonlinearity
 * or numerical issues in the matrix solver.
 */
double superposition_verify(const CircuitTopology *ckt, uint32_t target_node) {
    if (!ckt || target_node >= ckt->num_nodes) return -1.0;

    /* Method 1: All sources active simultaneously */
    uint32_t n = ckt->num_nodes;
    double *Y = build_nodal_admittance_matrix(ckt, &n);
    if (!Y) return -1.0;

    double *I_vec = (double *)calloc(n, sizeof(double));
    double *V_all = (double *)calloc(n, sizeof(double));
    if (!I_vec || !V_all) { free(Y); free(I_vec); free(V_all); return -1.0; }

    build_current_vector(ckt, I_vec, n);
    int ret = solve_nodal_voltages(Y, I_vec, n, V_all);
    free(Y); free(I_vec);

    if (ret != 0) { free(V_all); return -1.0; }
    double v_all = V_all[target_node];
    free(V_all);

    /* Method 2: Superposition sum */
    SuperpositionResult spr;
    ret = superposition_solve(ckt, target_node, &spr);
    if (ret != 0) return -1.0;

    double v_super = spr.total_response;
    free(spr.individual_contributions);

    return fabs(v_all - v_super);
}

/* ==========================================================================
 * L5: Source Deactivation Helper
 * ==========================================================================
 *
 * Given a circuit element array, deactivates all independent sources
 * to compute equivalent passive impedance.
 *
 * Deactivation:
 *   - Voltage source → short (set type to SHORT, value to 0)
 *   - Current source → open (set type to OPEN)
 *   - Dependent sources → remain active (their controlling variables
 *     are still present in the circuit)
 */
int deactivate_independent_sources(CircuitElement *elements, uint32_t count) {
    if (!elements || count == 0) return -1;

    for (uint32_t i = 0; i < count; i++) {
        switch (elements[i].type) {
            case ELEM_VOLTAGE_SOURCE:
                elements[i].type = ELEM_SHORT;
                elements[i].value = 0.0;
                break;
            case ELEM_CURRENT_SOURCE:
                elements[i].type = ELEM_OPEN;
                elements[i].value = 0.0;
                break;
            default:
                /* Passive elements and dependent sources stay as-is */
                break;
        }
    }
    return 0;
}
