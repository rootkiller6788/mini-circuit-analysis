/**
 * @file spice_core.c
 * @brief Top-level SPICE simulator engine implementation
 *
 * Orchestrates the complete simulation workflow:
 *   netlist → allocate workspace → DC OP → AC → TRAN
 *
 * This module manages the lifecycle of the simulator object,
 * coordinates between netlist parsing, matrix assembly,
 * numerical solving, and result reporting.
 *
 * Reference: Nagel (1975) SPICE2 architecture §5
 */

#include "spice_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════
 * Simulator Lifecycle
 * ═══════════════════════════════════════════════════════════════════════ */

void spice_simulator_init(spice_simulator_t *sim) {
    if (!sim) return;
    memset(sim, 0, sizeof(*sim));
    spice_netlist_init(&sim->netlist);
    spice_convergence_defaults(&sim->conv_params);
    sim->initialized = 0;
}

int spice_simulator_allocate_workspace(spice_simulator_t *sim) {
    if (!sim) return -1;

    sim->mna_size  = spice_netlist_mna_size(&sim->netlist);
    sim->num_nodes = sim->netlist.num_nodes;

    if (sim->mna_size <= 0) return -1;

    /* Allocate MNA workspace */
    size_t mat_bytes = (size_t)sim->mna_size * sim->mna_size * sizeof(double);
    size_t vec_bytes = (size_t)sim->mna_size * sizeof(double);

    sim->mna_matrix = calloc(1, mat_bytes);
    sim->mna_rhs    = calloc(1, vec_bytes);
    sim->pivot      = calloc((size_t)sim->mna_size, sizeof(int32_t));
    sim->node_voltages = calloc((size_t)sim->mna_size, sizeof(double));
    sim->branch_currents = calloc((size_t)sim->mna_size, sizeof(double));

    if (!sim->mna_matrix || !sim->mna_rhs || !sim->pivot ||
        !sim->node_voltages) {
        free(sim->mna_matrix);
        free(sim->mna_rhs);
        free(sim->pivot);
        free(sim->node_voltages);
        free(sim->branch_currents);
        return -1;
    }

    sim->initialized = 1;
    return 0;
}

int spice_simulator_load(spice_simulator_t *sim, const char *filename) {
    if (!sim || !filename) return -1;

    /* Parse the netlist */
    int ret = spice_netlist_parse_file(&sim->netlist, filename);
    if (ret != 0) {
        fprintf(stderr, "Error parsing netlist: %s\n", filename);
        return -1;
    }

    /* Transfer netlist options to convergence params */
    if (sim->netlist.has_options) {
        sim->conv_params.abstol = sim->netlist.options_abstol;
        sim->conv_params.vntol  = sim->netlist.options_vntol;
        sim->conv_params.reltol = sim->netlist.options_reltol;
        sim->conv_params.temp   = sim->netlist.options_temp;
        sim->conv_params.itl1   = sim->netlist.options_itl1;
        sim->conv_params.itl4   = sim->netlist.options_itl4;
    }

    /* Allocate workspace */
    ret = spice_simulator_allocate_workspace(sim);
    if (ret != 0) {
        fprintf(stderr, "Error allocating workspace\n");
        return -1;
    }

    return 0;
}

int spice_simulator_run(spice_simulator_t *sim) {
    if (!sim || !sim->initialized) return -1;

    /* Always compute DC operating point first if transient or AC is requested */
    int need_dc = (sim->netlist.analysis_type == SPICE_ANALYSIS_DC ||
                   sim->netlist.analysis_type == SPICE_ANALYSIS_AC ||
                   sim->netlist.analysis_type == SPICE_ANALYSIS_TRAN ||
                   sim->netlist.analysis_type == SPICE_ANALYSIS_TF ||
                   sim->netlist.analysis_type == SPICE_ANALYSIS_NOISE);

    if (need_dc) {
        int ret = spice_simulator_run_dc(sim);
        if (ret != 0 && sim->netlist.analysis_type == SPICE_ANALYSIS_DC) {
            fprintf(stderr, "DC analysis failed to converge\n");
            return -1;
        }
        /* For AC/TRAN, we try to continue even if DC fails partially */
    }

    switch (sim->netlist.analysis_type) {
    case SPICE_ANALYSIS_DC:
        /* DC already done above */
        break;
    case SPICE_ANALYSIS_AC:
        spice_simulator_run_ac(sim);
        break;
    case SPICE_ANALYSIS_TRAN:
        spice_simulator_run_tran(sim);
        break;
    case SPICE_ANALYSIS_TF:
        if (sim->dc_op) {
            sim->tf_result = calloc(1, sizeof(spice_tf_result_t));
            if (sim->tf_result)
                spice_tf_analysis(&sim->netlist, sim->dc_op, sim->tf_result);
        }
        break;
    default:
        break;
    }

    return 0;
}

int spice_simulator_run_dc(spice_simulator_t *sim) {
    if (!sim || !sim->initialized) return -1;

    int num_nodes = sim->netlist.num_nodes;
    int num_branches = sim->netlist.num_vsources + sim->netlist.num_inductors;

    /* Free previous DC result */
    if (sim->dc_op) {
        spice_dc_result_free(sim->dc_op);
        sim->dc_op = NULL;
    }

    sim->dc_op = spice_dc_result_alloc(num_nodes, num_branches);
    if (!sim->dc_op) return -1;

    int ret = spice_dc_analysis(&sim->netlist, &sim->conv_params, sim->dc_op);
    if (ret == 0) {
        /* Store converged node voltages */
        for (int32_t i = 0; i < sim->mna_size && i <= num_nodes; i++) {
            sim->node_voltages[i] = sim->dc_op->node_voltages[i];
        }
    }

    return ret;
}

int spice_simulator_run_ac(spice_simulator_t *sim) {
    if (!sim || !sim->initialized) return -1;

    /* Free previous AC result */
    if (sim->ac_result) {
        spice_ac_result_free(sim->ac_result);
        sim->ac_result = NULL;
    }

    int num_nodes = sim->netlist.num_nodes;
    int n_freqs = sim->netlist.ac_points_per_decade;
    if (n_freqs < 10) n_freqs = 100;

    sim->ac_result = spice_ac_result_alloc(n_freqs, num_nodes);
    if (!sim->ac_result) return -1;

    return spice_ac_analysis(&sim->netlist, sim->dc_op,
                              &sim->conv_params, sim->ac_result);
}

int spice_simulator_run_tran(spice_simulator_t *sim) {
    if (!sim || !sim->initialized) return -1;

    /* Free previous transient result */
    if (sim->tran_result) {
        spice_tran_result_free(sim->tran_result);
        sim->tran_result = NULL;
    }

    int num_nodes = sim->netlist.num_nodes;
    int num_branches = sim->netlist.num_vsources + sim->netlist.num_inductors;
    int max_steps = 100000;

    sim->tran_result = spice_tran_result_alloc(max_steps, num_nodes, num_branches);
    if (!sim->tran_result) return -1;

    return spice_transient_analysis(&sim->netlist, &sim->conv_params,
                                     sim->tran_result);
}

void spice_simulator_print_results(const spice_simulator_t *sim) {
    if (!sim) return;

    printf("\n");
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║       SPICE Simulation Results               ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    if (sim->dc_op && sim->dc_op->converged) {
        printf("── DC Operating Point ──────────────────────────\n");
        printf("  Status: Converged in %d iterations\n", sim->dc_op->iterations);
        printf("  Node Voltages:\n");
        printf("    %-4s %12s\n", "Node", "Voltage (V)");
        printf("    %-4s %12s\n", "────", "───────────");
        printf("    %-4s %12.6e\n", "GND", 0.0);
        for (int32_t i = 1; i <= sim->dc_op->num_nodes; i++) {
            printf("    %-4d %12.6e\n", i, sim->dc_op->node_voltages[i]);
        }

        if (sim->dc_op->num_branches > 0) {
            printf("  Branch Currents:\n");
            for (int32_t i = 0; i < sim->dc_op->num_branches; i++) {
                printf("    Branch %d: %12.6e A\n", i + 1,
                       sim->dc_op->branch_currents[i]);
            }
        }
        printf("  Total DC Power: %12.6e W\n", sim->dc_op->total_power);
    }

    if (sim->ac_result && sim->ac_result->converged) {
        printf("\n── AC Analysis ─────────────────────────────────\n");
        printf("  Status: Converged\n");
        printf("  Frequency points: %d\n", sim->ac_result->num_freqs);
        printf("  Range: %.3e Hz — %.3e Hz\n",
               sim->ac_result->frequencies[0],
               sim->ac_result->frequencies[sim->ac_result->num_freqs - 1]);
    }

    if (sim->tran_result && sim->tran_result->converged) {
        printf("\n── Transient Analysis ──────────────────────────\n");
        printf("  Status: Converged\n");
        printf("  Time steps: %d\n", sim->tran_result->num_steps);
        printf("  Failed steps: %d\n", sim->tran_result->failed_steps);
        if (sim->tran_result->num_steps > 0) {
            printf("  Final time: %.6e s\n",
                   sim->tran_result->time_points[sim->tran_result->num_steps - 1]);
        }
    }

    if (sim->tf_result && sim->tf_result->converged) {
        printf("\n── Transfer Function ───────────────────────────\n");
        printf("  Gain: %12.6e\n", sim->tf_result->transfer_gain);
        printf("  Input impedance:  %12.6e Ω\n", sim->tf_result->input_impedance);
        printf("  Output impedance: %12.6e Ω\n", sim->tf_result->output_impedance);
    }
}

void spice_simulator_destroy(spice_simulator_t *sim) {
    if (!sim) return;

    spice_netlist_destroy(&sim->netlist);

    if (sim->dc_op)   spice_dc_result_free(sim->dc_op);
    if (sim->ac_result) spice_ac_result_free(sim->ac_result);
    if (sim->tran_result) spice_tran_result_free(sim->tran_result);
    free(sim->tf_result);

    free(sim->mna_matrix);
    free(sim->mna_rhs);
    free(sim->pivot);
    free(sim->node_voltages);
    free(sim->branch_currents);

    sim->dc_op = NULL;
    sim->ac_result = NULL;
    sim->tran_result = NULL;
    sim->tf_result = NULL;
    sim->mna_matrix = NULL;
    sim->mna_rhs = NULL;
    sim->pivot = NULL;
    sim->node_voltages = NULL;
    sim->branch_currents = NULL;
    sim->initialized = 0;
}

int spice_get_node_voltage(const spice_simulator_t *sim,
                            const char *nodename, double *voltage) {
    if (!sim || !nodename || !voltage) return -1;
    if (!sim->dc_op || !sim->dc_op->converged) return -1;

    if (strcmp(nodename, "0") == 0 || strcmp(nodename, "GND") == 0 ||
        strcmp(nodename, "gnd") == 0) {
        *voltage = 0.0;
        return 0;
    }

    /* Search node table */
    for (int32_t i = 1; i <= sim->netlist.num_nodes; i++) {
        if (strcmp(sim->netlist.nodes[i].name, nodename) == 0) {
            *voltage = sim->dc_op->node_voltages[i];
            return 0;
        }
    }

    return -1;
}

int spice_get_source_current(const spice_simulator_t *sim,
                              const char *srcname, double *current) {
    if (!sim || !srcname || !current) return -1;
    if (!sim->dc_op || !sim->dc_op->converged) return -1;

    /* Find the voltage source and its branch index */
    int br = 0;
    for (int32_t i = 0; i < sim->netlist.num_components; i++) {
        spice_component_t *comp = sim->netlist.components[i];
        if (comp && (comp->type == SPICE_COMP_VSOURCE ||
                     comp->type == SPICE_COMP_INDUCTOR)) {
            if (strcmp(comp->name, srcname) == 0) {
                if (br < sim->dc_op->num_branches) {
                    *current = sim->dc_op->branch_currents[br];
                    return 0;
                }
            }
            br++;
        }
    }

    return -1;
}

void spice_print_node_voltages(const spice_simulator_t *sim) {
    if (!sim || !sim->dc_op || !sim->dc_op->converged) {
        printf("No DC operating point available.\n");
        return;
    }
    printf("Node voltages (DC):\n");
    printf("  GND = 0.0 V\n");
    for (int32_t i = 1; i <= sim->dc_op->num_nodes; i++) {
        printf("  V(%s) = %+.6e V\n",
               sim->netlist.nodes[i].name,
               sim->dc_op->node_voltages[i]);
    }
}

void spice_print_branch_currents(const spice_simulator_t *sim) {
    if (!sim || !sim->dc_op || !sim->dc_op->converged) {
        printf("No DC operating point available.\n");
        return;
    }
    int br = 0;
    printf("Branch currents (DC):\n");
    for (int32_t i = 0; i < sim->netlist.num_components; i++) {
        spice_component_t *comp = sim->netlist.components[i];
        if (comp && (comp->type == SPICE_COMP_VSOURCE ||
                     comp->type == SPICE_COMP_INDUCTOR)) {
            if (br < sim->dc_op->num_branches) {
                printf("  I(%s) = %+.6e A\n",
                       comp->name, sim->dc_op->branch_currents[br]);
            }
            br++;
        }
    }
}

int spice_build_mna_system(spice_simulator_t *sim) {
    if (!sim || !sim->initialized) return -1;

    /* Zero the MNA matrix and RHS */
    memset(sim->mna_matrix, 0,
           (size_t)sim->mna_size * sim->mna_size * sizeof(double));
    memset(sim->mna_rhs, 0,
           (size_t)sim->mna_size * sizeof(double));

    spice_dense_matrix_t G;
    G.nrows = sim->mna_size;
    G.ncols = sim->mna_size;
    G.ld    = sim->mna_size;
    G.data  = sim->mna_matrix;

    /* Wrap RHS for stamp API */
    spice_vector_t rhs_wrapper;
    rhs_wrapper.size = sim->mna_size;
    rhs_wrapper.data = sim->mna_rhs;

    int branch_idx = sim->num_nodes;

    for (int32_t i = 0; i < sim->netlist.num_components; i++) {
        spice_component_t *comp = sim->netlist.components[i];
        if (!comp) continue;

        switch (comp->type) {
        case SPICE_COMP_RESISTOR: {
            spice_resistor_t *r = (spice_resistor_t*)comp;
            spice_stamp_resistor(&G, &rhs_wrapper, r->base.nplus, r->base.nminus,
                                  r->resistance, sim->mna_size);
            break;
        }
        case SPICE_COMP_VSOURCE: {
            spice_vsource_t *vs = (spice_vsource_t*)comp;
            int br = branch_idx++;
            spice_stamp_vsource_dc(&G, &rhs_wrapper, vs->base.nplus, vs->base.nminus,
                                    vs->dc_value, br, sim->mna_size);
            break;
        }
        case SPICE_COMP_ISOURCE: {
            spice_isource_t *is = (spice_isource_t*)comp;
            spice_stamp_isource_dc(&G, &rhs_wrapper, is->base.nplus, is->base.nminus,
                                    is->dc_value, sim->mna_size);
            break;
        }
        default:
            break;
        }
    }

    return 0;
}

int spice_solve_mna_system(spice_simulator_t *sim) {
    if (!sim || !sim->initialized) return -1;

    int info = spice_dense_lu_factor(sim->mna_matrix, sim->pivot,
                                      sim->mna_size, sim->mna_size);
    if (info > 0) return info;

    spice_dense_lu_solve(sim->mna_matrix, sim->pivot, sim->mna_rhs,
                          sim->mna_size, sim->mna_size);

    /* Store results in node_voltages */
    memcpy(sim->node_voltages, sim->mna_rhs,
           (size_t)sim->mna_size * sizeof(double));

    return 0;
}

int spice_check_dc_convergence(const spice_simulator_t *sim,
                                const double *v_old, const double *v_new) {
    if (!sim || !v_old || !v_new) return 0;

    for (int32_t i = 0; i < sim->mna_size; i++) {
        double delta = fabs(v_new[i] - v_old[i]);
        double tol = sim->conv_params.reltol *
                     fmax(fabs(v_new[i]), fabs(v_old[i])) +
                     sim->conv_params.vntol;
        if (delta > tol) return 0;
    }
    return 1;
}
