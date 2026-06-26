/**
 * @file spice_core.h
 * @brief Top-level SPICE simulator engine
 *
 * Knowledge coverage:
 *   L2 (Core Concepts): Simulation workflow — parse → DC → AC → TRAN
 *   L6 (Canonical Problems): End-to-end circuit simulation pipeline
 *
 * The simulator orchestrates the three analysis phases in the
 * standard SPICE order:
 *   1. Parse netlist (build circuit topology)
 *   2. DC operating point (if needed for AC or as standalone)
 *   3. AC analysis (if requested)
 *   4. Transient analysis (if requested)
 *
 * Reference: Nagel (1975), "SPICE2: A Computer Program to Simulate
 *            Semiconductor Circuits", UCB/ERL M520
 *
 * Course alignment:
 *   Berkeley EE105 (Analog IC): Full simulation flow
 *   Michigan EECS 411 (Microwave): Frequency-domain simulation
 *   Georgia Tech ECE 6601 (Comm): Circuit-system co-simulation
 */

#ifndef SPICE_CORE_H
#define SPICE_CORE_H

#include <stdint.h>
#include "spice_netlist.h"
#include "spice_matrix.h"
#include "spice_analysis.h"
#include "spice_models.h"

/* ── L2: Simulator State ───────────────────────────────────────────── */

/**
 * @brief Top-level SPICE simulator object
 *
 * Owns the netlist, matrix workspace, solver state,
 * and analysis results for a single simulation run.
 */
typedef struct {
    spice_netlist_t             netlist;    /**< Parsed circuit netlist      */
    spice_convergence_params_t  conv_params;/**< Convergence/tolerance params  */
    spice_dc_result_t*          dc_op;      /**< DC operating point result   */
    spice_ac_result_t*          ac_result;  /**< AC analysis result          */
    spice_tran_result_t*        tran_result;/**< Transient analysis result    */
    spice_tf_result_t*          tf_result;  /**< Transfer function result     */
    double*                     mna_matrix; /**< Dense MNA workspace (n^2)    */
    double*                     mna_rhs;   /**< MNA RHS workspace (n)         */
    int32_t*                    pivot;      /**< LU pivot workspace (n)       */
    int32_t                     mna_size;   /**< MNA system dimension         */
    int32_t                     num_nodes;  /**< Number of non-ground nodes   */
    double*                     node_voltages; /**< Current node voltage vector */
    double*                     branch_currents; /**< Current branch currents  */
    int32_t                     initialized;/**< 1 if workspace allocated     */
} spice_simulator_t;

/* ── Simulator Lifecycle ───────────────────────────────────────────── */

/**
 * @brief Initialize an empty simulator
 *
 * @param sim Simulator to initialize
 *
 * Sets all pointers to NULL and parameters to defaults.
 */
void spice_simulator_init(spice_simulator_t *sim);

/**
 * @brief Load a SPICE netlist file into the simulator
 *
 * Parses the netlist file and allocates internal workspace.
 *
 * @param sim      Simulator
 * @param filename Path to netlist file (.cir)
 * @return 0 on success, negative on error
 */
int spice_simulator_load(spice_simulator_t *sim, const char *filename);

/**
 * @brief Allocate MNA workspace for the loaded netlist
 *
 * Called internally by spice_simulator_load.
 *
 * @param sim Simulator with loaded netlist
 * @return 0 on success
 */
int spice_simulator_allocate_workspace(spice_simulator_t *sim);

/**
 * @brief Run all requested analyses in proper order
 *
 * Sequence: DC → AC → TRAN (or whatever the netlist requests)
 *
 * @param sim Loaded simulator
 * @return 0 on success
 */
int spice_simulator_run(spice_simulator_t *sim);

/**
 * @brief Run only DC operating point
 *
 * @param sim Loaded simulator
 * @return 0 on success
 */
int spice_simulator_run_dc(spice_simulator_t *sim);

/**
 * @brief Run only AC analysis
 *
 * Requires DC operating point to be computed first.
 *
 * @param sim Loaded simulator with DC OP
 * @return 0 on success
 */
int spice_simulator_run_ac(spice_simulator_t *sim);

/**
 * @brief Run only transient analysis
 *
 * Uses DC operating point as initial condition.
 *
 * @param sim Loaded simulator
 * @return 0 on success
 */
int spice_simulator_run_tran(spice_simulator_t *sim);

/**
 * @brief Print simulation results to stdout
 *
 * @param sim Simulator with completed analyses
 */
void spice_simulator_print_results(const spice_simulator_t *sim);

/**
 * @brief Free all resources associated with the simulator
 *
 * @param sim Simulator to destroy
 */
void spice_simulator_destroy(spice_simulator_t *sim);

/**
 * @brief Get voltage at a named node
 *
 * Looks up the node by name and returns its DC voltage.
 * Requires DC analysis to have completed.
 *
 * @param sim      Simulator with DC results
 * @param nodename Node name string
 * @param voltage  Output: node voltage (referenced to GND)
 * @return 0 on success, -1 if node not found
 */
int spice_get_node_voltage(const spice_simulator_t *sim,
                            const char *nodename, double *voltage);

/**
 * @brief Get current through a named voltage source
 *
 * Requires DC analysis to have completed.
 *
 * @param sim      Simulator with DC results
 * @param srcname  Voltage source name (e.g., "V1")
 * @param current  Output: branch current
 * @return 0 on success, -1 if source not found
 */
int spice_get_source_current(const spice_simulator_t *sim,
                              const char *srcname, double *current);

/**
 * @brief Print voltage at every node
 *
 * @param sim Simulator with DC results
 */
void spice_print_node_voltages(const spice_simulator_t *sim);

/**
 * @brief Print branch currents for all voltage sources and inductors
 *
 * @param sim Simulator with DC results
 */
void spice_print_branch_currents(const spice_simulator_t *sim);

/**
 * @brief Build the full MNA system for the current Newton iterate
 *
 * Assembles G matrix and RHS vector from all components,
 * evaluating nonlinear devices at the current voltage estimate.
 *
 * @param sim Simulator with netlist and voltage estimate
 * @return 0 on success
 *
 * This is the core assembly routine called in each Newton iteration
 * and at each transient time step.
 */
int spice_build_mna_system(spice_simulator_t *sim);

/**
 * @brief Solve the assembled MNA system
 *
 * Performs LU factorization and forward/backward substitution
 * to obtain the next voltage estimate.
 *
 * @param sim Simulator with assembled MNA system
 * @return 0 on success, >0 if singular
 */
int spice_solve_mna_system(spice_simulator_t *sim);

/**
 * @brief Check DC convergence
 *
 * Tests whether all node voltage changes and branch current changes
 * are within the specified tolerances.
 *
 * @param sim        Simulator
 * @param v_old      Previous voltage vector
 * @param v_new      Current voltage vector
 * @return 1 if converged, 0 if not
 */
int spice_check_dc_convergence(const spice_simulator_t *sim,
                                const double *v_old, const double *v_new);

#endif /* SPICE_CORE_H */
