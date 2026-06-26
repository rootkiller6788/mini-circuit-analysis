/**
 * @file spice_analysis.h
 * @brief SPICE analysis types — DC, AC, transient, and solver interfaces
 *
 * Knowledge coverage:
 *   L2 (Core Concepts): DC operating point, AC small-signal, transient analysis
 *   L5 (Algorithms): Newton-Raphson iteration, numerical integration methods
 *   L6 (Canonical Problems): Bias point finding, frequency response, step response
 *
 * Analysis flow:
 *   1. DC: Find operating point via Newton-Raphson on nonlinear system
 *   2. AC: Linearize around DC OP, solve complex MNA at each frequency
 *   3. TRAN: Time-stepping with companion models + Newton at each step
 *
 * Reference: Nagel & Pederson, "SPICE — Simulation Program with Integrated
 *            Circuit Emphasis" (1973), UCB/ERL M382
 *
 * Course alignment:
 *   MIT 6.003 (Signal Processing): frequency response / Bode plots
 *   Stanford EE102A (Signal Processing): transient analysis, convolution
 *   ETH 227-0427 (Signal Processing): time/frequency domain analysis
 */

#ifndef SPICE_ANALYSIS_H
#define SPICE_ANALYSIS_H

#include <stdint.h>
#include <stddef.h>
#include "spice_matrix.h"
#include "spice_netlist.h"

/* ── L2: Convergence Parameters ───────────────────────────────────── */

typedef struct {
    double    abstol;       /**< Absolute current tolerance (A), default 1pA */
    double    vntol;        /**< Absolute voltage tolerance (V), default 1µV */
    double    reltol;       /**< Relative tolerance, default 0.001 (0.1%)    */
    double    chgtol;       /**< Charge tolerance for transient, default 1fC  */
    int32_t   itl1;         /**< DC iteration limit, default 100             */
    int32_t   itl2;         /**< DC transfer curve iteration limit           */
    int32_t   itl4;         /**< Transient analysis iteration limit          */
    double    temp;          /**< Circuit temperature in °C, default 27       */
    int32_t   max_time_steps; /**< Maximum transient time steps              */
    double    time_step_min;  /**< Minimum allowed time step (s)             */
    double    gmin;           /**< Minimum conductance (S) for convergence   */
} spice_convergence_params_t;

/**
 * @brief Initialize convergence parameters to SPICE defaults
 *
 * @param params Structure to initialize
 */
void spice_convergence_defaults(spice_convergence_params_t *params);

/* ── L2: Analysis Result Types ────────────────────────────────────── */

/**
 * @brief DC operating point result
 *
 * Contains node voltages and branch currents at the
 * converged DC solution.
 */
typedef struct {
    int32_t    num_nodes;       /**< Number of nodes (excluding GND)  */
    int32_t    num_branches;    /**< Number of branch currents        */
    double*    node_voltages;   /**< V[1..num_nodes], GND is 0        */
    double*    branch_currents; /**< Branch currents for V-sources, L */
    int32_t    converged;       /**< 1 if Newton converged             */
    int32_t    iterations;      /**< Number of Newton iterations      */
    double     total_power;     /**< Total DC power dissipation (W)   */
} spice_dc_result_t;

/**
 * @brief AC analysis result
 *
 * Frequency-dependent complex node voltages.
 * Used for Bode plots, impedance analysis, transfer functions.
 */
typedef struct {
    int32_t            num_freqs;          /**< Number of frequency points      */
    int32_t            num_nodes;          /**< Number of nodes                 */
    double*            frequencies;        /**< Frequency array (Hz)            */
    spice_complex_t**  node_voltages;      /**< V[freq_idx][node_idx]           */
    int32_t            converged;          /**< 1 if all points solved OK       */
} spice_ac_result_t;

/**
 * @brief Transient analysis result
 *
 * Time-domain waveforms for node voltages and branch currents.
 */
typedef struct {
    int32_t    num_steps;        /**< Number of time points saved       */
    int32_t    num_nodes;        /**< Number of nodes                  */
    int32_t    num_branches;     /**< Number of branch currents        */
    double*    time_points;      /**< Time array (seconds)             */
    double**   node_voltages;    /**< V[step][node]                    */
    double**   branch_currents;  /**< I[step][branch]                  */
    int32_t    converged;        /**< 1 if full simulation completed    */
    int32_t    failed_steps;     /**< Number of non-convergent steps   */
} spice_tran_result_t;

/**
 * @brief Transfer function result (.TF analysis)
 *
 * Finds DC small-signal transfer function V(out)/V(in) or V(out)/I(in),
 * plus input and output impedances.
 */
typedef struct {
    double    transfer_gain;      /**< V(out)/V(in) or V(out)/I(in)  */
    double    input_impedance;    /**< Input resistance (ohms)        */
    double    output_impedance;   /**< Output resistance (ohms)       */
    char      output_name[64];    /**< Name of output variable        */
    char      input_name[64];     /**< Name of input source           */
    int32_t   converged;
} spice_tf_result_t;

/* ── L2/L6: Analysis Executors ────────────────────────────────────── */

/**
 * @brief Perform DC operating point analysis
 *
 * Solves the nonlinear algebraic system F(V) = 0 using Newton-Raphson
 * iteration. Convergence is checked against abstol/vntol/reltol.
 *
 * @param netlist Parsed circuit netlist
 * @param params  Convergence parameters
 * @param result  Output DC result (caller-allocated)
 * @return 0 on success, -1 on non-convergence
 *
 * Algorithm: Newton-Raphson iteration
 *   V_{k+1} = V_k - J^{-1}(V_k) * F(V_k)
 * where J is the Jacobian (conductance + derivative stamps).
 *
 * Reference: Nagel (1975) §4.3, Kielkowski (1998) §3
 * Complexity: O(N^3) per iteration (dense LU)
 */
int spice_dc_analysis(const spice_netlist_t *netlist,
                       const spice_convergence_params_t *params,
                       spice_dc_result_t *result);

/**
 * @brief Perform AC small-signal analysis
 *
 * Sweeps frequency from fstart to fstop logarithmically.
 * At each frequency, linearizes all nonlinear devices around
 * the DC operating point and solves the complex MNA system.
 *
 * @param netlist Parsed circuit
 * @param dc_op   Previously computed DC operating point
 * @param params  Convergence parameters
 * @param result  Output AC result
 * @return 0 on success
 *
 * The complex MNA matrix is: G + jωC, where G is the DC conductance
 * and C is the capacitance matrix (including device capacitances).
 *
 * Reference: Vlach & Singhal, "Computer Methods for Circuit Analysis
 *            and Design" (1994) §6
 * Complexity: O(N_freq * N^3)
 */
int spice_ac_analysis(const spice_netlist_t *netlist,
                       const spice_dc_result_t *dc_op,
                       const spice_convergence_params_t *params,
                       spice_ac_result_t *result);

/**
 * @brief Perform transient analysis
 *
 * Time-domain simulation using numerical integration (trapezoidal rule
 * default, backward Euler for problematic steps). At each time step,
 * nonlinear devices are solved via Newton-Raphson (inner loop).
 *
 * Time step control: local truncation error (LTE) estimate determines
 * whether to accept/reject a step and size the next step.
 *
 * @param netlist Parsed circuit
 * @param params  Convergence parameters (includes time step limits)
 * @param result  Output transient result
 * @return 0 on success
 *
 * Integration methods:
 *   TRAP: y_{n+1} = y_n + (h/2)*(f_n + f_{n+1})   — A-stable, 2nd order
 *   BE:   y_{n+1} = y_n + h * f_{n+1}             — L-stable, 1st order
 *
 * Reference: Chua & Lin, "Computer-Aided Analysis of Electronic
 *            Circuits" (1975) §7-8
 * Complexity: O(N_steps * N_iter * N^3)
 */
int spice_transient_analysis(const spice_netlist_t *netlist,
                              const spice_convergence_params_t *params,
                              spice_tran_result_t *result);

/**
 * @brief Perform DC transfer function analysis (.TF)
 *
 * Computes small-signal transfer function V(out)/V(in) and
 * input/output impedances by solving the adjoint system.
 *
 * @param netlist Parsed circuit
 * @param dc_op   DC operating point
 * @param result  Output TF result
 * @return 0 on success
 *
 * Reference: Director & Rohrer, "The Generalized Adjoint Network
 *            and Network Sensitivities" (1969), IEEE CT-16
 */
int spice_tf_analysis(const spice_netlist_t *netlist,
                       const spice_dc_result_t *dc_op,
                       spice_tf_result_t *result);

/* ── Result Management ─────────────────────────────────────────────── */

/**
 * @brief Allocate a DC result structure
 *
 * @param num_nodes    Number of nodes
 * @param num_branches Number of branch currents
 * @return Allocated result or NULL
 */
spice_dc_result_t* spice_dc_result_alloc(int32_t num_nodes, int32_t num_branches);

/**
 * @brief Free a DC result
 */
void spice_dc_result_free(spice_dc_result_t *result);

/**
 * @brief Allocate an AC result structure
 *
 * @param num_freqs Number of frequency points
 * @param num_nodes Number of nodes
 * @return Allocated result or NULL
 */
spice_ac_result_t* spice_ac_result_alloc(int32_t num_freqs, int32_t num_nodes);

/**
 * @brief Free an AC result
 */
void spice_ac_result_free(spice_ac_result_t *result);

/**
 * @brief Allocate a transient result with initial capacity
 *
 * @param max_steps   Maximum number of time steps to store
 * @param num_nodes   Number of nodes
 * @param num_branches Number of branch currents
 * @return Allocated result or NULL
 */
spice_tran_result_t* spice_tran_result_alloc(int32_t max_steps,
                                              int32_t num_nodes,
                                              int32_t num_branches);

/**
 * @brief Free a transient result
 */
void spice_tran_result_free(spice_tran_result_t *result);

/**
 * @brief Write transient result to CSV file for plotting
 *
 * @param result   Transient result to export
 * @param filename Output CSV file path
 * @return 0 on success
 */
int spice_tran_export_csv(const spice_tran_result_t *result, const char *filename);

/**
 * @brief Write AC result to CSV file (frequency, magnitude, phase per node)
 *
 * @param result   AC result
 * @param filename Output CSV file path
 * @return 0 on success
 */
int spice_ac_export_csv(const spice_ac_result_t *result, const char *filename);

#endif /* SPICE_ANALYSIS_H */
