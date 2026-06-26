/**
 * @file circuit_mna.h
 * @brief Modified Nodal Analysis (MNA) for circuit simulation
 *
 * MNA is the standard formulation used in SPICE and all modern circuit
 * simulators. It combines nodal analysis (KCL at each node) with branch
 * constitutive relations to handle voltage sources and controlled sources
 * without the limitations of classical node analysis.
 *
 * The MNA system has the block structure:
 *   [ Y  B ] [ v_n ]   [ i_s ]
 *   [ C  D ] [ i_v ] = [ v_s ]
 *
 * where:
 *   Y = node admittance matrix (n x n)
 *   B = voltage source incidence (n x nv)
 *   C = voltage source KVL constraints (nv x n)
 *   D = zeros for independent sources (nv x nv)
 *   v_n = node voltage vector (unknown)
 *   i_v = voltage source current vector (unknown)
 *   i_s = equivalent current source vector (known)
 *   v_s = voltage source values (known)
 *
 * References:
 *   - C.W. Ho, A.E. Ruehli, P.A. Brennan, "The Modified Nodal Approach
 *     to Network Analysis", IEEE TCAS, 1975
 *   - L.W. Nagel, "SPICE2: A Computer Program to Simulate Semiconductor
 *     Circuits", UCB ERL-M520, 1975
 *   - Berkeley EE16A / MIT 6.002 / Tsinghua Circuit Principles
 *
 * Knowledge coverage:
 *   L5 (Algorithms): MNA stamping, LU decomposition, forward-back substitution
 *   L6 (Canonical):  DC operating point, AC frequency sweep, transient analysis
 *   L7 (Applications): SPICE-compatible netlist solver
 */

#ifndef CIRCUIT_MNA_H
#define CIRCUIT_MNA_H

#include "circuit_topology.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * MNA Matrix Dimensions
 * ========================================================================== */

/** Maximum MNA system size (nodes + voltage sources) */
#define CT_MNA_MAX_SIZE  (CT_MAX_NODES + CT_MAX_BRANCHES)

/* ==========================================================================
 * L1: MNA Data Structures
 * ========================================================================== */

/**
 * @brief MNA system matrix and right-hand-side vector.
 *
 * The MNA matrix is of size (n + nv) x (n + nv), where n = number of
 * non-ground nodes, nv = number of voltage source branches.
 * The RHS vector has the same dimension.
 */
typedef struct {
    int32_t  size;                       /**< MNA system size (n + nv) */
    int32_t  num_nodes;                  /**< Number of non-ground nodes (n) */
    int32_t  num_vsrc;                   /**< Number of voltage sources (nv) */

    /** MNA matrix G (size x size), stored row-major.
     * For DC: real-valued (conductance matrix).
     * For AC: complex-valued (admittance matrix at frequency omega).
     */
    double complex G[CT_MNA_MAX_SIZE][CT_MNA_MAX_SIZE];

    /** Right-hand-side vector */
    double complex rhs[CT_MNA_MAX_SIZE];

    /** Solution vector (unknowns) */
    double complex sol[CT_MNA_MAX_SIZE];

    /** Node-to-MNA-row mapping: mna_row[node_idx] = row in G */
    int32_t  node_map[CT_MAX_NODES];

    /** Voltage-source-to-MNA-row mapping: vsrc_map[br_idx] = row in G */
    int32_t  vsrc_map[CT_MAX_BRANCHES];

    /** Inverse: MNA-row to node/vsrc type and index */
    int8_t   row_type[CT_MNA_MAX_SIZE];  /**< 0=node voltage, 1=vsrc current */
    int32_t  row_index[CT_MNA_MAX_SIZE]; /**< Node ID or branch ID */

    /** Pivot array for LU decomposition (row permutation) */
    int32_t  pivot[CT_MNA_MAX_SIZE];

    /** Flag: LU decomposition has been computed */
    uint8_t  is_factored;
} ct_mna_system_t;

/**
 * @brief DC operating point result.
 */
typedef struct {
    double    node_voltages[CT_MAX_NODES];   /**< V[n] for each node (V[0]=0) */
    double    branch_currents[CT_MAX_BRANCHES]; /**< I[b] for each branch */
    double    branch_voltages[CT_MAX_BRANCHES]; /**< V[b] for each branch */
    double    branch_power[CT_MAX_BRANCHES];    /**< P[b] = V[b] * I[b] */
    double    total_power_dissipated;            /**< Sum of P[b] for passive elements */
    double    total_power_supplied;              /**< Sum of P[b] for active elements */
    uint8_t   converged;                         /**< Flag: solution converged */
    int32_t   iterations;                        /**< Newton iterations taken */
} ct_dc_result_t;

/**
 * @brief AC analysis result at a single frequency.
 */
typedef struct {
    double    frequency;                       /**< Analysis frequency [Hz] */
    double complex node_voltages[CT_MAX_NODES];   /**< Phasor node voltages */
    double complex branch_currents[CT_MAX_BRANCHES]; /**< Phasor branch currents */
    double    magnitude_v[CT_MAX_NODES];        /**< |V| at each node */
    double    phase_v[CT_MAX_NODES];            /**< phase(V) in degrees */
    double    input_impedance;                  /**< Z_in at the input port */
    double    output_impedance;                 /**< Z_out at the output port */
    double    voltage_gain;                     /**< |V_out/V_in| */
    double    phase_shift;                      /**< phase(V_out) - phase(V_in) */
} ct_ac_result_t;

/**
 * @brief Transient analysis configuration.
 */
typedef struct {
    double    t_start;          /**< Start time [s] */
    double    t_stop;           /**< Stop time [s] */
    double    t_step;           /**< Time step [s] */
    double    t_max_step;       /**< Maximum allowed time step [s] */
    double    integration_method; /**< 0=Backward Euler, 1=Trapezoidal, 2=Gear-2 */
    double    rel_tol;          /**< Relative tolerance for convergence */
    double    abs_tol;          /**< Absolute tolerance for convergence */
    int32_t   max_iterations;   /**< Maximum Newton iterations per time step */
    uint8_t   store_all_steps;  /**< Flag: store results at every step */
} ct_transient_config_t;

/* ==========================================================================
 * L5: MNA Stamping Functions
 * ========================================================================== */

/**
 * @brief Initialize an MNA system from a circuit topology.
 *
 * Allocates the MNA matrix with proper dimensions: size = n + nv.
 * Sets up node-to-row and vsrc-to-row mappings.
 *
 * @param circuit  Circuit topology
 * @param mna      MNA system to initialize
 * @return 0 on success, -1 on error
 */
int ct_mna_init(const ct_circuit_t *circuit, ct_mna_system_t *mna);

/**
 * @brief Stamp a linear resistor into the MNA matrix.
 *
 * For a resistor R between nodes n+ and n-, adds:
 *   G[n+][n+] += 1/R,  G[n+][n-] -= 1/R
 *   G[n-][n+] -= 1/R,  G[n-][n-] += 1/R
 *
 * @param mna       MNA system
 * @param node_p    Positive terminal node ID
 * @param node_n    Negative terminal node ID
 * @param resistance Resistance value [Ohm] (must be > 0)
 * @return 0 on success, -1 on invalid parameters
 */
int ct_mna_stamp_resistor(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                          double resistance);

/**
 * @brief Stamp a linear capacitor into the MNA matrix.
 *
 * For DC: capacitor is an open circuit (no stamp).
 * For AC: adds j*omega*C to the admittance matrix.
 *   G[n+][n+] += jwC,  G[n+][n-] -= jwC
 *   G[n-][n+] -= jwC,  G[n-][n-] += jwC
 *
 * For transient analysis, the capacitor is handled via companion models
 * (trapezoidal or backward Euler integration).
 *
 * @param mna          MNA system
 * @param node_p       Positive terminal node ID
 * @param node_n       Negative terminal node ID
 * @param capacitance  Capacitance value [Farad]
 * @param omega        Angular frequency [rad/s] (0 for DC)
 * @return 0 on success, -1 on error
 */
int ct_mna_stamp_capacitor(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                           double capacitance, double omega);

/**
 * @brief Stamp a linear inductor into the MNA matrix.
 *
 * For DC: inductor is a short circuit (handled as special branch).
 * For AC: adds 1/(j*omega*L) to the admittance matrix.
 *   G[n+][n+] += 1/(jwL),  G[n+][n-] -= 1/(jwL)
 *   G[n-][n+] -= 1/(jwL),  G[n-][n-] += 1/(jwL)
 *
 * @param mna         MNA system
 * @param node_p      Positive terminal node ID
 * @param node_n      Negative terminal node ID
 * @param inductance  Inductance value [Henry]
 * @param omega       Angular frequency [rad/s]
 * @return 0 on success, -1 on error
 */
int ct_mna_stamp_inductor(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                          double inductance, double omega);

/**
 * @brief Stamp an independent voltage source into the MNA matrix.
 *
 * A voltage source V between nodes n+ and n- introduces:
 *   - A new unknown: the current i_v through the source
 *   - One row: v(n+) - v(n-) = V (KVL constraint)
 *   - Columns in existing rows: +/- 1 for the i_v unknown
 *
 * @param mna      MNA system
 * @param node_p   Positive terminal node ID
 * @param node_n   Negative terminal node ID
 * @param vsrc_id  Voltage source branch ID
 * @param voltage  Source voltage value [V]
 * @return 0 on success, -1 on error
 */
int ct_mna_stamp_vsource(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                         int32_t vsrc_id, double voltage);

/**
 * @brief Stamp an independent current source into the MNA RHS.
 *
 * A current source I from node n+ to n- adds:
 *   rhs[n+] -= I
 *   rhs[n-] += I
 *
 * @param mna      MNA system
 * @param node_p   Current entering node (arrow tip)
 * @param node_n   Current leaving node (arrow tail)
 * @param current  Source current value [A]
 * @return 0 on success, -1 on error
 */
int ct_mna_stamp_isource(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                         double current);

/**
 * @brief Stamp a voltage-controlled current source (VCCS) into MNA.
 *
 * A VCCS: i_out = gm * v_control, where v_control is the voltage between
 * nodes nc+ and nc-. The stamp adds gm to G at (nout+, nc+) and (nout-, nc-)
 * and subtracts at cross terms.
 *
 * This is a transconductance stamp — fundamental for active circuit modeling
 * (MOSFET small-signal, BJT hybrid-pi, op-amp macro-models, gyrators).
 *
 * @param mna        MNA system
 * @param nc_p       Positive control node
 * @param nc_n       Negative control node
 * @param nout_p     Positive output node
 * @param nout_n     Negative output node
 * @param gm         Transconductance [Siemens]
 * @return 0 on success, -1 on error
 */
int ct_mna_stamp_vccs(ct_mna_system_t *mna, int32_t nc_p, int32_t nc_n,
                      int32_t nout_p, int32_t nout_n, double gm);

/**
 * @brief Stamp a current-controlled current source (CCCS) into MNA.
 *
 * A CCCS: i_out = beta * i_control, where i_control is the current through
 * a specified branch (typically a voltage source used as an ammeter).
 *
 * @param mna           MNA system
 * @param control_br    Branch ID of the controlling current
 * @param nout_p        Positive output node
 * @param nout_n        Negative output node
 * @param beta          Current gain [dimensionless]
 * @return 0 on success, -1 on error
 */
int ct_mna_stamp_cccs(ct_mna_system_t *mna, int32_t control_br,
                      int32_t nout_p, int32_t nout_n, double beta);

/**
 * @brief Stamp an ideal op-amp into the MNA matrix.
 *
 * An ideal op-amp enforces v(+) = v(-) with infinite input impedance
 * and zero output impedance. Implemented via a nullor (nullator + norator pair).
 * This adds one constraint equation: v(n+) - v(n-) = 0.
 *
 * @param mna       MNA system
 * @param node_p    Non-inverting input node ID
 * @param node_n    Inverting input node ID
 * @param node_out  Output node ID
 * @return 0 on success, -1 on error
 */
int ct_mna_stamp_opamp(ct_mna_system_t *mna, int32_t node_p, int32_t node_n,
                       int32_t node_out);

/* ==========================================================================
 * L5: MNA Solving Functions
 * ========================================================================== */

/**
 * @brief Perform LU decomposition of the MNA matrix with partial pivoting.
 *
 * Decomposes G = P * L * U, where P is a permutation matrix (stored in pivot[]),
 * L is unit lower triangular, and U is upper triangular.
 *
 * Complexity: O(size^3 / 3) floating-point operations.
 *
 * Implements the Crout algorithm variant for in-place LU factorization.
 * For ill-conditioned matrices (near-singular), pivot growth factor is
 * monitored and reported.
 *
 * @param mna       MNA system (G is factorized in-place)
 * @param singular  Output: 1 if matrix is numerically singular, 0 otherwise
 * @return 0 on success, -1 on error
 */
int ct_mna_lu_factor(ct_mna_system_t *mna, int *singular);

/**
 * @brief Solve MNA system after LU factorization using forward/back substitution.
 *
 * Solves G * x = rhs given L and U factors:
 *   1. Forward substitution: L * y = P * rhs
 *   2. Back substitution:   U * x = y
 *
 * Complexity: O(size^2).
 *
 * @param mna  MNA system (must have is_factored = 1)
 * @return 0 on success, -1 on error
 */
int ct_mna_solve(ct_mna_system_t *mna);

/**
 * @brief Build the complete MNA system from a circuit (all elements stamped).
 *
 * Iterates through all branches and applies appropriate stamps based on
 * element type. This is the main assembly function.
 *
 * @param circuit  Circuit topology
 * @param mna      MNA system to populate
 * @return 0 on success, -1 on error
 */
int ct_mna_assemble(const ct_circuit_t *circuit, ct_mna_system_t *mna);

/**
 * @brief Solve for DC operating point using MNA.
 *
 * For linear circuits: one assembly + one LU solve suffices.
 * For nonlinear circuits: Newton-Raphson iteration with linearized models.
 *
 * @param circuit  Circuit topology
 * @param result   Output DC operating point
 * @return 0 on success, -1 on failure to converge
 */
int ct_mna_dc_solve(const ct_circuit_t *circuit, ct_dc_result_t *result);

/**
 * @brief Perform AC analysis at a specified frequency.
 *
 * Builds the complex admittance matrix at the given frequency,
 * solves for phasor node voltages and branch currents.
 *
 * @param circuit   Circuit topology
 * @param freq_hz   Frequency in Hz
 * @param result    Output AC analysis result
 * @return 0 on success, -1 on error
 */
int ct_mna_ac_solve(const ct_circuit_t *circuit, double freq_hz,
                    ct_ac_result_t *result);

/**
 * @brief Extract branch voltages and currents from MNA solution.
 *
 * After solving the MNA system, compute:
 *   - Branch voltages: v_b = A^T * v_n (from node potentials)
 *   - Branch currents: from branch constitutive relations
 *
 * @param circuit         Circuit topology
 * @param mna             Solved MNA system
 * @param branch_voltages Output: branch voltage array [b]
 * @param branch_currents Output: branch current array [b]
 * @return 0 on success, -1 on error
 */
int ct_mna_extract_branch_results(const ct_circuit_t *circuit,
                                  const ct_mna_system_t *mna,
                                  double *branch_voltages,
                                  double *branch_currents);

#ifdef __cplusplus
}
#endif

#endif /* CIRCUIT_MNA_H */
