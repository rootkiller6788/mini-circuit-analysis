/**
 * @file circuit_analysis.h
 * @brief Classical circuit analysis methods
 *
 * Implements the fundamental analysis techniques predating MNA:
 * node-voltage method, mesh-current method, superposition, and
 * Thevenin/Norton equivalent computation using topological information.
 *
 * While MNA is the industrial standard, these classical methods provide
 * essential insights into circuit behavior and form the theoretical
 * foundation for understanding more advanced techniques.
 *
 * References:
 *   - W.H. Hayt, J.E. Kemmerly, S.M. Durbin, "Engineering Circuit Analysis"
 *   - C.A. Desoer, E.S. Kuh, "Basic Circuit Theory" (1969)
 *   - MIT 6.002 / Berkeley EE16A/B / Michigan EECS 215
 *
 * Knowledge coverage:
 *   L2 (Concepts):  Node-voltage, mesh-current, superposition
 *   L3 (Math):      Nodal admittance matrix, mesh impedance matrix
 *   L6 (Canonical): Thevenin/Norton equivalents, bridge circuits
 */

#ifndef CIRCUIT_ANALYSIS_H
#define CIRCUIT_ANALYSIS_H

#include "circuit_topology.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L2: Node-Voltage Method
 * ========================================================================== */

/**
 * @brief Build the node admittance matrix Y_n for nodal analysis.
 *
 * Y_n is an (n-1) x (n-1) matrix where:
 *   Y_n[i][i] = sum of admittances connected to node i (self-admittance)
 *   Y_n[i][j] = negative sum of admittances between nodes i and j (mutual)
 *
 * The node-voltage equation is: Y_n * v_n = i_n
 * where i_n is the vector of equivalent current sources entering each node.
 *
 * Limitation: Cannot directly handle floating voltage sources.
 * Use MNA for general circuits; node-voltage for purely passive/current-source networks.
 *
 * @param circuit  Circuit topology (DC domain)
 * @param Y_n      Output: node admittance matrix (n x n, row-major flat array)
 * @param n        Number of independent nodes (= num_nodes - 1)
 * @return 0 on success, -1 on error
 */
int ct_build_nodal_admittance(const ct_circuit_t *circuit,
                              double *Y_n, int32_t n);

/**
 * @brief Solve circuit using the node-voltage method.
 *
 * 1. Build Y_n and i_n
 * 2. Solve Y_n * v_n = i_n via Gaussian elimination
 * 3. Compute branch variables from node voltages
 *
 * @param circuit          Circuit topology
 * @param node_voltages    Output: node voltages (index 0 = ground = 0V)
 * @param branch_currents  Output: branch currents
 * @param branch_voltages  Output: branch voltages
 * @return 0 on success, -1 on error (e.g., singular matrix)
 */
int ct_node_voltage_solve(const ct_circuit_t *circuit,
                          double *node_voltages,
                          double *branch_currents,
                          double *branch_voltages);

/* ==========================================================================
 * L2: Mesh-Current Method
 * ========================================================================== */

/**
 * @brief Build the mesh impedance matrix Z_m for mesh analysis.
 *
 * Requires a planar circuit with identified meshes (window panes).
 * Z_m is an m x m matrix where m = number of independent meshes:
 *   Z_m[i][i] = sum of impedances in mesh i (self-impedance)
 *   Z_m[i][j] = +/- sum of impedances shared by meshes i and j
 *              (positive if mesh currents flow same direction through shared branch)
 *
 * Dual of node-voltage method. Preferred when b - n + 1 < n - 1
 * (fewer meshes than nodes).
 *
 * @param circuit  Circuit topology (must be planar)
 * @param Z_m      Output: mesh impedance matrix (m x m, row-major flat array)
 * @param m        Number of independent meshes
 * @return 0 on success, -1 if circuit is non-planar
 */
int ct_build_mesh_impedance(const ct_circuit_t *circuit,
                            double *Z_m, int32_t m);

/**
 * @brief Identify all meshes (window panes) of a planar circuit.
 *
 * Each mesh is a minimal cycle that does not enclose other branches.
 * For a planar embedding with f faces, there are f-1 interior faces
 * (the exterior face is not a mesh for analysis purposes).
 *
 * Uses the Left-Hand Rule (Maze Algorithm): traverse the planar embedding
 * keeping the "left hand on the wall" to trace each mesh boundary.
 *
 * @param circuit       Circuit topology (must be planar)
 * @param max_meshes    Maximum meshes to identify
 * @param mesh_branches Output: mesh_branches[mesh][branch] = orientation (+1/-1/0)
 * @param num_meshes    Output: number of meshes found
 * @return 0 on success, -1 on error
 */
int ct_identify_meshes(const ct_circuit_t *circuit, int32_t max_meshes,
                       int8_t *mesh_branches, int32_t *num_meshes);

/**
 * @brief Solve circuit using the mesh-current method.
 *
 * 1. Identify meshes of the planar circuit
 * 2. Build Z_m and v_m (mesh voltage source vector)
 * 3. Solve Z_m * i_m = v_m
 * 4. Compute branch variables from mesh currents
 *
 * @param circuit          Circuit topology (must be planar)
 * @param mesh_currents    Output: mesh currents [m]
 * @param branch_currents  Output: branch currents [b]
 * @param branch_voltages  Output: branch voltages [b]
 * @param num_meshes       Output: number of meshes
 * @return 0 on success, -1 on error
 */
int ct_mesh_current_solve(const ct_circuit_t *circuit,
                          double *mesh_currents,
                          double *branch_currents,
                          double *branch_voltages,
                          int32_t *num_meshes);

/* ==========================================================================
 * L2: Superposition Principle
 * ========================================================================== */

/**
 * @brief Apply the superposition principle for linear circuits.
 *
 * For a linear circuit with multiple independent sources:
 *   1. Deactivate all sources except one (V-sources → shorts, I-sources → opens)
 *   2. Solve the single-source circuit
 *   3. Repeat for each source
 *   4. Sum all partial solutions
 *
 * Superposition follows from the linearity of Maxwell's equations for
 * lumped circuits. It does NOT apply to power (which is quadratic).
 *
 * @param circuit           Circuit topology (must be linear)
 * @param node_voltages     Output: total node voltages (sum of contributions)
 * @param branch_currents   Output: total branch currents
 * @param num_sources_used  Output: number of sources superposed
 * @return 0 on success, -1 on error
 */
int ct_superposition_solve(const ct_circuit_t *circuit,
                           double *node_voltages,
                           double *branch_currents,
                           int32_t *num_sources_used);

/* ==========================================================================
 * L6: Thevenin and Norton Equivalents
 * ========================================================================== */

/**
 * @brief Compute the Thevenin equivalent circuit at a specified port.
 *
 * Thevenin's Theorem (Helmholtz, 1853; Thevenin, 1883):
 * Any linear two-terminal network can be replaced by an equivalent circuit
 * consisting of a voltage source V_th in series with a resistance R_th.
 *
 * Method:
 *   1. V_th = open-circuit voltage at the port
 *   2. R_th = V_th / I_sc (ratio of open-circuit voltage to short-circuit current)
 *      OR R_th = equivalent resistance with all independent sources deactivated
 *
 * @param circuit   Circuit topology
 * @param node_a    Port terminal A node ID
 * @param node_b    Port terminal B node ID
 * @param V_th      Output: Thevenin equivalent voltage [V]
 * @param R_th      Output: Thevenin equivalent resistance [Ohm]
 * @return 0 on success, -1 on error
 */
int ct_thevenin_equivalent(const ct_circuit_t *circuit,
                           int32_t node_a, int32_t node_b,
                           double *V_th, double *R_th);

/**
 * @brief Compute the Norton equivalent circuit at a specified port.
 *
 * Norton's Theorem (Mayer, 1926; Norton, 1926):
 * Any linear two-terminal network can be replaced by an equivalent circuit
 * consisting of a current source I_n in parallel with a resistance R_n.
 *
 * Relationship to Thevenin:
 *   I_n = V_th / R_th
 *   R_n = R_th
 *
 * @param circuit   Circuit topology
 * @param node_a    Port terminal A node ID
 * @param node_b    Port terminal B node ID
 * @param I_n       Output: Norton equivalent current [A]
 * @param R_n       Output: Norton equivalent resistance [Ohm]
 * @return 0 on success, -1 on error
 */
int ct_norton_equivalent(const ct_circuit_t *circuit,
                         int32_t node_a, int32_t node_b,
                         double *I_n, double *R_n);

/* ==========================================================================
 * L6: Maximum Power Transfer Theorem
 * ========================================================================== */

/**
 * @brief Compute the load resistance for maximum power transfer.
 *
 * Maximum Power Transfer Theorem (Jacobi, 1840):
 * Maximum power is delivered to a load when the load resistance equals
 * the Thevenin resistance of the source network:
 *   R_load = R_th  →  P_max = V_th^2 / (4 * R_th)
 *
 * For AC circuits (complex load):
 *   Z_load = Z_th*  (complex conjugate match)
 *   P_max = |V_th|^2 / (8 * Re{Z_th})
 *
 * @param V_th     Thevenin voltage magnitude
 * @param R_th     Thevenin resistance
 * @param P_max    Output: maximum deliverable power [W]
 * @param R_opt    Output: optimal load resistance [Ohm]
 * @return 0 on success, -1 on invalid parameters
 */
int ct_max_power_transfer(double V_th, double R_th,
                          double *P_max, double *R_opt);

/* ==========================================================================
 * L3: Y-Delta (Star-Mesh) Transformation
 * ========================================================================== */

/**
 * @brief Perform Y-to-Delta (star-to-mesh) transformation.
 *
 * Converts a three-terminal Y (star/T) network to an equivalent
 * Delta (pi/mesh) network. The transformation preserves terminal behavior
 * at the three external nodes.
 *
 * Y → Delta formulas (Kennelly, 1899):
 *   R_ab = (R_a*R_b + R_b*R_c + R_c*R_a) / R_c
 *   R_bc = (R_a*R_b + R_b*R_c + R_c*R_a) / R_a
 *   R_ca = (R_a*R_b + R_b*R_c + R_c*R_a) / R_b
 *
 * @param R_a, R_b, R_c  Y-connected resistances [Ohm]
 * @param R_ab, R_bc, R_ca Output: Delta-connected resistances [Ohm]
 * @return 0 on success, -1 on error
 */
int ct_y_to_delta(double R_a, double R_b, double R_c,
                  double *R_ab, double *R_bc, double *R_ca);

/**
 * @brief Perform Delta-to-Y (mesh-to-star) transformation.
 *
 * Delta → Y formulas:
 *   R_a = R_ab * R_ca / (R_ab + R_bc + R_ca)
 *   R_b = R_ab * R_bc / (R_ab + R_bc + R_ca)
 *   R_c = R_bc * R_ca / (R_ab + R_bc + R_ca)
 *
 * @param R_ab, R_bc, R_ca Delta-connected resistances [Ohm]
 * @param R_a, R_b, R_c    Output: Y-connected resistances [Ohm]
 * @return 0 on success, -1 on error
 */
int ct_delta_to_y(double R_ab, double R_bc, double R_ca,
                  double *R_a, double *R_b, double *R_c);

#ifdef __cplusplus
}
#endif

#endif /* CIRCUIT_ANALYSIS_H */
