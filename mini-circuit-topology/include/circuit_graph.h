/**
 * @file circuit_graph.h
 * @brief Graph-theoretic operations for circuit topology
 *
 * Implements graph algorithms specialized for electrical network analysis:
 * spanning trees, fundamental loops/cut-sets, graph traversal, connectivity.
 *
 * References:
 *   - S. Seshu, M.B. Reed, "Linear Graphs and Electrical Networks" (1961)
 *   - W.K. Chen, "Applied Graph Theory: Graphs and Electrical Networks" (1976)
 *   - UIUC ECE 451 / Georgia Tech ECE 6350 / ETH 227-0455
 *
 * Knowledge coverage:
 *   L2 (Concepts): Spanning tree, fundamental loops, fundamental cut-sets
 *   L3 (Math):     DFS/BFS traversal, union-find, connectivity algorithms
 *   L5 (Algorithms): Minimum spanning tree, graph search, cycle detection
 */

#ifndef CIRCUIT_GRAPH_H
#define CIRCUIT_GRAPH_H

#include "circuit_topology.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L2: Graph Search Algorithms
 * ========================================================================== */

/**
 * @brief Depth-first search (DFS) on the circuit graph.
 *
 * Visits all nodes reachable from start_node, recording discovery order.
 * Used as a building block for tree selection, connectivity checks,
 * and fundamental loop detection.
 *
 * Complexity: O(n + b) where n = number of nodes, b = number of branches.
 *
 * @param circuit       Circuit graph
 * @param start_node    Starting node ID
 * @param visited       Output: visited[node] = 1 if reached
 * @param parent        Output: parent[node] = parent in DFS tree (-1 for root)
 * @param discovery     Output: discovery[node] = order visited
 * @param visit_count   Output: number of nodes visited
 * @return 0 on success, -1 on error
 */
int ct_dfs(const ct_circuit_t *circuit, int32_t start_node,
           int8_t *visited, int32_t *parent, int32_t *discovery,
           int32_t *visit_count);

/**
 * @brief Breadth-first search (BFS) on the circuit graph.
 *
 * Visits nodes in order of increasing distance from start_node.
 * Used for finding shortest paths and level-based analysis.
 *
 * Complexity: O(n + b).
 *
 * @param circuit       Circuit graph
 * @param start_node    Starting node ID
 * @param visited       Output: visited[node] = 1 if reached
 * @param distance      Output: distance[node] = shortest path length from start
 * @param parent        Output: parent[node] = BFS tree parent (-1 for root)
 * @return 0 on success, -1 on error
 */
int ct_bfs(const ct_circuit_t *circuit, int32_t start_node,
           int8_t *visited, int32_t *distance, int32_t *parent);

/* ==========================================================================
 * L2: Tree and Co-Tree Operations
 * ========================================================================== */

/**
 * @brief Select a spanning tree prioritizing specific branch types.
 *
 * Implements a weighted DFS-based tree selection, giving priority to:
 *   1. Voltage source branches (highest priority)
 *   2. Capacitor branches
 *   3. Resistor branches
 *   4. Inductor branches
 *   5. Current source branches (lowest priority)
 *
 * This priority scheme minimizes the number of current variables in MNA.
 * Eliminates need for excessive zero-voltages source insertion.
 *
 * @param circuit     Circuit graph
 * @param tree        Output tree
 * @param priority    Priority array[b]: higher = prefer for tree
 * @return 0 on success, -1 if circuit is disconnected
 */
int ct_select_tree_weighted(const ct_circuit_t *circuit, ct_tree_t *tree,
                            const int *priority);

/**
 * @brief Enumerate ALL spanning trees of a circuit graph.
 *
 * Implements the Char ripple algorithm (Chen, 1976) for systematic
 * tree enumeration. Useful for topological formula evaluation.
 * WARNING: Number of trees grows exponentially with circuit size.
 * For a complete graph K_n, number of trees = n^(n-2) (Cayley's formula).
 *
 * @param circuit     Circuit graph
 * @param max_trees   Maximum trees to enumerate
 * @param trees       Output array of trees
 * @param num_found   Output: actual number of trees found
 * @return 0 on success, -1 if max_trees exceeded
 */
int ct_enumerate_trees(const ct_circuit_t *circuit, int32_t max_trees,
                       ct_tree_t *trees, int32_t *num_found);

/**
 * @brief Compute the tree product (sum of tree admittance products).
 *
 * For a resistive network, the determinant of the node admittance matrix
 * equals the sum of products of admittances of all trees:
 *   det(Y_n) = sum_{all trees T} prod_{b in T} Y_b
 *
 * This is the topological formula for network determinants (Maxwell, 1892).
 * Related to the Matrix-Tree Theorem (Kirchhoff, 1847).
 *
 * @param circuit     Circuit graph (all branches must have valid admittance)
 * @param tree_product Output: sum of tree products
 * @return 0 on success, -1 on error
 */
int ct_tree_product(const ct_circuit_t *circuit, double *tree_product);

/* ==========================================================================
 * L3: Cycle and Cut-Set Operations
 * ========================================================================== */

/**
 * @brief Find a fundamental cycle (loop) given a link branch.
 *
 * For a selected tree T and a link branch l not in T, the fundamental
 * cycle F(l) consists of l plus the unique path in T connecting the
 * endpoints of l.
 *
 * Each fundamental cycle forms one row of the fundamental loop matrix B.
 *
 * @param circuit    Circuit graph
 * @param tree       Selected spanning tree
 * @param link_id    ID of the link branch (must NOT be in tree)
 * @param cycle      Output: branch IDs in the fundamental cycle
 * @param cycle_len  Output: number of branches in the cycle
 * @return 0 on success, -1 on error
 */
int ct_fundamental_cycle(const ct_circuit_t *circuit, const ct_tree_t *tree,
                         int32_t link_id, int32_t *cycle, int32_t *cycle_len);

/**
 * @brief Find a fundamental cut-set given a tree branch.
 *
 * For a selected tree T and a tree branch t, the fundamental cut-set
 * C(t) is defined by removing t from T, splitting the graph into two
 * components. C(t) consists of t plus all links that connect the two
 * components.
 *
 * Each fundamental cut-set forms one row of the fundamental cut-set matrix Q.
 *
 * @param circuit    Circuit graph
 * @param tree       Selected spanning tree
 * @param tree_br_id ID of the tree branch (must be in tree)
 * @param cutset     Output: branch IDs in the fundamental cut-set
 * @param cutset_len Output: number of branches in the cut-set
 * @return 0 on success, -1 on error
 */
int ct_fundamental_cutset(const ct_circuit_t *circuit, const ct_tree_t *tree,
                          int32_t tree_br_id, int32_t *cutset,
                          int32_t *cutset_len);

/**
 * @brief Detect all independent cycles using Paton's algorithm.
 *
 * Finds a cycle basis: a set of b-n+1 independent cycles that span
 * the entire cycle space of the graph. Each branch can be assigned
 * a unique fundamental cycle.
 *
 * Complexity: O(n * b) using DFS-based back-edge detection.
 *
 * @param circuit     Circuit graph
 * @param max_cycles  Maximum number of cycles to detect
 * @param cycles      Output: cycles[cycle_idx][branch_idx] = orientation
 * @param num_cycles  Output: number of cycles found (= b-n+1 for connected graph)
 * @return 0 on success, -1 on error
 */
int ct_cycle_basis(const ct_circuit_t *circuit, int32_t max_cycles,
                   int8_t *cycles, int32_t *num_cycles);

/* ==========================================================================
 * L3: Graph Statistics and Properties
 * ========================================================================== */

/**
 * @brief Compute the graph density of the circuit.
 *
 * Graph density = 2*b / (n*(n-1)), the ratio of actual branches to
 * possible branches in a simple graph. Dense circuits (density > 0.5)
 * benefit from nodal analysis; sparse circuits from sparse matrix methods.
 *
 * @param circuit Circuit graph
 * @return Graph density [0, 1], or -1 on error
 */
double ct_graph_density(const ct_circuit_t *circuit);

/**
 * @brief Compute the cyclomatic number (nullity) of the circuit graph.
 *
 * Cyclomatic number mu = b - n + c, where c is the number of connected
 * components. This equals the number of independent loops (KVL equations).
 * For a connected graph: mu = b - n + 1.
 *
 * The cyclomatic number determines the dimension of the cycle space.
 * It appears in the Betti numbers of the graph: beta_1 = mu.
 *
 * @param circuit Circuit graph
 * @return Cyclomatic number, or -1 on error
 */
int ct_cyclomatic_number(const ct_circuit_t *circuit);

/**
 * @brief Compute the rank of the circuit graph.
 *
 * Graph rank r = n - c, where c is the number of connected components.
 * This equals the number of independent KCL equations (cut-sets).
 * For a connected graph: r = n - 1.
 *
 * The rank determines the dimension of the cut-set space.
 *
 * @param circuit Circuit graph
 * @return Graph rank, or -1 on error
 */
int ct_graph_rank(const ct_circuit_t *circuit);

/**
 * @brief Check if two branches are in series (topologically).
 *
 * Two branches are in series if they share exactly one common node,
 * and that node has degree exactly 2 (no other branches connected).
 *
 * @param circuit  Circuit graph
 * @param br1      First branch ID
 * @param br2      Second branch ID
 * @param common_node Output: the common node ID (if in series)
 * @return 1 if in series, 0 if not, -1 on error
 */
int ct_is_series(const ct_circuit_t *circuit, int32_t br1, int32_t br2,
                 int32_t *common_node);

/**
 * @brief Check if two branches are in parallel (topologically).
 *
 * Two branches are in parallel if they connect exactly the same pair of nodes
 * (both endpoints identical, regardless of orientation).
 *
 * @param circuit  Circuit graph
 * @param br1      First branch ID
 * @param br2      Second branch ID
 * @return 1 if in parallel, 0 if not, -1 on error
 */
int ct_is_parallel(const ct_circuit_t *circuit, int32_t br1, int32_t br2);

#ifdef __cplusplus
}
#endif

#endif /* CIRCUIT_GRAPH_H */
