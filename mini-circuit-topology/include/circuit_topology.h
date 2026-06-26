/**
 * @file circuit_topology.h
 * @brief Core data structures for electrical circuit topology analysis
 *
 * Defines the fundamental graph-theoretic representations used in circuit
 * analysis: nodes, branches, incidence matrices, and network graphs.
 *
 * References:
 *   - L.O. Chua, C.A. Desoer, E.S. Kuh, "Linear and Nonlinear Circuits" (1987)
 *   - N. Balabanian, T.A. Bickart, "Electrical Network Theory" (1969)
 *   - MIT 6.002 / Berkeley EE16A / Tsinghua Circuit Principles
 *
 * Knowledge coverage:
 *   L1 (Definitions): Node, Branch, Graph, Tree, Co-tree, Cut-set, Loop
 *   L2 (Concepts):   KCL, KVL, incidence matrix, fundamental loops/cutsets
 *   L3 (Math):       Matrix representations of graph topology
 */

#ifndef CIRCUIT_TOPOLOGY_H
#define CIRCUIT_TOPOLOGY_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Core Definitions -- Circuit Topology Primitives
 * ========================================================================== */

/** Maximum nodes/branches supported in a single circuit graph.
 *  Sized for educational circuits; use heap allocation for production. */
#define CT_MAX_NODES      64
#define CT_MAX_BRANCHES   128
#define CT_MAX_NAME_LEN   32

/**
 * @brief Element type enumeration for circuit branches.
 * Covers passive, active, and controlled sources per standard
 * SPICE-compatible classification.
 */
typedef enum {
    CT_ELEM_NONE       = 0,   /**< Empty / unused slot */
    CT_ELEM_RESISTOR   = 1,   /**< Linear resistor R [Ohm] */
    CT_ELEM_CAPACITOR  = 2,   /**< Linear capacitor C [Farad] */
    CT_ELEM_INDUCTOR   = 3,   /**< Linear inductor L [Henry] */
    CT_ELEM_VSOURCE    = 4,   /**< Independent voltage source V */
    CT_ELEM_ISOURCE    = 5,   /**< Independent current source I */
    CT_ELEM_VCVS       = 6,   /**< Voltage-controlled voltage source E */
    CT_ELEM_CCCS       = 7,   /**< Current-controlled current source F */
    CT_ELEM_VCCS       = 8,   /**< Voltage-controlled current source G */
    CT_ELEM_CCVS       = 9,   /**< Current-controlled voltage source H */
    CT_ELEM_DIODE      = 10,  /**< PN junction diode */
    CT_ELEM_BJT_NPN    = 11,  /**< NPN bipolar junction transistor */
    CT_ELEM_BJT_PNP    = 12,  /**< PNP bipolar junction transistor */
    CT_ELEM_MOS_NMOS   = 13,  /**< N-channel MOSFET */
    CT_ELEM_MOS_PMOS   = 14,  /**< P-channel MOSFET */
    CT_ELEM_OPAMP      = 15,  /**< Ideal operational amplifier */
    CT_ELEM_TRANSFORMER= 16,  /**< Ideal transformer */
    CT_ELEM_GYRATOR    = 17,  /**< Gyrator (non-reciprocal) */
    CT_ELEM_NULLATOR   = 18,  /**< Nullator (pathological element) */
    CT_ELEM_NORATOR    = 19,  /**< Norator (pathological element) */
    CT_ELEM_MUTUAL_IND = 20,  /**< Coupled inductors (mutual inductance M) */
    CT_ELEM_TRANSM_LINE= 21,  /**< Transmission line segment */
    CT_ELEM_VCO        = 22,  /**< Voltage-controlled oscillator */
    CT_ELEM_MEMRISTOR  = 23,  /**< Memristor (Chua, 1971) */
    CT_ELEM_COUNT
} ct_element_type_t;

/**
 * @brief Domain classification for circuit analysis mode.
 */
typedef enum {
    CT_DOMAIN_DC       = 0,   /**< DC steady-state analysis */
    CT_DOMAIN_AC       = 1,   /**< AC small-signal phasor analysis */
    CT_DOMAIN_TRANSIENT= 2,   /**< Transient (time-domain) analysis */
    CT_DOMAIN_S_PARAM  = 3,   /**< Scattering parameter analysis */
    CT_DOMAIN_NOISE    = 4,   /**< Noise analysis */
    CT_DOMAIN_SENS     = 5,   /**< Sensitivity analysis */
} ct_domain_t;

/**
 * @brief A circuit node -- a point where two or more branches connect.
 *
 * The ground/reference node is conventionally node-0.
 * Each node has a unique integer ID and an optional name for netlisting.
 */
typedef struct {
    int32_t  id;                        /**< Unique node identifier (0 = ground) */
    char     name[CT_MAX_NAME_LEN];     /**< Human-readable node name */
    uint32_t color;                     /**< Graph coloring for planar test */
    uint8_t  is_ground : 1;             /**< Flag: this is the reference node */
    uint8_t  is_terminal : 1;           /**< Flag: external terminal node */
    uint8_t  is_internal : 1;           /**< Flag: internal (non-terminal) node */
    uint8_t  reserved : 5;
    int32_t  partition;                 /**< Partition ID for hierarchical analysis */
} ct_node_t;

/**
 * @brief A circuit branch -- a two-terminal element with voltage and current.
 *
 * Each branch connects a "from" node (n+) to a "to" node (n-).
 * The reference direction convention is passive sign convention:
 * current enters n+ and leaves n-.
 */
typedef struct {
    int32_t  id;                        /**< Unique branch identifier */
    char     name[CT_MAX_NAME_LEN];     /**< Human-readable branch name */
    ct_element_type_t elem_type;        /**< Type of circuit element */
    int32_t  node_from;                 /**< Positive terminal node ID (n+) */
    int32_t  node_to;                   /**< Negative terminal node ID (n-) */
    double   value;                     /**< Primary parameter (R, C, L, gain, etc.) */
    double   value2;                    /**< Secondary parameter (e.g., mutual M, beta) */
    double   initial_v;                 /**< Initial voltage (for transient) */
    double   initial_i;                 /**< Initial current (for transient) */
    uint8_t  is_voltage_source : 1;     /**< Flag: branch is a voltage source */
    uint8_t  is_current_measured : 1;   /**< Flag: branch current is a variable */
    uint8_t  is_nonlinear : 1;          /**< Flag: element is nonlinear */
    uint8_t  is_time_varying : 1;       /**< Flag: element parameters vary with time */
    uint8_t  is_active : 1;             /**< Flag: element can deliver power */
    uint8_t  is_reciprocal : 1;         /**< Flag: element obeys reciprocity */
    uint8_t  reserved : 2;
} ct_branch_t;

/**
 * @brief Complete specification of a circuit's topology.
 */
typedef struct {
    int32_t  num_nodes;
    int32_t  num_branches;
    int32_t  num_vsrc;
    int32_t  num_isrc;
    int32_t  num_controlled;
    ct_node_t    nodes[CT_MAX_NODES];
    ct_branch_t  branches[CT_MAX_BRANCHES];
    int32_t  adj_matrix[CT_MAX_NODES][CT_MAX_NODES];
    ct_domain_t domain;
    double      frequency;
    double      temperature;
    char     title[128];
} ct_circuit_t;

/* ==========================================================================
 * L1: Graph-Theoretic Topology Structures
 * ========================================================================== */

/**
 * @brief Reduced incidence matrix A (size: (n-1) x b).
 *
 * A_ij = +1 if branch j leaves node i
 * A_ij = -1 if branch j enters node i
 * A_ij =  0 otherwise
 *
 * Encodes KCL: A * i_b = 0 (currents sum to zero at each node).
 */
typedef struct {
    int32_t  rows;
    int32_t  cols;
    int8_t   data[CT_MAX_NODES][CT_MAX_BRANCHES];
} ct_incidence_t;

/**
 * @brief Fundamental cut-set matrix Q (size: (n-1) x b).
 *
 * For a given tree, each fundamental cut-set is defined by a single tree branch.
 * Q_ij = 1 if branch j is in cut-set i with same orientation
 * Q_ij = -1 if branch j is in cut-set i with opposite orientation
 * Q_ij = 0 otherwise
 *
 * Encodes KCL: Q * i_b = 0
 */
typedef struct {
    int32_t  rows;
    int32_t  cols;
    int8_t   data[CT_MAX_NODES][CT_MAX_BRANCHES];
} ct_cutset_matrix_t;

/**
 * @brief Fundamental loop matrix B (size: (b-n+1) x b).
 *
 * For a given tree, each fundamental loop contains exactly one link (co-tree branch).
 * B_ij = 1 if branch j is in loop i with same orientation
 * B_ij = -1 if branch j is in loop i with opposite orientation
 * B_ij = 0 otherwise
 *
 * Encodes KVL: B * v_b = 0  (voltage sum around each loop is zero).
 */
typedef struct {
    int32_t  rows;
    int32_t  cols;
    int8_t   data[CT_MAX_BRANCHES][CT_MAX_BRANCHES];
} ct_loop_matrix_t;

/**
 * @brief Tree enumeration result for the circuit graph.
 *
 * A tree of a connected graph with n nodes is a connected subgraph containing
 * all n nodes and (n-1) branches without any loops.
 * The complement of a tree is called the co-tree (or link set).
 *
 * Core theorem: In any connected graph, number of links = b - n + 1.
 * Follows from Euler's formula for planar graphs: V - E + F = 2.
 */
typedef struct {
    int32_t  num_tree_branches;
    int32_t  num_links;
    int32_t  tree_branches[CT_MAX_BRANCHES];
    int32_t  link_branches[CT_MAX_BRANCHES];
    double   tree_weight;
    uint8_t  is_connected;
} ct_tree_t;

/* ==========================================================================
 * L2: Core Concept Functions -- Graph Construction and Query
 * ========================================================================== */

int ct_circuit_init(ct_circuit_t *circuit, const char *title);
int ct_add_node(ct_circuit_t *circuit, const char *name, uint8_t flags);
int ct_add_branch(ct_circuit_t *circuit, ct_element_type_t elem_type,
                  int32_t node_from, int32_t node_to,
                  double value, double value2, const char *name);
int ct_validate_connectivity(const ct_circuit_t *circuit, int32_t *error_node,
                             char *error_msg, size_t msg_len);
int ct_count_components(const ct_circuit_t *circuit);
int ct_node_degree(const ct_circuit_t *circuit, int32_t node_id);
int ct_is_planar(const ct_circuit_t *circuit);
int ct_dual_circuit(const ct_circuit_t *src, ct_circuit_t *dst);

/* L2: Incidence matrix construction */
int ct_build_incidence_matrix(const ct_circuit_t *circuit, ct_incidence_t *A);

/* L2: Tree selection -- builds a spanning tree (DFS-based) */
int ct_select_tree(const ct_circuit_t *circuit, ct_tree_t *tree,
                   int prefer_voltage_sources);

/* L2: Build fundamental cut-set matrix from a tree */
int ct_build_cutset_matrix(const ct_circuit_t *circuit,
                           const ct_tree_t *tree, ct_cutset_matrix_t *Q);

/* L2: Build fundamental loop matrix from a tree */
int ct_build_loop_matrix(const ct_circuit_t *circuit,
                         const ct_tree_t *tree, ct_loop_matrix_t *B);

/* L3: Matrix operations for topology */
int ct_incidence_to_adjacency(const ct_incidence_t *A, int32_t n_nodes,
                              int32_t adj_out[CT_MAX_NODES][CT_MAX_NODES]);
int ct_verify_kcl(const ct_circuit_t *circuit, const double *branch_currents,
                  double tolerance, int32_t *violated_node);
int ct_verify_kvl(const ct_circuit_t *circuit, const double *branch_voltages,
                  double tolerance, int32_t *violated_loop);

/* L4: Tellegen's Theorem verification */
int ct_verify_tellegen(const ct_circuit_t *circuit,
                       const double *branch_voltages,
                       const double *branch_currents,
                       double tolerance);

/* L4: Tellegen's Theorem — cross-circuit (strong form) */
int ct_verify_tellegen_cross(const ct_circuit_t *circuit,
                             const double *v1, const double *i1,
                             const double *v2, const double *i2,
                             double *sum_v1i2, double *sum_v2i1,
                             double tolerance);

/* L4: Power balance (energy conservation) */
int ct_verify_power_balance(const ct_circuit_t *circuit,
                            const double *v, const double *i,
                            double *p_supplied, double *p_dissipated,
                            double tolerance);

/* L4: Reciprocity check */
int ct_check_reciprocity(const ct_circuit_t *circuit,
                         double z12, double z21, double tolerance);

#ifdef __cplusplus
}
#endif

#endif /* CIRCUIT_TOPOLOGY_H */
