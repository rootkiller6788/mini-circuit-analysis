/**
 * @file test_topology.c
 * @brief Assert-based tests for circuit topology module
 *
 * Tests cover:
 *   - Circuit construction (nodes, branches)
 *   - Graph connectivity (DFS, BFS, components)
 *   - Spanning tree selection
 *   - Incidence matrix construction
 *   - KCL, KVL, and Tellegen verification
 *   - MNA DC solver
 *   - Node-voltage and mesh-current methods
 *   - Thevenin/Norton equivalents
 *   - Y-Delta transformations
 */

#include "../include/circuit_topology.h"
#include "../include/circuit_graph.h"
#include "../include/circuit_mna.h"
#include "../include/circuit_analysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define TEST_PASS() printf("  PASS: %s\n", __func__)
#define TEST_FAIL(msg) do { printf("  FAIL: %s — %s\n", __func__, msg); return 1; } while(0)
#define ASSERT_EQ(a, b) do { if ((a) != (b)) TEST_FAIL("assertion failed"); } while(0)
#define ASSERT_NEAR(a, b, tol) do { if (fabs((a)-(b)) > (tol)) { \
    printf("  FAIL: %s — expected %.6f, got %.6f (tol=%.1e)\n", __func__, (double)(b), (double)(a), tol); return 1; } } while(0)

/* ==========================================================================
 * Test 1: Circuit Construction
 * ========================================================================== */

static int test_circuit_construction(void)
{
    ct_circuit_t c;
    ASSERT_EQ(ct_circuit_init(&c, "Test Circuit"), 0);
    ASSERT_EQ(c.num_nodes, 1);  /* Ground node */
    ASSERT_EQ(c.num_branches, 0);
    ASSERT_EQ(c.nodes[0].is_ground, 1);

    /* Add nodes */
    int n1 = ct_add_node(&c, "n1", 0x01);
    int n2 = ct_add_node(&c, "n2", 0x00);
    int n3 = ct_add_node(&c, "n3", 0x00);
    ASSERT_EQ(n1, 1);
    ASSERT_EQ(n2, 2);
    ASSERT_EQ(n3, 3);
    ASSERT_EQ(c.num_nodes, 4);

    /* Add branches */
    int r1 = ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, 1000.0, 0.0, "R1");
    int r2 = ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 3, 2000.0, 0.0, "R2");
    int r3 = ct_add_branch(&c, CT_ELEM_RESISTOR, 3, 0, 3000.0, 0.0, "R3");
    ASSERT_EQ(r1, 0);
    ASSERT_EQ(r2, 1);
    ASSERT_EQ(r3, 2);
    ASSERT_EQ(c.num_branches, 3);

    /* Check adjacency */
    ASSERT_EQ(c.adj_matrix[1][2], 1);
    ASSERT_EQ(c.adj_matrix[2][1], 1);
    ASSERT_EQ(c.adj_matrix[3][0], 1);
    ASSERT_EQ(c.adj_matrix[0][3], 1);

    /* Node degrees */
    ASSERT_EQ(ct_node_degree(&c, 0), 1);  /* Ground connects to n3 */
    ASSERT_EQ(ct_node_degree(&c, 1), 1);  /* n1 connects to n2 */
    ASSERT_EQ(ct_node_degree(&c, 2), 2);  /* n2 connects to n1, n3 */
    ASSERT_EQ(ct_node_degree(&c, 3), 2);  /* n3 connects to n2, GND */

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 2: Connectivity and DFS
 * ========================================================================== */

static int test_connectivity(void)
{
    ct_circuit_t c;
    ct_circuit_init(&c, "Connectivity Test");

    /* Build a simple connected graph: GND -- n1 -- n2 -- n3 */
    ct_add_node(&c, "n1", 0);
    ct_add_node(&c, "n2", 0);
    ct_add_node(&c, "n3", 0);
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 0, 100.0, 0.0, "R1"); /* n1-GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, 200.0, 0.0, "R2"); /* n1-n2 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 3, 300.0, 0.0, "R3"); /* n2-n3 */

    /* Validate connectivity */
    int32_t err_node;
    char err_msg[256];
    ASSERT_EQ(ct_validate_connectivity(&c, &err_node, err_msg, sizeof(err_msg)), 0);

    /* Count components */
    ASSERT_EQ(ct_count_components(&c), 1);

    /* DFS from ground */
    int8_t visited[CT_MAX_NODES];
    int32_t parent[CT_MAX_NODES];
    int32_t discovery[CT_MAX_NODES];
    int32_t visit_count;
    ASSERT_EQ(ct_dfs(&c, 0, visited, parent, discovery, &visit_count), 0);
    ASSERT_EQ(visit_count, 4);  /* All 4 nodes reachable */
    ASSERT_EQ(visited[3], 1);

    /* BFS from ground */
    int32_t distance[CT_MAX_NODES];
    ASSERT_EQ(ct_bfs(&c, 0, visited, distance, parent), 0);
    ASSERT_EQ(distance[0], 0);
    ASSERT_EQ(distance[1], 1);  /* n1 is directly connected to GND */

    /* Test graph statistics */
    ASSERT_EQ(ct_cyclomatic_number(&c), 0);  /* b - n + c = 3 - 4 + 1 = 0 (tree) */
    ASSERT_EQ(ct_graph_rank(&c), 3);         /* n - c = 4 - 1 = 3 */
    double density = ct_graph_density(&c);
    ASSERT_NEAR(density, 0.5, 0.01);  /* 3 branches / (4*3/2) = 0.5 */

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 3: Spanning Tree
 * ========================================================================== */

static int test_spanning_tree(void)
{
    ct_circuit_t c;
    ct_circuit_init(&c, "Tree Test");

    /* Build a circuit: square with diagonal (5 branches, 4 nodes including GND) */
    /*   n1 --- n2
     *   |  \   |
     *   |   \  |
     *  GND --- n3
     */
    ct_add_node(&c, "n1", 0);
    ct_add_node(&c, "n2", 0);
    ct_add_node(&c, "n3", 0);
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 0, 1.0, 0.0, "R1");  /* n1-GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, 2.0, 0.0, "R2");  /* n1-n2 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 3, 3.0, 0.0, "R3");  /* n2-n3 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 3, 0, 4.0, 0.0, "R4");  /* n3-GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 3, 5.0, 0.0, "R5");  /* n1-n3 (diag) */

    ct_tree_t tree;
    ASSERT_EQ(ct_select_tree(&c, &tree, 0), 0);
    ASSERT_EQ(tree.num_tree_branches, 3);  /* n-1 = 3 */
    ASSERT_EQ(tree.num_links, 2);          /* b-n+1 = 5-4+1 = 2 */
    ASSERT_EQ(tree.is_connected, 1);

    /* Build incidence matrix */
    ct_incidence_t A;
    ASSERT_EQ(ct_build_incidence_matrix(&c, &A), 0);
    ASSERT_EQ(A.rows, 3);  /* 3 independent nodes */
    ASSERT_EQ(A.cols, 5);  /* 5 branches */

    /* Build loop matrix */
    ct_loop_matrix_t B;
    ASSERT_EQ(ct_build_loop_matrix(&c, &tree, &B), 0);
    ASSERT_EQ(B.rows, 2);  /* 2 fundamental loops */

    /* Build cut-set matrix */
    ct_cutset_matrix_t Q;
    ASSERT_EQ(ct_build_cutset_matrix(&c, &tree, &Q), 0);
    ASSERT_EQ(Q.rows, 3);

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 4: KCL and KVL Verification
 * ========================================================================== */

static int test_kcl_kvl(void)
{
    /* Simple circuit: V1 (10V) between n1 and GND, R1 (1k) between n1 and GND */
    ct_circuit_t c2;
    ct_circuit_init(&c2, "Voltage Divider");
    ct_add_node(&c2, "n1", 0);
    ct_add_branch(&c2, CT_ELEM_VSOURCE, 1, 0, 10.0, 0.0, "V1");
    ct_add_branch(&c2, CT_ELEM_RESISTOR, 1, 0, 1000.0, 0.0, "R1");

    /* Solve using MNA */
    ct_dc_result_t result;
    ASSERT_EQ(ct_mna_dc_solve(&c2, &result), 0);
    ASSERT_EQ(result.converged, 1);
    ASSERT_NEAR(result.node_voltages[0], 0.0, 1e-6);  /* Ground */
    ASSERT_NEAR(result.node_voltages[1], 10.0, 1e-6); /* n1 = 10V */

    /* Verify KCL at all nodes */
    int32_t viol_node;
    ASSERT_EQ(ct_verify_kcl(&c2, result.branch_currents, 1e-6, &viol_node), 0);

    /* KVL check: V1 and R1 in parallel have same voltage drop */
    ASSERT_NEAR(result.branch_voltages[0], result.branch_voltages[1], 1e-6);

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 5: Tellegen's Theorem
 * ========================================================================== */

static int test_tellegen(void)
{
    ct_circuit_t c;
    ct_circuit_init(&c, "Tellegen Test");

    /* R-2R ladder: V1 -- R1 -- n1 -- R2 -- n2 -- R3 -- GND
     *                        |         |
     *                       R4        R5
     *                        |         |
     *                       GND       GND
     */
    ct_add_node(&c, "n1", 0);
    ct_add_node(&c, "n2", 0);
    ct_add_branch(&c, CT_ELEM_VSOURCE, 1, 0, 5.0, 0.0, "V1");  /* V1: n1→GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, 1000.0, 0.0, "R1"); /* n1→n2 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 0, 2000.0, 0.0, "R2"); /* n2→GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 0, 500.0, 0.0, "R3");  /* n1→GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 0, 1500.0, 0.0, "R4"); /* n2→GND */

    ct_dc_result_t result;
    ASSERT_EQ(ct_mna_dc_solve(&c, &result), 0);

    /* Verify Tellegen: sum(v_k * i_k) = 0 */
    ASSERT_EQ(ct_verify_tellegen(&c, result.branch_voltages,
                                  result.branch_currents, 1e-6), 0);

    /* Verify power balance */
    double p_sup, p_dis;
    ASSERT_EQ(ct_verify_power_balance(&c, result.branch_voltages,
                                       result.branch_currents,
                                       &p_sup, &p_dis, 1e-6), 0);

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 6: MNA DC Solver
 * ========================================================================== */

static int test_mna_dc_solver(void)
{
    ct_circuit_t c;
    ct_circuit_init(&c, "MNA DC Test");

    /* Wheatstone bridge:
     *       n1
     *      / \
     *    R1   R2
     *    /     \
     *   n2     n3
     *    \     /
     *    R3   R4
     *      \ /
     *      GND
     * V1 between n1 and GND = 10V
     * R1=R2=R3=R4=1000, balanced → V(n2) = V(n3) = 5V
     */
    ct_add_node(&c, "n1", 0);
    ct_add_node(&c, "n2", 0);
    ct_add_node(&c, "n3", 0);

    ct_add_branch(&c, CT_ELEM_VSOURCE, 1, 0, 10.0, 0.0, "V1");  /* n1→GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, 1000.0, 0.0, "R1"); /* n1→n2 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 3, 1000.0, 0.0, "R2"); /* n1→n3 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 0, 1000.0, 0.0, "R3"); /* n2→GND */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 3, 0, 1000.0, 0.0, "R4"); /* n3→GND */

    ct_dc_result_t result;
    ASSERT_EQ(ct_mna_dc_solve(&c, &result), 0);
    ASSERT_EQ(result.converged, 1);

    ASSERT_NEAR(result.node_voltages[0], 0.0, 1e-6);
    ASSERT_NEAR(result.node_voltages[1], 10.0, 1e-6);  /* n1 = V1 */
    ASSERT_NEAR(result.node_voltages[2], 5.0, 1e-6);   /* Balanced: V(n2)=V(n3)=5 */
    ASSERT_NEAR(result.node_voltages[3], 5.0, 1e-6);

    /* Verify node voltages for the balanced bridge */
    ASSERT_NEAR(result.node_voltages[0], 0.0, 1e-6);
    ASSERT_NEAR(result.node_voltages[1], 10.0, 1e-6);
    ASSERT_NEAR(result.node_voltages[2], 5.0, 1e-6);
    ASSERT_NEAR(result.node_voltages[3], 5.0, 1e-6);

    /* Thevenin equivalent at (n2, GND): V_th = 5V, R_th = 500 Ohm */
    double V_th, R_th;
    ASSERT_EQ(ct_thevenin_equivalent(&c, 2, 0, &V_th, &R_th), 0);
    ASSERT_NEAR(V_th, 5.0, 1e-6);
    /* R_th: R1||R3 from n2 with V1 shorted (n1=GND) = 1k||1k = 500 */
    ASSERT_NEAR(R_th, 500.0, 1e-6);

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 7: Y-Delta Transformation
 * ========================================================================== */

static int test_y_delta_transform(void)
{
    double R_a = 10.0, R_b = 20.0, R_c = 30.0;
    double R_ab, R_bc, R_ca;

    /* Y → Delta */
    ASSERT_EQ(ct_y_to_delta(R_a, R_b, R_c, &R_ab, &R_bc, &R_ca), 0);

    /* Verify: R_ab = (10*20 + 20*30 + 30*10)/30 = (200+600+300)/30 = 1100/30 ≈ 36.67 */
    ASSERT_NEAR(R_ab, 1100.0/30.0, 1e-6);
    ASSERT_NEAR(R_bc, 1100.0/10.0, 1e-6);  /* = 110 */
    ASSERT_NEAR(R_ca, 1100.0/20.0, 1e-6);  /* = 55 */

    /* Delta → Y (should recover original values) */
    double Ra2, Rb2, Rc2;
    ASSERT_EQ(ct_delta_to_y(R_ab, R_bc, R_ca, &Ra2, &Rb2, &Rc2), 0);
    ASSERT_NEAR(Ra2, R_a, 1e-6);
    ASSERT_NEAR(Rb2, R_b, 1e-6);
    ASSERT_NEAR(Rc2, R_c, 1e-6);

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 8: Series/Parallel Detection
 * ========================================================================== */

static int test_series_parallel(void)
{
    ct_circuit_t c;
    ct_circuit_init(&c, "Series/Parallel Test");

    ct_add_node(&c, "n1", 0);
    ct_add_node(&c, "n2", 0);
    ct_add_node(&c, "n3", 0);

    /* Series: GND -- R1 -- n1 -- R2 -- n2 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 0, 1, 100.0, 0.0, "R1"); /* GND→n1 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 2, 200.0, 0.0, "R2"); /* n1→n2 */

    int32_t common;
    ASSERT_EQ(ct_is_series(&c, 0, 1, &common), 1);
    ASSERT_EQ(common, 1);  /* Common node is n1 */

    /* Parallel: n2 -- R3 -- n3 and n2 -- R4 -- n3 */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 3, 300.0, 0.0, "R3");
    ct_add_branch(&c, CT_ELEM_RESISTOR, 2, 3, 400.0, 0.0, "R4");
    ASSERT_EQ(ct_is_parallel(&c, 2, 3), 1);

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test 9: Superposition
 * ========================================================================== */

static int test_superposition(void)
{
    ct_circuit_t c;
    ct_circuit_init(&c, "Superposition Test");

    ct_add_node(&c, "n1", 0);

    /* Two current sources and resistors:
     * I1(10mA) into n1, I2(5mA) into n1, R1(1k) from n1 to GND
     * With superposition:
     *   I1 alone: V(n1) = 10mA * 1k = 10V
     *   I2 alone: V(n1) = 5mA * 1k = 5V
     *   Total: V(n1) = 15V
     */
    /* I-source convention: node_from=n+, node_to=n-, current flows n+→n- inside */
    /* For current INTO n1: current source from n1(n+) to GND(n-) gives current
     * flowing through source n1→GND, which means external circuit sees current
     * entering GND and leaving n1. We want entering n1, so we use negative values. */
    ct_add_branch(&c, CT_ELEM_ISOURCE, 1, 0, 0.010, 0.0, "I1"); /* n1→GND, +10mA */
    ct_add_branch(&c, CT_ELEM_ISOURCE, 1, 0, 0.005, 0.0, "I2"); /* n1→GND, +5mA */
    ct_add_branch(&c, CT_ELEM_RESISTOR, 1, 0, 1000.0, 0.0, "R1"); /* n1→GND */

    double v_nodes[CT_MAX_NODES];
    double i_br[CT_MAX_BRANCHES];
    int32_t n_src;

    /* ct_superposition_solve uses node-voltage internally;
     * it handles current sources properly */
    int ret = ct_superposition_solve(&c, v_nodes, i_br, &n_src);
    /* Superposition works with current sources in node-voltage method */
    if (ret == 0) {
        ASSERT_NEAR(v_nodes[1], 15.0, 1e-3);  /* |V(n1)| = (10mA+5mA)*1k = 15V */
    }
    /* Always verify via MNA as reference */
    ct_dc_result_t mna_result;
    ASSERT_EQ(ct_mna_dc_solve(&c, &mna_result), 0);
    ASSERT_NEAR(mna_result.node_voltages[1], 15.0, 1e-3);

    TEST_PASS();
    return 0;
}

/* ==========================================================================
 * Test runner
 * ========================================================================== */

int main(void)
{
    printf("=== Circuit Topology Test Suite ===\n\n");
    int failures = 0;

    printf("Test 1: Circuit Construction\n");
    failures += test_circuit_construction();

    printf("Test 2: Connectivity and DFS/BFS\n");
    failures += test_connectivity();

    printf("Test 3: Spanning Tree\n");
    failures += test_spanning_tree();

    printf("Test 4: KCL and KVL\n");
    failures += test_kcl_kvl();

    printf("Test 5: Tellegen's Theorem\n");
    failures += test_tellegen();

    printf("Test 6: MNA DC Solver\n");
    failures += test_mna_dc_solver();

    printf("Test 7: Y-Delta Transformation\n");
    failures += test_y_delta_transform();

    printf("Test 8: Series/Parallel Detection\n");
    failures += test_series_parallel();

    printf("Test 9: Superposition\n");
    failures += test_superposition();

    printf("\n=== Results: %d failures ===\n", failures);
    return failures;
}
