/**
 * @file circuit_graph.c
 * @brief Graph-theoretic operations for circuit topology — Implementation
 *
 * Implements DFS/BFS traversal, spanning tree selection, fundamental
 * cycle and cut-set identification, graph connectivity testing, and
 * graph-theoretic metrics (cyclomatic number, rank, density).
 *
 * All functions operate on the ct_circuit_t adjacency representation
 * for O(1) neighbor queries.
 */

#include "circuit_topology.h"
#include "circuit_graph.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* ==========================================================================
 * L2: Graph Construction Helpers
 * ========================================================================== */

int ct_circuit_init(ct_circuit_t *circuit, const char *title)
{
    if (!circuit) return -1;
    memset(circuit, 0, sizeof(ct_circuit_t));

    /* Node 0 is always ground */
    circuit->nodes[0].id = 0;
    circuit->nodes[0].is_ground = 1;
    snprintf(circuit->nodes[0].name, CT_MAX_NAME_LEN, "GND");

    if (title) {
        snprintf(circuit->title, sizeof(circuit->title), "%s", title);
    } else {
        snprintf(circuit->title, sizeof(circuit->title), "Untitled Circuit");
    }

    circuit->temperature = 300.0;  /* Room temperature in Kelvin */
    circuit->frequency = 0.0;
    circuit->domain = CT_DOMAIN_DC;
    circuit->num_nodes = 1;  /* Ground node accounts for 1 */

    return 0;
}

int ct_add_node(ct_circuit_t *circuit, const char *name, uint8_t flags)
{
    if (!circuit) return -1;
    if (circuit->num_nodes >= CT_MAX_NODES) return -1;

    int32_t id = circuit->num_nodes;
    ct_node_t *node = &circuit->nodes[id];
    memset(node, 0, sizeof(ct_node_t));

    node->id = id;
    if (name) {
        snprintf(node->name, CT_MAX_NAME_LEN, "%s", name);
    } else {
        snprintf(node->name, CT_MAX_NAME_LEN, "N%d", id);
    }
    node->is_terminal = (flags & 0x01) ? 1 : 0;
    node->is_internal = (flags & 0x02) ? 1 : 0;

    /* Initialize adjacency row/col for this node */
    for (int i = 0; i < CT_MAX_NODES; i++) {
        circuit->adj_matrix[id][i] = 0;
        circuit->adj_matrix[i][id] = 0;
    }

    circuit->num_nodes++;
    return id;
}

int ct_add_branch(ct_circuit_t *circuit, ct_element_type_t elem_type,
                  int32_t node_from, int32_t node_to,
                  double value, double value2, const char *name)
{
    if (!circuit) return -1;
    if (circuit->num_branches >= CT_MAX_BRANCHES) return -1;
    if (node_from < 0 || node_from >= circuit->num_nodes) return -1;
    if (node_to < 0 || node_to >= circuit->num_nodes) return -1;

    int32_t id = circuit->num_branches;
    ct_branch_t *br = &circuit->branches[id];
    memset(br, 0, sizeof(ct_branch_t));

    br->id = id;
    br->elem_type = elem_type;
    br->node_from = node_from;
    br->node_to = node_to;
    br->value = value;
    br->value2 = value2;

    if (name) {
        snprintf(br->name, CT_MAX_NAME_LEN, "%s", name);
    } else {
        snprintf(br->name, CT_MAX_NAME_LEN, "B%d", id);
    }

    /* Set flags based on element type */
    switch (elem_type) {
    case CT_ELEM_VSOURCE:
        br->is_voltage_source = 1;
        br->is_active = 1;
        circuit->num_vsrc++;
        break;
    case CT_ELEM_ISOURCE:
        br->is_active = 1;
        circuit->num_isrc++;
        break;
    case CT_ELEM_VCVS:
    case CT_ELEM_CCCS:
    case CT_ELEM_VCCS:
    case CT_ELEM_CCVS:
        br->is_active = 1;
        circuit->num_controlled++;
        break;
    case CT_ELEM_OPAMP:
        br->is_active = 1;
        br->is_reciprocal = 0;
        circuit->num_controlled++;
        break;
    case CT_ELEM_GYRATOR:
        br->is_active = 1;
        br->is_reciprocal = 0;
        break;
    case CT_ELEM_NULLATOR:
    case CT_ELEM_NORATOR:
        br->is_active = 1;
        br->is_reciprocal = 0;
        break;
    case CT_ELEM_DIODE:
    case CT_ELEM_BJT_NPN:
    case CT_ELEM_BJT_PNP:
    case CT_ELEM_MOS_NMOS:
    case CT_ELEM_MOS_PMOS:
        br->is_nonlinear = 1;
        break;
    default:
        /* Passive elements: R, C, L are reciprocal */
        br->is_reciprocal = 1;
        break;
    }

    /* Update adjacency matrix */
    circuit->adj_matrix[node_from][node_to]++;
    if (node_from != node_to) {
        circuit->adj_matrix[node_to][node_from]++;
    }

    circuit->num_branches++;
    return id;
}

int ct_node_degree(const ct_circuit_t *circuit, int32_t node_id)
{
    if (!circuit) return -1;
    if (node_id < 0 || node_id >= circuit->num_nodes) return -1;

    int degree = 0;
    for (int i = 0; i < circuit->num_nodes; i++) {
        degree += circuit->adj_matrix[node_id][i];
    }
    return degree;
}

/* ==========================================================================
 * L2: DFS — Depth-First Search
 * ========================================================================== */

int ct_dfs(const ct_circuit_t *circuit, int32_t start_node,
           int8_t *visited, int32_t *parent, int32_t *discovery,
           int32_t *visit_count)
{
    if (!circuit || !visited || !parent || !discovery || !visit_count)
        return -1;
    if (start_node < 0 || start_node >= circuit->num_nodes)
        return -1;

    /* Initialize arrays */
    int n = circuit->num_nodes;
    for (int i = 0; i < n; i++) {
        visited[i] = 0;
        parent[i] = -1;
        discovery[i] = -1;
    }

    /* Manual stack for iterative DFS (avoids recursion depth limits) */
    int32_t stack[CT_MAX_NODES];
    int32_t stack_top = 0;
    int32_t order = 0;

    visited[start_node] = 1;
    discovery[start_node] = order++;
    parent[start_node] = -1;
    stack[stack_top++] = start_node;

    while (stack_top > 0) {
        int32_t u = stack[--stack_top];

        /* Explore neighbors of u through adjacency matrix */
        for (int32_t v = 0; v < n; v++) {
            if (circuit->adj_matrix[u][v] > 0 && !visited[v]) {
                visited[v] = 1;
                parent[v] = u;
                discovery[v] = order++;
                stack[stack_top++] = v;
            }
        }
    }

    *visit_count = order;
    return 0;
}

/* ==========================================================================
 * L2: BFS — Breadth-First Search
 * ========================================================================== */

int ct_bfs(const ct_circuit_t *circuit, int32_t start_node,
           int8_t *visited, int32_t *distance, int32_t *parent)
{
    if (!circuit || !visited || !distance || !parent)
        return -1;
    if (start_node < 0 || start_node >= circuit->num_nodes)
        return -1;

    int n = circuit->num_nodes;
    for (int i = 0; i < n; i++) {
        visited[i] = 0;
        distance[i] = -1;
        parent[i] = -1;
    }

    /* Queue for BFS — simple ring buffer */
    int32_t queue[CT_MAX_NODES];
    int32_t q_head = 0, q_tail = 0;

    visited[start_node] = 1;
    distance[start_node] = 0;
    parent[start_node] = -1;
    queue[q_tail++] = start_node;

    while (q_head < q_tail) {
        int32_t u = queue[q_head++];

        for (int32_t v = 0; v < n; v++) {
            if (circuit->adj_matrix[u][v] > 0 && !visited[v]) {
                visited[v] = 1;
                distance[v] = distance[u] + 1;
                parent[v] = u;
                queue[q_tail++] = v;
            }
        }
    }

    return 0;
}

/* ==========================================================================
 * L2: Connectivity Validation
 * ========================================================================== */

int ct_validate_connectivity(const ct_circuit_t *circuit, int32_t *error_node,
                             char *error_msg, size_t msg_len)
{
    if (!circuit) return -1;

    int8_t visited[CT_MAX_NODES];
    int32_t parent[CT_MAX_NODES];
    int32_t discovery[CT_MAX_NODES];
    int32_t visit_count;

    /* Check each node has at least one connection (no floating nodes) */
    for (int32_t i = 1; i < circuit->num_nodes; i++) {
        int deg = ct_node_degree(circuit, i);
        if (deg == 0) {
            if (error_node) *error_node = i;
            if (error_msg) {
                snprintf(error_msg, msg_len,
                         "Floating node %d ('%s'): no connections",
                         i, circuit->nodes[i].name);
            }
            return -1;
        }
    }

    /* Check all non-ground nodes are reachable from ground */
    ct_dfs(circuit, 0, visited, parent, discovery, &visit_count);

    for (int32_t i = 1; i < circuit->num_nodes; i++) {
        if (!visited[i]) {
            if (error_node) *error_node = i;
            if (error_msg) {
                snprintf(error_msg, msg_len,
                         "Disconnected subgraph at node %d ('%s')",
                         i, circuit->nodes[i].name);
            }
            return -2;
        }
    }

    return 0;
}

int ct_count_components(const ct_circuit_t *circuit)
{
    if (!circuit) return -1;

    int8_t visited[CT_MAX_NODES];
    int32_t parent[CT_MAX_NODES];
    int32_t discovery[CT_MAX_NODES];
    int32_t visit_count;
    int components = 0;

    /* Initialize visited array */
    for (int i = 0; i < circuit->num_nodes; i++) {
        visited[i] = 0;
    }

    /* Run DFS from each unvisited node */
    for (int32_t i = 0; i < circuit->num_nodes; i++) {
        if (!visited[i]) {
            ct_dfs(circuit, i, visited, parent, discovery, &visit_count);
            components++;
        }
    }

    return components;
}

/* ==========================================================================
 * L2: Spanning Tree Selection
 * ========================================================================== */

int ct_select_tree(const ct_circuit_t *circuit, ct_tree_t *tree,
                   int prefer_voltage_sources)
{
    if (!circuit || !tree) return -1;

    /* Build priority array */
    int priority[CT_MAX_BRANCHES];
    for (int32_t i = 0; i < circuit->num_branches; i++) {
        ct_element_type_t et = circuit->branches[i].elem_type;
        if (prefer_voltage_sources) {
            if (et == CT_ELEM_VSOURCE) {
                priority[i] = 100;  /* Highest priority */
            } else if (et == CT_ELEM_CAPACITOR) {
                priority[i] = 50;
            } else if (et == CT_ELEM_RESISTOR) {
                priority[i] = 30;
            } else if (et == CT_ELEM_INDUCTOR) {
                priority[i] = 20;
            } else {
                priority[i] = 10;  /* Current sources last */
            }
        } else {
            priority[i] = 1;  /* All equal priority */
        }
    }

    return ct_select_tree_weighted(circuit, tree, priority);
}

int ct_select_tree_weighted(const ct_circuit_t *circuit, ct_tree_t *tree,
                            const int *priority)
{
    if (!circuit || !tree || !priority) return -1;

    int32_t n = circuit->num_nodes;
    int32_t b = circuit->num_branches;

    memset(tree, 0, sizeof(ct_tree_t));

    /* Union-Find data structure for cycle detection */
    int32_t parent[CT_MAX_NODES];
    int32_t rank[CT_MAX_NODES];
    for (int32_t i = 0; i < n; i++) {
        parent[i] = i;
        rank[i] = 0;
    }

    /* Find with path compression */
    /* Using a local helper via inline loop to avoid function call overhead */
    #define FIND(x) \
        ({ int32_t __x = (x); \
           while (parent[__x] != __x) { \
               parent[__x] = parent[parent[__x]]; \
               __x = parent[__x]; \
           } \
           __x; })

    /* Union by rank */
    #define UNION(x, y) \
        do { \
            int32_t __rx = FIND(x); \
            int32_t __ry = FIND(y); \
            if (__rx != __ry) { \
                if (rank[__rx] < rank[__ry]) { \
                    parent[__rx] = __ry; \
                } else if (rank[__rx] > rank[__ry]) { \
                    parent[__ry] = __rx; \
                } else { \
                    parent[__ry] = __rx; \
                    rank[__rx]++; \
                } \
            } \
        } while(0)

    /* Sort branches by priority (simple selection sort for clarity) */
    int32_t branch_order[CT_MAX_BRANCHES];
    for (int32_t i = 0; i < b; i++) {
        branch_order[i] = i;
    }

    /* Insertion sort by priority (stable, O(b^2) but b is typically small) */
    for (int32_t i = 1; i < b; i++) {
        int32_t key = branch_order[i];
        int key_prio = priority[key];
        int32_t j = i - 1;
        while (j >= 0 && priority[branch_order[j]] < key_prio) {
            branch_order[j + 1] = branch_order[j];
            j--;
        }
        branch_order[j + 1] = key;
    }

    /* Kruskal-like tree construction with priority order */
    int32_t tree_count = 0;
    int32_t link_count = 0;

    for (int32_t idx = 0; idx < b; idx++) {
        int32_t br_id = branch_order[idx];
        const ct_branch_t *br = &circuit->branches[br_id];
        int32_t u = br->node_from;
        int32_t v = br->node_to;

        if (FIND(u) != FIND(v)) {
            UNION(u, v);
            tree->tree_branches[tree_count++] = br_id;
        } else {
            tree->link_branches[link_count++] = br_id;
        }
    }

    tree->num_tree_branches = tree_count;
    tree->num_links = link_count;
    tree->is_connected = (tree_count == n - 1) ? 1 : 0;

    #undef FIND
    #undef UNION

    return tree->is_connected ? 0 : -1;
}

/* ==========================================================================
 * L3: Fundamental Cycle (Loop) Detection
 * ========================================================================== */

int ct_fundamental_cycle(const ct_circuit_t *circuit, const ct_tree_t *tree,
                         int32_t link_id, int32_t *cycle, int32_t *cycle_len)
{
    if (!circuit || !tree || !cycle || !cycle_len) return -1;

    const ct_branch_t *link = &circuit->branches[link_id];
    int32_t u = link->node_from;
    int32_t v = link->node_to;

    /* Use BFS on tree branches to find path from u to v */
    int8_t visited[CT_MAX_NODES];
    int32_t parent[CT_MAX_NODES];
    int32_t tree_adj[CT_MAX_NODES][CT_MAX_NODES];
    memset(tree_adj, 0, sizeof(tree_adj));

    /* Build tree adjacency (undirected) */
    for (int32_t i = 0; i < tree->num_tree_branches; i++) {
        int32_t br = tree->tree_branches[i];
        int32_t nu = circuit->branches[br].node_from;
        int32_t nv = circuit->branches[br].node_to;
        tree_adj[nu][nv] = 1;
        tree_adj[nv][nu] = 1;
    }

    /* BFS on tree */
    for (int32_t i = 0; i < circuit->num_nodes; i++) {
        visited[i] = 0;
        parent[i] = -1;
    }

    int32_t queue[CT_MAX_NODES];
    int32_t head = 0, tail = 0;
    visited[u] = 1;
    queue[tail++] = u;

    while (head < tail) {
        int32_t cur = queue[head++];
        if (cur == v) break;
        for (int32_t w = 0; w < circuit->num_nodes; w++) {
            if (tree_adj[cur][w] && !visited[w]) {
                visited[w] = 1;
                parent[w] = cur;
                queue[tail++] = w;
            }
        }
    }

    /* Trace path from v back to u */
    int32_t path[CT_MAX_NODES];
    int32_t path_len = 0;
    int32_t cur = v;
    while (cur != -1 && cur != u) {
        path[path_len++] = cur;
        cur = parent[cur];
    }
    path[path_len++] = u;

    /* Convert node path to branch IDs */
    *cycle_len = 0;
    for (int32_t i = 0; i < path_len - 1; i++) {
        int32_t na = path[i];
        int32_t nb = path[i + 1];
        /* Find the tree branch connecting na and nb */
        int found = 0;
        for (int32_t j = 0; j < tree->num_tree_branches; j++) {
            int32_t tbr = tree->tree_branches[j];
            if ((circuit->branches[tbr].node_from == na &&
                 circuit->branches[tbr].node_to == nb) ||
                (circuit->branches[tbr].node_from == nb &&
                 circuit->branches[tbr].node_to == na)) {
                cycle[(*cycle_len)++] = tbr;
                found = 1;
                break;
            }
        }
        if (!found && i == path_len - 2) {
            /* The last "tree branch" is actually the link itself */
            /* Handle edge case */
        }
    }
    /* Add the link branch to complete the cycle */
    cycle[(*cycle_len)++] = link_id;

    return 0;
}

/* ==========================================================================
 * L3: Graph Statistics
 * ========================================================================== */

double ct_graph_density(const ct_circuit_t *circuit)
{
    if (!circuit || circuit->num_nodes < 2) return -1.0;

    int32_t n = circuit->num_nodes;
    int32_t b = circuit->num_branches;

    /* Maximum possible branches in a simple graph: n*(n-1)/2 */
    /* For multigraph (parallel branches allowed): use this formula conservatively */
    double max_branches = (double)n * (double)(n - 1) / 2.0;
    if (max_branches < 1.0) return -1.0;

    return (double)b / max_branches;
}

int ct_cyclomatic_number(const ct_circuit_t *circuit)
{
    if (!circuit) return -1;

    int c = ct_count_components(circuit);
    if (c < 0) return -1;

    /* mu = b - n + c */
    return circuit->num_branches - circuit->num_nodes + c;
}

int ct_graph_rank(const ct_circuit_t *circuit)
{
    if (!circuit) return -1;

    int c = ct_count_components(circuit);
    if (c < 0) return -1;

    /* r = n - c */
    return circuit->num_nodes - c;
}

/* ==========================================================================
 * L3: Series/Parallel Detection
 * ========================================================================== */

int ct_is_series(const ct_circuit_t *circuit, int32_t br1, int32_t br2,
                 int32_t *common_node)
{
    if (!circuit) return -1;
    if (br1 < 0 || br1 >= circuit->num_branches) return -1;
    if (br2 < 0 || br2 >= circuit->num_branches) return -1;

    const ct_branch_t *b1 = &circuit->branches[br1];
    const ct_branch_t *b2 = &circuit->branches[br2];

    /* Find common node */
    int32_t common = -1;
    if (b1->node_from == b2->node_from || b1->node_from == b2->node_to)
        common = b1->node_from;
    else if (b1->node_to == b2->node_from || b1->node_to == b2->node_to)
        common = b1->node_to;

    if (common < 0) return 0;  /* No common node */

    /* Check common node has degree exactly 2 */
    int deg = ct_node_degree(circuit, common);
    if (deg != 2) return 0;

    if (common_node) *common_node = common;
    return 1;
}

int ct_is_parallel(const ct_circuit_t *circuit, int32_t br1, int32_t br2)
{
    if (!circuit) return -1;
    if (br1 < 0 || br1 >= circuit->num_branches) return -1;
    if (br2 < 0 || br2 >= circuit->num_branches) return -1;

    const ct_branch_t *b1 = &circuit->branches[br1];
    const ct_branch_t *b2 = &circuit->branches[br2];

    /* Same node pair (ignoring orientation) */
    if ((b1->node_from == b2->node_from && b1->node_to == b2->node_to) ||
        (b1->node_from == b2->node_to && b1->node_to == b2->node_from)) {
        return 1;
    }
    return 0;
}

/* ==========================================================================
 * L2: Incidence Matrix Construction
 * ========================================================================== */

int ct_build_incidence_matrix(const ct_circuit_t *circuit, ct_incidence_t *A)
{
    if (!circuit || !A) return -1;

    int32_t n = circuit->num_nodes;
    int32_t b = circuit->num_branches;

    memset(A, 0, sizeof(ct_incidence_t));
    A->rows = n - 1;  /* Exclude ground (node 0) */
    A->cols = b;

    /* Fill incidence matrix for non-ground nodes (rows 0..n-2 map to nodes 1..n-1) */
    for (int32_t j = 0; j < b; j++) {
        const ct_branch_t *br = &circuit->branches[j];
        int32_t from = br->node_from;
        int32_t to = br->node_to;

        if (from > 0) {
            A->data[from - 1][j] = 1;   /* Branch leaves node */
        }
        if (to > 0) {
            A->data[to - 1][j] = -1;    /* Branch enters node */
        }
    }

    return 0;
}

int ct_incidence_to_adjacency(const ct_incidence_t *A, int32_t n_nodes,
                              int32_t adj_out[CT_MAX_NODES][CT_MAX_NODES])
{
    if (!A || !adj_out) return -1;

    /* Clear output */
    for (int32_t i = 0; i < n_nodes; i++) {
        for (int32_t j = 0; j < n_nodes; j++) {
            adj_out[i][j] = 0;
        }
    }

    /* For each branch (column), find the two endpoints */
    for (int32_t j = 0; j < A->cols; j++) {
        int32_t from = -1, to = -1;
        for (int32_t i = 0; i < A->rows; i++) {
            if (A->data[i][j] == 1) from = i + 1;      /* node i+1 */
            if (A->data[i][j] == -1) to = i + 1;
        }
        /* Also check ground node (implicit) */
        int32_t sum = 0;
        for (int32_t i = 0; i < A->rows; i++) {
            sum += A->data[i][j];
        }
        if (sum == 1) to = 0;       /* Branch goes to ground */
        if (sum == -1) from = 0;    /* Branch comes from ground */

        if (from >= 0 && to >= 0) {
            adj_out[from][to]++;
            if (from != to) adj_out[to][from]++;
        }
    }

    return 0;
}

/* ==========================================================================
 * L2: Build Matrices from Tree
 * ========================================================================== */

int ct_build_cutset_matrix(const ct_circuit_t *circuit,
                           const ct_tree_t *tree, ct_cutset_matrix_t *Q)
{
    if (!circuit || !tree || !Q) return -1;

    int32_t n = circuit->num_nodes;
    int32_t b = circuit->num_branches;

    memset(Q, 0, sizeof(ct_cutset_matrix_t));
    Q->rows = n - 1;
    Q->cols = b;

    /* For each tree branch, build its fundamental cut-set */
    for (int32_t t = 0; t < tree->num_tree_branches; t++) {
        int32_t tree_br = tree->tree_branches[t];

        /* Mark the tree branch itself (orientation +1) */
        Q->data[t][tree_br] = 1;

        /* Find links in this cut-set: links that connect the two components
         * formed by removing the tree branch. */
        /* Simplified: mark links that form fundamental cycles with this tree branch */
        int32_t cycle[CT_MAX_BRANCHES];
        int32_t cycle_len;

        /* For each link, check if this tree branch is in its fundamental cycle */
        for (int32_t l = 0; l < tree->num_links; l++) {
            int32_t link_br = tree->link_branches[l];
            ct_fundamental_cycle(circuit, tree, link_br, cycle, &cycle_len);

            for (int32_t c = 0; c < cycle_len; c++) {
                if (cycle[c] == tree_br) {
                    /* This link is in the cut-set of this tree branch */
                    /* Orientation: same direction as tree branch */
                    (void)circuit; /* suppress unused warning in minimal builds */
                    Q->data[t][link_br] = 1;
                    break;
                }
            }
        }
    }

    return 0;
}

int ct_build_loop_matrix(const ct_circuit_t *circuit,
                         const ct_tree_t *tree, ct_loop_matrix_t *B)
{
    if (!circuit || !tree || !B) return -1;

    int32_t b = circuit->num_branches;
    int32_t num_links = tree->num_links;

    memset(B, 0, sizeof(ct_loop_matrix_t));
    B->rows = num_links;
    B->cols = b;

    /* For each link, its fundamental loop forms one row of B */
    for (int32_t l = 0; l < num_links; l++) {
        int32_t link_br = tree->link_branches[l];
        int32_t cycle[CT_MAX_BRANCHES];
        int32_t cycle_len;

        ct_fundamental_cycle(circuit, tree, link_br, cycle, &cycle_len);

        for (int32_t c = 0; c < cycle_len; c++) {
            int32_t br = cycle[c];
            if (br == link_br) {
                B->data[l][br] = 1;  /* Link defines positive orientation */
            } else {
                /* Tree branch orientation in this loop */
                B->data[l][br] = 1;
            }
        }
    }

    return 0;
}
