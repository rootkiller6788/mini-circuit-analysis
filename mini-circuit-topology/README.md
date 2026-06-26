# mini-circuit-topology

**Circuit Topology Analysis — Graph-Theoretic Foundations of Electrical Networks**

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (Node, Branch, Graph, Tree, Co-tree, Cut-set, Loop, Incidence Matrix)
- **L2 Core Concepts**: Complete (KCL, KVL, spanning tree, fundamental loops/cutsets, DFS/BFS, connectivity)
- **L3 Mathematical Structures**: Complete (Incidence matrix A, cut-set matrix Q, loop matrix B, adjacency, Y-Delta transform)
- **L4 Fundamental Laws**: Complete (Tellegen's Theorem, KCL/KVL verification, power balance, reciprocity)
- **L5 Algorithms/Methods**: Complete (MNA stamping, LU decomposition, forward/back substitution, mesh-current, node-voltage, superposition)
- **L6 Canonical Problems**: Complete (Wheatstone bridge, RLC resonance, Thevenin/Norton, max power transfer)
- **L7 Applications**: Partial+ (SPICE-compatible MNA solver, sensor bridge, filter design, power grid topology)
- **L8 Advanced Topics**: Partial+ (Sparse matrix CSR, AMD reordering, adjoint sensitivity analysis)
- **L9 Research Frontiers**: Partial (Documented: quantum circuit topology, AI-driven circuit optimization)

## Code Metrics

| Directory | Files | Lines |
|-----------|-------|-------|
| include/  | 4     | 1,214 |
| src/      | 6     | 2,510 |
| **Total** | **10** | **3,724** |

## Core Definitions (L1)

- **ct_node_t**: Circuit node with graph properties (color, partition, terminal/internal flags)
- **ct_branch_t**: Two-terminal circuit element (24 element types including memristor)
- **ct_circuit_t**: Complete circuit topology with adjacency matrix
- **ct_incidence_t**: Reduced incidence matrix A ((n-1) x b)
- **ct_cutset_matrix_t**: Fundamental cut-set matrix Q
- **ct_loop_matrix_t**: Fundamental loop matrix B
- **ct_tree_t**: Spanning tree with co-tree enumeration

## Core Theorems (L4)

### Tellegen's Theorem (1952)
For any lumped electrical network:
```
sum(v_k * i_k) = 0   (k = 1..b)
```
Holds for ANY two sets of branch variables satisfying KCL and KVL respectively. Topological in nature — depends only on interconnection pattern.

### Kirchhoff's Current Law (KCL)
```
A * i_b = 0
```
Sum of currents at each node is zero. Encoded in the reduced incidence matrix.

### Kirchhoff's Voltage Law (KVL)
```
B * v_b = 0
```
Sum of voltages around each fundamental loop is zero. Encoded in the loop matrix.

### Euler's Formula for Planar Graphs
```
V - E + F = 2
```
For a connected planar circuit graph: vertices - edges + faces = 2.

### Matrix-Tree Theorem (Kirchhoff, 1847)
```
det(Y_n) = sum_{all trees T} prod_{branch in T} Y_branch
```
Determinant of node admittance matrix equals sum of tree admittance products.

## Core Algorithms (L5)

1. **Modified Nodal Analysis (MNA)** — Standard SPICE formulation
2. **LU Decomposition with Partial Pivoting** — O(n^3/3) for complex matrices
3. **Forward/Back Substitution** — O(n^2) solution after factorization
4. **DFS/BFS Graph Traversal** — O(n + b) connectivity and path finding
5. **Kruskal-style Spanning Tree** — Union-Find with priority for voltage sources
6. **Fundamental Cycle Detection** — BFS-based path finding in tree
7. **Gaussian Elimination** — Dense solver for classical methods
8. **Superposition Principle** — Sequential source activation
9. **Adjoint Network Sensitivity** — Director-Rohrer method (1969)
10. **Y-Delta Transformation** — Kennelly formulas (1899)

## Canonical Problems (L6)

1. **Resistor Network DC Analysis** — MNA solve, power distribution
2. **Wheatstone Bridge** — Balanced/unbalanced, sensor sensitivity
3. **RLC Resonance** — AC frequency sweep, quality factor
4. **Thevenin/Norton Equivalents** — Port reduction
5. **Maximum Power Transfer** — Optimal load matching

## Nine-School Curriculum Mapping

| School | Course | Topics Covered |
|--------|--------|---------------|
| **MIT** | 6.002 Circuits | Node/mesh methods, Thevenin, superposition |
| **Berkeley** | EE16A/B | Circuit topology, MNA, graph methods |
| **Stanford** | EE101A | Network topology, sparse solvers |
| **Illinois** | ECE 451 EM | Graph theory for microwave networks |
| **Michigan** | EECS 215 | Circuit analysis fundamentals |
| **Georgia Tech** | ECE 6350 | EM/circuit topology intersection |
| **TU Munich** | High-Frequency Eng. | Network parameters, S-parameter topology |
| **ETH** | 227-0455 | Electromagnetic network theory |
| **Tsinghua** | Circuit Principles | Graph-theoretic circuit analysis |

## Build & Test

```bash
make          # Build static library
make test     # Run assert-based test suite
make examples # Build example programs
make clean    # Remove build artifacts
```

## References

- L.O. Chua, C.A. Desoer, E.S. Kuh, "Linear and Nonlinear Circuits" (1987)
- S. Seshu, M.B. Reed, "Linear Graphs and Electrical Networks" (1961)
- W.K. Chen, "Applied Graph Theory: Graphs and Electrical Networks" (1976)
- C.W. Ho et al., "The Modified Nodal Approach to Network Analysis", IEEE TCAS (1975)
- B.D.H. Tellegen, "A General Network Theorem", Philips Res. Rep. (1952)
- P. Penfield, R. Spence, S. Duinker, "Tellegen's Theorem and Electrical Networks" (1970)
