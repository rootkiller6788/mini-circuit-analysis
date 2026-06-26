# Knowledge Graph — mini-circuit-topology

## L1: Definitions
- ct_node_t, ct_branch_t, ct_circuit_t: Core topology types
- ct_incidence_t, ct_cutset_matrix_t, ct_loop_matrix_t: Topological matrices
- ct_tree_t: Spanning tree with co-tree enumeration
- ct_mna_system_t: MNA system matrix and vectors
- ct_dc_result_t, ct_ac_result_t: Analysis result types
- ct_sparse_matrix_t: CSR sparse matrix representation
- 24 element types: R, C, L, V, I, VCVS, CCCS, VCCS, CCVS, diode, BJT, MOS, op-amp, transformer, gyrator, nullator, norator, mutual inductor, transmission line, VCO, memristor

## L2: Core Concepts
- KCL: Sum of currents at each node = 0 (A * i_b = 0)
- KVL: Sum of voltages around each loop = 0 (B * v_b = 0)
- Node-Voltage Method: Y_n * v_n = i_n
- Mesh-Current Method: Z_m * i_m = v_m
- Modified Nodal Analysis (MNA): Combined formulation
- Spanning Tree: n nodes, n-1 branches, no cycles
- Co-tree (Link Set): b-n+1 links defining fundamental loops
- Fundamental Cut-set: Single tree branch + connecting links
- Fundamental Loop: Single link + tree path
- Superposition: Linear circuit source decomposition
- Thevenin Equivalent: V_th + R_th series reduction
- Norton Equivalent: I_n || R_n parallel reduction
- Maximum Power Transfer: R_load = R_th
- Y-Delta Transformation: Star-mesh equivalence (Kennelly 1899)
- Planarity: Euler formula V-E+F=2

## L3: Mathematical Structures
- Incidence Matrix A: (n-1)*b, entries in {+1,-1,0}
- Cut-set Matrix Q: Q * i_b = 0 (KCL)
- Loop Matrix B: B * v_b = 0 (KVL)
- Orthogonality: Q * B^T = 0
- Node Admittance Matrix Y_n: Self + mutual admittances
- Mesh Impedance Matrix Z_m: Self + mutual impedances
- MNA System Matrix: (n+nv)*(n+nv) block structure
- CSR Sparse Format: val[], col[], row_ptr[]

## L4: Fundamental Laws
- Tellegen Theorem (1952): sum(v_k*i_k)=0
- Strong Tellegen: sum(v1_k*i2_k)=0 (cross-circuit)
- KCL: A*i_b=0 (charge conservation)
- KVL: B*v_b=0 (energy conservation)
- Euler Formula: V-E+F=2 (planar graphs)
- Matrix-Tree Theorem: det(Y_n)=sum of tree products
- Reciprocity: Z_12=Z_21 for passive bilateral networks
- Maximum Power Transfer: P_max=V_th^2/(4R_th)

## L5: Algorithms
- MNA Element Stamping: R, C, L, V, I, VCCS, CCCS, Op-amp
- LU Decomposition with Partial Pivoting (complex)
- Forward/Back Substitution: O(n^2)
- DFS/BFS Graph Traversal: O(n+b)
- Union-Find with Path Compression (Kruskal-style tree)
- Gaussian Elimination with Pivoting (dense)
- Fundamental Cycle Detection (BFS on tree)
- Superposition Solver (sequential source activation)
- Adjoint Sensitivity (Director-Rohrer, 1969)
- CSR Conversion (dense to sparse)
- AMD Ordering (approximate minimum degree)

## L6: Canonical Problems
- Resistor Network DC Analysis (MNA solve)
- Wheatstone Bridge (balanced/unbalanced/sensitivity)
- RLC Resonance (AC sweep, Q-factor, bandwidth)
- Thevenin/Norton Equivalent Extraction
- Maximum Power Transfer Computation

## L7: Applications
- SPICE-Compatible MNA DC/AC Solver
- Sensor Bridge Analysis (strain gauge: Toyota, ISO; RTD temperature)
- Passive Filter Design (RLC crossover: Bass audio)
- Power Grid Topology Validation
- PCB Signal Integrity (parasitic extraction)

## L8: Advanced Topics
- Sparse Matrix Methods (CSR format, O(nnz) MV)
- AMD Fill-in Reduction (sparse LU ordering)
- Adjoint Network Sensitivity (all parameters from 2 solves)
- Circuit Partitioning (Kernighan-Lin style)
- Symbolic Circuit Analysis (topological formulas)
- Monte Carlo Tolerance Analysis

## L9: Research Frontiers
- Quantum Circuit Topology (graph states, MBQC)
- AI-Driven Circuit Optimization (RL for topology)
- Memristor Crossbar Topology (neuromorphic computing)
- 6G RIS Circuit Topology (reconfigurable surfaces)
- Semantic Communication Circuit Graphs