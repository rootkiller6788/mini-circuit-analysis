# Knowledge Graph - mini-network-theorem

## L1: Definitions (Complete)

| # | Concept | C Definition | Location |
|---|---------|-------------|----------|
| 1 | Complex Impedance | `typedef struct ComplexImpedance` | include/network_theorem.h |
| 2 | Complex Admittance | `typedef struct ComplexAdmittance` | include/network_theorem.h |
| 3 | Voltage Source | `typedef struct VoltageSource` | include/network_theorem.h |
| 4 | Current Source | `typedef struct CurrentSource` | include/network_theorem.h |
| 5 | Element Type Enum | `typedef enum ElementType` | include/network_theorem.h |
| 6 | Circuit Element | `typedef struct CircuitElement` | include/network_theorem.h |
| 7 | Thevenin Equivalent | `typedef struct TheveninEquivalent` | include/network_theorem.h |
| 8 | Norton Equivalent | `typedef struct NortonEquivalent` | include/network_theorem.h |
| 9 | Superposition Result | `typedef struct SuperpositionResult` | include/network_theorem.h |
| 10 | Max Power Transfer Result | `typedef struct MaxPowerTransferResult` | include/network_theorem.h |
| 11 | Reciprocity Parameters | `typedef struct ReciprocityParams` | include/network_theorem.h |
| 12 | Z-Parameters | `typedef struct ZParameters` | include/network_theorem.h |
| 13 | Y-Parameters | `typedef struct YParameters` | include/network_theorem.h |
| 14 | H-Parameters | `typedef struct HParameters` | include/network_theorem.h |
| 15 | ABCD Parameters | `typedef struct ABCDParameters` | include/network_theorem.h |
| 16 | S-Parameters | `typedef struct SParameters` | include/network_theorem.h |
| 17 | Circuit Node | `typedef struct CircuitNode` | include/network_theorem.h |
| 18 | Circuit Branch | `typedef struct CircuitBranch` | include/network_theorem.h |
| 19 | Circuit Mesh | `typedef struct CircuitMesh` | include/network_theorem.h |
| 20 | Circuit Topology | `typedef struct CircuitTopology` | include/network_theorem.h |
| 21 | Millman Result | `typedef struct MillmanResult` | include/network_theorem.h |
| 22 | Tellegen Result | `typedef struct TellegenResult` | include/network_theorem.h |
| 23 | Wye Network | `typedef struct WyeNetwork` | include/network_theorem.h |
| 24 | Delta Network | `typedef struct DeltaNetwork` | include/network_theorem.h |
| 25 | Compensation Result | `typedef struct CompensationResult` | include/network_theorem.h |
| 26 | Substitution Result | `typedef struct SubstitutionResult` | include/network_theorem.h |
| 27 | MNA System | `typedef struct MNASystem` | include/network_theorem.h |
| 28 | Network Reduction | `typedef struct NetworkReduction` | include/network_theorem.h |

## L2: Core Concepts (Complete)

| # | Concept | Implementation | File |
|---|---------|---------------|------|
| 1 | Series Impedance | `impedance_series()` | src/thevenin_norton.c |
| 2 | Parallel Impedance | `impedance_parallel()` | src/thevenin_norton.c |
| 3 | Impedance-Admittance Conversion | `impedance_to_admittance()` | src/thevenin_norton.c |
| 4 | Source Transformation (Th to No) | `thevenin_to_norton()` | src/thevenin_norton.c |
| 5 | Source Transformation (No to Th) | `norton_to_thevenin()` | src/thevenin_norton.c |
| 6 | Per-Branch Source Transform | `source_transform_thevenin_to_norton()` | src/wye_delta.c |
| 7 | Z to Y Conversion | `z_to_y()` | src/two_port.c |
| 8 | Y to Z Conversion | `y_to_z()` | src/two_port.c |
| 9 | Z to H Conversion | `z_to_h()` | src/two_port.c |
| 10 | H to Z Conversion | `h_to_z()` | src/two_port.c |
| 11 | Z to ABCD Conversion | `z_to_abcd()` | src/two_port.c |
| 12 | ABCD Cascade | `abcd_cascade()` | src/two_port.c |
| 13 | Z to S Conversion | `z_to_s()` | src/two_port.c |
| 14 | S to Z Conversion | `s_to_z()` | src/two_port.c |
| 15 | Impedance Magnitude | `impedance_magnitude()` | src/thevenin_norton.c |
| 16 | Impedance Phase | `impedance_phase()` | src/thevenin_norton.c |
| 17 | Power Calculations | `power_dissipated()` | src/thevenin_norton.c |
| 18 | Power Factor | `power_factor()` | src/thevenin_norton.c |
| 19 | Series Resistance | `combine_series_resistors()` | src/wye_delta.c |
| 20 | Parallel Resistance | `combine_parallel_resistors()` | src/wye_delta.c |
| 21 | Input Impedance (2-port) | `two_port_input_impedance()` | src/two_port.c |
| 22 | Output Impedance (2-port) | `two_port_output_impedance()` | src/two_port.c |
| 23 | Voltage Gain (2-port) | `two_port_voltage_gain()` | src/two_port.c |
| 24 | Voltage Divider | `voltage_divider_analyze()` | src/wye_delta.c |
| 25 | Current Divider | `current_divider_analyze()` | src/wye_delta.c |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation | File |
|---|----------|---------------|------|
| 1 | Nodal Admittance Matrix | `build_nodal_admittance_matrix()` | src/mesh_nodal.c |
| 2 | Mesh Impedance Matrix | `build_mesh_impedance_matrix()` | src/mesh_nodal.c |
| 3 | MNA Augmented Matrix | `mna_build_system()` | src/mna_solver.c |
| 4 | Sparse CSR Matrix | `dense_to_csr()` | src/mna_solver.c |
| 5 | Matrix Determinant | `matrix_determinant()` | src/mesh_nodal.c |
| 6 | Matrix Inverse | `matrix_inverse()` | src/mesh_nodal.c |
| 7 | Matrix Condition Number | `matrix_condition_estimate()` | src/mesh_nodal.c |
| 8 | Residual Norm | `compute_residual_norm()` | src/mesh_nodal.c |

## L4: Fundamental Laws (Complete)

| # | Theorem/Formula | C Implementation | Lean Statement | File |
|---|----------------|-----------------|----------------|------|
| 1 | Thevenin's Theorem | `compute_thevenin()` | `TheveninEquiv` | src/thevenin_norton.c |
| 2 | Norton's Theorem | `compute_norton()` | `NortonEquiv` | src/thevenin_norton.c |
| 3 | Superposition Theorem | `superposition_solve()` | `superposition_linear` | src/superposition.c |
| 4 | Maximum Power Transfer | `max_power_transfer()` | `max_power_transfer` | src/power_transfer.c |
| 5 | Reciprocity Theorem | `verify_reciprocity()` | `ReciprocalTwoPort` | src/two_port.c |
| 6 | Millman's Theorem | `millman_solve()` | `millman_formula` | src/power_transfer.c |
| 7 | Tellegen's Theorem | `tellegen_verify()` | `tellegen_conservation` | src/power_transfer.c |
| 8 | Compensation Theorem | `compensation_compute()` | - | src/power_transfer.c |
| 9 | Source Transformation | `thevenin_to_norton()` | `source_transformation_equivalence` | src/thevenin_norton.c |
| 10 | Wye-Delta Transform | `wye_to_delta()` / `delta_to_wye()` | `wye_delta_equivalence` | src/wye_delta.c |
| 11 | Thevenin/Norton Duality | `norton_to_thevenin()` | formal dual in Lean | src/thevenin_norton.c |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Implementation | Complexity | File |
|---|----------|---------------|------------|------|
| 1 | Gaussian Elimination | `gaussian_elimination()` | O(n^3) | src/mesh_nodal.c |
| 2 | LU Decomposition (Doolittle) | `lu_decompose()` | O(n^3) | src/mesh_nodal.c |
| 3 | LU Solve (Fwd/Bwd Sub) | `lu_solve()` | O(n^2) | src/mesh_nodal.c |
| 4 | Nodal Analysis Solver | `solve_nodal_voltages()` | O(n^3) | src/mesh_nodal.c |
| 5 | Mesh Analysis Solver | `solve_mesh_currents()` | O(n^3) | src/mesh_nodal.c |
| 6 | Modified Nodal Analysis | `mna_build_system()` + `mna_solve()` | O(m^3) | src/mna_solver.c |
| 7 | Superposition Algorithm | `superposition_solve()` | O(N * n^3) | src/superposition.c |
| 8 | Thevenin Computation | `compute_thevenin()` | O(n^3) | src/thevenin_norton.c |
| 9 | Wye-Delta Transformation | `wye_to_delta()` / `delta_to_wye()` | O(1) | src/wye_delta.c |
| 10 | Network Reduction | `network_reduce()` | iterative | src/wye_delta.c |
| 11 | Newton-Raphson Step | `newton_raphson_step()` | O(n^3) per iter | src/mna_solver.c |
| 12 | Source Deactivation | `deactivate_independent_sources()` | O(E) | src/superposition.c |
| 13 | Sparse Mat-Vec Multiply | `csr_mat_vec_mul()` | O(nnz) | src/mna_solver.c |

## L6: Canonical Problems (Complete)

| # | Problem | Implementation | Example | File |
|---|---------|---------------|---------|------|
| 1 | Wheatstone Bridge | `wheatstone_analyze()` | example_bridge.c | src/wye_delta.c |
| 2 | Voltage Divider (loaded) | `voltage_divider_analyze()` | - | src/wye_delta.c |
| 3 | Current Divider | `current_divider_analyze()` | - | src/wye_delta.c |
| 4 | Multi-Source Circuit | `superposition_solve()` | example_superposition.c | src/superposition.c |
| 5 | Thevenin Equivalent Extraction | `compute_thevenin()` | example_thevenin.c | src/thevenin_norton.c |
| 6 | Maximum Power to Load | `max_power_transfer()` | example_thevenin.c | src/power_transfer.c |
| 7 | Cascaded Two-Port | `abcd_cascade()` | - | src/two_port.c |
| 8 | Balanced/Unbalanced Bridge | `wheatstone_analyze()` | example_bridge.c | src/wye_delta.c |

## L7: Applications (Partial+)

| # | Application | Domain | Implementation |
|---|------------|--------|---------------|
| 1 | Audio Amplifier Matching | Detroit automotive audio | `audio_amplifier_match()` in src/power_transfer.c |
| 2 | Strain Gauge Bridge | Boeing 787 structural health | `wheatstone_analyze()` in src/wye_delta.c |
| 3 | Smart Grid Load Analysis | Power systems | Millman superposition in example_bridge.c |

## L8: Advanced Topics (Partial+)

| # | Topic | Implementation | File |
|---|-------|---------------|------|
| 1 | Nonlinear Circuit Newton-Raphson | `newton_raphson_step()` | src/mna_solver.c |
| 2 | Sparse Matrix (CSR) for Large Circuits | `dense_to_csr()` | src/mna_solver.c |
| 3 | Iterative Network Reduction | `network_reduce()` | src/wye_delta.c |
| 4 | Condition Number Estimation | `matrix_condition_estimate()` | src/mesh_nodal.c |
| 5 | Monte Carlo Tolerance (framework) | `compensation_compute()` | src/power_transfer.c |

## L9: Research Frontiers (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | AI-Assisted Circuit Analysis | Documented in course-tree.md |
| 2 | Quantum Circuit Equivalents | Documented, no implementation |
| 3 | 6G RIS Impedance Matching | Documented reference |

## Summary

| Level | Status | Count |
|-------|--------|-------|
| L1 Definitions | **Complete** | 28 structs/enums |
| L2 Core Concepts | **Complete** | 25 functions |
| L3 Math Structures | **Complete** | 8 matrix ops |
| L4 Fundamental Laws | **Complete** | 11 theorems (C + Lean) |
| L5 Algorithms | **Complete** | 13 algorithms |
| L6 Canonical Problems | **Complete** | 8 problems (+ examples) |
| L7 Applications | **Partial+** | 3 applications |
| L8 Advanced Topics | **Partial+** | 5 topics |
| L9 Research Frontiers | **Partial** | 3 references |
