# Course Tree - mini-network-theorem

## Prerequisites

```
Basic Circuit Theory
├── Ohm's Law (V = I * R)
├── Kirchhoff's Current Law (KCL)
├── Kirchhoff's Voltage Law (KVL)
└── Complex Numbers (j = sqrt(-1))
```

## Dependency Tree

```
Level 1: Definitions (L1)
├── ComplexImpedance, ComplexAdmittance
├── CircuitElement, CircuitTopology
└── TheveninEquivalent, NortonEquivalent, TwoPort params

Level 2: Core Concepts (L2) [depends on L1]
├── Impedance computation (series, parallel, R, L, C)
├── Source Transformation (Thevenin <-> Norton)
├── Two-Port Parameter Conversions (Z <-> Y <-> H <-> ABCD <-> S)
└── Power calculations

Level 3: Mathematical Structures (L3) [depends on L1, L2]
├── Nodal Admittance Matrix
├── Mesh Impedance Matrix
├── MNA Augmented System
└── Sparse CSR Representation

Level 4: Fundamental Laws (L4) [depends on L1-L3]
├── Thevenin & Norton Theorems
├── Superposition Theorem
├── Maximum Power Transfer Theorem
├── Reciprocity Theorem
├── Millman's Theorem
├── Tellegen's Theorem
├── Compensation Theorem
└── Wye-Delta Transformation

Level 5: Algorithms (L5) [depends on L3, L4]
├── Gaussian Elimination (O(n^3))
├── LU Decomposition (Doolittle)
├── Nodal/Mesh Solvers
├── MNA Solver
├── Superposition Iteration
├── Network Reduction (iterative)
└── Newton-Raphson (nonlinear)

Level 6: Canonical Problems (L6) [depends on L4, L5]
├── Wheatstone Bridge Analysis
├── Multi-Source Superposition
├── Loaded Voltage Divider
├── Current Divider
├── Thevenin Extraction
└── Cascaded Two-Port

Level 7: Applications (L7) [depends on L6]
├── Audio Amplifier Output Matching (Detroit)
├── Strain Gauge Bridge (Boeing 787)
└── Smart Grid Load Analysis

Level 8: Advanced Topics (L8) [depends on L5, L7]
├── Nonlinear Analysis (Newton-Raphson)
├── Sparse Matrix Methods
├── Iterative Network Reduction
├── Condition Number Estimation
└── Monte Carlo Tolerance (Compensation)

Level 9: Research Frontiers (L9) [depends on L8]
├── AI-Assisted Circuit Analysis
├── Quantum Circuit Equivalents
└── 6G RIS Impedance Matching
```

## Cross-Module Dependencies

```
This module (mini-network-theorem) is a dependency for:
├── 2. mini-analog-electronics (two-port models, biasing networks)
├── 3. mini-digital-electronics (Thevenin for logic gate models)
├── 5. mini-communication-principle (matching networks, S-parameters)
├── 7. mini-electromagnetic-wave (transmission line equivalent circuits)
├── 10. mini-power-electronics (efficiency, max power transfer)
├── 11. mini-wireless-mobile-comm (impedance matching, antenna two-port)
├── 18. mini-emc-signal-integrity (Thevenin for noise sources)
└── 19. mini-electronic-mfg-test (bridge measurements)
```
