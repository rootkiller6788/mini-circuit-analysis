# mini-spice-simulation

A compact SPICE-compatible circuit simulator implementing Modified Nodal Analysis (MNA), Newton-Raphson DC solve, AC small-signal frequency analysis, and transient time-domain simulation.

**Reference**: Nagel, L.W. (1975). "SPICE2: A Computer Program to Simulate Semiconductor Circuits." UCB/ERL M520.

---

## Module Status: COMPLETE ✅

- **L1-L6**: Complete
- **L7**: Partial+ (3 applications: IC design simulation, Bode plot CSV export, waveform CSV export)
- **L8**: Partial+ (4 advanced topics: adaptive time-step, Gmin stepping, junction limiting, condition number)
- **L9**: Partial (3 research topics documented: parallel SPICE, ML convergence, GPU acceleration)

### Line Count Verification

| Directory | Files | Lines |
|-----------|-------|-------|
| `include/` | 5 headers | 1,659 |
| `src/` | 5 C sources | 3,465+ |
| **Total** | | **5,124+** ≥ 3,000 ✅ |

---

## Quick Start

```bash
# Build library and run tests
make test

# Build examples
make examples

# Run specific example
./examples/example_divider_dc
./examples/example_rc_transient
./examples/example_rlc_ac
```

---

## Architecture

```
┌─────────────────────────────────────────────┐
│               spice_simulator_t             │  ← Top-level engine
├─────────────────────────────────────────────┤
│  spice_netlist  │  spice_matrix             │
│  (parser)       │  (CSR, dense, LU)         │
├─────────────────────────────────────────────┤
│  spice_models   │  spice_analysis           │
│  (R,C,L,D,M,Q)  │  (DC, AC, TRAN, TF)       │
└─────────────────────────────────────────────┘
```

### Simulation Flow

```
Netlist (.cir) → Parse → MNA Assembly → DC OP (Newton-Raphson)
                                           ↓
                              AC Analysis ←┘ → TRAN Analysis
                              (complex LU)     (companion models)
```

---

## Knowledge Coverage

### L1 — Core Definitions (19/19 Complete)
`spice_component_type_t`, `spice_node_t`, `spice_resistor_t`, `spice_capacitor_t`, `spice_inductor_t`, `spice_vsource_t`, `spice_isource_t`, `spice_diode_model_t`, `spice_mos_model_t`, `spice_bjt_model_t`, `spice_netlist_t`, `spice_csr_matrix_t`, `spice_dense_matrix_t`, `spice_lu_factor_t`, `spice_dc_result_t`, `spice_ac_result_t`, `spice_tran_result_t`, `spice_convergence_params_t`, `spice_simulator_t`

### L2 — Core Concepts (10/10 Complete)
Modified Nodal Analysis, DC Operating Point, AC Small-Signal, Transient Analysis, Companion Models (Trapezoidal), Netlist Parsing, Convergence Criteria, Time-Step Control, Gmin Stepping, Device Model Stamping

### L3 — Mathematical Structures (6/6 Complete)
CSR Sparse Matrix, Dense Column-Major Matrix, C99 Complex Numbers, Vector BLAS-1 Operations, Complex Matrix/Vector, LU Decomposition

### L4 — Fundamental Laws (8/8 Complete)
| Law | Formula | Verification |
|-----|---------|-------------|
| Ohm's Law | V = I·R | `test_kcl_verification()` |
| KCL | ΣI = 0 at node | `test_kcl_verification()` |
| KVL | ΣV = 0 in loop | `test_kcl_verification()` |
| Capacitor I-V | I = C·dV/dt | Trapezoidal companion |
| Inductor I-V | V = L·dI/dt | Trapezoidal companion |
| Shockley Equation | I = I_S·(e^{V/(n·V_T)}-1) | `test_shockley_equation_assertion()` |
| Power Conservation | P_src = ΣP_diss | `test_kcl_verification()` |
| Voltage Divider | V_out = V_in·R2/(R1+R2) | `test_kcl_verification()` |

### L5 — Algorithms (10/10 Complete)
Newton-Raphson, LU Decomposition (Partial Pivoting), Forward/Back Substitution, Complex LU, Trapezoidal Integration, Sparse MatVec, Condition Number Estimation, MOSFET Shichman-Hodges, BJT Ebers-Moll, Safe Exponential (Junction Limiting)

### L6 — Canonical Problems (8/8 Complete)
Resistive Divider (DC), Series Resistor Network (DC), RC Transient (τ = RC), RLC Bandpass (AC), RC Lowpass (AC), Diode Biasing, MOSFET OP, BJT Forward-Active

### L7 — Applications (3/4, Partial+)
IC Design Simulation, Bode Plot CSV Export, Transient Waveform CSV Export

### L8 — Advanced Topics (4/5, Partial+)
Adaptive Time-Step, Gmin Stepping, PN Junction Limiting, Condition Number Estimation

### L9 — Research Frontiers (3 documented, Partial)
Parallel SPICE, ML Convergence Prediction, GPU Acceleration

---

## Supported SPICE Syntax

```
* Comment lines
R<name> <n+> <n-> <value>
C<name> <n+> <n-> <value> [IC=<v0>]
L<name> <n+> <n-> <value> [IC=<i0>]
V<name> <n+> <n-> [DC] <value> [AC <mag> <phase>]
I<name> <n+> <n-> [DC] <value>
D<name> <n+> <n-> <model>
Q<name> <nc> <nb> <ne> <model>
M<name> <nd> <ng> <ns> <nb> <model> [W=<w> L=<l>]

.MODEL <name> <type> (param1=val1 ...)
.DC <src> <start> <stop> <step>
.AC DEC <n> <fstart> <fstop>
.TRAN <tstep> <tstop> [tstart] [tmax]
.OPTIONS [param=value ...]
.END
```

---

## File Listing

```
mini-spice-simulation/
├── Makefile
├── README.md                     ← This file
├── include/
│   ├── spice_netlist.h           ← Netlist parser types & API
│   ├── spice_matrix.h            ← CSR sparse, dense, LU, BLAS-1
│   ├── spice_analysis.h          ← DC/AC/TRAN/TF analysis & results
│   ├── spice_models.h            ← R/L/C/V/I/D/M/Q model stamps
│   └── spice_core.h              ← Top-level simulator engine
├── src/
│   ├── spice_netlist.c           ← Netlist parser implementation
│   ├── spice_matrix.c            ← Matrix ops, LU, BLAS
│   ├── spice_models.c            ← Device evaluation & MNA stamps
│   ├── spice_analysis.c          ← DC/AC/TRAN analysis engines
│   └── spice_core.c              ← Simulator lifecycle & orchestration
├── tests/
│   └── test_spice.c              ← 24 tests (all passing)
├── examples/
│   ├── example_divider_dc.c      ← DC voltage divider
│   ├── example_rc_transient.c    ← RC charging transient
│   └── example_rlc_ac.c          ← RLC bandpass AC sweep
└── docs/
    ├── knowledge-graph.md        ← L1-L9 itemized knowledge map
    ├── coverage-report.md        ← Per-level coverage assessment
    ├── gap-report.md             ← Identified gaps & priorities
    ├── course-alignment.md       ← Nine-school curriculum mapping
    └── course-tree.md            ← Prerequisites & dependents
```

---

## Test Results

```
╔══════════════════════════════════════════════╗
║  Results:  24 passed,   0 failed            ║
╚══════════════════════════════════════════════╝

── Netlist Parser ───── 7/7 PASS
── Matrix Operations ─── 4/4 PASS
── Device Models ─────── 4/4 PASS
── Circuit Analysis ──── 4/4 PASS
── Simulator Core ────── 3/3 PASS
── Math Assertions ───── 2/2 PASS
```

---

## Nine-School Curriculum Alignment

This module covers core intersections of SPICE/circuit simulation across MIT (6.003/6.630), Stanford (EE102A/EE359), Berkeley (EE16A/B, EE105, EE117), Illinois (ECE 310/459), Michigan (EECS 351/411), Georgia Tech (ECE 4270/6601), TU Munich (HF Engineering), ETH (227-0427/0455), and Tsinghua (信号与系统/电磁场).

See [docs/course-alignment.md](docs/course-alignment.md) for detailed mapping.
