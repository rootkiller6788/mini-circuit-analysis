# Knowledge Graph — mini-spice-simulation

## L1: Definitions

| # | Concept | C Struct / typedef | Lean Definition | Status |
|---|---------|-------------------|-----------------|--------|
| 1 | Component Types | `spice_component_type_t` (15 types) | — | ✅ |
| 2 | Circuit Node | `spice_node_t`, `spice_node_id` | — | ✅ |
| 3 | Resistor (Ohm's Law) | `spice_resistor_t` | — | ✅ |
| 4 | Capacitor (I=C·dV/dt) | `spice_capacitor_t` | — | ✅ |
| 5 | Inductor (V=L·dI/dt) | `spice_inductor_t` | — | ✅ |
| 6 | Voltage Source | `spice_vsource_t` | — | ✅ |
| 7 | Current Source | `spice_isource_t` | — | ✅ |
| 8 | Diode Model Parameters | `spice_diode_model_t` | — | ✅ |
| 9 | MOSFET Level 1 Parameters | `spice_mos_model_t` | — | ✅ |
| 10 | BJT Model Parameters | `spice_bjt_model_t` | — | ✅ |
| 11 | Parsed Netlist | `spice_netlist_t` | — | ✅ |
| 12 | CSR Sparse Matrix | `spice_csr_matrix_t` | — | ✅ |
| 13 | Dense Matrix | `spice_dense_matrix_t` | — | ✅ |
| 14 | LU Decomposition | `spice_lu_factor_t` | — | ✅ |
| 15 | DC Operating Point Result | `spice_dc_result_t` | — | ✅ |
| 16 | AC Analysis Result | `spice_ac_result_t` | — | ✅ |
| 17 | Transient Analysis Result | `spice_tran_result_t` | — | ✅ |
| 18 | Convergence Parameters | `spice_convergence_params_t` | — | ✅ |
| 19 | Simulator Object | `spice_simulator_t` | — | ✅ |

## L2: Core Concepts

| # | Concept | Implementation | Status |
|---|---------|---------------|--------|
| 1 | Modified Nodal Analysis (MNA) | `build_mna_dc()`, `build_mna_ac()` | ✅ |
| 2 | DC Operating Point | `spice_dc_analysis()` | ✅ |
| 3 | AC Small-Signal Analysis | `spice_ac_analysis()` | ✅ |
| 4 | Transient Analysis | `spice_transient_analysis()` | ✅ |
| 5 | Companion Models (Trapezoidal) | `spice_stamp_capacitor_tran()`, `spice_stamp_inductor_tran()` | ✅ |
| 6 | Netlist Parsing | `spice_netlist_parse_file()` | ✅ |
| 7 | Convergence Criteria | `spice_check_dc_convergence()` | ✅ |
| 8 | Time-Step Control | Adaptive dt in `spice_transient_analysis()` | ✅ |
| 9 | Gmin Stepping | Added in `build_mna_dc()` | ✅ |
| 10 | Device Model Stamping | `spice_stamp_resistor()`, `spice_stamp_diode_dc()`, etc. | ✅ |

## L3: Mathematical Structures

| # | Concept | Data Type / Operation | Status |
|---|---------|----------------------|--------|
| 1 | CSR Sparse Matrix | `spice_csr_matrix_t` + `spice_csr_matvec()` | ✅ |
| 2 | Dense Column-Major Matrix | `spice_dense_matrix_t` | ✅ |
| 3 | Complex Numbers (C99) | `spice_complex_t` = `double complex` | ✅ |
| 4 | Real Vector Operations | `spice_vector_dot()`, `spice_vector_norm2()`, `spice_vector_norm_inf()` | ✅ |
| 5 | BLAS-1: scale, axpy, copy, zero | `spice_vector_scale()`, `spice_vector_axpy()`, etc. | ✅ |
| 6 | Complex Matrix/Vector | `spice_complex_vector_t`, `spice_complex_lu_solve()` | ✅ |

## L4: Fundamental Laws

| # | Law / Theorem | Code Verification | Formal Statement (Lean) | Status |
|---|--------------|-------------------|------------------------|--------|
| 1 | Ohm's Law (V=I·R) | `spice_stamp_resistor()` + DC test | — | ✅ |
| 2 | KCL (Kirchhoff's Current Law) | `test_kcl_verification()` — ΣI=0 at node | — | ✅ |
| 3 | KVL (Kirchhoff's Voltage Law) | `test_kcl_verification()` — ΣV=0 around loop | — | ✅ |
| 4 | Capacitor I-V (I=C·dV/dt) | `spice_stamp_capacitor_tran()` companion model | — | ✅ |
| 5 | Inductor I-V (V=L·dI/dt) | `spice_stamp_inductor_tran()` companion model | — | ✅ |
| 6 | Shockley Diode Equation | `spice_diode_evaluate()` + `test_shockley_equation_assertion()` | — | ✅ |
| 7 | Power Conservation | `test_kcl_verification()` — P_src = ΣP_dissipated | — | ✅ |
| 8 | Voltage Divider Formula | `test_kcl_verification()` — Vout = Vin·R2/(R1+R2) | — | ✅ |

## L5: Algorithms / Methods

| # | Algorithm | Implementation | Complexity | Status |
|---|-----------|---------------|------------|--------|
| 1 | Newton-Raphson Iteration | `spice_dc_analysis()` inner loop | O(N³)/iter | ✅ |
| 2 | LU Decomposition (Dense) | `spice_dense_lu_factor()` | O(N³) | ✅ |
| 3 | Forward/Back Substitution | `spice_dense_lu_solve()` | O(N²) | ✅ |
| 4 | Complex LU Solve | `spice_complex_lu_solve()` | O(N³) | ✅ |
| 5 | Trapezoidal Integration | `spice_stamp_capacitor_tran()` | O(1) per device | ✅ |
| 6 | Sparse Matrix-Vector Multiply | `spice_csr_matvec()` | O(NNZ) | ✅ |
| 7 | Condition Number Estimation | `spice_condition_number()` | O(N³) | ✅ |
| 8 | MOSFET Level 1 Evaluation | `spice_mosfet_evaluate()` | O(1) per device | ✅ |
| 9 | BJT Ebers-Moll Evaluation | `spice_bjt_evaluate()` | O(1) per device | ✅ |
| 10 | Safe Exponential (PN Junction Limiting) | `safe_exp()` in spice_models.c | O(1) | ✅ |

## L6: Canonical Problems

| # | Problem | Example / Test | Status |
|---|---------|---------------|--------|
| 1 | Resistive Voltage Divider (DC) | `example_divider_dc.c`, `test_dc_analysis_resistive_divider` | ✅ |
| 2 | Series Resistor Network (DC) | `test_dc_analysis_series_resistors` | ✅ |
| 3 | RC Circuit Transient Response | `example_rc_transient.c`, `test_transient_analysis_rc` | ✅ |
| 4 | RLC Bandpass Filter (AC) | `example_rlc_ac.c` | ✅ |
| 5 | RC Lowpass Filter (AC) | `test_ac_analysis_rc` | ✅ |
| 6 | Diode DC Biasing | `test_diode_model` | ✅ |
| 7 | MOSFET Operating Point | `test_mosfet_model` | ✅ |
| 8 | BJT Forward-Active Region | `test_bjt_model` | ✅ |

## L7: Applications

| # | Application | Status |
|---|------------|--------|
| 1 | Circuit Simulation for IC Design (SPICE-compatible) | ✅ Core engine |
| 2 | Filter Frequency Response (Bode plots via CSV export) | ✅ `spice_ac_export_csv()` |
| 3 | Transient Waveform Generation (oscilloscope-style CSV) | ✅ `spice_tran_export_csv()` |
| 4 | Device Characterization (I-V curves) | Partial — diode/MOSFET/BJT evaluation |

## L8: Advanced Topics

| # | Topic | Implementation | Status |
|---|-------|---------------|--------|
| 1 | Adaptive Time-Step Control | `spice_transient_analysis()` dt adjustment | ✅ |
| 2 | Gmin Stepping for Convergence | `build_mna_dc()` | ✅ |
| 3 | PN Junction Limiting (Safe Exp) | `safe_exp()` | ✅ |
| 4 | Condition Number Estimation | `spice_condition_number()` | ✅ |

## L9: Research Frontiers

| # | Topic | Status |
|---|-------|--------|
| 1 | Parallel SPICE (multi-threaded matrix solve) | Documented, not implemented |
| 2 | ML-Assisted Convergence Prediction | Documented, not implemented |
| 3 | GPU-Accelerated Circuit Simulation | Documented, not implemented |

---

## Summary

| Level | Coverage | Rating |
|-------|----------|--------|
| L1   | 19/19 definitions | **Complete** |
| L2   | 10/10 core concepts | **Complete** |
| L3   | 6/6 math structures | **Complete** |
| L4   | 8/8 fundamental laws | **Complete** |
| L5   | 10/10 algorithms | **Complete** |
| L6   | 8/8 canonical problems | **Complete** |
| L7   | 3/4 applications | **Partial+** |
| L8   | 4/5 advanced topics | **Partial+** |
| L9   | 3 documented | **Partial** |

**Total Score**: 2×6 + 1×2 + 0×1 = **14/18** → COMPLETE when L1-L6 all Complete.
