# Knowledge Graph — mini-power-factor

## L1: Definitions (Complete ✅)

| # | Concept | C Implementation | Lean Formalization |
|---|---------|-----------------|-------------------|
| 1 | Real Power P [Watts] | `pf_compute_single_phase()` in power_factor.c | `SinglePhasePower.realPowerMilliWatt` |
| 2 | Reactive Power Q [VAR] | `pf_phasor_power()` in power_factor.c | `SinglePhasePower.reactivePowerMilliVAR` |
| 3 | Apparent Power S [VA] | `pf_single_phase_t.s_apparent` | `SinglePhasePower.apparentPowerMilliVA` |
| 4 | Power Factor PF = P/S | `pf_classify()` | `powerFactorFromPQ` |
| 5 | Displacement PF (cos φ₁) | `hp_decompose_pf()` in harmonic_power.c | — |
| 6 | Distortion PF | `hp_decompose_pf()` | `distortionPF` |
| 7 | Total Harmonic Distortion (THD) | `hp_compute_thd()` | `thdSimple` |
| 8 | Total Demand Distortion (TDD) | `hp_compute_tdd()` | — |
| 9 | Crest Factor | `pf_crest_factor()` | `crestFactorPermille` |
| 10 | Form Factor | `pf_form_factor()` | `formFactorPermille` |
| 11 | K-Factor (transformer derating) | `hp_k_factor()` | `kFactor` |
| 12 | IEEE 1159 Event Types | `pq_event_type_t` | — |
| 13 | Phasor (magnitude + angle) | `phasor_t` | — |
| 14 | Symmetrical Components (012) | `phasor_012_t` | `symmetricalComponentsBalanced` |
| 15 | dq0 Reference Frame | `dq0_t` | — |
| 16 | PWHD (Partial Weighted HD) | `hp_compute_pwhd()` | — |

## L2: Core Concepts (Complete ✅)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Power Triangle: S² = P² + Q² | `pf_compute_single_phase()` verifies | `powerTriangleHolds` theorem |
| 2 | Leading/Lagging PF | `pf_type_t`, detected by phase lag in `pf_compute_from_samples()` |
| 3 | PF Classification (Good/Fair/Poor/Bad) | `pf_classify()` |
| 4 | Reactive Power Compensation | `pf_required_reactive_comp()` |
| 5 | Phasor Diagram Interpretation | `phasor_t` + phasor arithmetic |
| 6 | ITIC (CBEMA) Voltage Tolerance | `pq_check_itic()`, `pq_itic_upper_envelope()`, `pq_itic_lower_envelope()` |
| 7 | Balanced vs Unbalanced 3-Phase | `pf_compute_three_phase()` with balanced flag |
| 8 | Power in Series Elements (same I) | `cp_series_complex_power()` |
| 9 | Power in Parallel Elements (same V) | `cp_parallel_complex_power()` |
| 10 | Distortion PF Decomposition | `hp_decompose_pf()`: PF_true = PF_disp × PF_dist |

## L3: Mathematical Structures (Complete ✅)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Complex Power S = P + jQ | `cp_complex_power()` in complex_power.c | `ComplexPower` structure |
| 2 | Phasor Rotation e^{jθ} | `cp_rotate()` | — |
| 3 | Fortescue Transform (ABC→012) | `phasor_abc_to_012()`, `pf_symmetrical_components()` |
| 4 | Inverse Fortescue (012→ABC) | `phasor_012_to_abc()` |
| 5 | Clarke Transform (ABC→αβ0) | `phasor_clarke_transform()` |
| 6 | Inverse Clarke (αβ0→ABC) | `phasor_inverse_clarke()` |
| 7 | Park Transform (ABC→dq0) | `phasor_park_transform()` |
| 8 | Inverse Park (dq0→ABC) | `phasor_inverse_park()` |
| 9 | DFT/FFT (Cooley-Tukey) | `hp_fft_radix2()` (radix-2 DIT) |
| 10 | Goertzel Algorithm (single bin) | `hp_goertzel_detect()` |
| 11 | Fourier-Based Harmonic Power | `hp_analyze_spectrum()`, `hp_fft_power()` |
| 12 | Impedance from Power (Z = V²/S*) | `cp_impedance_from_power()` |
| 13 | Admittance from Power (Y = S*/V²) | `cp_admittance_from_power()` |
| 14 | Complex Power from V/Z and I/Z | `cp_power_from_v_z()`, `cp_power_from_i_z()` |

## L4: Fundamental Laws (Complete ✅)

| # | Law/Theorem | C Verification | Lean Statement |
|---|------------|----------------|----------------|
| 1 | Conservation of Real Power | `pf_verify_power_balance()` | `realPowerBalance` |
| 2 | Boucherot's Theorem (ΣQ=0) | `pf_verify_boucherot()` | `boucherotZeroSum` |
| 3 | Steinmetz Complex Power Balance | `cp_verify_complex_power_balance()` | `complexPowerBalance_zero` |
| 4 | Parseval's Theorem for Power | `hp_verify_parseval_power()` | `parsevalAdditivity` |
| 5 | IEEE 1459 Power Decomposition | `pf_ieee1459_t`, `hp_decompose_pf()` | — |
| 6 | Tellegen's Theorem (via power balance) | `pf_verify_power_balance()` cross-checks | — |
| 7 | Maximum Power Transfer (conjugate match) | — | `conjugateMatch_maximizes` |
| 8 | IEEE 1366 Reliability (SAIFI/SAIDI/MAIFI) | `pq_compute_saifi()`, `pq_compute_saidi()`, `pq_compute_maifi()` | — |
| 9 | IEC 61000-4-15 Flicker (Pst/Plt) | `pq_compute_pst()`, `pq_compute_plt()` | — |

## L5: Algorithms and Methods (Complete ✅)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Single-Phase PF from RMS Values | `pf_compute_single_phase()` |
| 2 | Time-Domain PF from Samples | `pf_compute_from_samples()` with cross-correlation |
| 3 | Phase Lag Detection (correlation peak) | `pf_detect_phase_lag()` |
| 4 | Sliding Window RMS | `pf_sliding_rms()` with incremental update |
| 5 | Exponential Moving Average RMS | `pf_ema_rms()` |
| 6 | Capacitor Sizing (Single-Phase) | `pfc_size_capacitor()` |
| 7 | Capacitor Sizing (Three-Phase) | `pfc_size_capacitor_3phase()` |
| 8 | Automatic Step Bank Design | `pfc_design_step_bank()` |
| 9 | Detuning Reactor Evaluation | `pfc_needs_detuning()`, `pfc_resonant_frequency()` |
| 10 | Boost PFC Average Current Mode | `pfc_boost_init()`, `pfc_boost_control_step()` |
| 11 | THD from Harmonic Magnitudes | `hp_compute_thd()` |
| 12 | TDD per IEEE 519 | `hp_compute_tdd()` |
| 13 | PWHD (IEC 61000-3) | `hp_compute_pwhd()` |
| 14 | Radix-2 DIT FFT | `hp_fft_radix2()` |
| 15 | Goertzel Single-Bin DFT | `hp_goertzel_detect()` |
| 16 | K-Factor Computation | `hp_k_factor()` |
| 17 | Transformer Harmonic Derating | `hp_transformer_derating()` |
| 18 | IEEE 519 Compliance Check | `hp_check_ieee519()` |
| 19 | SRF-PLL Grid Synchronization | `srf_pll_init()`, `srf_pll_step()` |
| 20 | Capacitor Step Optimization | `pfc_optimize_steps()` with hysteresis |
| 21 | Optimal Shunt Compensation Admittance | `cp_compensation_admittance()` |
| 22 | kVAR Rating Calculation | `cp_kvar_required()` |
| 23 | Demand Interval Computation | `pf_demand_interval()` |

## L6: Canonical Problems (Complete ✅)

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | Industrial PF Correction (end-to-end) | `pfc_solve_industrial()` + example_industrial_pfc.c |
| 2 | Harmonic Spectrum Analysis (rectifier load) | `hp_analyze_spectrum()` + example_harmonic_analysis.c |
| 3 | Harmonic Resonance Risk Assessment | `pfc_harmonic_risk()`, `pfc_needs_detuning()` |
| 4 | Three-Phase Unbalance Analysis | example_three_phase_pq.c |
| 5 | Capacitor Step Switching with Hysteresis | `pfc_optimize_steps()` |
| 6 | Harmonic Source Identification | `hp_dominant_harmonics()` |

## L7: Applications (Partial+ ✅)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | Data Center PQ Assessment (80 PLUS) | `pq_assess_datacenter()` |
| 2 | EV Charger PQ Compliance (SAE J2894) | `pq_assess_ev_charger()` |
| 3 | Cost of Poor PQ (EPRI method) | `pq_cost_of_poor_quality()` |
| 4 | PQ Monitoring Report (IEC 61000-4-30) | `pq_generate_report()` |
| 5 | Three-Phase Capacitor Bank (delta/wye) | `cp_three_phase_pf_capacitor()` |

## L8: Advanced Topics (Partial+ ✅)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Monte Carlo PF Uncertainty | `test_l8_pf_uncertainty_monte_carlo()` in tests |
| 2 | Boost PFC Control (ACMC) | `pfc_boost_control_step()` with PI cascaded loops |
| 3 | SRF-PLL Under Grid Disturbance | `srf_pll_step()` with anti-windup |
| 4 | Harmonic Resonance Magnification | `pfc_harmonic_risk()` |
| 5 | Symmetrical Components for Fault Analysis | `pf_symmetrical_components()` |

## L9: Research Frontiers (Partial ✅)

| # | Topic | Documentation |
|---|-------|--------------|
| 1 | Wide-Bandgap PFC (SiC/GaN) | Documented in course-tree.md |
| 2 | AI-Based Power Quality Prediction | Documented in gap-report.md |
| 3 | Digital Twin for Power Quality | Documented in gap-report.md |
| 4 | Quantum Power Metrology | Documented in coverage-report.md |
