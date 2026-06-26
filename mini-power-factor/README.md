# mini-power-factor

**AC Power Factor Analysis — Complex Power, Harmonic Distortion, and Power Quality**

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (16 core definitions with C structs + Lean types)
- **L2 Core Concepts**: Complete (10 concepts: power triangle, PF classification, ITIC curve)
- **L3 Mathematical Structures**: Complete (14 structures: complex power, Fortescue, Clarke, Park, FFT)
- **L4 Fundamental Laws**: Complete (9 laws: power conservation, Boucherot, Steinmetz, Parseval, IEEE 1459/1366)
- **L5 Algorithms/Methods**: Complete (23 algorithms: PF calculation, PFC sizing, FFT, Goertzel, SRF-PLL)
- **L6 Canonical Problems**: Complete (6 end-to-end problems with industrial PFC, harmonics, 3-phase PQ)
- **L7 Applications**: Partial+ (5 applications: data center 80 PLUS, EV charger SAE J2894, EPRI PQ cost)
- **L8 Advanced Topics**: Partial+ (5 topics: Monte Carlo PF uncertainty, ACMC PFC, SRF-PLL under disturbance)
- **L9 Research Frontiers**: Partial (4 frontiers documented: SiC/GaN PFC, AI PQ prediction, digital twin, quantum metrology)

## Code Metrics

| Directory | Files | Lines |
|-----------|-------|-------|
| include/  | 6     | 2,428 |
| src/      | 7     | 3,482 |
| **Total** | **13** | **5,910** |

## Core Definitions (L1)

- **pf_single_phase_t**: Real power P [W], reactive Q [VAR], apparent S [VA], PF, phase angle φ
- **pf_three_phase_t**: Three-phase balanced/unbalanced power with per-phase decomposition
- **pf_ieee1459_t**: IEEE 1459-2010 non-sinusoidal power decomposition (P₁, Q₁, D_I, D_V, S_H, S_N)
- **hp_spectrum_t**: Complete harmonic spectrum with THD, TDD, distortion PF, K-factor
- **phasor_t**: Phasor in polar (magnitude, angle) + rectangular (Re, Im) form
- **phasor_abc_t**: Three-phase phasor set (ABC frame) with balance detection
- **phasor_012_t**: Symmetrical components (zero, positive, negative sequence)
- **dq0_t**: Park transform dq0 reference frame quantities
- **pq_event_t**: IEEE 1159-2019 power quality event with ITIC curve classification
- **pfc_capacitor_bank_t**: Automatic capacitor bank with step configuration and detuning
- **pfc_boost_state_t**: Active boost PFC controller state (ACMC)
- **srf_pll_t**: SRF-PLL state for three-phase grid synchronization

## Core Theorems (L4)

### Power Factor Definition
```
PF = P / |S| = cos(φ)  [sinusoidal]
PF_true = PF_displacement × PF_distortion  [non-sinusoidal, IEEE 1459]
```

### Power Triangle (Steinmetz, 1893)
```
S² = P² + Q²
S = V_rms × I_rms  [VA]
P = V_rms × I_rms × cos(φ)  [W]
Q = V_rms × I_rms × sin(φ)  [VAR]
```

### Complex Power
```
S = V_phasor × I_phasor* = P + jQ
|S| = V_rms × I_rms
φ = atan2(Q, P)
```

### Boucherot's Theorem (1896)
```
Σ Q_k = 0  (for all branches at same frequency)
```
Reactive power is conserved across a network with sinusoidal sources.

### Parseval's Theorem for Power
```
P_time = (1/N) Σ v[n]×i[n] = Σ P_h  (h=1..∞)
```
Total real power in time domain equals sum of per-harmonic powers.

### IEEE 1459-2010 Power Decomposition
```
S² = P₁² + Q₁² + D_I² + D_V² + S_H²
S_N² = D_I² + D_V² + S_H²  (non-fundamental apparent power)
```

### Fortescue Symmetrical Components (1918)
```
[V0]   [1   1   1 ] [Va]
[V1] = [1   a   a²] [Vb] × 1/3   where a = 1∠120°
[V2]   [1   a²  a ] [Vc]
```

### ITIC (CBEMA) Voltage Tolerance Curve
```
Upper envelope: 500% @ 1ms → 200% @ 10ms → 120% @ 0.5s → 110% steady
Lower envelope: 0% @ 20ms → 70% @ 0.5s → 80% @ 10s → 90% steady
```

### IEEE 1366-2012 Reliability Indices
```
SAIFI = Σ(customers interrupted) / total customers
SAIDI = Σ(customer-minutes) / total customers
MAIFI = Σ(momentary events) / total customers
```

## Core Algorithms (L5)

1. **Single-Phase PF**: From RMS + phase angle, and from time-domain samples
2. **Phase Lag Detection**: Cross-correlation peak finding between V and I
3. **Three-Phase Power**: Balanced (√3·V_LL·I_L·cosφ) and unbalanced
4. **Symmetrical Components**: Fortescue transform (ABC→012) and inverse
5. **Clarke Transform**: ABC→αβ0 stationary frame
6. **Park Transform**: ABC→dq0 rotating frame and inverse
7. **Radix-2 DIT FFT**: Cooley-Tukey (1965) in-place O(N log N)
8. **Goertzel Algorithm**: Single-bin DFT detector O(2N) per bin
9. **THD/TDD/PWHD**: Total harmonic/distortion/weighted distortion per IEEE 519
10. **Capacitor Sizing**: Single-phase and three-phase (delta/wye) for PF correction
11. **Automatic Step Bank**: Binary-weighted capacitor step design
12. **Detuning Reactor**: Resonance avoidance with anti-resonance detuning
13. **Boost PFC ACMC**: Average current mode control with PI cascaded loops
14. **SRF-PLL**: Grid synchronization with frequency/angle/magnitude estimation
15. **K-Factor**: Transformer harmonic derating per IEEE C57.110
16. **IEEE 519 Compliance**: TDD + individual harmonic limit checking
17. **ITIC Curve**: Voltage tolerance evaluation for equipment ride-through
18. **Flicker Pst/Plt**: Short-term and long-term per IEC 61000-4-15
19. **SAIFI/SAIDI/MAIFI**: Distribution reliability indices per IEEE 1366
20. **Sliding RMS / EMA RMS**: Real-time power quality monitoring windows
21. **Crest Factor / Form Factor**: Waveform shape characterization
22. **Harmonic Resonance Magnification**: Capacitor-system interaction at PCC
23. **Industrial PFC Solution**: End-to-end sizing + resonance + payback analysis

## Canonical Problems (L6)

1. **Industrial PF Correction** (`example_industrial_pfc.c`): 500kW motor load, PF 0.70→0.95, 480V/60Hz, with capacitor sizing, resonance check, detuning assessment, step bank design, and economic payback
2. **Harmonic Power Analysis** (`example_harmonic_analysis.c`): Rectifier load with IEEE 1459 decomposition, DFT spectrum, THD, distortion PF, K-factor transformer derating, IEEE 519 compliance, and Parseval verification
3. **Three-Phase PQ Analysis** (`example_three_phase_pq.c`): Motor bus with unbalanced voltage, Fortescue components, Clarke/Park transforms, SRF-PLL synchronization, PF correction, data center PQ, and EPRI cost analysis

## Reference Curriculum Mapping (L7)

| School | Key Courses | Topics Covered |
|--------|------------|----------------|
| **MIT** | 6.061 Electric Power, 6.131 Power Electronics | Complex power, PF correction, boost PFC |
| **Stanford** | EE251 Power Electronics, EE292 Smart Grid | PFC ACMC, PQ monitoring, demand response |
| **Berkeley** | EE137A Power Electronics | PF standards, harmonic compliance |
| **Illinois** | ECE 431 Electric Machinery, ECE 464 Power Electronics | dq0 transforms, active PFC |
| **Michigan** | EECS 418 Power Electronics | Automotive PFC (EV chargers) |
| **Georgia Tech** | ECE 4330 Power Electronics | PFC topologies, efficiency optimization |
| **TU Munich** | Power Electronics, Power Systems | IEC harmonic standards, unbalance |
| **ETH Zurich** | 227-0517 Power Electronics, 227-0526 Power Systems | High-freq PFC, flicker measurement |
| **清华** | 电力电子技术, 电力系统分析 | PF correction, sequence components |

## Build and Test

```bash
make          # Build static library (libpower_factor.a)
make test     # Run 47 assertion-based tests (all must pass)
make examples # Build 3 end-to-end example programs
make clean    # Remove build artifacts
```

## References

- Steinmetz, C.P. "Complex Quantities and Their Use in Electrical Engineering" (1893)
- Fortescue, C.L. "Method of Symmetrical Co-ordinates" (AIEE Trans., 1918)
- Park, R.H. "Two-Reaction Theory of Synchronous Machines" (AIEE Trans., 1929)
- Cooley, J.W. & Tukey, J.W. "An Algorithm for the Machine Calculation of Complex Fourier Series" (Math. Comp., 1965)
- IEEE Std 1459-2010: Definitions for the Measurement of Electric Power Quantities
- IEEE Std 519-2014: Recommended Practice for Harmonic Control in Electric Power Systems
- IEEE Std 1159-2019: Recommended Practice for Monitoring Electric Power Quality
- IEEE Std 1366-2012: Guide for Electric Power Distribution Reliability Indices
- IEEE Std C57.110-2018: Recommended Practice for Establishing Liquid-Filled and Dry-Type Transformer Capability When Supplying Nonsinusoidal Load Currents
- IEC 61000-3-2: Limits for Harmonic Current Emissions (Equipment ≤ 16A)
- IEC 61000-4-15: Flickermeter — Functional and Design Specifications
- IEC 61000-4-30: Power Quality Measurement Methods
- Erickson, R.W. & Maksimovic, D. "Fundamentals of Power Electronics" (2nd ed., 2001) — Ch.18 PFC
- Arrillaga, J. & Watson, N.R. "Power System Harmonics" (2nd ed., 2003)
- Bollen, M.H.J. "Understanding Power Quality Problems: Voltage Sags and Interruptions" (2000)
- SAE J2894: Power Quality Requirements for Plug-In Electric Vehicle Chargers
- Kaura, V. & Blasko, V. "Operation of a Phase Locked Loop System Under Distorted Utility Conditions" (IEEE Trans. IA, 1997)
