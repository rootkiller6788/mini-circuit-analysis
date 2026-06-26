# Knowledge Prerequisite Tree — mini-power-factor

```
mini-power-factor
│
├── L1: Definitions (prereq: basic AC circuit theory)
│   ├── Real/Reactive/Apparent Power
│   ├── Power Factor (PF = P/S)
│   ├── Phasor representation (magnitude + angle)
│   ├── THD, TDD, Crest Factor, Form Factor
│   └── IEEE 1159 event types
│
├── L2: Core Concepts (prereq: L1)
│   ├── Power Triangle: S² = P² + Q²
│   ├── Leading vs. Lagging PF
│   ├── PF Classification (Good/Fair/Poor/Bad)
│   ├── ITIC (CBEMA) voltage tolerance curve
│   ├── Balanced/Unbalanced three-phase
│   └── Distortion PF decomposition
│
├── L3: Mathematical Structures (prereq: L1, complex numbers)
│   ├── Complex Power: S = V × I* = P + jQ
│   ├── Fortescue Transform (ABC → 012)
│   ├── Clarke/Park Transforms (αβ0, dq0)
│   ├── DFT/FFT (Cooley-Tukey, Goertzel)
│   └── Impedance from Power: Z = V²/S*
│
├── L4: Fundamental Laws (prereq: L1-L3)
│   ├── Conservation of Real Power (ΣP_in = ΣP_out + losses)
│   ├── Boucherot's Theorem (ΣQ_k = 0 at same frequency)
│   ├── Steinmetz Complex Power Balance (ΣS_k = 0)
│   ├── Parseval's Theorem for Power (time domain = frequency domain)
│   ├── IEEE 1459 Power Decomposition
│   ├── IEEE 1366 Reliability Indices (SAIFI/SAIDI/MAIFI)
│   └── IEC 61000-4-15 Flicker (Pst/Plt)
│
├── L5: Algorithms and Methods (prereq: L1-L4)
│   ├── PF from RMS values / Time-domain samples
│   ├── Phase lag detection (cross-correlation)
│   ├── Sliding window / EMA RMS
│   ├── Capacitor sizing (single/three-phase)
│   ├── Automatic capacitor bank step design
│   ├── Detuning reactor resonance analysis
│   ├── Boost PFC average current mode control
│   ├── FFT radix-2 / Goertzel harmonic detection
│   ├── K-Factor / transformer derating
│   ├── IEEE 519 compliance checking
│   ├── SRF-PLL grid synchronization
│   └── Capacitor step switching optimization
│
├── L6: Canonical Problems (prereq: L1-L5)
│   ├── Industrial PF correction (sizing + resonance + payback)
│   ├── Harmonic spectrum analysis (rectifier load)
│   ├── Harmonic resonance risk (capacitor bank + system)
│   ├── Three-phase unbalance (motor derating)
│   └── Capacitor step switching (hysteresis control)
│
├── L7: Applications (prereq: L1-L6)
│   ├── Data center PQ (80 PLUS, Energy Star)
│   ├── EV charger PQ (SAE J2894, IEC 61851)
│   ├── Cost of poor PQ (EPRI methodology)
│   └── PQ monitoring (IEC 61000-4-30)
│
├── L8: Advanced Topics (prereq: L1-L7)
│   ├── Monte Carlo PF uncertainty estimation
│   ├── Boost PFC ACMC with anti-windup
│   ├── SRF-PLL under distorted grid conditions
│   ├── Harmonic resonance magnification factor
│   └── Symmetrical components for fault analysis
│
└── L9: Research Frontiers (documented, not implemented)
    ├── Wide-bandgap PFC (SiC/GaN) — higher f_sw, lower losses
    ├── AI/ML power quality prediction — LSTM, transformer models
    ├── Digital twin for PQ — real-time simulation + monitoring
    └── Quantum power metrology — Josephson voltage standards
```
