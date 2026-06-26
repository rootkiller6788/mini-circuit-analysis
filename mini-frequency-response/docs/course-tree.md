# Course Tree — mini-frequency-response

## Prerequisite Knowledge Dependencies

```
                      ┌─────────────────────┐
                      │  Complex Numbers    │
                      │  (L3 Math Struct)   │
                      └─────────┬───────────┘
                                │
              ┌─────────────────┼─────────────────┐
              ▼                 ▼                   ▼
    ┌─────────────────┐ ┌──────────────┐  ┌──────────────────┐
    │ Laplace Transform│ │ Polynomial   │  │ Frequency Domain  │
    │ H(s) = N(s)/D(s) │ │ Arithmetic   │  │ H(jω) = |H|∠φ    │
    │ (L4 Fundamental) │ │ (L3 Math)    │  │ (L2 Concept)      │
    └────────┬────────┘ └──────┬───────┘  └────────┬─────────┘
             │                 │                    │
    ┌────────┴────────┐        │           ┌────────┴─────────┐
    │ Pole-Zero       │        │           │ Frequency         │
    │ Analysis        │        │           │ Response Data     │
    │ (L5 Algorithm)  │        │           │ (L1 Definition)   │
    └────────┬────────┘        │           └────────┬─────────┘
             │                 │                    │
    ┌────────┴─────────────────┴────────────────────┴─────────┐
    │                                                          │
    │              BODE PLOT CONSTRUCTION (L5)                  │
    │              Exact + Asymptotic Methods                   │
    │                                                          │
    └──────┬───────────────────┬───────────────────┬───────────┘
           │                   │                   │
    ┌──────┴──────┐   ┌────────┴────────┐  ┌───────┴──────────┐
    │ Filter      │   │ Stability       │  │ Resonance       │
    │ Design      │   │ Analysis        │  │ Analysis        │
    │ (L5/L6)     │   │ (L4/L5/L6)      │  │ (L6)            │
    └──────┬──────┘   └────────┬────────┘  └───────┬──────────┘
           │                   │                    │
    ┌──────┴──────┐   ┌────────┴────────┐  ┌───────┴──────────┐
    │ Active RC   │   │ Routh-Hurwitz   │  │ Series RLC      │
    │ Sallen-Key  │   │ Nyquist         │  │ Parallel RLC    │
    │ MFB, Tow-Th │   │ Root Locus      │  │ Coupled Res.    │
    └─────────────┘   └─────────────────┘  └──────────────────┘
```

## Dependency Rules

1. **Complex numbers** are the mathematical foundation for ALL frequency-domain analysis.
2. **Laplace transform** bridges time-domain circuits to s-domain transfer functions.
3. **Transfer function H(s)** is the central object; everything else derives from it.
4. **Pole-zero analysis** reveals stability and frequency response character.
5. **Bode plots** are the primary visualization and design tool.
6. **Filter design**, **stability analysis**, and **resonance analysis** are three independent branches that all depend on transfer functions and Bode plots.

## Knowledge Dependencies (What You Must Know First)

| This Module Requires | From Module |
|---------------------|-------------|
| Complex algebra | Prerequisite math |
| Polynomial arithmetic | Prerequisite math |
| Laplace transform basics | mini-signal-system-theory |
| Circuit analysis basics (KVL, KCL) | mini-dc-ac-circuit |
| Impedance concept (R, L, C) | mini-circuit-topology |

## What Depends on This Module

| Depends On This | Module |
|----------------|--------|
| Analog filter design | mini-analog-electronics |
| Communication system frequency planning | mini-communication-principle |
| Digital filter design (prototype transformations) | mini-digital-signal-process |
| RF amplifier design (S-parameters, stability) | mini-wireless-mobile-comm |
| Oscillator design (resonance, crystal model) | mini-analog-electronics |
