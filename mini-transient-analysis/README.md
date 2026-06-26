# mini-transient-analysis

Transient circuit analysis module covering time-domain behavior of linear and nonlinear circuits with energy storage elements.

## Module Status: COMPLETE ✅

| Level | Name | Status |
|-------|------|--------|
| L1 | Definitions | Complete |
| L2 | Core Concepts | Complete |
| L3 | Mathematical Structures | Complete |
| L4 | Fundamental Laws | Complete |
| L5 | Algorithms/Methods | Complete |
| L6 | Canonical Problems | Complete |
| L7 | Applications | Complete (20+ applications) |
| L8 | Advanced Topics | Partial (5 topics) |
| L9 | Research Frontiers | Partial (documented) |

**Score: 16/18 — COMPLETE**

## Quick Build

```bash
make        # build library, tests, examples
make test   # run test suite (46 checks, 23 test groups)
make examples  # run example programs
```

## Core Definitions (L1)

| Concept | Symbol | Formula |
|---------|--------|---------|
| RC time constant | τ | R × C |
| RL time constant | τ | L / R |
| Damping ratio | ζ | (R/2)√(C/L) series, (1/2R)√(L/C) parallel |
| Natural frequency | ωₙ | 1/√(LC) |
| Damped frequency | ω_d | ωₙ√(1-ζ²) |
| Neper frequency | α | R/(2L) series, 1/(2RC) parallel |
| Quality factor | Q | (1/R)√(L/C) series, R√(C/L) parallel |
| Capacitor energy | E_C | ½CV² |
| Inductor energy | E_L | ½LI² |

## Core Theorems (L4)

1. **RC Step Response**: v_c(t) = Vs + (v0 - Vs)·e^(-t/τ)
2. **RL Step Response**: i_l(t) = Vs/R + (i0 - Vs/R)·e^(-Rt/L)
3. **Overdamped RLC**: v_c(t) = Vss + A₁e^(s₁t) + A₂e^(s₂t), s₁,₂ = -α ± √(α² - ωₙ²)
4. **Critically Damped RLC**: v_c(t) = Vss + (A₁ + A₂t)e^(-αt)
5. **Underdamped RLC**: v_c(t) = Vss + e^(-αt)(A₁cos(ω_d t) + A₂sin(ω_d t))
6. **Percent Overshoot**: PO = 100·e^(-πζ/√(1-ζ²))
7. **Settling Time (2%)**: t_s = 4/(ζωₙ)
8. **Peak Time**: t_p = π/(ωₙ√(1-ζ²))

## Core Algorithms (L5)

| Method | Order | Stability | Use Case |
|--------|-------|-----------|----------|
| Forward Euler | O(h) | Conditional | Simple, fast |
| Backward Euler | O(h) | A-stable | Stiff systems |
| Trapezoidal | O(h²) | A-stable | General purpose |
| RK4 | O(h⁴) | Conditional | High accuracy |
| RK45 (D-P) | O(h⁵) | Conditional | Adaptive step |
| Gear BDF2 | O(h²) | Stiff-stable | Very stiff |

## Canonical Problems (L6)

- RC step, impulse, ramp, pulse, sinusoidal response
- RL step, impulse, ramp response
- RLC series: all three damping regimes
- RLC parallel: step and natural response
- Sequential switching
- Cascaded RC (two stages)
- Transmission line bounce diagram
- Third-order dominant pole approximation
- Coupled oscillator transient

## Nine-School Course Mapping

| School | Key Course | Topics |
|--------|-----------|--------|
| MIT | 6.002 Circuits | RC/RL/RLC transients |
| Stanford | EE101A | First/second-order analysis |
| Berkeley | EE16A/B, EE105 | State-space, BJT/MOSFET switching |
| Illinois | ECE 210, 342 | Laplace, switching transients |
| Michigan | EECS 215, 411 | RLC, transmission lines |
| Georgia Tech | ECE 2040, 3040 | Transient response, device switching |
| TU Munich | Circuit Theory, PE | Transient analysis, converters |
| ETH Zurich | 227-0056, 227-0116 | Numerical methods, PE dynamics |
| Tsinghua | Circuit Theory, PE | Full transient, converter design |

## File Structure

```
mini-transient-analysis/
├── Makefile              # make test builds and runs tests
├── README.md             # This file
├── include/              # 6 header files
│   ├── transient_defs.h      # L1: Core type definitions
│   ├── first_order.h         # L2-L4: First-order RC/RL
│   ├── second_order.h        # L2-L4: Second-order RLC
│   ├── higher_order.h         # L3-L5: State-space, eigenvalues
│   ├── numerical_transient.h  # L5: ODE solver methods
│   └── switching_transient.h  # L6-L7: Device switching
├── src/                  # 7 source files + Lean
│   ├── transient_defs.c      # 556 lines
│   ├── first_order.c         # 709 lines
│   ├── second_order.c        # 837 lines
│   ├── higher_order.c         # 753 lines
│   ├── numerical_transient.c  # 814 lines
│   ├── switching_transient.c  # 640 lines
│   └── transient_formal.lean  # Lean 4 formalization
├── tests/                # test_transient.c (23 test groups)
├── examples/             # 3 runnable examples
│   ├── example_rc_step.c
│   ├── example_rlc_response.c
│   └── example_switching.c
├── docs/                 # 5 knowledge documents
│   ├── knowledge-graph.md
│   ├── coverage-report.md
│   ├── gap-report.md
│   ├── course-alignment.md
│   └── course-tree.md
├── benches/
└── demos/
```

## References

- Hayt, Kemmerly & Durbin, "Engineering Circuit Analysis", 9th Ed, 2019
- Nilsson & Riedel, "Electric Circuits", 11th Ed, 2019
- Ogata, "Modern Control Engineering", 5th Ed, 2010
- Mohan, Undeland & Robbins, "Power Electronics", 3rd Ed, 2003
- Erickson & Maksimovic, "Fundamentals of Power Electronics", 2nd Ed, 2001
- Press et al., "Numerical Recipes in C", 3rd Ed, 2007
- Hairer, Norsett & Wanner, "Solving Ordinary Differential Equations I", 1993
- Sedra & Smith, "Microelectronic Circuits", 8th Ed, 2020
