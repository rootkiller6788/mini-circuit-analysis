# Course Tree — mini-two-port-network

## Prerequisites (前置依赖)

```
DC Circuit Analysis (Kirchhoff's Laws, Ohm's Law)
├── AC Circuit Analysis (phasors, complex impedance)
│   ├── Matrix Algebra (linear equations, inversion)
│   │   └── **Two-Port Network Theory** ← THIS MODULE
│   └── Laplace Transform (s-domain analysis)
└── Semiconductor Physics (PN junction, transistor operation)
    └── Small-Signal Models (hybrid-π, T-model)
```

## This Module's Position in the Tree

```
mini-two-port-network
│
├── Provides foundation for:
│   ├── mini-analog-electronics (amplifier analysis)
│   ├── mini-communication-principle (RF chain analysis)
│   ├── mini-digital-signal-process (filter design)
│   ├── mini-electromagnetic-wave (transmission line theory)
│   ├── mini-wireless-mobile-comm (RF front-end design)
│   └── mini-emc-signal-integrity (signal path analysis)
│
├── Depends on:
│   ├── Complex number arithmetic (junior-level math)
│   ├── Linear algebra (matrix operations)
│   ├── Basic circuit theory (KVL, KCL, Ohm's law)
│   └── Frequency-domain analysis (phasors, impedance)
│
└── Parallel to:
    ├── mini-signal-system-theory (transfer function = two-port)
    └── mini-control-automation (feedback = interconnected two-ports)
```

## Learning Path

1. **Start here**: `two_port.h` — complex numbers and 2×2 matrices (L3)
2. **Then**: `z_params.h` — impedance parameters, T-network (L1, L2)
3. **Then**: `y_params.h` — admittance parameters, π-network, MOSFET model (L1, L2)
4. **Then**: `h_params.h` — BJT h-parameter model, CE/CB/CC (L1, L6)
5. **Then**: `g_params.h` — FET inverse hybrid, source follower (L1, L6)
6. **Then**: `abcd_params.h` — transmission parameters, cascade, TL (L1, L4)
7. **Then**: `s_params.h` — scattering parameters, stability metrics (L1, L4)
8. **Advanced**: `conversion.h` — all 30 conversions (L5)
9. **Advanced**: `interconnection.h` — five topologies + Bartlett (L5)
10. **Advanced**: `stability.h` — full stability analysis suite (L4)
11. **Advanced**: `network_synthesis.h` — Cauer, Foster, Darlington (L5)
12. **Capstone**: `filter_design.h` — complete filter design pipeline (L6)
