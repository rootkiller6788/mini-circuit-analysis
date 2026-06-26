# Course Tree — mini-spice-simulation

## Prerequisites (Concepts this module depends on)

```
Linear Algebra
  ├── Matrix arithmetic (LU decomposition)
  ├── Linear system solving (forward/back substitution)
  └── Eigenvalues (condition number estimation)

Calculus
  ├── Derivatives (Newton-Raphson, device linearization)
  ├── Numerical integration (Trapezoidal rule, companion models)
  └── Complex numbers (AC phasor analysis)

Circuit Theory
  ├── Ohm's Law (V = I * R)
  ├── Kirchhoff's Current Law (KCL: ΣI = 0 at node)
  ├── Kirchhoff's Voltage Law (KVL: ΣV = 0 around loop)
  ├── Capacitor: I = C * dV/dt
  ├── Inductor: V = L * dI/dt
  └── Nodal Analysis → Modified Nodal Analysis (MNA)

Semiconductor Physics
  ├── PN junction → Shockley diode equation
  ├── MOS capacitor → MOSFET Shichman-Hodges model
  └── Bipolar junction → Ebers-Moll model

Numerical Methods
  ├── Gaussian elimination
  ├── Partial pivoting
  ├── Newton-Raphson iteration
  └── Convergence criteria (absolute + relative tolerance)

Programming
  ├── C11 (struct, union, dynamic allocation)
  ├── Makefile build system
  └── CSV data export
```

## Dependents (Concepts that build on this module)

```
1. mini-circuit-analysis (parent module) → uses SPICE for circuit design
2. mini-analog-electronics → DC/AC/TRAN for amplifier design
3. mini-digital-electronics → TRAN for gate delay, DC for logic levels
4. mini-communication-principle → AC for filter design, TRAN for modulation
5. mini-wireless-mobile-comm → RF circuit simulation
6. mini-optical-fiber-comm → transimpedance amplifier simulation
7. mini-radar-remote-sensing → RF front-end simulation
8. mini-emc-signal-integrity → parasitic extraction, TRAN for crosstalk
```

## Learning Path

```
Step 1: Understand Ohm's Law, KCL, KVL
  └── Step 2: Learn nodal analysis (manually solve small circuits)
      └── Step 3: Study MNA formulation (adds V/I sources and inductors)
          └── Step 4: Implement linear DC solver (Gaussian elimination)
              └── Step 5: Add nonlinear devices (Newton-Raphson)
                  └── Step 6: Add dynamic elements (companion models, TRAN)
                      └── Step 7: Add frequency domain (complex MNA, AC)
                          └── Step 8: Device modeling (MOSFET, BJT)
                              └── Step 9: Production SPICE features (convergence, time-step)
```

## L9: Research Frontiers

Future extensions to this module include:
- **Parallel SPICE**: Multi-threaded LU decomposition and device evaluation
- **ML-Assisted Convergence**: Neural network predictor for initial guess
- **GPU Acceleration**: CUDA-based sparse matrix solve for large circuits
- **Quantum Circuit Simulation**: Extension to quantum-classical co-simulation
