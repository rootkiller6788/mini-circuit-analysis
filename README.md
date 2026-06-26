# 1. mini-circuit-analysis

**Circuit Analysis — from-scratch, zero-dependency C implementations of
electrical network theory, spanning DC/AC steady-state, transient dynamics,
frequency-domain methods, and SPICE-compatible simulation.**

---

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-circuit-topology](mini-circuit-topology/) | Graph-theoretic foundations, MNA formulation, KCL/KVL, Tellegen's theorem, sparse solvers | MIT 6.002, Stanford EE101A |
| [mini-dc-ac-circuit](mini-dc-ac-circuit/) | Ohm's law, nodal/mesh analysis, complex impedance, phasors, resonance, coupled inductors | MIT 6.002, Berkeley EE16A |
| [mini-frequency-response](mini-frequency-response/) | Bode plots, transfer functions, analog filter design (Butterworth, Chebyshev), stability criteria | MIT 6.003, Stanford EE102A |
| [mini-network-theorem](mini-network-theorem/) | Superposition, Thévenin, Norton, Millman, maximum power transfer, reciprocity, substitution | MIT 6.002, Berkeley EE16B |
| [mini-power-factor](mini-power-factor/) | Complex power, power triangle, harmonics & THD, PF correction, phasor operations, power quality | MIT 6.061, Stanford EE253 |
| [mini-spice-simulation](mini-spice-simulation/) | SPICE netlist parser, MNA matrix assembly, device models (diode, BJT, MOSFET), DC/AC/TRAN analysis | MIT 6.002, Berkeley EE105 |
| [mini-transient-analysis](mini-transient-analysis/) | First/second-order transients, state-space methods, numerical ODE solvers (Euler, RK4, BDF2, Gear), switching transients | MIT 6.002, Stanford EE101B |
| [mini-two-port-network](mini-two-port-network/) | Z, Y, H, G, ABCD, S parameters; all 30 directional conversions; interconnection; network synthesis; stability | MIT 6.003, Stanford EE101B |

---

## Design Philosophy

1. **Zero dependencies** — Every sub-module compiles with only a C99/C11 compiler and the standard
   library. No external math or linear-algebra packages are required; all solvers are hand-written
   for pedagogical transparency.

2. **Layered learning** — Each sub-module follows an L1→L9 progression: Definitions → Core Concepts →
   Mathematical Structures → Fundamental Laws → Algorithms → Canonical Problems → Applications →
   Advanced Topics → Research Frontiers. This mirrors a complete undergraduate-to-graduate curriculum.

3. **From equation to code** — Every function maps directly to a textbook equation. The governing
   physical law (Ohm, Maxwell, Faraday) is cited in the source, so readers can trace each line of C
   back to first principles.

4. **Inter-module composability** — Shared type definitions (`circuit_elements.h`, `transient_defs.h`,
   `spice_core.h`) allow sub-modules to compose. The topology module's MNA matrix feeds directly into
   the SPICE simulator and the DC/AC solver, forming a miniature yet complete EDA pipeline.

---

## Building

Build all sub-modules at once:

```bash
# From this directory
for d in mini-*/; do (cd "$d" && make); done
```

Or build individual modules:

```bash
cd mini-circuit-topology && make && make test
cd mini-dc-ac-circuit && make && make test
cd mini-frequency-response && make && make test
cd mini-network-theorem && make && make test
cd mini-power-factor && make && make test
cd mini-spice-simulation && make && make test
cd mini-transient-analysis && make && make test
cd mini-two-port-network && make && make test
```

---

## Project Structure

```
1. mini-circuit-analysis/
├── mini-circuit-topology/          # Graph theory, MNA, KCL/KVL, Tellegen's theorem, sparse LU
├── mini-dc-ac-circuit/             # DC/AC steady-state, Ohm's law, nodal/mesh, complex power
├── mini-frequency-response/        # Bode plots, transfer functions, filter design, stability
├── mini-network-theorem/           # Superposition, Thévenin/Norton, maximum power transfer
├── mini-power-factor/              # Complex power, harmonics, THD, PF correction, power quality
├── mini-spice-simulation/          # Netlist parsing, MNA stamping, device models, multi-mode analysis
├── mini-transient-analysis/        # First/second-order, state-space, numerical ODE, switching transients
├── mini-two-port-network/          # Z/Y/H/G/ABCD/S parameters, 30 conversions, synthesis, stability
├── .gitignore                      # Build artifact and IDE exclusion rules
├── README.md                       # This file (English)
└── README-CN.md                    # Chinese version
```

---

## License

MIT — See each sub-module's source headers for details.
