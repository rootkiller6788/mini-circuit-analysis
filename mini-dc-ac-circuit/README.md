# mini-dc-ac-circuit

## Module Status: COMPLETE

- L1: Complete - Core definitions
- L2: Complete - Core concepts  
- L3: Complete - Mathematical structures
- L4: Complete - Fundamental laws/theorems
- L5: Complete - Algorithms/methods
- L6: Complete - Canonical problems
- L7: Complete - Applications
- L8: Partial - Advanced topics
- L9: Partial - Research frontiers

## Code Metrics
- include/ + src/ total lines: 3028 (C + H)
- Lean 4 formalization: 355 lines
- Headers: 6 files, Sources: 6 files
- Tests: 4 files, Examples: 3 files

## Core Theorems (L4)
1. Ohm's Law (1827): V = I*R
2. Joule's Law (1841): P = I^2*R
3. Kirchhoff's Laws (1845): KCL + KVL
4. Thevenin's Theorem (1883)
5. Norton's Theorem (1926)
6. Maximum Power Transfer (Jacobi 1840)
7. Tellegen's Theorem (1952)
8. Reciprocity Theorem
9. Friis Noise Formula (1944)

## Key Formulas
- Ohm's Law: V = I*R
- Voltage Divider: Vout = Vin*R2/(R1+R2)
- RLC Impedance: Z = R + j(wL - 1/wC)
- Resonance: w0 = 1/sqrt(LC), Q = w0*L/R
- Max Power: RL = Rth, Pmax = Vth^2/(4*Rth)
- RC Tau: tau = R*C
- 555 Astable: f = 1.44/((R1+2R2)C)

## Build
```
make        # Build and test
make test   # Run tests
make clean  # Clean
```

## Reference Texts
- Hayt, Kemmerly & Durbin: Engineering Circuit Analysis
- Sedra & Smith: Microelectronic Circuits
- Nilsson & Riedel: Electric Circuits
