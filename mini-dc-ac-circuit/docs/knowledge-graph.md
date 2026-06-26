# Knowledge Graph ? mini-dc-ac-circuit

## L1 Definitions
- Voltage, Current, Resistance (SI units)
- Capacitor: i = C*dv/dt (Maxwell 1861)
- Inductor: v = L*di/dt (Faraday 1831)
- Impedance Z = R + jX, Admittance Y = G + jB
- Phasor: magnitude-angle representation (Steinmetz 1893)
- Circuit topology: Node, Branch, Loop
- Element types: R, C, L, Vsrc, Isrc, dependent sources
- Multi-terminal: Op-Amp, Transformer, Gyrator

## L2 Core Concepts
- Ohm's Law relationship V=IR
- Voltage and Current divider rules
- Series and parallel equivalent elements
- Complex power S = P + jQ
- Power factor and power factor correction
- Time constant tau = RC, tau = L/R
- Resonance and quality factor
- Frequency response and Bode plots

## L3 Mathematical Structures
- Complex numbers for AC analysis
- Matrix algebra for circuit equations
- Gaussian elimination, LU decomposition
- Eigenvalues and eigenvectors
- SVD for 2x2 matrices
- Differential equations for transient analysis
- Laplace transform concepts

## L4 Fundamental Laws
- Ohm's Law (1827)
- Joule's Law (1841)
- Kirchhoff's Current Law (1845)
- Kirchhoff's Voltage Law (1845)
- Thevenin's Theorem (1883)
- Norton's Theorem (1926)
- Maximum Power Transfer (Jacobi 1840)
- Tellegen's Theorem (1952)
- Reciprocity Theorem
- Friis Noise Formula (1944)

## L5 Algorithms/Methods
- Modified Nodal Analysis (MNA)
- Mesh Analysis
- Gaussian elimination with partial pivoting
- LU decomposition (Doolittle)
- Cholesky decomposition
- QR decomposition (Gram-Schmidt)
- Newton-Raphson nonlinear DC solver
- Monte Carlo tolerance analysis
- Euler and RK4 transient integration
- Power iteration for dominant eigenvalue

## L6 Canonical Problems
- Wheatstone Bridge
- RLC series/parallel resonance
- RC/RL/RLC transient response
- Delta-Wye transformation
- Power factor correction
- Impedance matching (L, Pi, T networks)
- Filter design (RC, RLC, Sallen-Key)
- Bode plot generation

## L7 Applications
- Power supply design (rectifier, regulator, SMPS)
- 555 timer circuits
- Operational amplifier circuits
- RF impedance matching and VSWR
- ADC/DAC quantization analysis
- Temperature sensing (NTC, PT100, thermocouple)
- Strain gauge bridge
- LED current limiting
- Motor driver and solenoid circuits
- Battery charging (Li-Ion, NiMH)
- Automotive electrical (ISO 7637, CAN, LIN)

## L8 Advanced Topics
- Nonlinear device models (diode, BJT, MOSFET)
- Newton-Raphson DC solver
- Monte Carlo tolerance/sensitivity analysis
- SVD for circuit matrices
- Matrix exponential and companion matrix
- Phase noise, oscillator analysis
- S-parameter conversions
- Nonlinear capacitance (varactor)
- Saturable inductance

## L9 Research Frontiers
- Memristor device modeling (documented)
- Quantum circuit elements (documented)
- Neuromorphic circuit design (documented)
