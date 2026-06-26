# Knowledge Graph — mini-two-port-network

## L1: Definitions (Complete)
- Complex number (rectangular & polar form)
- 2×2 complex matrix (universal two-port container)
- Z-parameters (open-circuit impedance): V = Z·I
- Y-parameters (short-circuit admittance): I = Y·V
- H-parameters (hybrid): V1=h11·I1+h12·V2, I2=h21·I1+h22·V2
- G-parameters (inverse hybrid): I1=g11·V1+g12·I2, V2=g21·V1+g22·I2
- ABCD-parameters (transmission/cascade): V1=A·V2+B·(-I2)
- S-parameters (scattering): b = S·a (power waves)
- Input/output impedance, voltage/current/power gain
- Reflection coefficient Γ = (Z-Z0)/(Z+Z0)
- VSWR = (1+|Γ|)/(1-|Γ|)
- Return loss RL = -20·log10(|Γ|)
- T-network and π-network equivalent circuits
- Image impedance and image transfer constant
- Rollett stability factor K
- μ stability factor (Edwards-Sinsky)
- Stability circles on Smith chart
- Noise figure NF and noise parameters (NFmin, Rn, Γopt)

## L2: Core Concepts (Complete)
- Reciprocity (z12=z21, y12=y21, s12=s21, AD-BC=1)
- Symmetry (z11=z22, etc.)
- Passivity (Re(Z)≥0, positive semi-definite Hermitian part)
- Losslessness (pure imaginary Z/Y, unitary S)
- Unilateral vs bilateral networks
- Miller's theorem and Miller effect
- Conjugate matching for max power transfer
- Impedance transformation through two-ports
- Load-pull effect (Γin depends on ΓL via s12)
- Even/odd mode decomposition (Bartlett's theorem)
- Brune's test for valid interconnections
- Image parameter filter theory
- Stability classification (unconditional/conditional/unstable)

## L3: Mathematical Structures (Complete)
- Complex number arithmetic (add, sub, mul, div, mag, arg, conj, exp, polar, sqrt)
- 2×2 matrix operations (add, sub, mul, det, trace, inv, transpose, Hermitian)
- Matrix inversion and singular case detection
- Euler's formula: e^(jθ) = cos(θ) + j·sin(θ)
- Hyperbolic functions for transmission lines: cosh(γl), sinh(γl)
- Chebyshev polynomials T_N(ω) — defined piecewise for |ω|≤1 and |ω|>1
- Bessel polynomials (reverse coefficients for Thomson filter)
- Continued fraction expansion (Cauer synthesis)
- Positive-real function theory
- Jacobi elliptic functions (documented, elliptic filter approximation implemented)

## L4: Fundamental Laws (Complete)
- Maximum Power Transfer Theorem (Jacobi, 1840): ZL = ZS*
- Reciprocity Theorem (Lorentz): response symmetric for source/load swap
- Miller's Theorem (1920): Z1 = Zf/(1-Av), Z2 = Zf·Av/(Av-1)
- Bartlett's Bisection Theorem (1927): z11=(Zoc+Zsc)/2, z12=(Zoc-Zsc)/2
- Foster's Reactance Theorem (1924): dX/dω > 0 for lossless one-ports
- Rollett Stability Criterion (1962): K > 1 and |Δ| < 1
- Edwards-Sinsky μ-criterion (1992): μ > 1 ↔ unconditional stability
- Darlington Synthesis Theorem (1939): any PR impedance = lossless 2-port + 1Ω
- Kennelly Δ-Y Transform (1899): T ↔ π conversion generalized to complex Z
- Nyquist Stability Criterion: applied to two-port loop gain

## L5: Algorithms/Methods (Complete)
- 30-directional parameter conversion (Z↔Y↔H↔G↔ABCD↔S)
- T-network extraction from Z-parameters
- π-network extraction from Y-parameters
- T↔π conversion (generalized Δ-Y)
- Cauer ladder synthesis (continued fraction)
- L-network matching design (single-frequency)
- π-network matching design (loaded Q)
- T-network matching design
- Foster-I and Foster-II canonical form evaluation
- Darlington synthesis step (element extraction)
- Butterworth g-value computation: gk = 2·sin((2k-1)π/(2N))
- Chebyshev g-value computation (with ripple)
- LC ladder ABCD construction (cascade alternating series/shunt)
- Sallen-Key biquad design (ω0, Q → R1, R2, C1, C2)
- MFB biquad design
- Coupled-line filter even/odd impedance computation
- Stability circle center and radius computation
- Conjugate match optimization (ΓMS, ΓML)
- Noise figure computation from noise parameters
- De-embedding (remove fixture effects)

## L6: Canonical Problems (Complete)
- BJT CE amplifier analysis (H-parameters → Av, Zin, Zout)
- BJT CB and CC configuration conversion
- FET CS amplifier analysis (G-parameters → intrinsic gain)
- Source follower output impedance
- Butterworth filter design and frequency response
- Chebyshev-I/II filter design
- Bessel filter design (linear phase)
- Elliptic filter magnitude response
- Propagation constant and image impedance of cascade
- Stability improvement: series R or shunt G addition
- Oscillation margin computation
- Quarter-wave transformer: Zin = Z0²/ZL
- Half-wave impedance repeater: Zin = ZL
- Transmission line ABCD (lossy and lossless)

## L7: Applications (Partial — 4 implementations)
- GPS L1 (1.575 GHz) receiver chain cascade analysis
- 2.4 GHz WiFi LNA matching network design
- 100 MHz LC filter design for RF front-end
- BJT amplifier analysis for audio/RF

## L8: Advanced Topics (Partial)
- Noise figure optimization (NFmin, Rn, Γopt trade-off with gain)
- Broadband matching (π-network with selectable Q)
- Numerical stability improvement via resistive loading
- Group delay analysis for Bessel vs Butterworth comparison

## L9: Research Frontiers (Documented)
- mmWave CMOS two-port modeling (gm·ro degradation in sub-10nm)
- RIS (Reconfigurable Intelligent Surface) as tunable two-port arrays
- Sub-THz transistor modeling beyond quasi-static approximation
