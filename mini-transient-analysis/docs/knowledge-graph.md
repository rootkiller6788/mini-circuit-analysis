# Knowledge Graph - mini-transient-analysis

## L1: Definitions (Complete)
- Time constant (tau = RC, tau = L/R)
- Damping classes: overdamped, critically damped, underdamped, undamped
- Transient response types: natural, forced, complete, step, impulse, ramp
- Transient performance metrics: rise time, settling time, overshoot, peak time
- System order classification: first, second, third, higher-order
- Numerical solver types: Euler, RK4, RK45, trapezoidal, BDF2
- BJT switching regions: cutoff, forward active, saturation
- MOSFET regions: cutoff, linear, saturation, subthreshold
- Component parasitic models: ESR, ESL, DCR, EPR

## L2: Core Concepts (Complete)
- First-order RC/RL time constant tau
- Second-order characteristic parameters: zeta, omega_n, omega_d, alpha
- Damping ratio determines transient behavior
- Complete response = natural + forced (superposition)
- Natural response decays according to system poles
- Forced response follows input waveform shape
- Energy storage: E_c = 0.5*C*V^2, E_l = 0.5*L*I^2
- Continuity conditions: v_c(0+)=v_c(0-), i_l(0+)=i_l(0-)

## L3: Mathematical Structures (Complete)
- First-order ODE: tau*dy/dt + y = f(t)
- Second-order ODE: d2y/dt2 + 2*zeta*omega_n*dy/dt + omega_n^2*y = f(t)
- Characteristic equation: s^2 + 2*zeta*omega_n*s + omega_n^2 = 0
- Pole classification by discriminant
- State-space representation: dx/dt = A*x + B*u
- Complex frequency: s = sigma + j*omega
- Transfer function to state-space conversion
- Eigenvalue analysis (QR simplified)

## L4: Fundamental Laws (Complete)
- Ohm's Law: v = R*i
- Faraday's Law: v = L*di/dt (inductor voltage)
- Maxwell's displacement current: i = C*dv/dt (capacitor current)
- KVL in time domain for series RLC
- KCL in time domain for parallel RLC
- Energy conservation in RLC circuits
- Superposition principle for linear circuits
- Initial condition continuity laws
- Routh-Hurwitz stability criterion
- Time constant properties (63.2% rule, 5*tau rule)

## L5: Algorithms/Methods (Complete)
- Forward Euler (explicit, O(h))
- Backward Euler (implicit, A-stable, O(h))
- Trapezoidal rule (Crank-Nicolson, O(h^2), A-stable)
- RK2 Midpoint method
- RK2 Heun method
- RK3 Classical method
- RK4 Classical method (O(h^4))
- RK45 Dormand-Prince adaptive step-size control
- Adams-Bashforth 2-step (explicit)
- Adams-Moulton 2-step (implicit)
- BDF2/BDF3 Gear's methods for stiff systems
- Predictor-corrector (PECE)
- Newton-Raphson for implicit methods
- Stability region analysis
- Stiffness detection and classification
- Richardson extrapolation
- Step-size control algorithms
- Order verification via refinement

## L6: Canonical Problems (Complete)
- RC step response (charging/discharging)
- RL step response
- RC impulse response
- RL impulse response
- RC ramp response
- RC pulse response (single and train)
- RC sinusoidal complete response
- RLC series step response (all three damping regimes)
- RLC parallel step response
- RLC series/parallel natural response
- RLC series/parallel impulse response
- RLC sinusoidal complete response
- Sequential switching response
- Cascaded RC response
- Transmission line bounce diagram
- Third-order step response
- Coupled oscillator transient

## L7: Applications (Complete, 20+ applications)
- 555 timer astable/monostable timing
- RC relaxation oscillator
- Power-On Reset (POR) timeout
- Watchdog timer timeout
- Debounce circuit settling time
- Buck converter: startup inrush, load step, ripple, soft-start
- Boost converter: startup, ripple, RHP zero
- MOSFET gate driver: delay, rise/fall time, Miller plateau
- Diode reverse recovery: Qrr, Irr, softness factor
- Snubber circuit design (RC snubber)
- Relay/solenoid coil transient
- DC motor startup and PWM ripple
- H-bridge dead-time calculation
- HBM/CDM ESD transients
- TVS clamping voltage
- Power sequencing delay
- Inrush current limiter
- Transmission line: reflection, propagation delay, ringing
- Transformer inrush current
- Capacitor lifetime derating

## L8: Advanced Topics (Partial, 5 topics)
- PLL transient: natural frequency, damping, lock time
- Time-varying RLC parameters
- Stiff ODE systems and Gear's method
- State-space eigenvalue sensitivity
- Nonlinear circuit element models

## L9: Research Frontiers (Documented)
- Model order reduction (Pade, modal)
- Hankel singular values for system reduction
- Coupled oscillator transient analysis
- Distributed RC transmission line modeling
