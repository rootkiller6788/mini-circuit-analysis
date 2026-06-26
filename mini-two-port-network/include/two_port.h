/**
 * @file two_port.h
 * @brief Core Two-Port Network Definitions
 *
 * Two-port network theory is the foundation of modern circuit analysis,
 * providing a unified framework for analyzing linear time-invariant (LTI)
 * circuits with two electrical ports (four terminals).
 *
 * Reference textbooks:
 * - "Microelectronic Circuits" Sedra & Smith (2020) - Ch. 9: Frequency Response
 * - "Microwave Engineering" Pozar (2012) - Ch. 4: Microwave Network Analysis
 * - "Network Analysis and Synthesis" Kuo (2006)
 *
 * Course mapping (L1-L9):
 *   L1: Z/Y/H/G/ABCD/S parameter definitions
 *   L2: Reciprocity, symmetry, passivity, losslessness concepts
 *   L3: Complex 2x2 matrix algebra, eigen-decomposition
 *   L4: Reciprocity theorem, Maximum power transfer theorem
 *   L5: Parameter conversion algorithms
 *   L6: Transistor amplifier, filter design, matching network
 *   L7: RF front-end design
 *   L8: Noise analysis, broadband matching
 *   L9: mmWave integration (documented)
 *
 * Schools covered:
 *   MIT 6.002/6.003, Stanford EE101/EE114, Berkeley EE105,
 *   Illinois ECE 451, Michigan EECS 411, Georgia Tech ECE 6350,
 *   TU Munich HF Engineering, ETH 227-0455, Tsinghua 电路原理
 */

#ifndef TWO_PORT_H
#define TWO_PORT_H

#include <stddef.h>
#include <complex.h>
#include <math.h>

/* M_PI is not standard in C99; define if not provided */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * L1: Core Definitions — Complex Numbers and 2x2 Matrices
 * ============================================================================ */

/**
 * @brief Complex number representation in rectangular form.
 *
 * We use double precision for engineering accuracy (IEEE 754 double).
 * The rectangular form (Re + j*Im) is used for Z, Y, H, G, ABCD parameters.
 * For S-parameters, we also support polar form through helper functions.
 *
 * Course: MIT 6.002 — Phasors and complex impedance
 *         Stanford EE101 — Complex exponential representation
 */
typedef struct {
    double real;   /**< Real part, in appropriate units (Ω, S, dimensionless) */
    double imag;   /**< Imaginary part, in appropriate units */
} complex_t;

/**
 * @brief 2x2 complex matrix — the universal two-port parameter container.
 *
 * Convention:
 *   [m11  m12]   where m11 = input parameter, m22 = output parameter,
 *   [m21  m22]   m12 = reverse transfer, m21 = forward transfer.
 *
 * All six parameter types (Z, Y, H, G, ABCD, S) use this structure.
 * The physical meaning of each element depends on the parameter type.
 *
 * Course: Berkeley EE105 — Small-signal two-port models
 */
typedef struct {
    complex_t m11;  /**< Row 1, Col 1 — input-side parameter */
    complex_t m12;  /**< Row 1, Col 2 — reverse transfer parameter */
    complex_t m21;  /**< Row 2, Col 1 — forward transfer parameter */
    complex_t m22;  /**< Row 2, Col 2 — output-side parameter */
} matrix2x2_t;

/* ============================================================================
 * L2: Core Concepts — Parameter Types
 * ============================================================================ */

/**
 * @brief Enumeration of six fundamental two-port parameter types.
 *
 * Z-parameters (Impedance): V1 = z11*I1 + z12*I2, V2 = z21*I1 + z22*I2
 *   Open-circuit impedance parameters. Units: Ω.
 *   Used for: series-connected networks, T-equivalent circuits.
 *
 * Y-parameters (Admittance): I1 = y11*V1 + y12*V2, I2 = y21*V1 + y22*V2
 *   Short-circuit admittance parameters. Units: S (siemens).
 *   Used for: parallel-connected networks, π-equivalent circuits.
 *
 * H-parameters (Hybrid): V1 = h11*I1 + h12*V2, I2 = h21*I1 + h22*V2
 *   Mixed units: h11 in Ω, h22 in S, h12/h21 dimensionless.
 *   Used for: BJT small-signal models (h_ie, h_fe, h_re, h_oe).
 *   Course: Berkeley EE105 — BJT h-parameter model.
 *
 * G-parameters (Inverse Hybrid): I1 = g11*V1 + g12*I2, V2 = g21*V1 + g22*I2
 *   Mixed units: g11 in S, g22 in Ω, g12/g21 dimensionless.
 *   Used for: FET small-signal models.
 *
 * ABCD-parameters (Transmission/Cascade): V1 = A*V2 - B*I2, I1 = C*V2 - D*I2
 *   Mixed units: A/D dimensionless, B in Ω, C in S.
 *   Used for: cascade analysis, transmission lines, filters.
 *   Note: The sign convention uses -I2 for output current flowing out.
 *   Course: TU Munich — HF Engineering (Vierpoltheorie)
 *
 * S-parameters (Scattering): b1 = s11*a1 + s12*a2, b2 = s21*a1 + s22*a2
 *   Dimensionless complex ratios of power waves.
 *   Used for: RF/microwave circuits, network analyzers (VNA).
 *   Course: ETH 227-0455 — Microwave Circuit Design
 */
typedef enum {
    PARAM_Z,      /**< Impedance (open-circuit) parameters */
    PARAM_Y,      /**< Admittance (short-circuit) parameters */
    PARAM_H,      /**< Hybrid parameters (BJT model) */
    PARAM_G,      /**< Inverse hybrid parameters (FET model) */
    PARAM_ABCD,   /**< Transmission (cascade) parameters */
    PARAM_S       /**< Scattering parameters (RF/microwave) */
} param_type_t;

/**
 * @brief Two-port network configuration types.
 *
 * Series-Series:   Used with Z-parameters (additive Z).
 * Parallel-Parallel: Used with Y-parameters (additive Y).
 * Series-Parallel: Used with H-parameters (additive H).
 * Parallel-Series: Used with G-parameters (additive G).
 * Cascade:         Used with ABCD-parameters (multiplicative ABCD).
 *
 * Course: Illinois ECE 451 — Network interconnection rules
 */
typedef enum {
    CONFIG_SERIES_SERIES,       /**< Input series, output series → add Z */
    CONFIG_PARALLEL_PARALLEL,   /**< Input parallel, output parallel → add Y */
    CONFIG_SERIES_PARALLEL,     /**< Input series, output parallel → add H */
    CONFIG_PARALLEL_SERIES,     /**< Input parallel, output series → add G */
    CONFIG_CASCADE              /**< Cascade → multiply ABCD */
} config_type_t;

/**
 * @brief Stability classification for active two-ports.
 *
 * Based on Rollett's stability factor K and the auxiliary condition.
 *
 * Course: Georgia Tech ECE 6350 — Amplifier stability
 * Reference: Rollett, "Stability and Power-Gain Invariants of Linear Twoports", IRE, 1962
 */
typedef enum {
    STABLE_UNCONDITIONAL,   /**< |Γin| < 1 and |Γout| < 1 for all passive terminations */
    STABLE_CONDITIONAL,     /**< Stable only for certain load/source impedances */
    UNSTABLE                /**< Potential for oscillation */
} stability_t;

/* ============================================================================
 * Complex number arithmetic (L3: Mathematical Structures)
 * ============================================================================ */

/**
 * @brief Create a complex number from real and imaginary parts.
 * @param real Real part
 * @param imag Imaginary part
 * @return complex_t with given components
 *
 * Mathematical: z = a + jb  where j = √(-1)
 * IEEE convention: real part = a, imaginary part = b
 */
complex_t complex_make(double real, double imag);

/**
 * @brief Complex addition: (a+jb) + (c+jd) = (a+c) + j(b+d)
 * Complexity: O(1)
 */
complex_t complex_add(complex_t a, complex_t b);

/**
 * @brief Complex subtraction: (a+jb) - (c+jd) = (a-c) + j(b-d)
 * Complexity: O(1)
 */
complex_t complex_sub(complex_t a, complex_t b);

/**
 * @brief Complex multiplication: (a+jb)(c+jd) = (ac-bd) + j(ad+bc)
 * Complexity: O(1), 4 multiplications + 2 additions
 */
complex_t complex_mul(complex_t a, complex_t b);

/**
 * @brief Complex division: (a+jb)/(c+jd) = (ac+bd)/(c²+d²) + j(bc-ad)/(c²+d²)
 * @return complex_t result, or {NAN, NAN} if denominator is zero
 * Complexity: O(1), 6 multiplications + 2 divisions + 5 additions
 */
complex_t complex_div(complex_t a, complex_t b);

/**
 * @brief Complex magnitude: |a+jb| = sqrt(a² + b²)
 * Complexity: O(1)
 */
double complex_mag(complex_t z);

/**
 * @brief Complex phase angle: arg(a+jb) = atan2(b, a) in radians
 * Returns value in [-π, π].
 * Complexity: O(1)
 */
double complex_arg(complex_t z);

/**
 * @brief Complex conjugate: conj(a+jb) = a - jb
 * Complexity: O(1)
 */
complex_t complex_conj(complex_t z);

/**
 * @brief Complex exponential: e^(a+jb) = e^a * (cos(b) + j*sin(b))
 * Used for phasor calculations.
 * Course: MIT 6.003 — Euler's formula: e^(jθ) = cos(θ) + j sin(θ)
 * Complexity: O(1)
 */
complex_t complex_exp(complex_t z);

/**
 * @brief Create complex number from polar form: r * e^(jθ)
 * @param r magnitude (non-negative)
 * @param theta angle in radians
 * @return r*cos(theta) + j*r*sin(theta)
 * Complexity: O(1)
 */
complex_t complex_polar(double r, double theta);

/**
 * @brief Complex square root using principal branch.
 * sqrt(a+jb) = sqrt((|z|+a)/2) + j*sign(b)*sqrt((|z|-a)/2)
 * Complexity: O(1)
 */
complex_t complex_sqrt(complex_t z);

/* ============================================================================
 * 2x2 Matrix Operations (L3: Mathematical Structures)
 * ============================================================================ */

/**
 * @brief Create a 2x2 matrix from four complex elements.
 */
matrix2x2_t matrix2x2_make(complex_t m11, complex_t m12,
                            complex_t m21, complex_t m22);

/**
 * @brief Matrix addition: C = A + B, element-wise.
 * Used for: adding Z-parameters of series-series connected networks.
 * Complexity: O(1), 4 complex additions
 */
matrix2x2_t matrix2x2_add(matrix2x2_t a, matrix2x2_t b);

/**
 * @brief Matrix subtraction: C = A - B, element-wise.
 * Complexity: O(1), 4 complex subtractions
 */
matrix2x2_t matrix2x2_sub(matrix2x2_t a, matrix2x2_t b);

/**
 * @brief Matrix multiplication: C = A × B.
 *
 * [c11 c12]   [a11 a12]   [b11 b12]
 * [c21 c22] = [a21 a22] × [b21 b22]
 *
 * c11 = a11*b11 + a12*b21, c12 = a11*b12 + a12*b22
 * c21 = a21*b11 + a22*b21, c22 = a21*b12 + a22*b22
 *
 * Used for: cascading ABCD-parameter networks.
 * Complexity: O(1), 8 complex multiplications + 4 complex additions
 *
 * Course: TU Munich HF Engineering — Cascade theory
 */
matrix2x2_t matrix2x2_mul(matrix2x2_t a, matrix2x2_t b);

/**
 * @brief Matrix scalar multiplication: B = k * A.
 * Complexity: O(1), 4 complex multiplications
 */
matrix2x2_t matrix2x2_scale(complex_t k, matrix2x2_t a);

/**
 * @brief Matrix determinant: det(A) = a11*a22 - a12*a21.
 *
 * Used for: parameter conversion formulas, stability analysis.
 * Example: Δ = z11*z22 - z12*z21 (determinant of Z-matrix)
 * Complexity: O(1), 2 complex multiplications + 1 complex subtraction
 */
complex_t matrix2x2_det(matrix2x2_t a);

/**
 * @brief Matrix trace: tr(A) = a11 + a22.
 * Complexity: O(1)
 */
complex_t matrix2x2_trace(matrix2x2_t a);

/**
 * @brief Matrix inverse: A^(-1) = (1/det(A)) * [[a22, -a12], [-a21, a11]].
 *
 * @return Inverse matrix, or all-NAN if determinant is zero (singular).
 * Used for: converting Z↔Y (Z = Y^(-1) and Y = Z^(-1)).
 *
 * Mathematical basis:
 * For a 2x2 matrix, the inverse formula is exact:
 *   A^(-1) = [1/(a11*a22 - a12*a21)] * [ a22  -a12 ]
 *                                       [ -a21  a11 ]
 *
 * Course: Berkeley EE16A — Matrix inversion for circuit analysis
 * Complexity: O(1)
 */
matrix2x2_t matrix2x2_inv(matrix2x2_t a);

/**
 * @brief Matrix transpose: A^T = [[a11, a21], [a12, a22]].
 * Complexity: O(1)
 */
matrix2x2_t matrix2x2_transpose(matrix2x2_t a);

/**
 * @brief Matrix Hermitian conjugate: A^H = [[a11*, a21*], [a12*, a22*]].
 * Used for: S-parameter analysis, power calculations.
 * Complexity: O(1)
 */
matrix2x2_t matrix2x2_hermitian(matrix2x2_t a);

/**
 * @brief Check if matrix is symmetric: A = A^T.
 * Symmetric two-ports satisfy z12 = z21 (reciprocity in Z).
 * Complexity: O(1)
 */
int matrix2x2_is_symmetric(matrix2x2_t a, double tolerance);

/**
 * @brief Check if matrix is reciprocal.
 * For Z: z12 = z21; For Y: y12 = y21; For S: s12 = s21.
 * For ABCD: det(ABCD) = AD - BC = 1 (reciprocal network condition).
 * Complexity: O(1) or O(1) with determinant.
 */
int matrix2x2_is_reciprocal_z(matrix2x2_t z, double tolerance);
int matrix2x2_is_reciprocal_y(matrix2x2_t y, double tolerance);
int matrix2x2_is_reciprocal_abcd(matrix2x2_t abcd, double tolerance);

/* ============================================================================
 * L2: Core Concepts — Two-Port Properties
 * ============================================================================ */

/**
 * @brief Compute input impedance of a two-port terminated with load ZL.
 *
 * Formula: Zin = z11 - (z12*z21)/(z22 + ZL)
 *
 * This is the fundamental formula for impedance transformation through
 * a two-port network. Derived from the Z-parameter equations:
 *   V1 = z11*I1 + z12*I2
 *   V2 = z21*I1 + z22*I2
 * with termination: V2 = -ZL*I2 (negative sign due to current direction).
 *
 * @param z Z-parameter matrix
 * @param zl Load impedance (real for resistive, complex for reactive)
 * @return Input impedance at port 1
 *
 * Course: Stanford EE101 — Impedance transformation and matching
 * Ref: Sedra & Smith §9.4 — Miller's theorem and input impedance
 */
complex_t two_port_zin_from_zparams(matrix2x2_t z, complex_t zl);

/**
 * @brief Compute output impedance of a two-port driven by source ZS.
 *
 * Formula: Zout = z22 - (z12*z21)/(z11 + ZS)
 *
 * Application: Determining the Thevenin equivalent at the output port.
 * Also used for output matching in amplifier design.
 *
 * @param z Z-parameter matrix
 * @param zs Source impedance
 * @return Output impedance at port 2
 */
complex_t two_port_zout_from_zparams(matrix2x2_t z, complex_t zs);

/**
 * @brief Compute forward voltage gain Av = V2/V1 with load ZL.
 *
 * Formula using Z-parameters: Av = (z21 * ZL) / ((z11 * (z22 + ZL)) - (z12 * z21))
 *
 * Or equivalently using ABCD: Av = ZL / (A*ZL + B)
 *
 * This is the loaded voltage gain, different from the open-circuit gain.
 *
 * @param z Z-parameter matrix
 * @param zl Load impedance
 * @return Complex voltage gain (magnitude + phase shift)
 */
complex_t two_port_voltage_gain_z(matrix2x2_t z, complex_t zl);

/**
 * @brief Compute forward current gain Ai = I2/I1 with load ZL.
 *
 * Formula using Z-parameters: Ai = -z21 / (z22 + ZL)
 *
 * @param z Z-parameter matrix
 * @param zl Load impedance
 * @return Complex current gain
 */
complex_t two_port_current_gain_z(matrix2x2_t z, complex_t zl);

/**
 * @brief Compute power gain Gp = P_load / P_input.
 *
 * Power gain considers actual power delivered to load vs power into input.
 * Gp = |Av|^2 * Re(ZL) / Re(Zin) or equivalently
 * Gp = |s21|^2 * (1-|ΓL|^2) / (|1-s22*ΓL|^2 * (1-|Γin|^2))
 *
 * @param z Z-parameter matrix
 * @param zl Load impedance
 * @return Power gain (dimensionless, ≥ 0)
 */
double two_port_power_gain_z(matrix2x2_t z, complex_t zl);

/**
 * @brief Compute available power gain Ga = P_avail_load / P_avail_source.
 *
 * Available power gain is the maximum power gain achievable when both
 * ports are conjugately matched. Important for amplifier design.
 *
 * @param s S-parameter matrix
 * @return Available power gain
 */
double two_port_available_gain_s(matrix2x2_t s);

/**
 * @brief Compute transducer power gain Gt = P_load / P_avail_source.
 *
 * Transducer gain includes the effect of both source and load mismatches.
 * Gt = |s21|^2 * (1-|ΓS|^2)(1-|ΓL|^2) / |(1-s11*ΓS)(1-s22*ΓL)-s12*s21*ΓS*ΓL|^2
 *
 * @param s S-parameter matrix
 * @param gs Source reflection coefficient
 * @param gl Load reflection coefficient
 * @return Transducer power gain
 */
double two_port_transducer_gain_s(matrix2x2_t s, complex_t gs, complex_t gl);

/* ============================================================================
 * L2: Reflection Coefficient and VSWR
 * ============================================================================ */

/**
 * @brief Compute reflection coefficient Γ from impedance Z and reference Z0.
 *
 * Formula: Γ = (Z - Z0) / (Z + Z0)
 *
 * This is the fundamental bridge between impedance-based and wave-based
 * analysis. |Γ| = 1 for complete reflection (open/short), Γ = 0 for match.
 *
 * @param z Impedance
 * @param z0 Reference impedance (typically 50Ω for RF)
 * @return Complex reflection coefficient
 *
 * Course: Georgia Tech ECE 6350 — Reflection and transmission
 * Ref: Pozar §2.3 — The reflection coefficient
 */
complex_t reflection_coefficient(complex_t z, double z0);

/**
 * @brief Compute VSWR from reflection coefficient magnitude.
 *
 * Formula: VSWR = (1 + |Γ|) / (1 - |Γ|)
 *
 * VSWR (Voltage Standing Wave Ratio) ranges from 1 (perfect match) to ∞
 * (complete reflection). In practice, VSWR < 2 (return loss > 10 dB) is
 * considered acceptable for most applications.
 *
 * @param gamma Reflection coefficient
 * @return VSWR (≥ 1.0), or INFINITY if |Γ| = 1
 *
 * Course: ETH 227-0455 — Standing waves and VSWR
 */
double vswr_from_gamma(complex_t gamma);

/**
 * @brief Compute return loss from reflection coefficient.
 *
 * Formula: RL(dB) = -20 * log10(|Γ|)
 *
 * Return loss quantifies how much power is reflected. Higher RL = better match.
 * RL > 10 dB → |Γ| < 0.316 → VSWR < 1.92
 * RL > 20 dB → |Γ| < 0.1   → VSWR < 1.22
 *
 * @param gamma Reflection coefficient
 * @return Return loss in dB (≥ 0), or INFINITY for perfect match
 */
double return_loss_db(complex_t gamma);

/**
 * @brief Compute impedance from reflection coefficient and reference Z0.
 *
 * Inverse formula: Z = Z0 * (1 + Γ) / (1 - Γ)
 *
 * Used to convert VNA measurements (S-parameters) back to impedances.
 *
 * @param gamma Reflection coefficient
 * @param z0 Reference impedance
 * @return Impedance
 */
complex_t impedance_from_gamma(complex_t gamma, double z0);

/**
 * @brief Check if a passive termination is valid (non-negative real part).
 *
 * A passive impedance must satisfy Re(Z) ≥ 0. This is the physical realizability
 * condition — no negative resistance without active elements.
 *
 * @param z Impedance
 * @return 1 if passive (Re(Z) ≥ 0), 0 if active (Re(Z) < 0)
 */
int is_passive_impedance(complex_t z);

/**
 * @brief Check if a two-port network is reciprocal.
 *
 * For reciprocal networks: z12 = z21, y12 = y21, s12 = s21.
 * All passive networks (R, L, C only) are reciprocal.
 * Networks containing dependent sources may be non-reciprocal.
 *
 * @param m Parameter matrix
 * @param type Parameter type
 * @param tolerance Numerical tolerance (e.g., 1e-9)
 * @return 1 if reciprocal, 0 otherwise
 *
 * Course: Illinois ECE 451 — Reciprocity theorem
 */
int two_port_is_reciprocal(matrix2x2_t m, param_type_t type, double tolerance);

/**
 * @brief Check if a two-port network is lossless.
 *
 * A lossless network has zero real power dissipation.
 * For S-parameters: S^H * S = I (unitary condition).
 * For Z/Y: Re(Z) = 0 or Re(Y) = 0 (pure imaginary).
 *
 * @param m Parameter matrix
 * @param type Parameter type
 * @param tolerance Numerical tolerance
 * @return 1 if lossless, 0 otherwise
 */
int two_port_is_lossless(matrix2x2_t m, param_type_t type, double tolerance);

/**
 * @brief Check if a two-port network is symmetrical.
 *
 * A symmetrical network has identical input and output ports:
 * z11 = z22, y11 = y22, s11 = s22, A = D (for reciprocal ABCD).
 *
 * @param m Parameter matrix
 * @param type Parameter type
 * @param tolerance Numerical tolerance
 * @return 1 if symmetrical, 0 otherwise
 */
int two_port_is_symmetrical(matrix2x2_t m, param_type_t type, double tolerance);

/* ============================================================================
 * L2: Miller's Theorem
 * ============================================================================ */

/**
 * @brief Apply Miller's theorem to split an impedance.
 *
 * Miller's theorem: An impedance Z connected between two nodes with
 * voltage gain Av can be split into:
 *   Z1 = Z / (1 - Av)  (at input)
 *   Z2 = Z * Av / (Av - 1)  (at output)
 *
 * This is the cornerstone of high-frequency amplifier analysis,
 * explaining the Miller effect in common-emitter amplifiers.
 *
 * @param z Bridging impedance
 * @param av Voltage gain (Vout/Vin, complex for AC analysis)
 * @param z1 Output: equivalent input-side impedance
 * @param z2 Output: equivalent output-side impedance
 *
 * Course: Berkeley EE105 — Miller effect in BJT amplifiers
 * Ref: Sedra & Smith §9.4
 */
void miller_split_impedance(complex_t z, complex_t av,
                             complex_t *z1, complex_t *z2);

/**
 * @brief Compute Miller input capacitance from feedback capacitance.
 *
 * Cin_miller = Cf * (1 + |Av|)  for inverting amplifier
 * Cin_miller = Cf * (1 - Av)     for non-inverting (Av < 1 reduces cap)
 *
 * This formula explains why the input capacitance of a CE amplifier
 * is dominated by Cμ (collector-base capacitance) multiplied by gain.
 *
 * @param cf Feedback capacitance in Farads
 * @param av Voltage gain (real for mid-band)
 * @return Miller capacitance at input
 *
 * Course: Stanford EE114 — Frequency response of amplifiers
 */
double miller_input_capacitance(double cf, double av);

/* ============================================================================
 * L4: Maximum Power Transfer Theorem
 * ============================================================================ */

/**
 * @brief Find the conjugate match impedance for maximum power transfer.
 *
 * Theorem (Jacobi, 1840): Maximum power is transferred from source to load
 * when the load impedance equals the complex conjugate of the source impedance:
 *   ZL = ZS*
 *
 * This implies: RL = RS and XL = -XS (reactive cancellation).
 *
 * @param zsource Source impedance
 * @return Optimal load impedance for maximum power transfer
 *
 * Course: MIT 6.002 — Maximum power transfer
 * Ref: Desoer & Kuh, "Basic Circuit Theory", Ch. 9
 */
complex_t conjugate_match(complex_t zsource);

/**
 * @brief Compute maximum available power from a source.
 *
 * Pavail = |Vs|^2 / (8 * Re(ZS))
 *
 * This is the theoretical maximum power a source can deliver when
 * conjugately matched. Represents the "power of the source" independent
 * of loading.
 *
 * @param vs_mag Magnitude of Thevenin source voltage
 * @param zs Source impedance
 * @return Maximum available power in watts
 */
double maximum_available_power(double vs_mag, complex_t zs);

/* ============================================================================
 * L5: Numerical Utilities
 * ============================================================================ */

/**
 * @brief Check if two complex numbers are approximately equal.
 * @param a First complex number
 * @param b Second complex number
 * @param tolerance Absolute tolerance (e.g., 1e-9)
 * @return 1 if |a-b| < tolerance, 0 otherwise
 */
int complex_approx_equal(complex_t a, complex_t b, double tolerance);

/**
 * @brief Check if a complex number is zero within tolerance.
 */
int complex_is_zero(complex_t z, double tolerance);

/**
 * @brief Convert complex to string representation (for debugging).
 * Buffer must be at least 64 bytes.
 */
void complex_to_string(complex_t z, char *buf, size_t buf_size);

/**
 * @brief Print a 2x2 matrix to stdout (for debugging).
 */
void matrix2x2_print(matrix2x2_t m, const char *label);

#endif /* TWO_PORT_H */
