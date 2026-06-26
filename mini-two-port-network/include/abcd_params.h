/**
 * @file abcd_params.h
 * @brief ABCD (Transmission) Parameter Analysis
 *
 * ABCD-parameters (also called T-parameters or chain parameters) relate
 * input voltage/current to output voltage/current:
 *
 *   V1 = A * V2 + B * (-I2)    (note: -I2 convention, current out)
 *   I1 = C * V2 + D * (-I2)
 *
 * Matrix form: [[V1], [I1]] = [[A, B], [C, D]] * [[V2], [-I2]]
 *
 * Component units:
 *   A = V1/V2 | I2=0  (reverse voltage ratio, open output) — dimensionless
 *   B = V1/(-I2) | V2=0  (transfer impedance, shorted output) — Ω
 *   C = I1/V2 | I2=0  (transfer admittance, open output) — S
 *   D = I1/(-I2) | V2=0  (reverse current ratio, shorted output) — dimensionless
 *
 * The -I2 convention means current flowing OUT of port 2 (toward load),
 * making ABCD parameters cascade-multiplicative:
 *   ABCD_total = ABCD1 × ABCD2 × ... × ABCDn
 *
 * ABCD parameters are natural for:
 *   - Cascade analysis (the ONLY parameters with this property)
 *   - Transmission line analysis
 *   - Filter design (individual filter sections multiply)
 *   - Microwave circuit design
 *
 * Reciprocal condition: det(ABCD) = AD - BC = 1
 * Symmetrical condition: A = D (for reciprocal networks)
 *
 * Reference: Pozar, "Microwave Engineering" (2012), Ch. 4
 * Course: TU Munich HF Engineering — Chain matrix theory
 *         Georgia Tech ECE 6350 — ABCD parameters
 */

#ifndef ABCD_PARAMS_H
#define ABCD_PARAMS_H

#include "two_port.h"

/* ============================================================================
 * L1: ABCD-parameter Creation Functions
 * ============================================================================ */

/**
 * @brief Create ABCD-parameter matrix from four elements.
 *
 * @param a Reverse voltage ratio, open output (dimensionless)
 * @param b Transfer impedance, shorted output (Ω)
 * @param c Transfer admittance, open output (S)
 * @param d Reverse current ratio, shorted output (dimensionless)
 * @return ABCD-parameter matrix
 */
matrix2x2_t abcd_create(complex_t a, complex_t b, complex_t c, complex_t d);

/**
 * @brief Create ABCD parameters for a series impedance Z.
 *
 *   ABCD = [[1, Z], [0, 1]]
 *
 * This is one of the three fundamental ABCD building blocks.
 * A series impedance in a transmission path.
 *
 * @param z Series impedance (Ω)
 * @return ABCD matrix for series Z
 *
 * Course: Michigan EECS 411 — Series element ABCD
 */
matrix2x2_t abcd_series_z(complex_t z);

/**
 * @brief Create ABCD parameters for a shunt admittance Y.
 *
 *   ABCD = [[1, 0], [Y, 1]]
 *
 * Second fundamental ABCD building block.
 * A shunt admittance from the line to ground.
 *
 * @param y Shunt admittance (S)
 * @return ABCD matrix for shunt Y
 */
matrix2x2_t abcd_shunt_y(complex_t y);

/**
 * @brief Create ABCD parameters for an ideal transformer.
 *
 *   ABCD = [[1/n, 0], [0, n]]
 *
 * where n = turns ratio (secondary/primary).
 * V1 = (1/n)*V2, I1 = n*(-I2).
 *
 * @param n Turns ratio Ns/Np
 * @return ABCD matrix for ideal transformer
 */
matrix2x2_t abcd_ideal_transformer(double n);

/**
 * @brief Create ABCD parameters for a transmission line section.
 *
 *   ABCD = [[cosh(γl), Z0*sinh(γl)],
 *           [sinh(γl)/Z0, cosh(γl)]]
 *
 * where γ = α + jβ is the complex propagation constant,
 * α = attenuation constant (Np/m), β = phase constant (rad/m),
 * l = line length (m), Z0 = characteristic impedance (Ω).
 *
 * For lossless line (α = 0): cosh(jβl) = cos(βl), sinh(jβl) = j*sin(βl).
 *
 * @param z0 Characteristic impedance (Ω)
 * @param alpha Attenuation constant (Np/m)
 * @param beta Phase constant (rad/m), = 2π/λ
 * @param length Line length (m)
 * @return ABCD matrix for transmission line
 *
 * Course: ETH 227-0455 — Transmission line ABCD
 * Ref: Pozar §2.1 — Transmission line theory
 */
matrix2x2_t abcd_transmission_line(double z0, double alpha,
                                    double beta, double length);

/**
 * @brief Create ABCD parameters for an L-section matching network.
 *
 * L-network (C series + L shunt):
 *
 *        C
 *   P1 —||——┬——— P2
 *           │
 *           L
 *           │
 *          GND
 *
 *   ABCD = [[1, 1/(jωC)], [0, 1]] × [[1, 0], [1/(jωL), 1]]
 *
 * Supports 4 L-network types: series-shunt and shunt-series, with L or C.
 *
 * @param z_series Series element impedance (Ω), e.g., 1/(jωC) for C
 * @param y_shunt Shunt element admittance (S), e.g., 1/(jωL) for L
 * @return ABCD matrix for L-network
 *
 * Course: Stanford EE359 — Impedance matching networks
 * Ref: Bowick, "RF Circuit Design", Ch. 4
 */
matrix2x2_t abcd_l_network(complex_t z_series, complex_t y_shunt);

/**
 * @brief Create ABCD parameters for a π-network.
 *
 *        Yb
 *   P1 —┬———┬——— P2
 *       Ya  Yc
 *       │   │
 *      GND GND
 *
 * ABCD = ABCD_shunt_Y(Ya) × ABCD_series_Z(1/Yb) × ABCD_shunt_Y(Yc)
 *        [[1,0],[Ya,1]] × [[1,1/Yb],[0,1]] × [[1,0],[Yc,1]]
 *
 * Full expansion:
 *   A = 1 + Yc/Yb
 *   B = 1/Yb
 *   C = Ya + Yc + Ya*Yc/Yb
 *   D = 1 + Ya/Yb
 *
 * π-networks are the standard topology for:
 *   - Power amplifier output matching (harmonic suppression)
 *   - Band-pass filters at RF
 *   - Bias tee circuits (DC + RF combining)
 *
 * @param ya First shunt admittance (S)
 * @param yb Series admittance (S)
 * @param yc Second shunt admittance (S)
 * @return ABCD matrix for π-network
 *
 * Course: Georgia Tech ECE 6350 — π-network matching
 */
matrix2x2_t abcd_pi_network(complex_t ya, complex_t yb, complex_t yc);

/**
 * @brief Create ABCD parameters for a T-network.
 *
 * ABCD = ABCD_series_Z(Za) × ABCD_shunt_Y(1/Zc) × ABCD_series_Z(Zb)
 *
 * @param za First series impedance (Ω)
 * @param zb Second series impedance (Ω)
 * @param zc Shunt impedance (Ω)
 * @return ABCD matrix for T-network
 */
matrix2x2_t abcd_t_network(complex_t za, complex_t zb, complex_t zc);

/* ============================================================================
 * L3: ABCD-parameter Analysis Functions
 * ============================================================================ */

/**
 * @brief Compute the voltage transfer function from ABCD parameters.
 *
 * With output terminated by load impedance ZL:
 *   H(s) = V2/V1 = ZL / (A*ZL + B)
 *
 * This is the loaded voltage gain expressed in ABCD parameters.
 * For cascade analysis: H_total = Vn/V1.
 *
 * @param abcd ABCD-parameter matrix
 * @param zl Load impedance (Ω)
 * @return Complex voltage transfer function
 */
complex_t abcd_voltage_transfer(matrix2x2_t abcd, complex_t zl);

/**
 * @brief Compute the input impedance from ABCD parameters.
 *
 * Zin = V1/I1 = (A*ZL + B) / (C*ZL + D)
 *
 * This is the general impedance transformation formula through any two-port.
 * It encompasses all special cases:
 *   - Quarter-wave transformer: Zin = Z0²/ZL
 *   - Half-wave line: Zin = ZL
 *   - Series Z: Zin = Z + ZL
 *
 * @param abcd ABCD-parameter matrix
 * @param zl Load impedance (Ω)
 * @return Input impedance (Ω)
 *
 * Course: ETH 227-0455 — Impedance transformation through ABCD
 */
complex_t abcd_input_impedance(matrix2x2_t abcd, complex_t zl);

/**
 * @brief Compute the image impedance at port 1.
 *
 * Zi1 = sqrt(A*B / (C*D))
 *
 * The image impedance is the impedance that, when terminated at port 2,
 * results in the same impedance seen at port 1. For symmetrical networks
 * (A=D), Zi1 = Zi2 = sqrt(B/C).
 *
 * Image impedance concept is fundamental to filter design (constant-k
 * and m-derived filters), where it determines passband and stopband.
 *
 * @param abcd ABCD-parameter matrix
 * @return Image impedance at port 1 (Ω)
 *
 * Course: Illinois ECE 451 — Image parameter filter design
 * Ref: Zverev, "Handbook of Filter Synthesis", Ch. 2
 */
complex_t abcd_image_impedance_1(matrix2x2_t abcd);

/**
 * @brief Compute the image impedance at port 2.
 *
 * Zi2 = sqrt(D*B / (C*A))
 *
 * @param abcd ABCD-parameter matrix
 * @return Image impedance at port 2 (Ω)
 */
complex_t abcd_image_impedance_2(matrix2x2_t abcd);

/**
 * @brief Compute the image transfer constant θ.
 *
 * θ = acosh(sqrt(A*D)) = ln(sqrt(A*D) + sqrt(A*D - 1))
 *
 * The image transfer constant has:
 *   Re(θ) = image attenuation constant α (Np) — represents loss
 *   Im(θ) = image phase constant β (rad) — represents phase shift per section
 *
 * For a filter section: passband if θ is pure imaginary, stopband if pure real.
 *
 * @param abcd ABCD-parameter matrix
 * @return Image transfer constant θ (complex)
 */
complex_t abcd_image_transfer_constant(matrix2x2_t abcd);

/**
 * @brief Compute the propagation constant of a cascade.
 *
 * For n identical sections with ABCD each: γ = θ (of one section).
 * For a general cascade: cosh(θ) = (A + D) / 2.
 *
 * This formula derives from the eigenvalue analysis of the ABCD matrix
 * and is central to filter theory.
 *
 * @param abcd ABCD-parameter matrix
 * @return Propagation constant (complex)
 */
complex_t abcd_propagation_constant(matrix2x2_t abcd);

/**
 * @brief Check if an ABCD matrix represents a reciprocal network.
 *
 * For reciprocal networks: det(ABCD) = AD - BC = 1.
 * This is the fundamental reciprocal condition in ABCD form.
 *
 * @param abcd ABCD-parameter matrix
 * @param tolerance Numerical tolerance
 * @return 1 if reciprocal, 0 otherwise
 */
int abcd_is_reciprocal(matrix2x2_t abcd, double tolerance);

/**
 * @brief Compute the insertion loss of the two-port in a Z0 system.
 *
 * IL(dB) = -20*log10(|S21|) but computed from ABCD.
 * S21 = 2 / (A + B/Z0 + C*Z0 + D)
 *
 * @param abcd ABCD-parameter matrix
 * @param z0 Reference impedance (Ω), typically 50
 * @return Insertion loss in dB
 */
double abcd_insertion_loss_db(matrix2x2_t abcd, double z0);

/* ============================================================================
 * L4: Cascade Analysis
 * ============================================================================ */

/**
 * @brief Cascade multiple ABCD-parameter networks.
 *
 * ABCD_total = ABCD[0] × ABCD[1] × ... × ABCD[n-1]
 *
 * This is the unique property of ABCD parameters — matrix multiplication
 * corresponds to physical cascade connection.
 *
 * @param stages Array of n ABCD matrices
 * @param n Number of stages
 * @return Combined ABCD matrix
 *
 * Course: TU Munich HF Engineering — Kettenschaltung (cascade)
 */
matrix2x2_t abcd_cascade_n(const matrix2x2_t *stages, int n);

/**
 * @brief De-embed a two-port network from a cascaded measurement.
 *
 * Given ABCD_total = ABCD_left × ABCD_dut × ABCD_right,
 * extract ABCD_dut = ABCD_left⁻¹ × ABCD_total × ABCD_right⁻¹.
 *
 * This is the mathematical basis of VNA calibration and de-embedding,
 * removing the effects of cables, connectors, and fixtures from measurements.
 *
 * @param abcd_total Total measured ABCD
 * @param abcd_left Left fixture ABCD
 * @param abcd_right Right fixture ABCD
 * @return De-embedded DUT ABCD
 *
 * Course: ETH 227-0455 — TRL calibration and de-embedding
 */
matrix2x2_t abcd_deembed(matrix2x2_t abcd_total,
                          matrix2x2_t abcd_left,
                          matrix2x2_t abcd_right);

#endif /* ABCD_PARAMS_H */
