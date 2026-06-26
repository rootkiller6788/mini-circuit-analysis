/**
 * @file interconnection.h
 * @brief Two-Port Network Interconnection Analysis
 *
 * When two-port networks are connected together, the combined network's
 * parameters depend on the interconnection topology. Each topology has a
 * natural parameter type for simple combination:
 *
 *   Series-Series:      Z_total = Z1 + Z2        (additive Z)
 *   Parallel-Parallel:  Y_total = Y1 + Y2        (additive Y)
 *   Series-Parallel:    H_total = H1 + H2        (additive H)
 *   Parallel-Series:    G_total = G1 + G2        (additive G)
 *   Cascade:            ABCD_total = ABCD1 × ABCD2  (multiplicative ABCD)
 *
 * These rules are valid only when the interconnection does not violate
 * the port conditions (Brune's test — each port still has equal currents
 * flowing in and out).
 *
 * Reference: Kuo, "Network Analysis and Synthesis" (2006), Ch. 14
 * Course: Illinois ECE 451 — Network interconnections
 *         TU Munich HF Engineering — Zusammenschaltung von Vierpolen
 */

#ifndef INTERCONNECTION_H
#define INTERCONNECTION_H

#include "two_port.h"

/* ============================================================================
 * L5: Interconnection Functions for All Five Configurations
 * ============================================================================ */

/**
 * @brief Connect two two-port networks in series-series configuration.
 *
 * Input ports are in series, output ports are in series.
 * Natural representation: Z-parameters are additive.
 *
 *   Z_total = Z1 + Z2  (matrix addition)
 *
 * Validity condition (Brune's test): The port currents must remain
 * equal and opposite at each port after interconnection. This is
 * satisfied if the two networks share a common ground or if the
 * series connection does not create current loops.
 *
 * Applications:
 *   - Adding a feedback network to an amplifier
 *   - Cascading filter sections (though cascade is usually better)
 *   - Modeling parasitic elements in series with device
 *
 * @param z1 First network Z-parameters
 * @param z2 Second network Z-parameters
 * @return Combined Z-parameter matrix
 */
matrix2x2_t interconnect_series_series(matrix2x2_t z1, matrix2x2_t z2);

/**
 * @brief Connect two two-port networks in parallel-parallel configuration.
 *
 * Input ports are in parallel, output ports are in parallel.
 * Natural representation: Y-parameters are additive.
 *
 *   Y_total = Y1 + Y2  (admittance matrix addition)
 *
 * Validity condition: Port voltages must be identical for both networks,
 * which is automatically satisfied for parallel connection.
 *
 * Applications:
 *   - Connecting components in parallel (e.g., multiple capacitors)
 *   - Admittance-based circuit analysis
 *   - Modeling shunt parasitics
 *
 * @param y1 First network Y-parameters
 * @param y2 Second network Y-parameters
 * @return Combined Y-parameter matrix
 */
matrix2x2_t interconnect_parallel_parallel(matrix2x2_t y1, matrix2x2_t y2);

/**
 * @brief Connect two two-port networks in series-parallel configuration.
 *
 * Input ports in series, output ports in parallel.
 * Natural representation: H-parameters are additive.
 *
 *   H_total = H1 + H2  (hybrid matrix addition)
 *
 * This configuration is common in feedback amplifier analysis where
 * feedback is applied in series at the input and shunt at the output
 * (series-shunt feedback → voltage-series feedback).
 *
 * Applications:
 *   - BJT amplifier with emitter degeneration
 *   - Series-shunt feedback amplifiers
 *   - Modeling the effect of source/load on H-parameters
 *
 * @param h1 First network H-parameters
 * @param h2 Second network H-parameters
 * @return Combined H-parameter matrix
 *
 * Course: Berkeley EE105 — Feedback amplifier topologies
 */
matrix2x2_t interconnect_series_parallel(matrix2x2_t h1, matrix2x2_t h2);

/**
 * @brief Connect two two-port networks in parallel-series configuration.
 *
 * Input ports in parallel, output ports in series.
 * Natural representation: G-parameters are additive.
 *
 *   G_total = G1 + G2
 *
 * This configuration is used in shunt-series feedback amplifiers
 * (current-shunt feedback).
 *
 * @param g1 First network G-parameters
 * @param g2 Second network G-parameters
 * @return Combined G-parameter matrix
 */
matrix2x2_t interconnect_parallel_series(matrix2x2_t g1, matrix2x2_t g2);

/**
 * @brief Cascade two two-port networks.
 *
 * Output of network 1 connects to input of network 2.
 * Natural representation: ABCD-parameters are multiplicative.
 *
 *   ABCD_total = ABCD1 × ABCD2  (matrix multiplication, order matters!)
 *
 * This is the most common interconnection in practice:
 *   - Filter stages in cascade
 *   - Amplifier stage + matching network + transmission line
 *   - Complete RF chain analysis (PA → filter → antenna)
 *
 * IMPORTANT: The order is ABCD1 × ABCD2 (first network, then second).
 * Matrix multiplication is non-commutative — order matters!
 *
 * @param abcd1 First network ABCD-parameters (closest to source)
 * @param abcd2 Second network ABCD-parameters (closest to load)
 * @return Combined ABCD-parameter matrix
 *
 * Course: TU Munich HF Engineering — Kettenschaltung
 */
matrix2x2_t interconnect_cascade(matrix2x2_t abcd1, matrix2x2_t abcd2);

/* ============================================================================
 * L5: General Interconnection with Automatic Conversion
 * ============================================================================ */

/**
 * @brief Generic interconnection of two networks with automatic conversion.
 *
 * Given two networks in arbitrary parameter types, this function:
 * 1. Detects the natural parameter type for the configuration
 * 2. Converts both networks to that type
 * 3. Performs the interconnection (addition or multiplication)
 * 4. Converts result back to the requested output type
 *
 * @param net1 First network parameter matrix
 * @param type1 First network parameter type
 * @param net2 Second network parameter matrix
 * @param type2 Second network parameter type
 * @param config Interconnection configuration type
 * @param out_type Desired output parameter type
 * @param z0 Reference impedance for S-parameter conversions
 * @return Combined network in requested parameter type
 */
matrix2x2_t interconnect_general(matrix2x2_t net1, param_type_t type1,
                                  matrix2x2_t net2, param_type_t type2,
                                  config_type_t config,
                                  param_type_t out_type, double z0);

/**
 * @brief Check if a given interconnection configuration is valid
 * (Brune's test approximation).
 *
 * In practice, this checks if the port conditions would be violated.
 * Returns 1 if the interconnection is likely valid, 0 if likely invalid.
 *
 * The Brune test checks: after interconnection, do the port currents
 * still obey I1_in = -I1_out and I2_in = -I2_out at each port?
 *
 * @param net1 First network (any parameter type)
 * @param type1 Type of first network
 * @param net2 Second network (any parameter type)
 * @param type2 Type of second network
 * @param config Interconnection configuration
 * @return 1 if port conditions are maintained, 0 if violated
 */
int interconnect_is_valid(matrix2x2_t net1, param_type_t type1,
                           matrix2x2_t net2, param_type_t type2,
                           config_type_t config);

/* ============================================================================
 * L4: Bartlett's Bisection Theorem (Symmetry-based simplification)
 * ============================================================================ */

/**
 * @brief Apply Bartlett's bisection theorem to a symmetric network.
 *
 * A symmetric two-port network can be analyzed by splitting it into two
 * identical halves and attaching either an open-circuit (even-mode) or
 * short-circuit (odd-mode) at the symmetry plane.
 *
 * Zoc = open-circuit impedance of each half
 * Zsc = short-circuit impedance of each half
 *
 * Then: z11 = z22 = (Zoc + Zsc)/2
 *       z12 = z21 = (Zoc - Zsc)/2
 *
 * This yields the complete Z-parameters from two simpler one-port
 * measurements. The theorem is widely used in filter design and
 * differential circuit analysis.
 *
 * @param zoc Open-circuit impedance of one half (Ω)
 * @param zsc Short-circuit impedance of one half (Ω)
 * @return Complete Z-parameter matrix for symmetric two-port
 *
 * Course: Illinois ECE 451 — Bartlett's theorem
 * Ref: Bartlett, "An Extension of a Property of Artificial Lines", 1927
 */
matrix2x2_t bartlett_bisection(complex_t zoc, complex_t zsc);

/**
 * @brief Even/odd mode analysis of a symmetric network.
 *
 * For a symmetric two-port:
 *   Even mode: V1 = V2 (in-phase excitation)
 *     Zeven = z11 + z12 = z22 + z21
 *   Odd mode: V1 = -V2 (anti-phase excitation)
 *     Zodd = z11 - z12 = z22 - z21
 *
 * These are the eigenvalues of the Z-matrix and correspond to
 * differential and common-mode analysis in differential circuits.
 *
 * @param z Symmetric Z-parameter matrix
 * @param zeven Output: even-mode impedance (common mode)
 * @param zodd Output: odd-mode impedance (differential mode)
 *
 * Course: Stanford EE214 — Differential amplifier analysis
 */
void symmetric_even_odd_impedance(matrix2x2_t z, complex_t *zeven,
                                   complex_t *zodd);

#endif /* INTERCONNECTION_H */
