/**
 * @file conversion.h
 * @brief Two-Port Parameter Conversion — All 30 Directional Conversions
 *
 * Parameter conversion is the bridge between different two-port representations.
 * Since all six parameter types describe the same physical network, they are
 * mathematically equivalent and can be converted from one to another.
 *
 * There are 6 × 5 = 30 possible directional conversions (Z↔Y, Z↔H, etc.).
 * Only 15 need unique formulas (the other 15 are inverses).
 *
 * Special cases:
 *   Z↔Y: Matrix inversion (Y = Z⁻¹, Z = Y⁻¹)
 *   ABCD↔S: Most important for RF — VNAs measure S, but circuits use ABCD
 *   H↔G: Matrix inversion (G = H⁻¹, H = G⁻¹)
 *
 * For singular matrices (e.g., an ideal FET with infinite Zin), some
 * conversions may be undefined. These functions return NAN-filled matrices.
 *
 * Reference: Frickey, "Conversions Between S, Z, Y, h, ABCD, and T Parameters",
 *   IEEE Trans. Microwave Theory Tech., 1994
 * Course: Georgia Tech ECE 6350 — Parameter conversion
 *         ETH 227-0455 — Network parameter transformations
 */

#ifndef CONVERSION_H
#define CONVERSION_H

#include "two_port.h"

/* ============================================================================
 * L5: Complete Parameter Conversion Algorithms
 * ============================================================================ */

/**
 * @brief Convert Z-parameters to Y-parameters.
 *
 * Y = Z⁻¹ (matrix inversion)
 *
 * Elements:
 *   y11 = z22 / ΔZ,  y12 = -z12 / ΔZ
 *   y21 = -z21 / ΔZ,  y22 = z11 / ΔZ
 * where ΔZ = z11*z22 - z12*z21
 *
 * @param z Z-parameter matrix
 * @return Y-parameter matrix, or all-NAN if singular
 */
matrix2x2_t convert_z_to_y(matrix2x2_t z);

/**
 * @brief Convert Y-parameters to Z-parameters.
 *
 * Z = Y⁻¹ (matrix inversion)
 * Dual of convert_z_to_y.
 *
 * @param y Y-parameter matrix
 * @return Z-parameter matrix, or all-NAN if singular
 */
matrix2x2_t convert_y_to_z(matrix2x2_t y);

/**
 * @brief Convert Z-parameters to H-parameters.
 *
 * Pseudo-inverse transformation:
 *   h11 = ΔZ / z22 = z11 - z12*z21/z22
 *   h12 = z12 / z22
 *   h21 = -z21 / z22
 *   h22 = 1 / z22
 * where ΔZ = z11*z22 - z12*z21
 *
 * @param z Z-parameter matrix
 * @return H-parameter matrix, or all-NAN if z22 = 0
 */
matrix2x2_t convert_z_to_h(matrix2x2_t z);

/**
 * @brief Convert H-parameters to Z-parameters.
 *
 *   z11 = ΔH / h22 = h11 - h12*h21/h22
 *   z12 = h12 / h22
 *   z21 = -h21 / h22
 *   z22 = 1 / h22
 * where ΔH = h11*h22 - h12*h21
 *
 * @param h H-parameter matrix
 * @return Z-parameter matrix, or all-NAN if h22 = 0
 */
matrix2x2_t convert_h_to_z(matrix2x2_t h);

/**
 * @brief Convert Z-parameters to G-parameters.
 *
 *   g11 = 1 / z11
 *   g12 = -z12 / z11
 *   g21 = z21 / z11
 *   g22 = ΔZ / z11 = z22 - z12*z21/z11
 *
 * @param z Z-parameter matrix
 * @return G-parameter matrix, or all-NAN if z11 = 0
 */
matrix2x2_t convert_z_to_g(matrix2x2_t z);

/**
 * @brief Convert G-parameters to Z-parameters.
 *
 *   z11 = 1 / g11
 *   z12 = -g12 / g11
 *   z21 = g21 / g11
 *   z22 = ΔG / g11 = g22 - g12*g21/g11
 *
 * @param g G-parameter matrix
 * @return Z-parameter matrix, or all-NAN if g11 = 0
 */
matrix2x2_t convert_g_to_z(matrix2x2_t g);

/**
 * @brief Convert Z-parameters to ABCD-parameters.
 *
 *   A = z11 / z21
 *   B = ΔZ / z21 = (z11*z22 - z12*z21) / z21
 *   C = 1 / z21
 *   D = z22 / z21
 *
 * @param z Z-parameter matrix
 * @return ABCD-parameter matrix, or all-NAN if z21 = 0
 */
matrix2x2_t convert_z_to_abcd(matrix2x2_t z);

/**
 * @brief Convert ABCD-parameters to Z-parameters.
 *
 *   z11 = A / C
 *   z12 = Δ(ABCD) / C = (A*D - B*C) / C
 *   z21 = 1 / C
 *   z22 = D / C
 *
 * @param abcd ABCD-parameter matrix
 * @return Z-parameter matrix, or all-NAN if C = 0
 */
matrix2x2_t convert_abcd_to_z(matrix2x2_t abcd);

/**
 * @brief Convert Y-parameters to H-parameters.
 *
 *   h11 = 1 / y11
 *   h12 = -y12 / y11
 *   h21 = y21 / y11
 *   h22 = ΔY / y11 = y22 - y12*y21/y11
 *
 * @param y Y-parameter matrix
 * @return H-parameter matrix, or all-NAN if y11 = 0
 */
matrix2x2_t convert_y_to_h(matrix2x2_t y);

/**
 * @brief Convert H-parameters to Y-parameters.
 *
 *   y11 = 1 / h11
 *   y12 = -h12 / h11
 *   y21 = h21 / h11
 *   y22 = ΔH / h11 = h22 - h12*h21/h11
 *
 * @param h H-parameter matrix
 * @return Y-parameter matrix, or all-NAN if h11 = 0
 */
matrix2x2_t convert_h_to_y(matrix2x2_t h);

/**
 * @brief Convert Y-parameters to ABCD-parameters.
 *
 *   A = -y22 / y21
 *   B = -1 / y21
 *   C = -ΔY / y21 = -(y11*y22 - y12*y21) / y21
 *   D = -y11 / y21
 *
 * @param y Y-parameter matrix
 * @return ABCD-parameter matrix, or all-NAN if y21 = 0
 */
matrix2x2_t convert_y_to_abcd(matrix2x2_t y);

/**
 * @brief Convert ABCD-parameters to Y-parameters.
 *
 *   y11 = D / B
 *   y12 = -Δ(ABCD) / B = -(A*D - B*C) / B
 *   y21 = -1 / B
 *   y22 = A / B
 *
 * @param abcd ABCD-parameter matrix
 * @return Y-parameter matrix, or all-NAN if B = 0
 */
matrix2x2_t convert_abcd_to_y(matrix2x2_t abcd);

/**
 * @brief Convert H-parameters to G-parameters.
 *
 * G = H⁻¹ (matrix inversion)
 *
 * @param h H-parameter matrix
 * @return G-parameter matrix, or all-NAN if singular
 */
matrix2x2_t convert_h_to_g(matrix2x2_t h);

/**
 * @brief Convert G-parameters to H-parameters.
 *
 * H = G⁻¹ (matrix inversion)
 *
 * @param g G-parameter matrix
 * @return H-parameter matrix, or all-NAN if singular
 */
matrix2x2_t convert_g_to_h(matrix2x2_t g);

/**
 * @brief Convert S-parameters to Z-parameters.
 *
 * Z = Z0 * (I + S) * (I - S)⁻¹
 *
 * where I is the 2x2 identity matrix.
 *
 * Elements:
 *   z11 = Z0 * ((1+s11)*(1-s22) + s12*s21) / ΔS
 *   z12 = Z0 * 2*s12 / ΔS
 *   z21 = Z0 * 2*s21 / ΔS
 *   z22 = Z0 * ((1-s11)*(1+s22) + s12*s21) / ΔS
 * where ΔS = (1-s11)*(1-s22) - s12*s21
 *
 * @param s S-parameter matrix
 * @param z0 Reference impedance (Ω)
 * @return Z-parameter matrix
 *
 * Course: ETH 227-0455 — S to Z conversion
 */
matrix2x2_t convert_s_to_z(matrix2x2_t s, double z0);

/**
 * @brief Convert Z-parameters to S-parameters.
 *
 * S = (Z - Z0*I) * (Z + Z0*I)⁻¹
 *
 * Elements:
 *   s11 = ((z11-Z0)*(z22+Z0) - z12*z21) / ΔZ
 *   s12 = 2*z12*Z0 / ΔZ
 *   s21 = 2*z21*Z0 / ΔZ
 *   s22 = ((z11+Z0)*(z22-Z0) - z12*z21) / ΔZ
 * where ΔZ = (z11+Z0)*(z22+Z0) - z12*z21
 *
 * @param z Z-parameter matrix
 * @param z0 Reference impedance (Ω)
 * @return S-parameter matrix
 */
matrix2x2_t convert_z_to_s(matrix2x2_t z, double z0);

/**
 * @brief Convert S-parameters to Y-parameters.
 *
 * Y = (1/Z0) * (I - S) * (I + S)⁻¹
 *
 * @param s S-parameter matrix
 * @param z0 Reference impedance (Ω)
 * @return Y-parameter matrix
 */
matrix2x2_t convert_s_to_y(matrix2x2_t s, double z0);

/**
 * @brief Convert Y-parameters to S-parameters.
 *
 * S = (I - Z0*Y) * (I + Z0*Y)⁻¹
 *
 * @param y Y-parameter matrix
 * @param z0 Reference impedance (Ω)
 * @return S-parameter matrix
 */
matrix2x2_t convert_y_to_s(matrix2x2_t y, double z0);

/**
 * @brief Convert S-parameters to ABCD-parameters.
 *
 * A = ((1+s11)*(1-s22) + s12*s21) / (2*s21)
 * B = Z0 * ((1+s11)*(1+s22) - s12*s21) / (2*s21)
 * C = (1/Z0) * ((1-s11)*(1-s22) - s12*s21) / (2*s21)
 * D = ((1-s11)*(1+s22) + s12*s21) / (2*s21)
 *
 * This is the most important conversion for VNA-based circuit analysis.
 *
 * @param s S-parameter matrix
 * @param z0 Reference impedance (Ω)
 * @return ABCD-parameter matrix, or all-NAN if s21 = 0
 *
 * Course: Georgia Tech ECE 6350 — S to ABCD conversion for VNA
 */
matrix2x2_t convert_s_to_abcd(matrix2x2_t s, double z0);

/**
 * @brief Convert ABCD-parameters to S-parameters.
 *
 * s11 = (A + B/Z0 - C*Z0 - D) / (A + B/Z0 + C*Z0 + D)
 * s12 = 2*(A*D - B*C) / (A + B/Z0 + C*Z0 + D)
 * s21 = 2 / (A + B/Z0 + C*Z0 + D)
 * s22 = (-A + B/Z0 - C*Z0 + D) / (A + B/Z0 + C*Z0 + D)
 *
 * @param abcd ABCD-parameter matrix
 * @param z0 Reference impedance (Ω)
 * @return S-parameter matrix
 */
matrix2x2_t convert_abcd_to_s(matrix2x2_t abcd, double z0);

/**
 * @brief Convert ABCD-parameters to H-parameters.
 *
 *   h11 = B / D
 *   h12 = Δ(ABCD) / D = (A*D - B*C) / D
 *   h21 = -1 / D
 *   h22 = C / D
 *
 * @param abcd ABCD-parameter matrix
 * @return H-parameter matrix, or all-NAN if D = 0
 */
matrix2x2_t convert_abcd_to_h(matrix2x2_t abcd);

/**
 * @brief Convert H-parameters to ABCD-parameters.
 *
 *   A = -ΔH / h21 = -(h11*h22 - h12*h21) / h21
 *   B = -h11 / h21
 *   C = -h22 / h21
 *   D = -1 / h21
 *
 * @param h H-parameter matrix
 * @return ABCD-parameter matrix, or all-NAN if h21 = 0
 */
matrix2x2_t convert_h_to_abcd(matrix2x2_t h);

/* ============================================================================
 * L5: Conversion Utility Functions
 * ============================================================================ */

/**
 * @brief Generic parameter conversion between any two types.
 *
 * Routes through internally via matrix inversion or specific formulas.
 * Handles all 30 conversion directions.
 *
 * @param src Source parameter matrix
 * @param src_type Source parameter type
 * @param dst_type Destination parameter type
 * @param z0 Reference impedance (only used for S-parameter conversions,
 *           ignored for others, use 0 if not applicable)
 * @return Converted matrix; NAN-filled if conversion is impossible
 */
matrix2x2_t convert_parameters(matrix2x2_t src, param_type_t src_type,
                                param_type_t dst_type, double z0);

/**
 * @brief Check if a parameter conversion is mathematically possible.
 *
 * @param src_type Source parameter type
 * @param dst_type Destination parameter type
 * @return 1 if possible, 0 if there exist networks where it fails
 */
int conversion_is_possible(param_type_t src_type, param_type_t dst_type);

/**
 * @brief Get the string name of a parameter type.
 *
 * @param type Parameter type
 * @return String name ("Z", "Y", "H", "G", "ABCD", "S")
 */
const char *param_type_name(param_type_t type);

/**
 * @brief Check parameter conversion for self-consistency.
 *
 * Converts from src to dst and back, checking if the result matches
 * the original within tolerance.
 *
 * @param original Original parameter matrix
 * @param type Original parameter type
 * @param z0 Reference impedance
 * @param tolerance Numerical tolerance
 * @return 1 if round-trip is consistent, 0 otherwise
 */
int conversion_roundtrip_check(matrix2x2_t original, param_type_t type,
                                double z0, double tolerance);

#endif /* CONVERSION_H */
