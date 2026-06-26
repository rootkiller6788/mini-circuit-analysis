/**
 * @file two_port.c
 * @brief Two-Port Network Parameter Conversions and Analysis
 *
 * Implements conversions between all six two-port parameter representations:
 * Z, Y, H, G, ABCD (transmission), and S (scattering) parameters.
 *
 * Each conversion function implements the mathematical relationship
 * between parameter sets, handling singularities (e.g., det=0 cases).
 *
 * Knowledge Coverage:
 *   L2 - Core Concepts: Two-port network representations
 *   L3 - Mathematical Structures: Parameter matrix transformations
 *   L4 - Fundamental Laws: Reciprocity theorem verification
 *   L6 - Canonical Problems: Cascaded network analysis
 *
 * Reference:
 *   Pozar "Microwave Engineering" (2012) Ch 4 — S-parameters
 *   Hayt et al. "Engineering Circuit Analysis" (2019) Ch 17
 *   Carson "High-Frequency Amplifiers" (1975) — Two-port theory
 *
 * Course Mapping:
 *   MIT 6.002: Two-port networks
 *   Berkeley EE105: Small-signal two-port models
 *   Stanford EE359: S-parameters in wireless systems
 */

#include "network_theorem.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L2: Z → Y Conversion
 * ==========================================================================
 *
 * Y = Z⁻¹, i.e.:
 *   [y11 y12] = 1/det(Z) * [ z22  -z12]
 *   [y21 y22]               [-z21   z11]
 *
 * where det(Z) = z11*z22 - z12*z21
 *
 * Returns -1 if det(Z) = 0 (singular, no Y-parameter representation).
 */
int z_to_y(const ZParameters *z, YParameters *y) {
    if (!z || !y) return -1;

    /* Compute determinant of Z matrix */
    double a = z->z11.real, b = z->z11.imag;
    double c = z->z12.real, d = z->z12.imag;
    double e = z->z21.real, f = z->z21.imag;
    double g = z->z22.real, h = z->z22.imag;

    /* det = z11*z22 - z12*z21 (complex multiplication) */
    double det_real = (a*g - b*h) - (c*e - d*f);
    double det_imag = (a*h + b*g) - (c*f + d*e);

    double det_mag = det_real*det_real + det_imag*det_imag;
    if (det_mag < 1e-30) return -1;

    /* y11 =  z22 / det */
    double denom = det_mag;
    y->y11.real = (g*det_real + h*det_imag) / denom;
    y->y11.imag = (h*det_real - g*det_imag) / denom;

    /* y12 = -z12 / det */
    y->y12.real = (-c*det_real - d*det_imag) / denom;
    y->y12.imag = (-d*det_real + c*det_imag) / denom;

    /* y21 = -z21 / det */
    y->y21.real = (-e*det_real - f*det_imag) / denom;
    y->y21.imag = (-f*det_real + e*det_imag) / denom;

    /* y22 =  z11 / det */
    y->y22.real = (a*det_real + b*det_imag) / denom;
    y->y22.imag = (b*det_real - a*det_imag) / denom;

    return 0;
}

/* ==========================================================================
 * L2: Y → Z Conversion
 * ==========================================================================
 *
 * Z = Y⁻¹, symmetric to Z→Y conversion.
 */
int y_to_z(const YParameters *y, ZParameters *z) {
    if (!y || !z) return -1;

    /* Treat Y as "Z" and compute its inverse to get Z */
    ZParameters z_temp;
    z_temp.z11.real = y->y11.real; z_temp.z11.imag = y->y11.imag;
    z_temp.z12.real = y->y12.real; z_temp.z12.imag = y->y12.imag;
    z_temp.z21.real = y->y21.real; z_temp.z21.imag = y->y21.imag;
    z_temp.z22.real = y->y22.real; z_temp.z22.imag = y->y22.imag;

    YParameters y_temp;
    if (z_to_y(&z_temp, &y_temp) != 0) return -1;

    z->z11.real = y_temp.y11.real; z->z11.imag = y_temp.y11.imag;
    z->z12.real = y_temp.y12.real; z->z12.imag = y_temp.y12.imag;
    z->z21.real = y_temp.y21.real; z->z21.imag = y_temp.y21.imag;
    z->z22.real = y_temp.y22.real; z->z22.imag = y_temp.y22.imag;
    return 0;
}

/* ==========================================================================
 * L2: Z → H Conversion
 * ==========================================================================
 *
 * H-parameters from Z-parameters:
 *   h11 = det(Z) / z22
 *   h12 = z12 / z22
 *   h21 = -z21 / z22
 *   h22 = 1 / z22
 *
 * Requires z22 ≠ 0.
 */
int z_to_h(const ZParameters *z, HParameters *h) {
    if (!z || !h) return -1;

    double a = z->z22.real, b = z->z22.imag;
    double denom = a*a + b*b;
    if (denom < 1e-30) return -1;

    /* det(Z) = z11*z22 - z12*z21 */
    double det_real = (z->z11.real*z->z22.real - z->z11.imag*z->z22.imag)
                    - (z->z12.real*z->z21.real - z->z12.imag*z->z21.imag);
    double det_imag = (z->z11.real*z->z22.imag + z->z11.imag*z->z22.real)
                    - (z->z12.real*z->z21.imag + z->z12.imag*z->z21.real);

    /* h11 = det(Z) / z22 */
    h->h11.real = (det_real*a + det_imag*b) / denom;
    h->h11.imag = (det_imag*a - det_real*b) / denom;

    /* h12 = z12 / z22 */
    h->h12_real = (z->z12.real*a + z->z12.imag*b) / denom;
    h->h12_imag = (z->z12.imag*a - z->z12.real*b) / denom;

    /* h21 = -z21 / z22 */
    h->h21_real = -(z->z21.real*a + z->z21.imag*b) / denom;
    h->h21_imag = -(z->z21.imag*a - z->z21.real*b) / denom;

    /* h22 = 1 / z22 */
    h->h22.real =  a / denom;
    h->h22.imag = -b / denom;

    return 0;
}

/* ==========================================================================
 * L2: H → Z Conversion
 * ==========================================================================
 *
 * Z-parameters from H-parameters:
 *   z11 = det(H) / h22
 *   z12 = h12 / h22
 *   z21 = -h21 / h22
 *   z22 = 1 / h22
 *
 * where det(H) = h11*h22 - h12*h21
 * Requires h22 ≠ 0.
 */
int h_to_z(const HParameters *h, ZParameters *z) {
    if (!h || !z) return -1;

    double a = h->h22.real, b = h->h22.imag;
    double denom = a*a + b*b;
    if (denom < 1e-30) return -1;

    /* det(H) = h11*h22 - h12*h21 */
    /* h11*h22 */
    double p1_real = h->h11.real*a - h->h11.imag*b;
    double p1_imag = h->h11.real*b + h->h11.imag*a;
    /* h12*h21 */
    double p2_real = h->h12_real*h->h21_real - h->h12_imag*h->h21_imag;
    double p2_imag = h->h12_real*h->h21_imag + h->h12_imag*h->h21_real;
    double det_real = p1_real - p2_real;
    double det_imag = p1_imag - p2_imag;

    /* z11 = det(H) / h22 */
    z->z11.real = (det_real*a + det_imag*b) / denom;
    z->z11.imag = (det_imag*a - det_real*b) / denom;

    /* z12 = h12 / h22 */
    z->z12.real = (h->h12_real*a + h->h12_imag*b) / denom;
    z->z12.imag = (h->h12_imag*a - h->h12_real*b) / denom;

    /* z21 = -h21 / h22 */
    z->z21.real = -(h->h21_real*a + h->h21_imag*b) / denom;
    z->z21.imag = -(h->h21_imag*a - h->h21_real*b) / denom;

    /* z22 = 1 / h22 */
    z->z22.real =  a / denom;
    z->z22.imag = -b / denom;

    return 0;
}

/* ==========================================================================
 * L2: Z → ABCD (Transmission) Conversion
 * ==========================================================================
 *
 * ABCD from Z-parameters:
 *   A = z11 / z21
 *   B = det(Z) / z21
 *   C = 1 / z21
 *   D = z22 / z21
 *
 * Requires z21 ≠ 0 (otherwise network cannot be cascaded in this direction).
 */
int z_to_abcd(const ZParameters *z, ABCDParameters *abcd) {
    if (!z || !abcd) return -1;

    double a = z->z21.real, b = z->z21.imag;
    double denom = a*a + b*b;
    if (denom < 1e-30) return -1;

    /* A = z11 / z21 */
    abcd->A_real = (z->z11.real*a + z->z11.imag*b) / denom;
    abcd->A_imag = (z->z11.imag*a - z->z11.real*b) / denom;

    /* det(Z) = z11*z22 - z12*z21 */
    double det_real = (z->z11.real*z->z22.real - z->z11.imag*z->z22.imag)
                    - (z->z12.real*z->z21.real - z->z12.imag*z->z21.imag);
    double det_imag = (z->z11.real*z->z22.imag + z->z11.imag*z->z22.real)
                    - (z->z12.real*z->z21.imag + z->z12.imag*z->z21.real);

    /* B = det(Z) / z21 */
    abcd->B_real = (det_real*a + det_imag*b) / denom;
    abcd->B_imag = (det_imag*a - det_real*b) / denom;

    /* C = 1 / z21 */
    abcd->C_real =  a / denom;
    abcd->C_imag = -b / denom;

    /* D = z22 / z21 */
    abcd->D_real = (z->z22.real*a + z->z22.imag*b) / denom;
    abcd->D_imag = (z->z22.imag*a - z->z22.real*b) / denom;

    return 0;
}

/* ==========================================================================
 * L2: Cascade Two ABCD Networks
 * ==========================================================================
 *
 * For cascaded two-ports (output of first → input of second):
 *   [ABCD_total] = [ABCD_1] * [ABCD_2]
 *
 * This is the key advantage of ABCD parameters — cascading is simple
 * matrix multiplication.
 *
 *              ┌─────┐   ┌─────┐
 *     V1,I1 ──→│ABCD1│──→│ABCD2│──→ V2',I2'
 *              └─────┘   └─────┘
 */
void abcd_cascade(const ABCDParameters *a, const ABCDParameters *b,
                  ABCDParameters *result) {
    if (!a || !b || !result) return;

    /* Matrix multiplication: C = A * B
     * [A1 B1] * [A2 B2] = [A1*A2+B1*C2   A1*B2+B1*D2]
     * [C1 D1]   [C2 D2]   [C1*A2+D1*C2   C1*B2+D1*D2]
     */
    double A1r = a->A_real, A1i = a->A_imag;
    double B1r = a->B_real, B1i = a->B_imag;
    double C1r = a->C_real, C1i = a->C_imag;
    double D1r = a->D_real, D1i = a->D_imag;

    double A2r = b->A_real, A2i = b->A_imag;
    double B2r = b->B_real, B2i = b->B_imag;
    double C2r = b->C_real, C2i = b->C_imag;
    double D2r = b->D_real, D2i = b->D_imag;

    /* Row 1, Col 1: A1*A2 + B1*C2 */
    result->A_real = (A1r*A2r - A1i*A2i) + (B1r*C2r - B1i*C2i);
    result->A_imag = (A1r*A2i + A1i*A2r) + (B1r*C2i + B1i*C2r);

    /* Row 1, Col 2: A1*B2 + B1*D2 */
    result->B_real = (A1r*B2r - A1i*B2i) + (B1r*D2r - B1i*D2i);
    result->B_imag = (A1r*B2i + A1i*B2r) + (B1r*D2i + B1i*D2r);

    /* Row 2, Col 1: C1*A2 + D1*C2 */
    result->C_real = (C1r*A2r - C1i*A2i) + (D1r*C2r - D1i*C2i);
    result->C_imag = (C1r*A2i + C1i*A2r) + (D1r*C2i + D1i*C2r);

    /* Row 2, Col 2: C1*B2 + D1*D2 */
    result->D_real = (C1r*B2r - C1i*B2i) + (D1r*D2r - D1i*D2i);
    result->D_imag = (C1r*B2i + C1i*B2r) + (D1r*D2i + D1i*D2r);
}

/* ==========================================================================
 * L2: Z → S (Scattering) Parameter Conversion
 * ==========================================================================
 *
 * S-parameters normalized to characteristic impedance Z0:
 *   S = (Z - Z0*I) * (Z + Z0*I)⁻¹
 *
 * For each element:
 *   s11 = ((z11-Z0)*(z22+Z0) - z12*z21) / Δ
 *   s12 = (2*z12*Z0) / Δ
 *   s21 = (2*z21*Z0) / Δ
 *   s22 = ((z11+Z0)*(z22-Z0) - z12*z21) / Δ
 *
 * where Δ = (z11+Z0)*(z22+Z0) - z12*z21
 *
 * Reference: Pozar "Microwave Engineering" (2012) §4.3
 */
int z_to_s(const ZParameters *z, double Z0, SParameters *s) {
    if (!z || !s || Z0 <= 0.0) return -1;

    double z11r = z->z11.real, z11i = z->z11.imag;
    double z12r = z->z12.real, z12i = z->z12.imag;
    double z21r = z->z21.real, z21i = z->z21.imag;
    double z22r = z->z22.real, z22i = z->z22.imag;

    /* Compute Δ = (z11+Z0)*(z22+Z0) - z12*z21 */
    double a_real = z11r + Z0, a_imag = z11i;
    double b_real = z22r + Z0, b_imag = z22i;
    /* (z11+Z0)*(z22+Z0) */
    double prod1_real = a_real*b_real - a_imag*b_imag;
    double prod1_imag = a_real*b_imag + a_imag*b_real;
    /* z12*z21 */
    double prod2_real = z12r*z21r - z12i*z21i;
    double prod2_imag = z12r*z21i + z12i*z21r;
    /* Δ = prod1 - prod2 */
    double D_real = prod1_real - prod2_real;
    double D_imag = prod1_imag - prod2_imag;
    double D_mag = D_real*D_real + D_imag*D_imag;

    if (D_mag < 1e-30) return -1;

    /* s11 = ((z11-Z0)*(z22+Z0) - z12*z21) / Δ */
    double num11_real = ((z11r-Z0)*b_real - z11i*b_imag) - prod2_real;
    double num11_imag = ((z11r-Z0)*b_imag + z11i*b_real) - prod2_imag;
    s->s11.real = (num11_real*D_real + num11_imag*D_imag) / D_mag;
    s->s11.imag = (num11_imag*D_real - num11_real*D_imag) / D_mag;

    /* s12 = 2*z12*Z0 / Δ */
    s->s12.real = (2.0*Z0) * (z12r*D_real + z12i*D_imag) / D_mag;
    s->s12.imag = (2.0*Z0) * (z12i*D_real - z12r*D_imag) / D_mag;

    /* s21 = 2*z21*Z0 / Δ */
    s->s21.real = (2.0*Z0) * (z21r*D_real + z21i*D_imag) / D_mag;
    s->s21.imag = (2.0*Z0) * (z21i*D_real - z21r*D_imag) / D_mag;

    /* s22 = ((z11+Z0)*(z22-Z0) - z12*z21) / Δ */
    double c_real = z22r - Z0, c_imag = z22i;
    double num22_real = (a_real*c_real - a_imag*c_imag) - prod2_real;
    double num22_imag = (a_real*c_imag + a_imag*c_real) - prod2_imag;
    s->s22.real = (num22_real*D_real + num22_imag*D_imag) / D_mag;
    s->s22.imag = (num22_imag*D_real - num22_real*D_imag) / D_mag;

    s->Z0 = Z0;
    return 0;
}

/* ==========================================================================
 * L2: S → Z Parameter Conversion
 * ==========================================================================
 *
 * Inverse of Z→S: Z = Z0 * (I + S) * (I - S)⁻¹
 *
 * z11 = Z0 * ((1+s11)*(1-s22) + s12*s21) / Δ'
 * z12 = Z0 * (2*s12) / Δ'
 * z21 = Z0 * (2*s21) / Δ'
 * z22 = Z0 * ((1-s11)*(1+s22) + s12*s21) / Δ'
 *
 * where Δ' = (1-s11)*(1-s22) - s12*s21
 */
int s_to_z(const SParameters *s, double Z0, ZParameters *z) {
    if (!s || !z || Z0 <= 0.0) return -1;

    double s11r = s->s11.real, s11i = s->s11.imag;
    double s12r = s->s12.real, s12i = s->s12.imag;
    double s21r = s->s21.real, s21i = s->s21.imag;
    double s22r = s->s22.real, s22i = s->s22.imag;

    /* Δ' = (1-s11)*(1-s22) - s12*s21 */
    double a_real = 1.0 - s11r, a_imag = -s11i;
    double b_real = 1.0 - s22r, b_imag = -s22i;
    double prod1_real = a_real*b_real - a_imag*b_imag;
    double prod1_imag = a_real*b_imag + a_imag*b_real;
    double prod2_real = s12r*s21r - s12i*s21i;
    double prod2_imag = s12r*s21i + s12i*s21r;
    double D_real = prod1_real - prod2_real;
    double D_imag = prod1_imag - prod2_imag;
    double D_mag = D_real*D_real + D_imag*D_imag;

    if (D_mag < 1e-30) return -1;

    /* z11 = Z0 * ((1+s11)*(1-s22) + s12*s21) / Δ' */
    double c_real = 1.0 + s11r, c_imag = s11i;
    double num11_real = (c_real*b_real - c_imag*b_imag) + prod2_real;
    double num11_imag = (c_real*b_imag + c_imag*b_real) + prod2_imag;
    z->z11.real = Z0 * (num11_real*D_real + num11_imag*D_imag) / D_mag;
    z->z11.imag = Z0 * (num11_imag*D_real - num11_real*D_imag) / D_mag;

    /* z12 = Z0 * 2*s12 / Δ' */
    z->z12.real = Z0 * 2.0 * (s12r*D_real + s12i*D_imag) / D_mag;
    z->z12.imag = Z0 * 2.0 * (s12i*D_real - s12r*D_imag) / D_mag;

    /* z21 = Z0 * 2*s21 / Δ' */
    z->z21.real = Z0 * 2.0 * (s21r*D_real + s21i*D_imag) / D_mag;
    z->z21.imag = Z0 * 2.0 * (s21i*D_real - s21r*D_imag) / D_mag;

    /* z22 = Z0 * ((1-s11)*(1+s22) + s12*s21) / Δ' */
    double d_real = 1.0 + s22r, d_imag = s22i;
    double num22_real = (a_real*d_real - a_imag*d_imag) + prod2_real;
    double num22_imag = (a_real*d_imag + a_imag*d_real) + prod2_imag;
    z->z22.real = Z0 * (num22_real*D_real + num22_imag*D_imag) / D_mag;
    z->z22.imag = Z0 * (num22_imag*D_real - num22_real*D_imag) / D_mag;

    return 0;
}

/* ==========================================================================
 * L4: Reciprocity Theorem Verification
 * ==========================================================================
 *
 * A two-port network is reciprocal iff Z_12 = Z_21.
 * Equivalent conditions for other parameter sets:
 *   Y-parameters: y12 = y21
 *   H-parameters: h12 = -h21
 *   ABCD: det(ABCD) = AD - BC = 1
 *   S-parameters: s12 = s21 (for symmetric normalization)
 *
 * Reciprocity holds for networks containing only passive elements
 * (R, L, C) and transformers. It is violated by:
 *   - Dependent sources (transistors, op-amps)
 *   - Gyrators (active circulator components)
 *   - Isolators (ferrite-based non-reciprocal devices)
 */
int verify_reciprocity(const ZParameters *zp, ReciprocityParams *result) {
    if (!zp || !result) return -1;

    result->z11 = zp->z11.real;
    result->z12 = zp->z12.real;
    result->z21 = zp->z21.real;
    result->z22 = zp->z22.real;

    /* Compute reciprocity error: |z12 - z21| */
    double d_real = zp->z12.real - zp->z21.real;
    double d_imag = zp->z12.imag - zp->z21.imag;
    result->reciprocity_error = sqrt(d_real*d_real + d_imag*d_imag);

    /* Also store Y and H parameters for completeness */
    YParameters yp;
    HParameters hp;

    if (z_to_y(zp, &yp) == 0) {
        result->y12 = yp.y12.real;
        result->y21 = yp.y21.real;
    } else {
        result->y12 = 0.0;
        result->y21 = 0.0;
    }

    if (z_to_h(zp, &hp) == 0) {
        result->h12 = hp.h12_real;
        result->h21 = hp.h21_real;
        result->h11 = hp.h11.real;
        result->h22 = hp.h22.real;
    } else {
        result->h12 = 0.0;
        result->h21 = 0.0;
        result->h11 = 0.0;
        result->h22 = 0.0;
    }

    /* Determine reciprocity with tolerance 1e-9 */
    result->is_reciprocal = (result->reciprocity_error < 1e-9) ? 1 : 0;

    return 0;
}

/* ==========================================================================
 * L6: Two-Port Network Analyzer
 * ==========================================================================
 *
 * Computes all standard metrics for a two-port network given
 * source and load terminations:
 *   - Input impedance Z_in
 *   - Output impedance Z_out
 *   - Voltage gain Av = V2/V1
 *   - Current gain Ai = I2/I1
 *   - Power gain Gp = P_load / P_input
 *   - Available power gain Ga
 *
 * Using Z-parameters and terminations Z_s (source) and Z_L (load).
 */

/**
 * Compute input impedance of a terminated two-port:
 *   Z_in = z11 - (z12*z21) / (z22 + Z_L)
 */
ComplexImpedance two_port_input_impedance(const ZParameters *zp,
                                           ComplexImpedance ZL) {
    ComplexImpedance zin;
    double denom_real = zp->z22.real + ZL.real;
    double denom_imag = zp->z22.imag + ZL.imag;
    double denom = denom_real*denom_real + denom_imag*denom_imag;

    if (denom < 1e-30) {
        zin.real = zp->z11.real;
        zin.imag = zp->z11.imag;
        return zin;
    }

    /* z12*z21 product */
    double p_real = zp->z12.real*zp->z21.real - zp->z12.imag*zp->z21.imag;
    double p_imag = zp->z12.real*zp->z21.imag + zp->z12.imag*zp->z21.real;

    /* (z12*z21) / (z22 + ZL) */
    double q_real = (p_real*denom_real + p_imag*denom_imag) / denom;
    double q_imag = (p_imag*denom_real - p_real*denom_imag) / denom;

    zin.real = zp->z11.real - q_real;
    zin.imag = zp->z11.imag - q_imag;
    return zin;
}

/**
 * Compute output impedance looking back into port 2:
 *   Z_out = z22 - (z12*z21) / (z11 + Z_s)
 */
ComplexImpedance two_port_output_impedance(const ZParameters *zp,
                                            ComplexImpedance ZS) {
    ComplexImpedance zout;
    double denom_real = zp->z11.real + ZS.real;
    double denom_imag = zp->z11.imag + ZS.imag;
    double denom = denom_real*denom_real + denom_imag*denom_imag;

    if (denom < 1e-30) {
        zout.real = zp->z22.real;
        zout.imag = zp->z22.imag;
        return zout;
    }

    double p_real = zp->z12.real*zp->z21.real - zp->z12.imag*zp->z21.imag;
    double p_imag = zp->z12.real*zp->z21.imag + zp->z12.imag*zp->z21.real;

    double q_real = (p_real*denom_real + p_imag*denom_imag) / denom;
    double q_imag = (p_imag*denom_real - p_real*denom_imag) / denom;

    zout.real = zp->z22.real - q_real;
    zout.imag = zp->z22.imag - q_imag;
    return zout;
}

/**
 * Compute voltage gain Av = V2/V1 with load Z_L:
 *   Av = (z21 * Z_L) / (z11*(z22+Z_L) - z12*z21)
 */
double two_port_voltage_gain(const ZParameters *zp, double ZL) {
    double denom_real = zp->z22.real + ZL;
    double denom_imag = zp->z22.imag;

    /* z11*(z22+ZL) */
    double t1_real = zp->z11.real*denom_real - zp->z11.imag*denom_imag;
    double t1_imag = zp->z11.real*denom_imag + zp->z11.imag*denom_real;

    /* z12*z21 */
    double t2_real = zp->z12.real*zp->z21.real - zp->z12.imag*zp->z21.imag;
    double t2_imag = zp->z12.real*zp->z21.imag + zp->z12.imag*zp->z21.real;

    /* Denominator: t1 - t2 */
    double D_real = t1_real - t2_real;
    double D_imag = t1_imag - t2_imag;
    double D_mag = D_real*D_real + D_imag*D_imag;

    if (D_mag < 1e-30) return 0.0;

    /* Numerator: z21 * ZL */
    double N_real = zp->z21.real * ZL;
    double N_imag = zp->z21.imag * ZL;

    double Av_real = (N_real*D_real + N_imag*D_imag) / D_mag;
    double Av_imag = (N_imag*D_real - N_real*D_imag) / D_mag;

    return sqrt(Av_real*Av_real + Av_imag*Av_imag);
}
