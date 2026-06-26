/**
 * @file y_params.c
 * @brief Y-Parameter Implementation — Admittance-based Two-Port Analysis
 *
 * Y-parameters are the dual of Z-parameters, ideal for parallel-connected
 * networks and nodal analysis. They are the natural representation for
 * MOSFET and BJT high-frequency small-signal models.
 */

#include "../include/y_params.h"
#include <float.h>

/* ============================================================================
 * L1: Y-parameter Creation
 * ============================================================================ */

matrix2x2_t yparams_create(complex_t y11, complex_t y12,
                            complex_t y21, complex_t y22) {
    return matrix2x2_make(y11, y12, y21, y22);
}

/**
 * Shunt admittance Yp to ground.
 *
 * Physics: Both port voltages equal Vp (across Yp).
 * Y = [[Yp, -Yp], [-Yp, Yp]]
 *
 * The negative off-diagonal signs come from the current direction convention:
 * I1 flows in, but the shunt draws current that reduces I2.
 *
 * Verification: I1 = Yp*V1 - Yp*V2 = Yp*(V1-V2) ✓
 *              I2 = -Yp*V1 + Yp*V2 = -Yp*(V1-V2) = -I1 ✓
 *
 * Course: Georgia Tech ECE 6350 — Shunt element in Y-matrix
 */
matrix2x2_t yparams_shunt_admittance(complex_t yp) {
    complex_t neg_yp = complex_make(-yp.real, -yp.imag);
    return matrix2x2_make(yp, neg_yp, neg_yp, yp);
}

/**
 * Series admittance Ys between ports.
 *
 * Y = [[Ys, -Ys], [-Ys, Ys]]
 *
 * This has the same structure as the shunt case because a series
 * admittance between ports creates the same Y-matrix pattern.
 * The physical circuit is different but the two-port behavior
 * is equivalent when both ports share a common ground.
 */
matrix2x2_t yparams_series_admittance(complex_t ys) {
    complex_t neg_ys = complex_make(-ys.real, -ys.imag);
    return matrix2x2_make(ys, neg_ys, neg_ys, ys);
}

/**
 * π-network Y-parameters.
 *
 *   y11 = Ya + Yb (admittance at port 1 with port 2 shorted: Ya to gnd, Yb to short)
 *   y12 = -Yb (current at port 1 due to V2: through Yb, opposite direction)
 *   y21 = -Yb (current at port 2 due to V1; reciprocal)
 *   y22 = Yb + Yc (admittance at port 2 with port 1 shorted)
 *
 * The π-network is THE universal Y-parameter equivalent circuit.
 * Any reciprocal network can be represented as a π.
 *
 * Course: TU Munich HF Engineering — π-Ersatzschaltbild
 */
matrix2x2_t yparams_pi_network(complex_t ya, complex_t yb, complex_t yc) {
    complex_t y11 = complex_add(ya, yb);
    complex_t neg_yb = complex_make(-yb.real, -yb.imag);
    complex_t y22 = complex_add(yb, yc);
    return matrix2x2_make(y11, neg_yb, neg_yb, y22);
}

/**
 * T-network Y-parameters (computed from Z-parameter inverse).
 *
 * ΔZ = Za*Zb + Zb*Zc + Zc*Za
 * y11 = (Zb+Zc)/ΔZ, y12 = -Zc/ΔZ, y21 = -Zc/ΔZ, y22 = (Za+Zc)/ΔZ
 */
matrix2x2_t yparams_t_network(complex_t za, complex_t zb, complex_t zc) {
    complex_t delta_z = complex_add(
        complex_mul(za, zb),
        complex_add(
            complex_mul(zb, zc),
            complex_mul(zc, za)
        )
    );

    if (complex_is_zero(delta_z, 1e-30)) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }

    complex_t inv_delta = complex_div(complex_make(1.0, 0.0), delta_z);
    complex_t zb_zc = complex_add(zb, zc);
    complex_t za_zc = complex_add(za, zc);
    complex_t neg_zc = complex_make(-zc.real, -zc.imag);

    return matrix2x2_make(
        complex_mul(zb_zc, inv_delta),
        complex_mul(neg_zc, inv_delta),
        complex_mul(neg_zc, inv_delta),
        complex_mul(za_zc, inv_delta)
    );
}

/**
 * MOSFET high-frequency Y-parameters.
 *
 * The MOSFET small-signal model at high frequencies includes
 * capacitive coupling between all terminals:
 *
 *   Cgs: Gate-to-source (channel + overlap), typically 0.1-1 pF
 *   Cgd: Gate-to-drain (Miller capacitance), typically 0.01-0.1 pF
 *   Cds: Drain-to-source (junction), typically 0.01-0.1 pF
 *   gm:  Transconductance, typically 1-100 mS
 *   go:  Output conductance = 1/ro, typically 10-100 µS
 *
 * Y-parameters:
 *   y11 = jω(Cgs+Cgd)       [gate input admittance, no DC path]
 *   y12 = -jω*Cgd           [reverse through Cgd, negative due to direction]
 *   y21 = gm - jω*Cgd       [forward: gm dominates at LF, Cgd feedthrough at HF]
 *   y22 = go + jω(Cds+Cgd)  [output admittance]
 *
 * The y12 pole through Cgd creates the Miller effect.
 * The y21 term gm - jω*Cgd shows that Cgd creates a right-half-plane zero
 * at ωz = gm/Cgd (typically in the GHz range).
 *
 * Course: Stanford EE214 — High-frequency MOSFET models
 * Ref: Razavi §6.3 — MOSFET Y-parameters
 */
matrix2x2_t yparams_mosfet_hf(double gm, double go,
                               double cgs, double cgd, double cds,
                               double omega) {
    complex_t y11 = complex_make(go * 0.001, omega * (cgs + cgd));
    /* Add a tiny real part for numerical stability; MOSFET gate is ideally open at DC */
    y11.real = 0.0;  /* Override: purely capacitive at gate */

    complex_t y12 = complex_make(0.0, -omega * cgd);
    complex_t y21 = complex_make(gm, -omega * cgd);
    complex_t y22 = complex_make(go, omega * (cds + cgd));

    return matrix2x2_make(y11, y12, y21, y22);
}

/**
 * BJT high-frequency Y-parameters.
 *
 * From the hybrid-π model:
 *   y11 = gπ + jω(Cπ + Cμ)     [input: rπ || (Cπ + Cμ)]
 *   y12 = -jω*Cμ              [reverse via base-collector capacitance]
 *   y21 = gm - jω*Cμ           [forward transconductance minus Cμ feedforward]
 *   y22 = go + jω*Cμ           [output: ro || Cμ]
 *
 * where gπ = 1/rπ, go = 1/ro.
 *
 * Important: y11 has a real part (gπ) even at DC, unlike MOSFET.
 * This is because BJT is a current-controlled device (base current required).
 *
 * Course: Berkeley EE105 — BJT frequency response
 */
matrix2x2_t yparams_bjt_hf(double rpi, double gm, double ro,
                            double cpi, double cmu, double omega) {
    double gpi = 1.0 / rpi;
    double go_v = 1.0 / ro;

    complex_t y11 = complex_make(gpi, omega * (cpi + cmu));
    complex_t y12 = complex_make(0.0, -omega * cmu);
    complex_t y21 = complex_make(gm, -omega * cmu);
    complex_t y22 = complex_make(go_v, omega * cmu);

    return matrix2x2_make(y11, y12, y21, y22);
}

/* ============================================================================
 * L3: Y-parameter Analysis
 * ============================================================================ */

/**
 * Short-circuit current gain: Ai_sc = y21/y11.
 *
 * For a BJT at low frequencies:
 *   y11 ≈ gπ = 1/rπ, y21 ≈ gm = β/rπ
 *   Ai_sc = gm/gπ = gm*rπ = β ✓
 *
 * For a MOSFET at low frequencies:
 *   y11 ≈ jω*Cgs (→ 0 at DC), y21 ≈ gm
 *   Ai_sc → ∞ at DC (ideal gate draws no current)
 *   The current gain is not a meaningful metric for MOSFETs at DC.
 */
complex_t yparams_short_circuit_current_gain(matrix2x2_t y) {
    return complex_div(y.m21, y.m11);
}

/**
 * Open-circuit voltage gain: Av_oc = -y21/y22.
 *
 * For a CS MOSFET: Av = -gm/(go + jω*Cds).
 * At DC: Av = -gm/go = -gm*ro (intrinsic gain, typically 20-100).
 * At high frequencies: |Av| rolls off.
 *
 * This is the maximum voltage gain achievable from a single transistor stage.
 * Cascoding increases the effective ro, boosting intrinsic gain.
 */
complex_t yparams_open_circuit_voltage_gain(matrix2x2_t y) {
    complex_t neg_y21 = complex_make(-y.m21.real, -y.m21.imag);
    return complex_div(neg_y21, y.m22);
}

/**
 * Input admittance: Yin = y11 - y12*y21/(y22 + YL).
 *
 * The Miller effect manifests through the y12*y21 term.
 * For a CS amplifier: y12 = -jω*Cgd, y21 ≈ gm.
 * y12*y21 = -jω*Cgd*gm (negative imaginary → positive real + positive imaginary after division)
 *
 * The real part of y12*y21/(y22+YL) represents an effective conductance
 * at the input. The imaginary part adds to Cgs.
 */
complex_t yparams_input_admittance(matrix2x2_t y, complex_t yl) {
    complex_t denom = complex_add(y.m22, yl);
    complex_t term = complex_div(complex_mul(y.m12, y.m21), denom);
    return complex_sub(y.m11, term);
}

complex_t yparams_output_admittance(matrix2x2_t y, complex_t ys) {
    complex_t denom = complex_add(y.m11, ys);
    complex_t term = complex_div(complex_mul(y.m12, y.m21), denom);
    return complex_sub(y.m22, term);
}

/**
 * Unity-gain frequency fT for MOSFET.
 *
 * At fT, |iout/iin| = 1 meaning |y21/y11| = 1.
 *
 * y11 = jω*(Cgs+Cgd), y21 = gm (approximately, neglecting Cgd feedforward)
 * |gm / (jω*(Cgs+Cgd))| = 1
 * ωT = gm / (Cgs + Cgd)
 * fT = gm / (2π * (Cgs + Cgd))
 *
 * This is the most important figure of merit for RF transistors.
 * Modern CMOS: fT = 200-400 GHz for 5nm technology.
 *
 * Ref: Lee, "The Design of CMOS RFICs", Eq. 3.15
 */
double yparams_f_t_mosfet(double gm, double cgs, double cgd) {
    double c_total = cgs + cgd;
    if (c_total < 1e-30) return INFINITY;
    return gm / (2.0 * M_PI * c_total);
}

/**
 * Maximum Available Gain from Y-parameters.
 *
 * For unconditionally stable two-ports: MAG = |y21/y12|
 * For conditionally stable: MSG = |y21/y12|
 *
 * In practice, MSG is widely used as a device figure of merit
 * independent of matching network design.
 */
double yparams_max_available_gain(matrix2x2_t y) {
    double mag_y12 = complex_mag(y.m12);
    if (mag_y12 < 1e-30) return INFINITY;  /* Unilateral → infinite gain possible */
    return complex_mag(y.m21) / mag_y12;
}

/* ============================================================================
 * L4: Admittance Analysis
 * ============================================================================ */

complex_t yparams_driving_point_admittance(matrix2x2_t y, complex_t yl) {
    return yparams_input_admittance(y, yl);
}

/**
 * Passivity test for Y-parameters.
 *
 * A Y-matrix is passive if the Hermitian part is positive semi-definite:
 *   Re(y11) ≥ 0, Re(y22) ≥ 0
 *   |y12 + y21*|² ≤ 4 * Re(y11) * Re(y22)
 *
 * The second condition ensures that the network dissipates energy
 * for all possible excitations. Violation indicates potential
 * for oscillation (negative conductance).
 *
 * Reference: Desoer & Kuh, "Basic Circuit Theory", Ch. 10
 * Course: Illinois ECE 451 — Passivity and activity criteria
 */
int yparams_is_passive(matrix2x2_t y, double tolerance) {
    /* Check diagonal real parts */
    if (y.m11.real < -tolerance || y.m22.real < -tolerance) {
        return 0;
    }

    /* Check the off-diagonal energy condition */
    complex_t y12_plus_y21_conj = complex_add(y.m12, complex_conj(y.m21));
    double lhs = complex_mag(y12_plus_y21_conj) * complex_mag(y12_plus_y21_conj);
    double rhs = 4.0 * y.m11.real * y.m22.real;

    return lhs <= rhs + tolerance;
}
