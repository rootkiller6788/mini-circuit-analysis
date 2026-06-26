/**
 * @file z_params.c
 * @brief Z-Parameter Implementation — Creation, Analysis, and Equivalent Circuits
 *
 * Implements all Z-parameter functions: network creation (R, L, C, T, π,
 * transformer, BJT, FET), analysis (gains, impedances, stability), and
 * T-network extraction.
 */

#include "../include/z_params.h"
#include <float.h>

/* ============================================================================
 * L1: Z-parameter Creation
 * ============================================================================ */

matrix2x2_t zparams_create(complex_t z11, complex_t z12,
                            complex_t z21, complex_t z22) {
    return matrix2x2_make(z11, z12, z21, z22);
}

/**
 * Series impedance Zs between two ports with common ground.
 *
 * Physics: The same current (I1 = -I2) flows through Zs.
 *   V1 - V2 = I1 * Zs, and I1 + I2 = 0 (Kirchhoff's current law at ground).
 *   Therefore: V1 = Zs*I1 + Zs*I2, V2 = Zs*I1 + Zs*I2.
 *   Z = [[Zs, Zs], [Zs, Zs]]
 *
 * Application: Modeling a bond wire or PCB trace between two stages.
 */
matrix2x2_t zparams_series_impedance(complex_t zs) {
    return matrix2x2_make(zs, zs, zs, zs);
}

/**
 * Shunt impedance Zp from line to ground.
 *
 * Physics: The same voltage V1 = V2 appears across Zp.
 *   I1 = V1/Zp, I2 = V2/Zp.
 *   Z = [[Zp, 0], [0, Zp]]? No, let me reconsider.
 *
 * For a shunt element: V1 = V2 (since the ports are directly connected).
 * The Z-matrix: with I2 = 0, V1 = z11*I1 → z11 is the impedance at port 1
 * with port 2 open. Since port 2 is open but connected to port 1 through
 * a wire, all current flows through Zp. So z11 = Zp.
 *
 * With I1 = 0, V1 = z12*I2. The voltage at both ports is I2*Zp, so z12 = Zp.
 * Similarly z21 = Zp and z22 = Zp.
 * Therefore Z = [[Zp, Zp], [Zp, Zp]].
 *
 * This is the same as series! That's correct — a shunt element between two
 * connected ports creates the same Z-matrix.
 */
matrix2x2_t zparams_shunt_impedance(complex_t zp) {
    return matrix2x2_make(zp, zp, zp, zp);
}

/**
 * T-network Z-parameters.
 *
 * Circuit:
 *   P1 — Za — node — Zb — P2
 *                |
 *               Zc
 *                |
 *               GND
 *
 * With I2 = 0 (port 2 open):
 *   V1 = I1 * (Za + Zc)  → z11 = Za + Zc
 *   V2 = I1 * Zc         → z21 = Zc
 *
 * With I1 = 0 (port 1 open):
 *   V1 = I2 * Zc         → z12 = Zc
 *   V2 = I2 * (Zb + Zc)  → z22 = Zb + Zc
 *
 * The T-network is always reciprocal (z12 = z21 = Zc).
 *
 * Course: TU Munich HF Engineering — T-Ersatzschaltbild
 */
matrix2x2_t zparams_t_network(complex_t za, complex_t zb, complex_t zc) {
    complex_t z11 = complex_add(za, zc);
    complex_t z22 = complex_add(zb, zc);
    return matrix2x2_make(z11, zc, zc, z22);
}

/**
 * π-network Z-parameters.
 *
 * Given three admittances Ya, Yb, Yc (π topology):
 * The total admittance at port 1 (V2=0): Y11 = Ya + Yb
 * The transfer admittance (V1=0): Y12 = -Yb
 * The Z-parameters are then Z = Y⁻¹:
 *
 *   ΔY = (Ya+Yb)*(Yb+Yc) - Yb² = Ya*Yb + Yb*Yc + Ya*Yc
 *   z11 = (Yb+Yc) / ΔY
 *   z12 = Yb / ΔY = z21
 *   z22 = (Ya+Yb) / ΔY
 */
matrix2x2_t zparams_pi_network(complex_t ya, complex_t yb, complex_t yc) {
    complex_t ya_yb = complex_add(ya, yb);
    complex_t yb_yc = complex_add(yb, yc);
    complex_t delta_y = complex_sub(
        complex_mul(ya_yb, yb_yc),
        complex_mul(yb, yb)
    );

    if (complex_is_zero(delta_y, 1e-30)) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }

    complex_t inv_delta = complex_div(complex_make(1.0, 0.0), delta_y);
    return matrix2x2_make(
        complex_mul(yb_yc, inv_delta),
        complex_mul(yb, inv_delta),
        complex_mul(yb, inv_delta),
        complex_mul(ya_yb, inv_delta)
    );
}

/**
 * Practical transformer Z-parameters.
 *
 * Mutual inductance model:
 *   V1 = jωL1*I1 + jωM*I2
 *   V2 = jωM*I1 + jωL2*I2
 *
 * where M = k*sqrt(L1*L2), L1 = Ll1 + Lm, L2 = Ll2 + n²*Lm.
 *
 * Adding winding resistances:
 *   z11 = Rw1 + jω(Ll1 + Lm)
 *   z12 = jωM
 *   z21 = jωM (reciprocal)
 *   z22 = Rw2 + jω(Ll2 + n²*Lm)
 *
 * For a 1:1 ideal transformer (k=1, Ll1=Ll2=0, Rw1=Rw2=0):
 *   z11 = jωLm, z12 = jωLm, z21 = jωLm, z22 = jωLm
 *   (approaches the ideal as Lm → ∞)
 *
 * Course: MIT 6.002 — Transformer equivalent circuit
 * Ref: Erickson & Maksimovic §13.2
 */
matrix2x2_t zparams_transformer(double rw1, double rw2,
                                 double l_lk1, double l_lk2,
                                 double lm, double k, double n,
                                 double omega) {
    double l1 = l_lk1 + lm;
    double l2 = l_lk2 + n * n * lm;
    double m = k * sqrt(l1 * l2);

    complex_t z11 = complex_make(rw1, omega * l1);
    complex_t z12 = complex_make(0.0, omega * m);
    complex_t z21 = complex_make(0.0, omega * m);
    complex_t z22 = complex_make(rw2, omega * l2);

    return matrix2x2_make(z11, z12, z21, z22);
}

/**
 * Common-emitter BJT Z-parameters (low-frequency).
 *
 * From the hybrid-π model:
 *   z11 = rπ         (input: base-emitter diode resistance)
 *   z12 = 0          (no intrinsic feedback at LF, ignoring Early effect on rπ)
 *   z21 = -gm * rπ * ro = -β * ro  (forward transimpedance)
 *   z22 = ro         (output: collector-emitter resistance)
 *
 * β = gm * rπ. Typical values: β=100, ro=50kΩ → z21 = -5 MΩ!
 * This huge forward transimpedance reflects the high gain of a BJT.
 *
 * Course: Berkeley EE105 — BJT small-signal model
 */
matrix2x2_t zparams_bjt_ce(double rpi, double gm, double ro) {
    /* z21 = -gm * rpi * ro = -beta * ro */
    double forward_real = -gm * rpi * ro;
    return matrix2x2_make(
        complex_make(rpi, 0.0),       /* z11 = rπ */
        complex_make(0.0, 0.0),       /* z12 = 0 */
        complex_make(forward_real, 0.0), /* z21 = -β*ro */
        complex_make(ro, 0.0)         /* z22 = ro */
    );
}

/**
 * Common-source FET Z-parameters (low-frequency).
 *
 * At DC: gate is open → z11 ≈ ∞.
 * We use a finite large number (1e12) for numerical practicality.
 *
 *   z21 = -gm * ro (forward transimpedance, like BJT but without β factor)
 *   z22 = ro (output impedance)
 *
 * MOSFET intrinsic gain: |z21/z11| = gm*ro (since z11→large, this is
 * approximately the open-circuit voltage gain limit).
 */
matrix2x2_t zparams_fet_cs(double gm, double ro) {
    return matrix2x2_make(
        complex_make(1e12, 0.0),       /* z11 ≈ ∞ (gate open at DC) */
        complex_make(0.0, 0.0),        /* z12 = 0 */
        complex_make(-gm * ro, 0.0),   /* z21 = -gm*ro */
        complex_make(ro, 0.0)          /* z22 = ro */
    );
}

/* ============================================================================
 * L3: Z-parameter Analysis
 * ============================================================================ */

matrix2x2_t zparams_resistor(double r) {
    complex_t zr = complex_make(r, 0.0);
    return matrix2x2_make(zr, zr, zr, zr);
}

/**
 * Inductor impedance: ZL = jωL.
 * Pure imaginary → energy stored, not dissipated.
 * Phase: +90° (current lags voltage).
 */
matrix2x2_t zparams_inductor(double l, double omega) {
    complex_t zl = complex_make(0.0, omega * l);
    return matrix2x2_make(zl, zl, zl, zl);
}

/**
 * Capacitor impedance: ZC = 1/(jωC) = -j/(ωC).
 * Pure imaginary → energy stored, not dissipated.
 * Phase: -90° (current leads voltage).
 */
matrix2x2_t zparams_capacitor(double c, double omega) {
    if (omega < 1e-30 || c < 1e-30) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }
    complex_t zc = complex_make(0.0, -1.0 / (omega * c));
    return matrix2x2_make(zc, zc, zc, zc);
}

/**
 * Series RLC impedance:
 *   Z = R + jωL + 1/(jωC) = R + j(ωL - 1/(ωC))
 *
 * Resonance: ω0 = 1/√(LC) → Im(Z) = 0, Z = R (purely resistive).
 * Below ω0: capacitive (negative reactance).
 * Above ω0: inductive (positive reactance).
 *
 * Quality factor at resonance: Q = ω0*L/R = 1/(ω0*C*R).
 * The Q determines bandwidth: BW = ω0/Q.
 */
matrix2x2_t zparams_series_rlc(double r, double l, double c, double omega) {
    double x;
    if (omega < 1e-30 || c < 1e-30) {
        x = omega * l;  /* Approximate as inductor-only if C is tiny */
    } else {
        x = omega * l - 1.0 / (omega * c);
    }
    complex_t z = complex_make(r, x);
    return matrix2x2_make(z, z, z, z);
}

/**
 * Parallel RLC impedance:
 *   Z = 1 / (1/R + 1/(jωL) + jωC)
 *     = 1 / (1/R + j(ωC - 1/(ωL)))
 *
 * At resonance (ω0 = 1/√(LC)): Z = R (maximum, pure real).
 * Off resonance: |Z| drops.
 *
 * The parallel RLC is the tank circuit used in oscillators and
 * narrowband amplifiers. Its Q determines frequency selectivity.
 */
matrix2x2_t zparams_parallel_rlc(double r, double l, double c, double omega) {
    if (r < 1e-30) {
        complex_t zero = complex_make(0.0, 0.0);
        return matrix2x2_make(zero, zero, zero, zero);
    }

    double bl;
    if (omega < 1e-30 || l < 1e-30) {
        bl = INFINITY;  /* DC short through inductor */
        complex_t zero = complex_make(0.0, 0.0);
        return matrix2x2_make(zero, zero, zero, zero);
    } else {
        bl = 1.0 / (omega * l);
    }

    double bc = omega * c;
    double g = 1.0 / r;
    double susceptance = bc - bl;

    /* Z = 1/(G + jB) = (G - jB)/(G² + B²) */
    double denom = g * g + susceptance * susceptance;
    complex_t z = complex_make(g / denom, -susceptance / denom);
    return matrix2x2_make(z, z, z, z);
}

complex_t zparams_open_circuit_gain(matrix2x2_t z) {
    return complex_div(z.m21, z.m11);
}

complex_t zparams_short_circuit_gain(matrix2x2_t z) {
    complex_t neg_z21 = complex_make(-z.m21.real, -z.m21.imag);
    return complex_div(neg_z21, z.m22);
}

/**
 * Extract T-network elements from Z-parameters.
 *
 * For a reciprocal network (z12 = z21):
 *   Za = z11 - z12  (input series arm)
 *   Zb = z22 - z12  (output series arm)
 *   Zc = z12        (shunt arm)
 *
 * Verification:
 *   z11 = Za + Zc = (z11 - z12) + z12 = z11 ✓
 *   z12 = Zc = z12 ✓
 *   z22 = Zb + Zc = (z22 - z12) + z12 = z22 ✓
 *
 * The T-equivalent exists for any reciprocal two-port.
 * For non-reciprocal networks, a T with dependent sources is needed.
 *
 * Course: Illinois ECE 451 — Equivalent networks
 */
int zparams_to_t_network(matrix2x2_t z, complex_t *za, complex_t *zb,
                         complex_t *zc) {
    /* Check reciprocity */
    complex_t diff = complex_sub(z.m12, z.m21);
    if (!complex_approx_equal(diff, complex_make(0, 0), 1e-9)) {
        return -1;  /* non-reciprocal */
    }

    *za = complex_sub(z.m11, z.m12);
    *zb = complex_sub(z.m22, z.m12);
    *zc = z.m12;
    return 0;
}

/**
 * Input impedance with terminated output: Zin = z11 - z12*z21/(z22+ZL)
 * See two_port.c:two_port_zin_from_zparams for derivation.
 */
complex_t zparams_input_impedance(matrix2x2_t z, complex_t zl) {
    return two_port_zin_from_zparams(z, zl);
}

complex_t zparams_output_impedance(matrix2x2_t z, complex_t zs) {
    return two_port_zout_from_zparams(z, zs);
}

complex_t zparams_voltage_gain(matrix2x2_t z, complex_t zl) {
    return two_port_voltage_gain_z(z, zl);
}

/**
 * Reverse voltage gain: V1/V2 when driving from port 2.
 * Formula: Av_rev = z12 * ZS / (z22 * z11 - z12 * z21 + z22 * ZS)
 *
 * By symmetry with the forward gain, swap (z11↔z22) and the load ZL is
 * replaced by the source termination ZS at port 1.
 */
complex_t zparams_reverse_voltage_gain(matrix2x2_t z, complex_t zs) {
    complex_t num = complex_mul(z.m12, zs);
    complex_t z11_zs = complex_add(z.m11, zs);
    complex_t denom = complex_sub(
        complex_mul(z.m22, z11_zs),
        complex_mul(z.m12, z.m21)
    );
    return complex_div(num, denom);
}

/**
 * Rollett K-factor adapted for Z-parameters.
 *
 * K_z = (2*Re(z11)*Re(z22) - Re(z12*z21)) / |z12*z21|
 *
 * This is derived from the general Rollett stability condition
 * applied to the Z-matrix representation. It detects potential
 * oscillation in active two-ports.
 *
 * K > 1 indicates unconditional stability.
 *
 * Reference: Rollett (1962), generalized to impedance parameters
 */
double zparams_rollett_k(matrix2x2_t z) {
    double re_z11 = z.m11.real;
    double re_z22 = z.m22.real;
    double re_z12z21 = (z.m12.real * z.m21.real) -
                       (z.m12.imag * z.m21.imag);
    double mag_z12z21 = complex_mag(complex_mul(z.m12, z.m21));

    if (mag_z12z21 < 1e-30) {
        /* Zero feedback → unconditionally stable */
        return INFINITY;
    }

    return (2.0 * re_z11 * re_z22 - re_z12z21) / mag_z12z21;
}
