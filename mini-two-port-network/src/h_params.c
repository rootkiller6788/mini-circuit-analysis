/**
 * @file h_params.c
 * @brief H-Parameter Implementation — Hybrid Two-Port Analysis
 *
 * H-parameters are the classic BJT small-signal model, fundamental to
 * transistor amplifier design. This module covers CE, CB, CC configurations
 * and their frequency-dependent behavior.
 */

#include "../include/h_params.h"
#include <float.h>

/* ============================================================================
 * L1: H-parameter Creation
 * ============================================================================ */

matrix2x2_t hparams_create(complex_t h11, complex_t h12,
                            complex_t h21, complex_t h22) {
    return matrix2x2_make(h11, h12, h21, h22);
}

/**
 * Common-emitter BJT H-parameters.
 *
 * Mapping from hybrid-π to h-parameters:
 *   h_ie = rbb' + rπ     (input impedance, base to emitter)
 *   h_re = rπ / (rπ + ro) ≈ 0  (reverse voltage transfer, typically ~10⁻⁴)
 *   h_fe = β = Ic/Ib      (forward current gain, DC current gain)
 *   h_oe = 1/ro + gm*h_re ≈ 1/ro  (output admittance, typically ~10⁻⁵ S)
 *
 * The approximation h_re ≈ 0 is excellent for most hand calculations.
 * More precisely, h_re comes from the Early effect voltage feedback:
 * An increase in Vce slightly increases Ic, which slightly changes Vbe.
 *
 * Course: Berkeley EE105 — Transistor h-parameter model
 * Ref: Sedra & Smith §7.3.2
 */
matrix2x2_t hparams_bjt_ce(double rbb, double rpi, double beta, double ro) {
    double h_ie = rbb + rpi;
    double h_re;
    if (ro > 0 && rpi > 0) {
        h_re = rpi / (rpi + ro);
    } else {
        h_re = 0.0;
    }
    double h_fe = beta;
    double h_oe = 1.0 / ro;

    return matrix2x2_make(
        complex_make(h_ie, 0.0),    /* h11 = h_ie in Ω */
        complex_make(h_re, 0.0),    /* h12 = h_re dimensionless */
        complex_make(h_fe, 0.0),    /* h21 = h_fe dimensionless */
        complex_make(h_oe, 0.0)     /* h22 = h_oe in S */
    );
}

/**
 * CE → CB conversion formulas:
 *
 *   h_ib = h_ie / (1 + h_fe)         (CB input: base-grounded, input at emitter)
 *   h_rb = (h_ie*h_oe/(1+h_fe)) - h_re  (CB reverse voltage, very small)
 *   h_fb = -h_fe / (1 + h_fe)        (CB current gain ≈ -α)
 *   h_ob = h_oe / (1 + h_fe)         (CB output admittance, very small)
 *
 * For β = 100: α = 100/101 ≈ 0.99, so h_fb ≈ -0.99.
 * CB input impedance: h_ib = rπ/β ≈ 1/gm (typically 25Ω at Ic=1mA).
 * This low input impedance makes CB useful for current buffers and
 * wideband amplifiers (no Miller effect).
 *
 * Reference: Millman & Halkias, "Integrated Electronics", Ch. 8
 * Course: MIT 6.002 — Common-base configuration
 */
matrix2x2_t hparams_ce_to_cb(matrix2x2_t h_ce) {
    double hie = h_ce.m11.real;
    double hre = h_ce.m12.real;
    double hfe = h_ce.m21.real;
    double hoe = h_ce.m22.real;

    double factor = 1.0 + hfe;

    double hib = hie / factor;
    double hrb = (hie * hoe) / factor - hre;
    double hfb = -hfe / factor;
    double hob = hoe / factor;

    return matrix2x2_make(
        complex_make(hib, 0.0),
        complex_make(hrb, 0.0),
        complex_make(hfb, 0.0),
        complex_make(hob, 0.0)
    );
}

/**
 * CE → CC (emitter follower) conversion:
 *
 *   h_ic = h_ie                     (CC input: base to collector, high)
 *   h_rc = 1 - h_re ≈ 1            (CC voltage gain, just under unity)
 *   h_fc = -(1 + h_fe)              (CC current gain, ≈ -β)
 *   h_oc = h_oe                     (CC output admittance, same as CE)
 *
 * The emitter follower provides:
 *   - Very high input impedance (≈ hfe*RE for RE in emitter)
 *   - Voltage gain slightly less than 1
 *   - Low output impedance (≈ 1/gm)
 *   - Current gain ≈ β
 *
 * Used as: voltage buffers, impedance transformers, output stages.
 *
 * Course: Stanford EE101 — Emitter follower
 */
matrix2x2_t hparams_ce_to_cc(matrix2x2_t h_ce) {
    double hie = h_ce.m11.real;
    double hre = h_ce.m12.real;
    double hfe = h_ce.m21.real;
    double hoe = h_ce.m22.real;

    double hic = hie;
    double hrc = 1.0 - hre;
    double hfc = -(1.0 + hfe);
    double hoc = hoe;

    return matrix2x2_make(
        complex_make(hic, 0.0),
        complex_make(hrc, 0.0),
        complex_make(hfc, 0.0),
        complex_make(hoc, 0.0)
    );
}

/**
 * High-frequency BJT H-parameters.
 *
 * At high frequencies, the h-parameters become complex due to
 * transistor capacitances:
 *
 *   h11(ω) = rbb' + rπ / (1 + jω*rπ*(Cπ+Cμ))  [Zπ(ω) = rπ || (Cπ+Cμ)]
 *   h12(ω) = jω*Cμ / (gπ + jω*Cπ)           [Through Cμ]
 *   h21(ω) = hfe / (1 + jω/ωβ)              [Beta rolloff]
 *   h22(ω) = hoe + jω*Cμ                     [Output admittance]
 *
 * where ωβ = 1/(rπ*Cπ) ≈ 2π*fβ.
 */
matrix2x2_t hparams_bjt_ce_hf(double rbb, double rpi, double beta, double ro,
                               double cpi, double cmu, double omega) {
    /* Zπ(ω) = rπ || (Cπ+Cμ) = rπ / (1 + jω*rπ*(Cπ+Cμ)) */
    double gpi = 1.0 / rpi;
    double c_total = cpi + cmu;
    double tau_pi = rpi * c_total;  /* RC time constant of base node */

    double denom_mag = 1.0 + (omega * tau_pi) * (omega * tau_pi);
    double real_zpi = rpi / denom_mag;
    double imag_zpi = -omega * tau_pi * rpi / denom_mag;

    complex_t h11 = complex_add(
        complex_make(rbb, 0.0),
        complex_make(real_zpi, imag_zpi)
    );

    /* h12(ω): voltage feedback through Cμ */
    /* h12 ≈ jω*Cμ * Zπ(ω) for small hre */
    complex_t y_mu = complex_make(0.0, omega * cmu);
    complex_t zpi = complex_make(real_zpi, imag_zpi);
    complex_t h12 = complex_mul(y_mu, zpi);

    /* h21(ω) = hfe / (1 + jω/ωβ) */
    double omega_beta = gpi / cpi;  /* ωβ = 1/(rπ*Cπ) */
    double denom_h21_mag = 1.0 + (omega / omega_beta) * (omega / omega_beta);
    double real_h21 = beta / denom_h21_mag;
    double imag_h21 = -beta * (omega / omega_beta) / denom_h21_mag;
    complex_t h21 = complex_make(real_h21, imag_h21);

    /* h22(ω) = hoe + jω*Cμ */
    complex_t h22 = complex_make(1.0 / ro, omega * cmu);

    return matrix2x2_make(h11, h12, h21, h22);
}

/* ============================================================================
 * L3: H-parameter Analysis
 * ============================================================================ */

/**
 * Voltage gain: Av = -h21 / (h11*h22 - h12*h21 + h11*YL).
 *
 * For typical BJT (h12 ≈ 0, h22 small):
 * Av ≈ -h21 / (h11 * YL) = -hfe * ZL / hie = -gm * ZL
 */
complex_t hparams_voltage_gain(matrix2x2_t h, complex_t yl) {
    complex_t delta_h = complex_sub(
        complex_mul(h.m11, h.m22),
        complex_mul(h.m12, h.m21)
    );
    complex_t denom = complex_add(delta_h, complex_mul(h.m11, yl));
    complex_t neg_h21 = complex_make(-h.m21.real, -h.m21.imag);
    return complex_div(neg_h21, denom);
}

/**
 * Current gain: Ai = h21 / (1 + h22*ZL).
 *
 * For typical loads (|ZL| << 1/h22): Ai ≈ h21 = hfe.
 */
complex_t hparams_current_gain(matrix2x2_t h, complex_t zl) {
    complex_t denom = complex_add(
        complex_make(1.0, 0.0),
        complex_mul(h.m22, zl)
    );
    return complex_div(h.m21, denom);
}

complex_t hparams_input_impedance(matrix2x2_t h, complex_t yl) {
    complex_t denom = complex_add(h.m22, yl);
    complex_t term = complex_div(complex_mul(h.m12, h.m21), denom);
    return complex_sub(h.m11, term);
}

complex_t hparams_output_admittance(matrix2x2_t h, complex_t zs) {
    complex_t denom = complex_add(h.m11, zs);
    complex_t term = complex_div(complex_mul(h.m12, h.m21), denom);
    return complex_sub(h.m22, term);
}

/**
 * Power gain from H-parameters:
 *   Gp = |h21|² * Re(YL) / (Re(Yin) * |1 + h22*ZL|²)
 */
double hparams_power_gain(matrix2x2_t h, complex_t yl) {
    complex_t yin = hparams_input_impedance(h, yl);
    /* Yin is actually Zin in h-parameter convention — the input is impedance */
    /* For power: Zin = h11 - h12*h21/(h22+YL), Re(Zin) is the real part */
    double re_zin = yin.real;  /* h11 is Ω, so yin is actually Zin */
    double re_yl = yl.real;
    complex_t zl = complex_div(complex_make(1.0, 0.0), yl);

    double h21_mag_sq = complex_mag(h.m21) * complex_mag(h.m21);
    complex_t denom_comp = complex_add(
        complex_make(1.0, 0.0),
        complex_mul(h.m22, zl)
    );
    double denom_mag_sq = complex_mag(denom_comp) * complex_mag(denom_comp);

    if (re_zin < 1e-30 || denom_mag_sq < 1e-30) return 0.0;
    return (h21_mag_sq * re_yl) / (re_zin * denom_mag_sq);
}

/**
 * Beta cutoff frequency: fβ = 1/(2π * rπ * Cπ).
 *
 * At fβ, the current gain drops by 3 dB from its DC value.
 * fT = β * fβ ≈ 1/(2π*(Cπ/β)) = gm/(2π*Cπ) approximately.
 *
 * For a typical BJT: rπ = 2.5kΩ (β=100 at Ic=1mA), Cπ = 10pF
 * fβ = 1/(2π*2500*1e-11) = 6.37 MHz.
 * fT = 100 * 6.37 MHz = 637 MHz.
 */
double hparams_beta_cutoff_freq(double rpi, double cpi) {
    if (rpi < 1e-30 || cpi < 1e-30) return INFINITY;
    return 1.0 / (2.0 * M_PI * rpi * cpi);
}

/**
 * hfe magnitude at frequency f: |hfe(f)| = β_DC / sqrt(1 + (f/fβ)²).
 *
 * Single-pole magnitude response (6 dB/octave rolloff beyond fβ).
 */
double hparams_hfe_magnitude(double beta_dc, double f_beta, double freq) {
    if (f_beta < 1e-30) return beta_dc;  /* DC case */
    double ratio = freq / f_beta;
    return beta_dc / sqrt(1.0 + ratio * ratio);
}
