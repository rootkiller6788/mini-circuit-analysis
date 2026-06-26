/**
 * @file s_params.c
 * @brief S-Parameter Implementation — Scattering Matrix Analysis
 *
 * S-parameters are the de facto standard for RF and microwave two-port
 * characterization. They describe the network in terms of incident and
 * reflected power waves, making them directly measurable with VNAs.
 */

#include "../include/s_params.h"
#include "../include/abcd_params.h"
#include <float.h>

/* ============================================================================
 * L1: S-parameter Creation
 * ============================================================================ */

matrix2x2_t sparams_create(complex_t s11, complex_t s12,
                            complex_t s21, complex_t s22) {
    return matrix2x2_make(s11, s12, s21, s22);
}

/**
 * S-parameters for a series impedance Z in a Z0 system.
 *
 * Derivation using reflection/transmission theory:
 *
 * A series Z in a Z0 transmission line:
 *   s11 = s22 = Z/(2*Z0 + Z)   [reflection, symmetric]
 *   s21 = s12 = 2*Z0/(2*Z0+Z)  [transmission]
 *
 * Check: Z = 0 (short) → s11 = 0/(2Z0) = 0, s21 = 2Z0/(2Z0) = 1 ✓
 *        Z = ∞ (open) → s11 = 1, s21 = 0 ✓
 *        Z = Z0 → s11 = Z0/(3Z0) = 1/3, s21 = 2Z0/(3Z0) = 2/3
 *
 * Power conservation: |s11|² + |s21|² = |Z|² + 4Z0²... wait let me check.
 * For passive Z (real): |s11|²+|s21|² = Z²/(2Z0+Z)² + 4Z0²/(2Z0+Z)²
 * = (Z²+4Z0²)/(2Z0+Z)² ≠ 1 in general. This is correct — a series resistor
 * dissipates power, so the network is not lossless.
 *
 * For reactive Z (Z = jX): s11 = jX/(2Z0+jX), s21 = 2Z0/(2Z0+jX).
 * Then |s11|²+|s21|² = X²/(4Z0²+X²) + 4Z0²/(4Z0²+X²) = 1 ✓ (lossless).
 *
 * Course: Georgia Tech ECE 6350 — Series element S-parameters
 */
matrix2x2_t sparams_series_z(complex_t z, double z0) {
    complex_t two_z0 = complex_make(2.0 * z0, 0.0);
    complex_t denom = complex_add(z, two_z0);
    complex_t s11 = complex_div(z, denom);
    complex_t s21 = complex_div(two_z0, denom);
    return matrix2x2_make(s11, s21, s21, s11);
}

/**
 * S-parameters for a shunt admittance Y in a Z0 system.
 *
 * Derivation: The shunt Y at a node with two Z0 lines on each side.
 * The total admittance at the node: Y_total = Y + 2/Z0.
 *
 *   s11 = s22 = (-Y) / (Y + 2/Z0) = -Y*Z0 / (2 + Y*Z0)
 *   s21 = s12 = (2/Z0) / (Y + 2/Z0) = 2 / (2 + Y*Z0)
 *
 * Check: Y = 0 (open circuit shunt = no element) → s11 = 0, s21 = 1 ✓
 *        Y = ∞ (short to ground) → s11 → -1, s21 → 0 ✓
 */
matrix2x2_t sparams_shunt_y(complex_t y, double z0) {
    complex_t yz0 = complex_make(y.real * z0, y.imag * z0);
    complex_t two = complex_make(2.0, 0.0);
    complex_t denom = complex_add(two, yz0);

    complex_t neg_yz0 = complex_make(-yz0.real, -yz0.imag);
    complex_t s11 = complex_div(neg_yz0, denom);
    complex_t s21 = complex_div(two, denom);

    return matrix2x2_make(s11, s21, s21, s11);
}

/**
 * S-parameters for a transmission line section.
 *
 * If Z0 = Z0_sys (matched): s11 = s22 = 0, s21 = s12 = e^(-γl).
 * The line is perfectly matched — no reflections.
 *
 * If Z0 ≠ Z0_sys: reflections occur at both ends.
 * This is handled by the impedance-mismatch boundary conditions
 * at the interface between the system impedance and the line.
 *
 * For the matched case (Z0_line = Z0_sys):
 *   s21 = e^(-αl) * e^(-jβl) = e^(-αl) * (cos(βl) - j*sin(βl))
 *
 * Attenuation: |s21| = e^(-αl) = -α*l*8.686 dB.
 * Phase shift: arg(s21) = -βl radians.
 */
matrix2x2_t sparams_transmission_line(double z0, double alpha,
                                       double beta, double length,
                                       double z0_sys) {
    /* If line Z0 equals system Z0 → perfect match */
    if (fabs(z0 - z0_sys) < 1e-9) {
        double mag = exp(-alpha * length);
        double phase = -beta * length;
        complex_t s21 = complex_make(mag * cos(phase), mag * sin(phase));
        complex_t s11 = complex_make(0.0, 0.0);
        return matrix2x2_make(s11, s21, s21, s11);
    }

    /* For mismatched line: use ABCD conversion */
    matrix2x2_t abcd = abcd_transmission_line(z0, alpha, beta, length);
    double z0r = z0_sys;
    complex_t a = abcd.m11;
    complex_t b = abcd.m12;
    complex_t c = abcd.m21;
    complex_t d = abcd.m22;

    complex_t b_div_z0 = complex_make(b.real / z0r, b.imag / z0r);
    complex_t c_mul_z0 = complex_make(c.real * z0r, c.imag * z0r);

    complex_t denom = complex_add(
        complex_add(a, b_div_z0),
        complex_add(c_mul_z0, d)
    );

    /* s11 = (A + B/Z0 - C*Z0 - D) / denom */
    complex_t s11_num = complex_sub(
        complex_add(a, b_div_z0),
        complex_add(c_mul_z0, d)
    );

    /* s22 = (-A + B/Z0 - C*Z0 + D) / denom */
    complex_t neg_a = complex_make(-a.real, -a.imag);
    complex_t s22_num = complex_add(
        complex_add(neg_a, b_div_z0),
        complex_sub(d, c_mul_z0)
    );

    /* s21 = 2*(AD-BC) / denom */
    complex_t det = complex_sub(complex_mul(a, d), complex_mul(b, c));
    complex_t s21_num = complex_make(2.0 * det.real, 2.0 * det.imag);

    /* s12 = 2 / denom */
    complex_t s12_num = complex_make(2.0, 0.0);

    complex_t s11 = complex_div(s11_num, denom);
    complex_t s12 = complex_div(s12_num, denom);
    complex_t s21 = complex_div(s21_num, denom);
    complex_t s22 = complex_div(s22_num, denom);

    return matrix2x2_make(s11, s12, s21, s22);
}

/**
 * Ideal attenuator: |s21| = 10^(-ATT/20).
 *
 *   s11 = s22 = 0 (matched)
 *   s12 = s21 = 10^(-ATT/20) ∠0°
 *
 * For ATT = 3 dB: |s21| = 10^(-0.15) = 0.7079 → power ratio = 0.5.
 * For ATT = 10 dB: |s21| = 0.3162 → power ratio = 0.1.
 * For ATT = 20 dB: |s21| = 0.1 → power ratio = 0.01.
 */
matrix2x2_t sparams_attenuator(double att_db) {
    double s21_mag = pow(10.0, -att_db / 20.0);
    complex_t s21 = complex_make(s21_mag, 0.0);
    complex_t s11 = complex_make(0.0, 0.0);
    return matrix2x2_make(s11, s21, s21, s11);
}

/**
 * Amplifier S-parameters with finite isolation and reflections.
 *
 * Real amplifiers have:
 *   |s21| > 1 (or 0 dB → gain)
 *   |s11| > 0 (input mismatch)
 *   |s22| > 0 (output mismatch)
 *   |s12| > 0 (finite reverse isolation)
 */
matrix2x2_t sparams_amplifier(complex_t gain, complex_t s11_input,
                               complex_t s22_input, complex_t s12_input) {
    return matrix2x2_make(s11_input, s12_input, gain, s22_input);
}

/* ============================================================================
 * L3: S-parameter Analysis
 * ============================================================================ */

double sparams_s21_db(matrix2x2_t s) {
    double mag = complex_mag(s.m21);
    if (mag < 1e-30) return -INFINITY;
    return 20.0 * log10(mag);
}

double sparams_input_return_loss_db(matrix2x2_t s) {
    double mag = complex_mag(s.m11);
    if (mag < 1e-30) return INFINITY;
    return -20.0 * log10(mag);
}

double sparams_output_return_loss_db(matrix2x2_t s) {
    double mag = complex_mag(s.m22);
    if (mag < 1e-30) return INFINITY;
    return -20.0 * log10(mag);
}

double sparams_isolation_db(matrix2x2_t s) {
    double mag = complex_mag(s.m12);
    if (mag < 1e-30) return INFINITY;
    return -20.0 * log10(mag);
}

/**
 * Maximum Stable Gain: MSG = |s21| / |s12|.
 *
 * MSG is the maximum gain for which K = 1 can be achieved by
 * adding loss (resistive loading). It's a widely-used figure of
 * merit for potentially unstable devices.
 *
 * When K > 1: MAG = MSG * (K - sqrt(K²-1)).
 * When K < 1: MSG is the gain with K stabilized to exactly 1.
 */
double sparams_msg(matrix2x2_t s) {
    double mag_s12 = complex_mag(s.m12);
    if (mag_s12 < 1e-30) return INFINITY;
    return complex_mag(s.m21) / mag_s12;
}

/**
 * Maximum Available Gain: MAG = |s21/s12| * (K - sqrt(K²-1)).
 *
 * Valid only when K > 1. For K ≤ 1, returns -1 to indicate
 * that MAG is not defined (use MSG instead).
 */
double sparams_mag(matrix2x2_t s) {
    double k = sparams_rollett_k(s);
    double msg = sparams_msg(s);

    if (k <= 1.0) return -1.0;  /* MAG undefined for K ≤ 1 */
    if (msg > 1e12) return INFINITY;

    return msg * (k - sqrt(k * k - 1.0));
}

/**
 * Γin = s11 + (s12*s21*ΓL) / (1 - s22*ΓL).
 *
 * This shows how the input reflection depends on the load.
 * The load-pull effect: changing the output match changes the
 * input impedance, which is critical for power amplifier design.
 *
 * For well-matched devices (small s11, s22): Γin ≈ s12*s21*ΓL.
 */
complex_t sparams_gamma_in(matrix2x2_t s, complex_t gl) {
    complex_t one = complex_make(1.0, 0.0);
    complex_t denom = complex_sub(one, complex_mul(s.m22, gl));
    complex_t term = complex_div(
        complex_mul(complex_mul(s.m12, s.m21), gl),
        denom
    );
    return complex_add(s.m11, term);
}

complex_t sparams_gamma_out(matrix2x2_t s, complex_t gs) {
    complex_t one = complex_make(1.0, 0.0);
    complex_t denom = complex_sub(one, complex_mul(s.m11, gs));
    complex_t term = complex_div(
        complex_mul(complex_mul(s.m12, s.m21), gs),
        denom
    );
    return complex_add(s.m22, term);
}

/**
 * Rollett K-factor: the most important stability metric.
 *
 * K = (1 - |s11|² - |s22|² + |Δ|²) / (2 * |s12*s21|)
 *
 * where Δ = s11*s22 - s12*s21.
 *
 * Physical meaning: K > 1 means the device is potentially
 * unconditionally stable (together with the auxiliary condition).
 *
 * K → ∞: The device is perfectly unilateral (s12 = 0), unconditionally stable.
 * K = 1: Marginal stability.
 * K < 1: Potentially unstable for some source/load terminations.
 *
 * Reference: Rollett (1962)
 */
double sparams_rollett_k(matrix2x2_t s) {
    double s11_mag_sq = complex_mag(s.m11) * complex_mag(s.m11);
    double s22_mag_sq = complex_mag(s.m22) * complex_mag(s.m22);
    complex_t delta = matrix2x2_det(s);
    double delta_mag_sq = complex_mag(delta) * complex_mag(delta);
    double s12s21_mag = complex_mag(s.m12) * complex_mag(s.m21);

    if (s12s21_mag < 1e-30) return INFINITY;  /* Unilateral → unconditionally stable */

    return (1.0 - s11_mag_sq - s22_mag_sq + delta_mag_sq) / (2.0 * s12s21_mag);
}

double sparams_delta_mag(matrix2x2_t s) {
    return complex_mag(matrix2x2_det(s));
}

/**
 * μ-factor (Edwards-Sinsky, 1992).
 *
 * μ = (1 - |s11|²) / (|s22 - conj(s11)*Δ| + |s12*s21|)
 *
 * The μ-factor has the advantage of being a single criterion:
 * μ > 1 ↔ unconditionally stable.
 *
 * Geometric meaning: μ is the minimum distance from the center
 * of the Smith chart to the unstable region. μ > 1 means the
 * entire Smith chart center is in the stable region.
 */
double sparams_mu_factor(matrix2x2_t s) {
    double s11_mag_sq = complex_mag(s.m11) * complex_mag(s.m11);
    complex_t delta = matrix2x2_det(s);
    complex_t s11_conj = complex_conj(s.m11);
    complex_t term = complex_sub(s.m22, complex_mul(s11_conj, delta));
    double denom = complex_mag(term) + complex_mag(s.m12) * complex_mag(s.m21);

    if (denom < 1e-30) return INFINITY;
    return (1.0 - s11_mag_sq) / denom;
}

/**
 * Source stability circle.
 *
 * Cs = conj(s11 - Δ*s22*) / (|s11|² - |Δ|²)
 * Rs = |s12*s21| / | |s11|² - |Δ|² |
 *
 * If |Cs| - Rs > 1: circle entirely outside Smith chart → unconditionally stable.
 * If Rs + |Cs| < 1: circle entirely inside → potentially unstable.
 */
void sparams_source_stability_circle(matrix2x2_t s, complex_t *cs, double *rs) {
    complex_t delta = matrix2x2_det(s);
    complex_t s22_conj = complex_conj(s.m22);
    complex_t num = complex_conj(complex_sub(s.m11, complex_mul(delta, s22_conj)));
    double s11_mag_sq = complex_mag(s.m11) * complex_mag(s.m11);
    double delta_mag_sq = complex_mag(delta) * complex_mag(delta);
    double denom = s11_mag_sq - delta_mag_sq;

    if (fabs(denom) < 1e-30) {
        *cs = complex_make(INFINITY, INFINITY);
        *rs = INFINITY;
        return;
    }

    *cs = complex_make(num.real / denom, num.imag / denom);
    double s12s21_mag = complex_mag(s.m12) * complex_mag(s.m21);
    *rs = s12s21_mag / fabs(denom);
}

void sparams_load_stability_circle(matrix2x2_t s, complex_t *cl, double *rl) {
    complex_t delta = matrix2x2_det(s);
    complex_t s11_conj = complex_conj(s.m11);
    complex_t num = complex_conj(complex_sub(s.m22, complex_mul(delta, s11_conj)));
    double s22_mag_sq = complex_mag(s.m22) * complex_mag(s.m22);
    double delta_mag_sq = complex_mag(delta) * complex_mag(delta);
    double denom = s22_mag_sq - delta_mag_sq;

    if (fabs(denom) < 1e-30) {
        *cl = complex_make(INFINITY, INFINITY);
        *rl = INFINITY;
        return;
    }

    *cl = complex_make(num.real / denom, num.imag / denom);
    double s12s21_mag = complex_mag(s.m12) * complex_mag(s.m21);
    *rl = s12s21_mag / fabs(denom);
}

/**
 * Simultaneous conjugate match for unconditionally stable two-port.
 *
 * ΓMS = (B1 ± sqrt(B1² - 4|C1|²)) / (2*C1)
 * B1 = 1 + |s11|² - |s22|² - |Δ|²
 * C1 = s11 - Δ*conj(s22)
 *
 * Use minus sign for |ΓMS| < 1 (physically realizable).
 */
int sparams_conjugate_match(matrix2x2_t s, complex_t *gms, complex_t *gml) {
    double k = sparams_rollett_k(s);
    double delta_mag = sparams_delta_mag(s);
    if (k <= 1.0) return -1;  /* Not unconditionally stable */

    double s11_mag_sq = complex_mag(s.m11) * complex_mag(s.m11);
    double s22_mag_sq = complex_mag(s.m22) * complex_mag(s.m22);
    complex_t delta = matrix2x2_det(s);

    /* B1 and C1 for ΓMS */
    double b1 = 1.0 + s11_mag_sq - s22_mag_sq - delta_mag * delta_mag;
    complex_t c1 = complex_sub(s.m11, complex_mul(delta, complex_conj(s.m22)));

    /* B2 and C2 for ΓML */
    double b2 = 1.0 + s22_mag_sq - s11_mag_sq - delta_mag * delta_mag;
    complex_t c2 = complex_sub(s.m22, complex_mul(delta, complex_conj(s.m11)));

    double c1_mag = complex_mag(c1);
    double c2_mag = complex_mag(c2);

    if (c1_mag < 1e-30 || c2_mag < 1e-30) return -1;

    /* ΓMS */
    double sqrt_term1 = sqrt(b1 * b1 - 4.0 * c1_mag * c1_mag);
    if (b1 < 0) sqrt_term1 = -sqrt_term1;
    *gms = complex_make(
        (b1 - sqrt_term1) * c1.real / (2.0 * c1_mag * c1_mag),
        -(b1 - sqrt_term1) * c1.imag / (2.0 * c1_mag * c1_mag)
    );

    /* ΓML */
    double sqrt_term2 = sqrt(b2 * b2 - 4.0 * c2_mag * c2_mag);
    if (b2 < 0) sqrt_term2 = -sqrt_term2;
    *gml = complex_make(
        (b2 - sqrt_term2) * c2.real / (2.0 * c2_mag * c2_mag),
        -(b2 - sqrt_term2) * c2.imag / (2.0 * c2_mag * c2_mag)
    );

    return 0;
}

/**
 * Noise figure of a two-port.
 *
 * NF = NFmin + 4*Rn/Z0 * |ΓS - Γopt|² / ((1-|ΓS|²)*|1+Γopt|²)
 *
 * The noise parameters (NFmin, Rn, Γopt) characterize the noise
 * behavior of the device. They are measured or supplied by the
 * manufacturer at each frequency.
 *
 * For a low-noise amplifier: choose ΓS close to Γopt.
 * For a power amplifier: choose ΓS for maximum gain (ΓMS).
 * The design tradeoff: noise vs. gain matching.
 *
 * Course: Stanford EE359 — LNA design
 */
double sparams_noise_figure(double nfmin, double rn,
                             complex_t gopt, complex_t gs) {
    complex_t diff = complex_sub(gs, gopt);
    double num = 4.0 * rn * complex_mag(diff) * complex_mag(diff);
    complex_t one = complex_make(1.0, 0.0);
    double denom = (1.0 - complex_mag(gs) * complex_mag(gs)) *
                   complex_mag(complex_add(one, gopt)) *
                   complex_mag(complex_add(one, gopt));

    if (denom < 1e-30) return INFINITY;
    return nfmin + num / denom;
}

/* ============================================================================
 * L2: S-parameter to Impedance
 * ============================================================================ */

complex_t sparams_input_impedance(matrix2x2_t s, double z0, complex_t gl) {
    complex_t gin = sparams_gamma_in(s, gl);
    return impedance_from_gamma(gin, z0);
}

double sparams_transmission_phase_deg(matrix2x2_t s) {
    return complex_arg(s.m21) * 180.0 / M_PI;
}

/**
 * Group delay: τg = -Δφ/Δω.
 *
 * Computed as the negative phase slope between two closely-spaced
 * frequency samples. Group delay measures the time delay experienced
 * by the signal envelope.
 *
 * Constant group delay = linear phase = no dispersion.
 *
 * Applications:
 *   - Filter characterization (passband group delay flatness)
 *   - Amplifier linearity (AM-PM distortion)
 *   - Antenna design (phase center)
 */
double sparams_group_delay(matrix2x2_t s1, matrix2x2_t s2,
                            double omega1, double omega2) {
    double phi1 = complex_arg(s1.m21);
    double phi2 = complex_arg(s2.m21);
    double omega_step = omega2 - omega1;
    if (fabs(omega_step) < 1e-30) return 0.0;
    return -(phi2 - phi1) / omega_step;
}
