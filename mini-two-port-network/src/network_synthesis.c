/**
 * @file network_synthesis.c
 * @brief Network Synthesis from Two-Port Parameters
 *
 * Implements network synthesis methods: T/π equivalent extraction,
 * T↔π conversion, Cauer ladder synthesis, L/π/T matching network
 * design, Foster canonical forms, and Darlington synthesis.
 */

#include "../include/network_synthesis.h"
#include <float.h>
#include <stdlib.h>

/* ============================================================================
 * L5: Equivalent Circuit Synthesis
 * ============================================================================ */

int synthesize_t_network(matrix2x2_t z, complex_t *za, complex_t *zb,
                         complex_t *zc) {
    /* Check reciprocity: z12 ≈ z21 */
    complex_t diff = complex_sub(z.m12, z.m21);
    if (!complex_approx_equal(diff, complex_make(0, 0), 1e-9)) {
        return -1;
    }

    *za = complex_sub(z.m11, z.m12);
    *zb = complex_sub(z.m22, z.m12);
    *zc = z.m12;
    return 0;
}

int synthesize_pi_network(matrix2x2_t y, complex_t *ya, complex_t *yb,
                          complex_t *yc) {
    /* Check reciprocity: y12 ≈ y21 */
    complex_t diff = complex_sub(y.m12, y.m21);
    if (!complex_approx_equal(diff, complex_make(0, 0), 1e-9)) {
        return -1;
    }

    /* ya = y11 + y12 (y12 is typically negative in π representation) */
    *ya = complex_add(y.m11, y.m12);
    /* yb = -y12 (series admittance) */
    *yb = complex_make(-y.m12.real, -y.m12.imag);
    /* yc = y22 + y12 */
    *yc = complex_add(y.m22, y.m12);
    return 0;
}

/**
 * T ↔ π conversion (generalized Kennelly Δ-Y transform).
 *
 * T → π (impedances to admittances):
 *   The T-network with impedances Za, Zb, Zc is equivalent to
 *   a π-network with admittances Ya, Yb, Yc where:
 *
 *   Ya = Zb / Δ, Yb = Zc / Δ, Yc = Za / Δ
 *   where Δ = Za*Zb + Zb*Zc + Zc*Za
 *
 * π → T (admittances to impedances):
 *   Za = Yc / Δ, Zb = Ya / Δ, Zc = Yb / Δ
 *   where Δ = Ya*Yb + Yb*Yc + Yc*Ya
 *
 * This is the complex-impedance generalization of the resistive
 * Δ-Y (or Π-T) transform widely used in three-phase power systems
 * and resistor network simplification.
 *
 * Reference: Kennelly, Electrical World, 1899
 * Course: Berkeley EE16B — Δ-Y transform
 */
void synthesize_t_pi_convert(complex_t za_in, complex_t zb_in,
                              complex_t zc_in, int to_pi,
                              complex_t *elem1_out, complex_t *elem2_out,
                              complex_t *elem3_out) {
    complex_t delta = complex_add(
        complex_mul(za_in, zb_in),
        complex_add(
            complex_mul(zb_in, zc_in),
            complex_mul(zc_in, za_in)
        )
    );

    if (complex_is_zero(delta, 1e-30)) {
        *elem1_out = complex_make(NAN, NAN);
        *elem2_out = complex_make(NAN, NAN);
        *elem3_out = complex_make(NAN, NAN);
        return;
    }

    complex_t inv_delta = complex_div(complex_make(1.0, 0.0), delta);

    if (to_pi) {
        /* T → π: Ya = Zb/Δ, Yb = Zc/Δ, Yc = Za/Δ */
        *elem1_out = complex_mul(zb_in, inv_delta);
        *elem2_out = complex_mul(zc_in, inv_delta);
        *elem3_out = complex_mul(za_in, inv_delta);
    } else {
        /* π → T: Za = Yc/Δ, Zb = Ya/Δ, Zc = Yb/Δ */
        *elem1_out = complex_mul(zc_in, inv_delta);
        *elem2_out = complex_mul(za_in, inv_delta);
        *elem3_out = complex_mul(zb_in, inv_delta);
    }
}

/**
 * Cauer ladder synthesis — one iteration.
 *
 * This function determines one series-L or shunt-C extraction
 * from a given impedance at a single frequency.
 *
 * For a full synthesis, this would be called repeatedly with
 * decreasing network complexity.
 *
 * At frequency ω:
 *   If Im(Zin) > 0: extract series L = Im(Zin)/ω
 *     Z_remaining = Zin - jωL = Re(Zin)
 *   If Im(Zin) < 0: extract shunt C = -1/(ω*Im(Zin))
 *     Y_remaining = 1/Zin - jωC (admittance subtraction)
 *
 * This process corresponds to continued fraction expansion.
 *
 * Reference: Cauer, "Die Verwirklichung von Wechselstromwiderständen
 *   vorgeschriebener Frequenzabhängigkeit", 1926
 * Course: Illinois ECE 451 — Cauer synthesis
 */
int synthesize_cauer_ladder(complex_t zin, double omega, int n_stages,
                             double *series_l, double *shunt_c) {
    complex_t z_rem = zin;

    for (int i = 0; i < n_stages; i++) {
        if (complex_is_zero(z_rem, 1e-30)) {
            /* Remaining impedance is zero → short to ground */
            series_l[i] = 0.0;
            shunt_c[i] = 0.0;
            continue;
        }

        if (z_rem.imag > 0) {
            /* Inductive → extract series L */
            series_l[i] = z_rem.imag / omega;
            shunt_c[i] = 0.0;
            z_rem = complex_make(z_rem.real, 0.0);
        } else if (z_rem.imag < 0) {
            /* Capacitive → extract shunt C */
            series_l[i] = 0.0;
            double y_mag = complex_mag(z_rem);
            if (y_mag < 1e-30) return -1;
            shunt_c[i] = -1.0 / (omega * z_rem.imag);
            /* Convert to admittance, subtract jωC, convert back */
            complex_t y_rem = complex_div(complex_make(1.0, 0.0), z_rem);
            y_rem = complex_sub(y_rem, complex_make(0.0, omega * shunt_c[i]));
            if (complex_is_zero(y_rem, 1e-30)) {
                z_rem = complex_make(INFINITY, 0.0);
            } else {
                z_rem = complex_div(complex_make(1.0, 0.0), y_rem);
            }
        } else {
            /* Purely resistive → done, remaining stages are zero */
            series_l[i] = 0.0;
            shunt_c[i] = 0.0;
        }
    }
    return 0;
}

/**
 * L-network matching design (single-frequency).
 *
 * Goal: Match ZS (source) to ZL (load) at frequency ω.
 *
 * Algorithm:
 * 1. Normalize impedances
 * 2. Determine if ZL is inside or outside the 1+jx circle
 * 3. Choose series-shunt or shunt-series topology
 * 4. Compute the required reactance/susceptance
 *
 * The L-network is the simplest matching network (2 elements).
 * It provides a perfect match at a single frequency with a specific
 * bandwidth determined by the impedance ratio.
 *
 * Course: Stanford EE359 — L-network matching
 */
int synthesize_l_match(complex_t zs, complex_t zl, double omega,
                        double *l_value, double *c_value) {
    /* For simplicity: assume ZS = R0 + j0 (resistive source) and
       ZL is complex. Design L-network to absorb ZL reactance. */

    double r0 = zs.real;
    if (r0 < 1e-30) return -1;

    double rl = zl.real;
    double xl = zl.imag;

    /* Case: RL > R0 (step down impedance) */
    if (rl > r0) {
        /* Shunt C at load, then series L */
        /* Q = sqrt(RL/R0 - 1) */
        double q = sqrt(rl / r0 - 1.0);
        if (q < 0) return -1;

        double xc_shunt = rl / q;
        double c_val = 1.0 / (omega * xc_shunt);

        /* The shunt C absorbs XL partially */
        /* Series L adds remaining reactance */
        double xl_series = q * r0 - xl;
        if (xl_series < 0) {
            /* Need series C instead */
            *l_value = 0.0;
            *c_value = -1.0 / (omega * xl_series);
        } else {
            *l_value = xl_series / omega;
            *c_value = c_val;
        }
        return 0;
    } else {
        /* RL < R0 (step up impedance) */
        /* Series L, then shunt C */
        double q = sqrt(r0 / rl - 1.0);
        if (q < 0) return -1;

        double xl_series = q * rl;
        /* Absorb XL into series L */
        xl_series -= xl;
        if (xl_series < 0) {
            *l_value = 0.0;
            *c_value = 0.0;
            return -1;
        }
        *l_value = xl_series / omega;

        double xc_shunt = r0 / q;
        *c_value = 1.0 / (omega * xc_shunt);
        return 0;
    }
}

/**
 * π-network matching design.
 *
 * A π-network provides an extra degree of freedom: the loaded Q
 * can be chosen independently of the impedance transformation ratio.
 *
 * Higher Q → narrower bandwidth, better harmonic rejection.
 * Lower Q → wider bandwidth, lower component sensitivity.
 *
 * Design steps:
 * 1. Choose Q (typically 3-20)
 * 2. Design virtual resistance R_v = R_high / (Q²+1)
 * 3. Design L-network from R_low to R_v
 * 4. Design L-network from R_v to R_high
 * 5. Combine into π topology
 */
int synthesize_pi_match(complex_t zs, complex_t zl, double omega,
                         double q_loaded,
                         complex_t *ya, complex_t *yb, complex_t *yc) {
    (void)omega;  /* Frequency is used by caller to convert susceptance to L/C values */
    double r1 = zs.real;
    double r2 = zl.real;

    if (r1 < 1e-30 || r2 < 1e-30) return -1;
    if (q_loaded < 0.5) return -1;

    /* Virtual resistance: the intermediate impedance */
    double r_high = (r1 > r2) ? r1 : r2;
    double r_low = (r1 < r2) ? r1 : r2;

    /* R_virtual = R_high / (Q² + 1) */
    double r_virtual = r_high / (q_loaded * q_loaded + 1.0);

    if (r_virtual < r_low) {
        /* Need to increase Q */
        r_virtual = r_low;
    }

    /* Compute the shunt reactances at each end */
    double q1 = sqrt(r1 / r_virtual - 1.0);
    double q2 = sqrt(r2 / r_virtual - 1.0);

    if (q1 < 0 || q2 < 0) return -1;

    /* Input shunt: jB1 = ±j*Q1/R1 */
    /* Output shunt: jB2 = ±j*Q2/R2 */
    double b1 = q1 / r1;
    double b2 = q2 / r2;

    /* Series reactance: jX = j*(Q1+Q2)*R_virtual */
    double x_series = (q1 + q2) * r_virtual;

    /* Determine signs: typically use capacitive shunts and inductive series
       for low-pass π-network */
    *ya = complex_make(0.0, -b1);  /* Capacitive shunt (negative susceptance) */
    *yc = complex_make(0.0, -b2);  /* Wait — for capacitor, Y = jωC, so imag positive */
    /* Let's use the convention: capacitive shunt = positive susceptance */
    *ya = complex_make(0.0, b1);
    *yc = complex_make(0.0, b2);
    *yb = complex_make(0.0, -1.0 / x_series);  /* Inductive series = negative susceptance */

    return 0;
}

/**
 * T-network matching design (dual of π-network).
 *
 * Same approach but working in the impedance domain.
 */
int synthesize_t_match(complex_t zs, complex_t zl, double omega,
                        double q_loaded,
                        complex_t *za, complex_t *zb, complex_t *zc) {
    /* Use admittance-domain dual of π-network */
    /* Convert to admittances, design π, convert back to T */
    complex_t ys, yl;
    ys = complex_div(complex_make(1.0, 0.0), zs);
    yl = complex_div(complex_make(1.0, 0.0), zl);

    complex_t ya, yb, yc;
    int ret = synthesize_pi_match(ys, yl, omega, q_loaded, &ya, &yb, &yc);
    if (ret != 0) return ret;

    /* Convert π to T */
    synthesize_t_pi_convert(ya, yb, yc, 0, za, zb, zc);
    return 0;
}

/* ============================================================================
 * L4: Foster's Reactance Theorem
 * ============================================================================ */

int foster_reactance_check(complex_t z1, complex_t z2) {
    /* Check purely reactive: Re(Z) ≈ 0 */
    if (fabs(z1.real) > 1e-9 || fabs(z2.real) > 1e-9) {
        return -1;  /* Not purely reactive */
    }

    /* Check monotonic increase: X(ω2) > X(ω1) for ω2 > ω1 */
    if (z2.imag > z1.imag) {
        return 1;  /* Foster condition satisfied */
    }
    return 0;  /* Foster condition violated */
}

void foster_pole_zero(double l, double c, int is_parallel,
                       double *pole_freq, double *zero_freq) {
    double lc = l * c;
    if (lc < 1e-30) {
        *pole_freq = INFINITY;
        *zero_freq = 0.0;
        return;
    }
    double omega_res = 1.0 / sqrt(lc);

    if (is_parallel) {
        /* Parallel LC: pole at resonance (impedance → ∞) */
        *pole_freq = omega_res;
        *zero_freq = 0.0;  /* Zero at DC (L shorts) */
    } else {
        /* Series LC: zero at resonance (impedance → 0) */
        *zero_freq = omega_res;
        *pole_freq = 0.0;  /* Pole at DC (C opens) */
    }
}

/**
 * Foster-I canonical form evaluation at frequency ω.
 *
 * Z(s) = Σ Ki*s/(s²+ωi²) + K∞*s + K0/s
 *
 * At s = jω:
 *   Z(jω) = Σ Ki*jω/(ωi²-ω²) + jω*K∞ + K0/(jω)
 *         = j * [ Σ Ki*ω/(ωi²-ω²) + ω*K∞ - K0/ω ]
 */
void synthesize_foster_one(const double *residues, const double *pole_freqs,
                            int n_resonators, double k_inf, double k0,
                            complex_t *zin_out, double omega) {
    double x_total = 0.0;

    /* Sum over parallel LC resonators in series */
    for (int i = 0; i < n_resonators; i++) {
        double wi = pole_freqs[i];
        double denom = wi * wi - omega * omega;
        if (fabs(denom) < 1e-30) {
            /* At resonance → pole (impedance → ∞) */
            *zin_out = complex_make(INFINITY, 0.0);
            return;
        }
        x_total += residues[i] * omega / denom;
    }

    /* K∞ * ω: series inductor contribution */
    x_total += k_inf * omega;

    /* K0/ω: series capacitor contribution (negative sign from 1/(jωC)) */
    if (omega > 1e-30) {
        x_total -= k0 / omega;
    }

    *zin_out = complex_make(0.0, x_total);
}

/**
 * Foster-II canonical form (dual): admittance synthesis.
 *
 * Y(s) = Σ Ki*s/(s²+ωi²) + K∞*s + K0/s
 *
 * Evaluated similarly but in the admittance domain.
 */
void synthesize_foster_two(const double *residues, const double *pole_freqs,
                            int n_resonators, double k_inf, double k0,
                            complex_t *yin_out, double omega) {
    double b_total = 0.0;

    for (int i = 0; i < n_resonators; i++) {
        double wi = pole_freqs[i];
        double denom = wi * wi - omega * omega;
        if (fabs(denom) < 1e-30) {
            *yin_out = complex_make(INFINITY, 0.0);
            return;
        }
        b_total += residues[i] * omega / denom;
    }

    b_total += k_inf * omega;

    if (omega > 1e-30) {
        b_total -= k0 / omega;
    }

    *yin_out = complex_make(0.0, b_total);
}

/* ============================================================================
 * L6: Darlington Synthesis
 * ============================================================================ */

/**
 * Darlington synthesis step: extract one reactive element.
 *
 * For type = 0 (series L extraction):
 *   Z_remaining(s) = Z(s) - sL
 *   Choose L such that Z_remaining has a zero on the jω axis.
 *
 * For type = 1 (shunt C extraction):
 *   Y_remaining(s) = Y(s) - sC
 *   Similarly.
 *
 * This is the fundamental operation in Darlington's method for
 * synthesizing a lossless two-port terminated in a resistor.
 *
 * Reference: Darlington, J. Math. Phys., 1939
 */
int synthesize_darlington_step(complex_t zin, double omega,
                                int extracted_type,
                                double *element_value,
                                complex_t *z_remaining) {
    if (extracted_type == 0) {
        /* Extract series L: Z_rem = Z_in - jωL */
        if (zin.imag <= 0) {
            /* Cannot extract L from capacitive impedance */
            return -1;
        }
        *element_value = zin.imag / omega;
        *z_remaining = complex_make(zin.real, 0.0);
    } else {
        /* Extract shunt C: Y_rem = Y_in - jωC */
        complex_t yin = complex_div(complex_make(1.0, 0.0), zin);
        if (yin.imag >= 0) {
            /* Cannot extract C from inductive admittance */
            return -1;
        }
        *element_value = -yin.imag / omega;
        complex_t y_rem = complex_make(yin.real, 0.0);
        if (complex_is_zero(y_rem, 1e-30)) {
            *z_remaining = complex_make(INFINITY, 0.0);
        } else {
            *z_remaining = complex_div(complex_make(1.0, 0.0), y_rem);
        }
    }
    return 0;
}

/**
 * Butterworth termination impedance.
 *
 * For an N-th order Butterworth filter prototype terminated in 1Ω:
 * Zin(s) is a positive-real function whose magnitude squared at s=jω
 * gives the Butterworth response.
 *
 * For N=1: Zin(s) = (s+1)/(s+1) = 1 (trivial)
 * For N=2: Zin(s) = (s²+√2s+1)/(s²+√2s+1) for matched case
 *
 * More precisely, for the doubly-terminated case, the input reflection
 * coefficient is designed to be:
 *   |s11(jω)|² = ω^(2N) / (1 + ω^(2N))
 *
 * This function computes Zin(jω) for a given frequency.
 */
complex_t synthesize_butterworth_z(double omega, int order) {
    /* For a Butterworth filter with N sections, the input impedance
       at frequency ω can be computed from the reflection coefficient. */

    /* |s11|² = ω^(2N) / (1 + ω^(2N)) */
    double omega_n = pow(fabs(omega), order);
    double s11_mag_sq = (omega_n * omega_n) / (1.0 + omega_n * omega_n);
    double s11_mag = sqrt(s11_mag_sq);

    /* Phase: For the doubly-terminated Butterworth, arg(s11) = N*π/2 for ω > 1 */
    double s11_phase;
    if (omega > 1.0) {
        s11_phase = M_PI * order / 2.0;
    } else {
        s11_phase = 0.0;
    }

    complex_t s11 = complex_make(s11_mag * cos(s11_phase), s11_mag * sin(s11_phase));

    /* Zin = Z0 * (1 + s11) / (1 - s11) */
    complex_t one = complex_make(1.0, 0.0);
    complex_t num = complex_add(one, s11);
    complex_t denom = complex_sub(one, s11);

    return complex_div(num, denom);  /* Normalized to Z0 = 1 */
}
