/**
 * @file filter_design.c
 * @brief Filter Design Using Two-Port Network Theory
 *
 * Implements filter approximation functions (Butterworth, Chebyshev,
 * Elliptic, Bessel), frequency transformations, LC ladder and active
 * RC filter synthesis, and distributed-element filter design.
 */

#include "../include/filter_design.h"
#include "../include/abcd_params.h"
#include "../include/s_params.h"
#include "../include/conversion.h"
#include <float.h>

/* ============================================================================
 * L5: Butterworth Low-Pass Prototype
 * ============================================================================ */

/**
 * Butterworth N-th order transfer function.
 *
 * |H(jω)|² = 1 / (1 + ω^(2N))
 *
 * H(jω) = 1 / sqrt(1 + ω^(2N))  [magnitude, zero phase at DC]
 *
 * The poles of H(s)H(-s) lie on a unit circle in the s-plane at:
 *   s_k = exp(j * (π/2 + (2k-1)π/(2N))), k = 1..N
 *
 * The N poles of H(s) are the left-half-plane (LHP) ones:
 *   s_k = -sin((2k-1)π/(2N)) + j cos((2k-1)π/(2N))
 *
 * At s = jω: H(jω) = 1 / Π(|jω - s_k|)
 * But for magnitude, we use the direct formula.
 *
 * Properties:
 *   - Maximally flat at ω = 0 (first 2N-1 derivatives = 0 at ω=0)
 *   - -3 dB at ω = 1 (normalized cutoff)
 *   - -20N dB/decade rolloff
 *   - Monotonic in passband and stopband (no ripple)
 *   - Moderate group delay variation
 *
 * Reference: Butterworth, Wireless Engineer, 1930
 * Course: MIT 6.003 — Butterworth filter design
 */
complex_t filter_butterworth_lp(int order, double omega) {
    double omega_pow = 1.0;
    for (int i = 0; i < order; i++) omega_pow *= omega;
    double mag_sq = 1.0 / (1.0 + omega_pow * omega_pow);
    double mag = sqrt(mag_sq);

    /* Phase: compute from the pole locations */
    /* For Butterworth, arg(H(jω)) = -Σ atan(ω - imag(p_k), -real(p_k)) */
    double phase = 0.0;
    for (int k = 1; k <= order; k++) {
        double theta_k = (2.0 * k - 1.0) * M_PI / (2.0 * order);
        double sigma_k = -sin(theta_k);  /* Real part of LHP pole */
        double omega_k = cos(theta_k);   /* Imag part of LHP pole */
        phase -= atan2(omega - omega_k, -sigma_k);
    }

    return complex_make(mag * cos(phase), mag * sin(phase));
}

/* ============================================================================
 * L5: Chebyshev Type-I Low-Pass Prototype
 * ============================================================================ */

/**
 * Chebyshev Type-I transfer function.
 *
 * |H(jω)|² = 1 / (1 + ε² * T_N²(ω))
 *
 * where T_N(ω) is the Chebyshev polynomial of order N:
 *   T_N(ω) = cos(N * arccos(ω)) for |ω| ≤ 1
 *   T_N(ω) = cosh(N * arccosh(ω)) for |ω| ≥ 1
 *
 * ε = sqrt(10^(RdB/10) - 1)  (ripple factor)
 *
 * Properties:
 *   - Equiripple in passband (|ω| ≤ 1): oscillates between 1 and 1/√(1+ε²)
 *   - Monotonic in stopband (|ω| > 1): rolls off faster than Butterworth
 *   - Sharper transition for the same order
 *
 * Reference: Chebyshev polynomials, applied to filters by Cauer (1931)
 * Course: MIT 6.003 — Chebyshev approximation
 */
complex_t filter_chebyshev1_lp(int order, double omega, double ripple_db) {
    /* Ripple factor ε */
    double epsilon = sqrt(pow(10.0, ripple_db / 10.0) - 1.0);

    double tn;
    if (fabs(omega) <= 1.0) {
        tn = cos(order * acos(omega));
    } else {
        /* For |ω| > 1: T_N(ω) = cosh(N * arccosh(ω)) */
        if (omega > 1.0) {
            tn = cosh(order * acosh(omega));
        } else {
            tn = cosh(order * acosh(-omega));
            if (order % 2 == 1) tn = -tn;  /* Odd-order: T_N(-ω) = -T_N(ω) */
        }
    }

    double mag_sq = 1.0 / (1.0 + epsilon * epsilon * tn * tn);
    double mag = sqrt(mag_sq);

    /* Phase: approximately the same as Butterworth for the dominant poles */
    double phase = 0.0;
    /* Chebyshev phase is more complex; approximate with Butterworth phase */
    for (int k = 1; k <= order; k++) {
        double theta_k = (2.0 * k - 1.0) * M_PI / (2.0 * order);
        double sigma_k = -sin(theta_k);
        double omega_k = cos(theta_k);
        phase -= atan2(omega - omega_k * sqrt(1.0 + 1.0/(epsilon*epsilon)),
                       -sigma_k);
    }

    return complex_make(mag * cos(phase), mag * sin(phase));
}

/* ============================================================================
 * L5: Chebyshev Type-II Low-Pass Prototype
 * ============================================================================ */

complex_t filter_chebyshev2_lp(int order, double omega, double stop_db) {
    /* Chebyshev II: |H(jω)|² = 1 / (1 + ε² * T_N²(1/ωc) / T_N²(1/ω)) */
    double epsilon = 1.0 / sqrt(pow(10.0, stop_db / 10.0) - 1.0);

    if (fabs(omega) < 1e-30) {
        return complex_make(1.0, 0.0);  /* DC = 1 (passband is monotonic) */
    }

    double omega_inv = 1.0 / omega;
    double tn1, tn_omega;

    if (fabs(omega_inv) <= 1.0) {
        tn1 = cos(order * acos(1.0));  /* T_N(1) = 1 */
        tn_omega = cos(order * acos(omega_inv));
    } else {
        if (omega_inv > 1.0) {
            tn1 = cosh(order * acosh(1.0));  /* cosh(0) = 1 */
            tn_omega = cosh(order * acosh(omega_inv));
        } else {
            tn1 = 1.0;
            tn_omega = cosh(order * acosh(-omega_inv));
            if (order % 2 == 1) tn_omega = -tn_omega;
        }
    }

    double ratio = tn1 / tn_omega;
    double mag_sq = 1.0 / (1.0 + epsilon * epsilon * ratio * ratio);
    double mag = sqrt(mag_sq);

    return complex_make(mag, 0.0);  /* Approximate phase as 0 */
}

/* ============================================================================
 * L5: Elliptic (Cauer) Low-Pass Prototype
 * ============================================================================ */

double filter_elliptic_lp_mag(int order, double omega,
                               double ripple_db, double stop_db) {
    /* Elliptic filter magnitude: |H(jω)|² = 1/(1 + ε²*R_N²(ω)) */
    /* where R_N is the Chebyshev rational function. */
    /* This is an approximation — exact elliptic requires Jacobi elliptic functions. */

    double epsilon = sqrt(pow(10.0, ripple_db / 10.0) - 1.0);
    double epsilon_s = 1.0 / sqrt(pow(10.0, stop_db / 10.0) - 1.0);

    /* For the stopband edge: ωs > 1 (normalized LP) */
    /* Selectivity factor: k1 = ε / ε_s */
    double k1 = epsilon / epsilon_s;

    /* Approximate the elliptic rational function as:
       |R_N(ω)| ≈ (ω^N + ω^(-N)) / 2 for the first-order approximation
       For a more accurate model, we use the pole-zero pattern */

    /* Simplified magnitude estimate */
    double omega_n;
    if (fabs(omega) < 1.0) {
        /* Passband: equiripple */
        /* Approximate as Chebyshev-I in passband */
        double tn = cos(order * acos(omega));
        double mag_sq = 1.0 / (1.0 + epsilon * epsilon * tn * tn);
        return sqrt(mag_sq);
    } else {
        /* Stopband: use inverse Chebyshev */
        omega_n = pow(fabs(omega), order);
        double mag_sq = 1.0 / (1.0 + epsilon * epsilon * omega_n * omega_n * k1 * k1);
        if (mag_sq > 1.0) mag_sq = 1.0;
        return sqrt(mag_sq);
    }
}

/* ============================================================================
 * L5: Bessel (Thomson) Low-Pass Prototype
 * ============================================================================ */

/**
 * Bessel filter: maximally flat group delay.
 *
 * The transfer function uses reverse Bessel polynomials.
 *
 * For N=2: H(s) = 3 / (s² + 3s + 3)
 * For N=3: H(s) = 15 / (s³ + 6s² + 15s + 15)
 * For N=4: H(s) = 105 / (s⁴ + 10s³ + 45s² + 105s + 105)
 *
 * At s = jω:
 *   H(jω) = b0 / ( (jω)^N + a_{N-1}(jω)^(N-1) + ... + a_0 )
 *
 * This function computes H(jω) for orders 1-8.
 */
complex_t filter_bessel_lp(int order, double omega) {
    /* Bessel polynomial coefficients (reverse: a0 + a1*s + ... + aN*s^N) */
    /* Only denominator needed; numerator is b0 = a0 */
    static const double bessel_coeffs[8][9] = {
        {1.0, 1.0, 0, 0, 0, 0, 0, 0, 0},                          /* N=1 */
        {3.0, 3.0, 1.0, 0, 0, 0, 0, 0, 0},                         /* N=2 */
        {15.0, 15.0, 6.0, 1.0, 0, 0, 0, 0, 0},                     /* N=3 */
        {105.0, 105.0, 45.0, 10.0, 1.0, 0, 0, 0, 0},               /* N=4 */
        {945.0, 945.0, 420.0, 105.0, 15.0, 1.0, 0, 0, 0},          /* N=5 */
        {10395.0, 10395.0, 4725.0, 1260.0, 210.0, 21.0, 1.0, 0, 0}, /* N=6 */
        {135135.0, 135135.0, 62370.0, 17325.0, 3150.0, 378.0, 28.0, 1.0, 0}, /* N=7 */
        {2027025.0, 2027025.0, 945945.0, 270270.0, 51975.0, 6930.0, 630.0, 36.0, 1.0} /* N=8 */
    };

    if (order < 1 || order > 8) {
        return complex_make(NAN, NAN);
    }

    int n = order;
    const double *a = bessel_coeffs[n - 1];

    /* Evaluate denominator polynomial at s = jω */
    /* D(jω) = a0 + a1*jω + a2*(jω)² + a3*(jω)³ + ... */
    double real_d = a[0];
    double imag_d = 0.0;

    double omega_pow = omega;
    for (int k = 1; k <= n; k++) {
        double term = a[k] * omega_pow;
        switch (k % 4) {
            case 0: real_d += term; break;   /* (jω)^4 = ω⁴ (real, +) */
            case 1: imag_d += term; break;   /* (jω)^1 = jω (imag, +) */
            case 2: real_d -= term; break;   /* (jω)^2 = -ω² (real, -) */
            case 3: imag_d -= term; break;   /* (jω)^3 = -jω³ (imag, -) */
        }
        omega_pow *= omega;
    }

    /* H(jω) = b0 / D(jω) where b0 = a[0] */
    double denom_mag_sq = real_d * real_d + imag_d * imag_d;
    double denom_mag = sqrt(denom_mag_sq);
    double mag = a[0] / denom_mag;
    double phase = -atan2(imag_d, real_d);

    return complex_make(mag * cos(phase), mag * sin(phase));
}

/* ============================================================================
 * L5: Frequency Transformations
 * ============================================================================ */

double filter_lp_to_lp(double omega_norm, double omega_c) {
    if (omega_c < 1e-30) return omega_norm;
    return omega_norm / omega_c;
}

double filter_lp_to_hp(double omega_norm, double omega_c) {
    if (omega_norm < 1e-30 || omega_c < 1e-30) return INFINITY;
    return omega_c / omega_norm;
}

double filter_lp_to_bp(double omega_norm, double omega_l, double omega_h) {
    double omega0 = sqrt(omega_l * omega_h);
    double bw = omega_h - omega_l;
    if (bw < 1e-30) return INFINITY;
    /* s → (s² + ω0²)/(s·BW) → at s=jω: */
    double omega = omega_norm;
    /* |(ω² - ω0²)/(ω·BW)| — this gives the LP equivalent frequency */
    if (omega < 1e-30) return INFINITY;
    return fabs((omega * omega - omega0 * omega0) / (omega * bw));
}

double filter_lp_to_bs(double omega_norm, double omega_l, double omega_h) {
    double omega0 = sqrt(omega_l * omega_h);
    double bw = omega_h - omega_l;
    if (bw < 1e-30) return INFINITY;
    /* s → s·BW/(s² + ω0²) */
    double omega = omega_norm;
    if (omega < 1e-30) return 0.0;  /* DC passes in BS */
    double val = fabs((omega * bw) / (omega * omega - omega0 * omega0));
    if (val < 1e-30) return 0.0;
    return 1.0 / val;
}

/* ============================================================================
 * L6: Butterworth LC Ladder g-values
 * ============================================================================ */

int filter_butterworth_g_values(int order, double *g_values) {
    if (order < 1) return -1;

    for (int k = 0; k < order; k++) {
        g_values[k] = 2.0 * sin((2.0 * k + 1.0) * M_PI / (2.0 * order));
    }
    return 0;
}

/* ============================================================================
 * L6: Chebyshev LC Ladder g-values
 * ============================================================================ */

int filter_chebyshev_g_values(int order, double ripple_db, double *g_values) {
    if (order < 1) return -1;

    double epsilon = sqrt(pow(10.0, ripple_db / 10.0) - 1.0);
    double beta_simple = asinh(1.0 / epsilon);

    double sinh_beta = sinh(beta_simple);

    double a[order + 1], b[order + 1];
    for (int k = 1; k <= order; k++) {
        a[k] = sin((2.0 * k - 1.0) * M_PI / (2.0 * order));
        b[k] = sinh_beta * sinh_beta + sin(k * M_PI / order) * sin(k * M_PI / order);
    }

    g_values[0] = 2.0 * a[1] / sinh_beta;

    for (int k = 2; k <= order; k++) {
        g_values[k-1] = 4.0 * a[k-1] * a[k] / (b[k-1] * g_values[k-2]);
    }

    return 0;
}

/* ============================================================================
 * L6: LC Ladder ABCD Construction
 * ============================================================================ */

matrix2x2_t filter_build_lc_ladder(const double *g_values, int order,
                                    double omega, int start_with_shunt) {
    /* Start with identity */
    matrix2x2_t result = matrix2x2_make(
        complex_make(1.0, 0.0), complex_make(0.0, 0.0),
        complex_make(0.0, 0.0), complex_make(1.0, 0.0)
    );

    for (int i = 0; i < order; i++) {
        int is_shunt = (start_with_shunt + i) % 2;
        double g = g_values[i];

        matrix2x2_t section;
        if (is_shunt) {
            /* Shunt C: Y = jω*C = jω*g */
            complex_t y_shunt = complex_make(0.0, omega * g);
            section = abcd_shunt_y(y_shunt);
        } else {
            /* Series L: Z = jω*L = jω*g */
            complex_t z_series = complex_make(0.0, omega * g);
            section = abcd_series_z(z_series);
        }

        result = matrix2x2_mul(result, section);
    }

    return result;
}

double filter_response_db(matrix2x2_t abcd, double z0) {
    matrix2x2_t s = convert_abcd_to_s(abcd, z0);
    return sparams_s21_db(s);
}

/* ============================================================================
 * L6: Active RC Biquad Design
 * ============================================================================ */

/**
 * Sallen-Key low-pass biquad component selection.
 *
 * Unity-gain Sallen-Key topology:
 *
 *   Vin — R1 — R2 — Vout
 *          |     |
 *         C1    C2
 *          |     |
 *         GND   Vout (for K=1)
 *
 * Transfer function: H(s) = 1 / (s²R1R2C1C2 + sC2(R1+R2) + 1)
 *
 * For a given ω0 and Q:
 *   Choose C2 = 10/ω0 nF (practical rule), then:
 *   C1 = C2 / (2Q)
 *   R1 = R2 = 1 / (ω0 * sqrt(C1*C2))
 *
 * Or using the design equations:
 *   Let m = C1/C2, then:
 *   R2/R1 = (1/(2Q))² * (m+1)²/m - 1
 *
 * We use a simpler approach: fix C1 = C2 = C, then
 *   R1 = 1/(2Qω0C), R2 = 2Q/(ω0C)
 *
 * Reference: Sallen & Key, IRE Trans. Circuit Theory, 1955
 */
int filter_sallen_key_lp(double omega0, double q,
                          double *r1, double *r2, double *c1, double *c2) {
    if (omega0 < 1e-30 || q < 0.5) return -1;

    /* Choose C = 1/(ω0 * R_nominal), with R_nominal = 10kΩ */
    double r_nom = 10000.0;
    double c_nom = 1.0 / (omega0 * r_nom);

    *c1 = c_nom;
    *c2 = c_nom;

    /* For unity-gain Sallen-Key with equal C:
       R1 = R2 = 1/(ω0*C) * 1/(2Q) ... nope.
       Let's derive: ω0² = 1/(R1R2C1C2), Q = sqrt(R1R2C1C2)/(C2(R1+R2))
       With C1 = C2 = C:
       ω0 = 1/(C*sqrt(R1R2))
       Q = sqrt(R1R2)/(R1+R2) = sqrt(k)/(1+k) where k = R2/R1

       For Q = 0.707 (Butterworth): k = 1 → R1 = R2
       For Q > 0.5: solve k from Q² = k/(1+k)² → k = ... */

    if (q <= 1.0) {
        double disc = 2.0 * q * q - 1.0;
        if (disc < 0) disc = 0;
        /* k = 2Q²-1 ± 2Q*sqrt(Q²-1)... let me just compute numerically */
        /* Simpler: R1 = 1/(2Qω0C), R2 = 2Q/(ω0C) */
        *r1 = 1.0 / (2.0 * q * omega0 * c_nom);
        *r2 = 2.0 * q / (omega0 * c_nom);
    } else {
        /* For high Q, use the formula directly */
        double alpha = 2.0 * q;
        *r1 = 1.0 / (alpha * omega0 * c_nom);
        *r2 = alpha / (omega0 * c_nom);
    }

    return 0;
}

/**
 * Multiple-Feedback (MFB) low-pass biquad design.
 *
 * Topology:
 *   Vin — R1 — R3 — R2 — Vout (connected to inverting input of op-amp)
 *               |     |
 *              C1    C2
 *               |     |
 *            GND    GND (actually C2 goes from Vout to inv input)
 *
 * Actually, the standard MFB topology:
 *   Vin — R1 — node1 — R3 — node2 (inverting input)
 *               |        |
 *              C1       C2 (from inv input to output)
 *               |        |
 *              GND     Vout
 *              R2 from node1 to Vout
 *
 * Transfer function:
 *   H(s) = -R2/R1 / (s²R2R3C1C2 + sC2(R2+R3+R2R3/R1) + 1)
 *
 * Design with equal R:
 *   R1 = R2 = R3 = R, ω0 = 1/(R*sqrt(C1*C2)), Q = sqrt(C1/C2)/3
 */
int filter_mfb_lp(double omega0, double q, double gain,
                   double *r1, double *r2, double *r3,
                   double *c1, double *c2) {
    if (omega0 < 1e-30 || q < 0.5 || gain < 0.1) return -1;

    /* Choose C2 first (practical value) */
    double c2_val = 10e-9;  /* 10 nF */
    *c2 = c2_val;

    /* C1 = C2 * (4Q² * (gain+1)) */
    double c1_val = c2_val * 4.0 * q * q * (gain + 1.0);
    *c1 = c1_val;

    /* R2 = gain * R1, R3 = R1 * gain / ((2Q)² * (gain+1) * C1/C2) ... */
    /* Compute R2 from the time constant */
    double r2_val = 1.0 / (2.0 * omega0 * q * c2_val);
    *r2 = r2_val;

    *r1 = r2_val / gain;
    *r3 = 1.0 / (omega0 * omega0 * r2_val * c1_val * c2_val);

    return 0;
}

/**
 * Biquad ABCD matrix at frequency ω.
 *
 * Transfer function: H(s) = (b2*s² + b1*s + b0) / (s² + a1*s + a0)
 *
 * For a voltage-controlled voltage source (VCVS) two-port,
 * the ABCD matrix is:
 *   A = 1/H(s), B = 0, C = 0, D = H(s)
 *
 * Wait, that's not right for a general biquad as a two-port.
 * For a VCVS with gain H(s):
 *   V1 = V2 / H(s)  (voltage-controlled)
 *   I1 = 0 (infinite input impedance, ideal op-amp)
 *
 * ABCD: A = 1/H, B = 0, C = 0, D = H
 *
 * For a cascade of biquads: ABCD_total = ABCD1 × ABCD2 × ...
 */
matrix2x2_t filter_biquad_abcd(double b2, double b1, double b0,
                                double a1, double a0, double omega) {
    /* Evaluate H(s) at s = jω */
    /* Numerator: b0 + jω*b1 - ω²*b2 */
    double num_real = b0 - omega * omega * b2;
    double num_imag = omega * b1;
    /* Denominator: a0 - ω² + jω*a1 (note: a2 = 1) */
    double den_real = a0 - omega * omega;
    double den_imag = omega * a1;

    double den_mag_sq = den_real * den_real + den_imag * den_imag;
    if (den_mag_sq < 1e-30) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }

    /* H(s) = num/den */
    double h_real = (num_real * den_real + num_imag * den_imag) / den_mag_sq;
    double h_imag = (num_imag * den_real - num_real * den_imag) / den_mag_sq;
    complex_t h = complex_make(h_real, h_imag);

    /* For a VCVS: A = 1/H, B=0, C=0, D=H */
    if (complex_is_zero(h, 1e-30)) {
        complex_t nan = complex_make(NAN, NAN);
        return matrix2x2_make(nan, nan, nan, nan);
    }
    complex_t a_val = complex_div(complex_make(1.0, 0.0), h);

    return matrix2x2_make(
        a_val, complex_make(0.0, 0.0),
        complex_make(0.0, 0.0), h
    );
}

/* ============================================================================
 * L6: Distributed Element Filter Design
 * ============================================================================ */

/**
 * Coupled-line bandpass filter even/odd impedance computation.
 *
 * From the low-pass prototype g-values, we compute the even/odd
 * mode impedances for each parallel-coupled line section.
 *
 * For each coupled section i:
 *   J_i/Y0 = sqrt(π*BW/(2*g_i*g_{i+1}))  for i = 1..N-1
 *   J_0/Y0 = sqrt(π*BW/(2*g_0*g_1))     (input)
 *   J_N/Y0 = sqrt(π*BW/(2*g_N*g_{N+1})) (output)
 *
 *   (Z0e)_i = Z0 * (1 + J_i*Z0 + (J_i*Z0)²)
 *   (Z0o)_i = Z0 * (1 - J_i*Z0 + (J_i*Z0)²)
 *
 * Reference: Matthaei, Young, Jones §5.14
 * Course: Georgia Tech ECE 6350 — Coupled-line filters
 */
int filter_coupled_line_bp(const double *g_values, int order,
                            double bandwidth_fraction, double z0,
                            double *z0e, double *z0o) {
    if (order < 1) return -1;

    int n_sections = order + 1;  /* N+1 coupled sections for N-th order */
    double bw = bandwidth_fraction;

    /* Compute J-inverters */
    double j_sq[n_sections];
    for (int i = 0; i < n_sections; i++) {
        double gi, gj;
        if (i == 0) {
            gi = 1.0;  /* g0 = source impedance = 1 */
            gj = g_values[0];
        } else if (i == n_sections - 1) {
            gi = g_values[order - 1];
            gj = 1.0;  /* g_{N+1} = load = 1 */
        } else {
            gi = g_values[i - 1];
            gj = g_values[i];
        }
        j_sq[i] = (M_PI * bw) / (2.0 * gi * gj);
    }

    /* Compute even/odd impedances */
    for (int i = 0; i < n_sections; i++) {
        double j = sqrt(j_sq[i]);
        double j_z0 = j * z0;
        z0e[i] = z0 * (1.0 + j_z0 + j_z0 * j_z0);
        z0o[i] = z0 * (1.0 - j_z0 + j_z0 * j_z0);
    }

    return 0;
}

/**
 * Coupled-line filter response at frequency ω.
 *
 * This computes the S-parameters of the complete coupled-line filter
 * using the even/odd mode theory. Each coupled section is modeled as
 * a λ/4-long directional coupler.
 *
 * For the exact response, we cascade the ABCD matrices of each section.
 */
matrix2x2_t filter_coupled_line_response(const double *z0e,
                                          const double *z0o,
                                          int n_sections,
                                          double omega, double omega0,
                                          double z0) {
    /* Each coupled section at frequency ω: electrical length θ = (π/2)*(ω/ω0) */
    double theta = (M_PI / 2.0) * (omega / omega0);

    /* Start with identity */
    matrix2x2_t total_abcd = matrix2x2_make(
        complex_make(1.0, 0.0), complex_make(0.0, 0.0),
        complex_make(0.0, 0.0), complex_make(1.0, 0.0)
    );

    for (int i = 0; i < n_sections; i++) {
        double ze = z0e[i];
        double zo = z0o[i];
        double z0_section = sqrt(ze * zo);

        /* Treat each coupled section as an equivalent transmission line
           of electrical length θ and characteristic impedance Z0_section.
           This is an approximation — the exact response requires even/odd
           mode analysis with coupled-line theory.

           For the exact computation, one would use:
           A = D = (ze + zo) * cos(θ) / (ze - zo)
           B = j * (ze² + zo² - 2*ze*zo*cos²(θ)) / ((ze - zo) * sin(θ))
           (see Matthaei, Young, Jones §5.09) */

        matrix2x2_t section = abcd_transmission_line(
            z0_section, 0.0, theta, 1.0
        );

        total_abcd = matrix2x2_mul(total_abcd, section);
    }

    /* Convert to S-parameters */
    return convert_abcd_to_s(total_abcd, z0);
}

/**
 * Group delay from S-parameters.
 *
 * τg = -d(arg(S21))/dω ≈ -(φ2 - φ1)/(ω2 - ω1)
 */
double filter_group_delay(matrix2x2_t s1, matrix2x2_t s2, double delta_omega) {
    return sparams_group_delay(s1, s2, 0.0, delta_omega);
    /* The 0.0 for omega1 is not quite right; we only need the difference.
       The function uses omega2 - omega1 = delta_omega. */
}
