/**
 * @file complex_power.c
 * @brief Complex power computations, impedance analysis, and conservation laws
 *
 * Complex power S = P + jQ is the cornerstone of AC power theory,
 * introduced by C.P. Steinmetz in 1893. It transforms AC circuit analysis
 * from differential equations to algebraic complex-number operations.
 *
 * Key relationships:
 *   S = V_phasor × conj(I_phasor) = V²/Z* = I²×Z
 *   P = Re(S) = V×I×cos(φ) — real (active) power [Watts]
 *   Q = Im(S) = V×I×sin(φ) — reactive power [VAR]
 *   |S| = V×I — apparent power [VA]
 *   PF = P/|S| = cos(φ) — power factor
 *
 * Knowledge points:
 *   L2: Complex power S from V, I phasors
 *   L2: Power from impedance/admittance
 *   L3: Impedance/admittance from power measurements
 *   L3: Series and parallel complex power decomposition
 *   L4: Steinmetz complex power balance
 *   L5: Optimal shunt compensation admittance
 *   L6: Capacitor sizing for PF correction
 *   L7: Three-phase capacitor bank configuration
 *
 * References:
 *   - Steinmetz (1893) AIEE Proceedings
 *   - IEEE Std 1459-2010
 *   - Erickson & Maksimovic (2001) Ch.18
 */

#include "complex_power.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==========================================================================
 * L2: Complex Power from Voltage and Current Phasors
 * ========================================================================== */

double complex cp_complex_power(double complex v_phasor, double complex i_phasor)
{
    /* S = V × I*  — the fundamental complex power equation */
    return v_phasor * conj(i_phasor);
}

/* ==========================================================================
 * L3: Complex Power from Voltage and Impedance
 * ========================================================================== */

double complex cp_power_from_v_z(double v_magnitude, double complex z_impedance)
{
    /* S = V² / Z*
     *
     * Derivation:
     *   I = V / Z
     *   S = V × I* = V × (V/Z)* = V² / Z*
     *
     * For Z = R + jX:
     *   S = V² / (R - jX) = V² (R + jX) / (R² + X²)
     *   P = V² R / (R² + X²)
     *   Q = V² X / (R² + X²)
     */
    if (cabs(z_impedance) < CP_EPS) return 0.0 + I * 0.0;

    double complex z_conj = conj(z_impedance);
    return (v_magnitude * v_magnitude) / z_conj;
}

/* ==========================================================================
 * L3: Complex Power from Current and Impedance
 * ========================================================================== */

double complex cp_power_from_i_z(double i_magnitude, double complex z_impedance)
{
    /* S = I² × Z
     *
     * Derivation:
     *   V = I × Z
     *   S = V × I* = (I×Z) × I* = |I|² × Z
     *
     * For Z = R + jX:
     *   P = I² × R
     *   Q = I² × X
     */
    return (i_magnitude * i_magnitude) * z_impedance;
}

/* ==========================================================================
 * L3: Load Impedance from Power Measurements
 * ========================================================================== */

double complex cp_impedance_from_power(double v_rms, double p, double q)
{
    /* Z = V² / S* = V² / (P - jQ)
     *
     * Alternative: Z = V² × (P + jQ) / (P² + Q²)
     *   R = V² × P / (P² + Q²)
     *   X = V² × Q / (P² + Q²)
     *
     * This is fundamental to electrical metrology: determining
     * the equivalent circuit of a load from V, P, Q measurements.
     */
    double s_sq = p * p + q * q;
    if (s_sq < CP_EPS) return 0.0 + I * 0.0;

    double v_sq = v_rms * v_rms;
    double r = v_sq * p / s_sq;
    double x = v_sq * q / s_sq;

    return r + I * x;
}

/* ==========================================================================
 * L3: Admittance from Power Measurements
 * ========================================================================== */

double complex cp_admittance_from_power(double v_rms, double p, double q)
{
    /* Y = 1/Z = S* / V² = (P - jQ) / V²
     *
     * G = P / V²  (conductance)
     * B = -Q / V² (susceptance)
     *
     * Note: B > 0 for capacitive (leading), B < 0 for inductive (lagging).
     */
    double v_sq = v_rms * v_rms;
    if (v_sq < CP_EPS) return 0.0 + I * 0.0;

    double g = p / v_sq;
    double b = -q / v_sq;

    return g + I * b;
}

/* ==========================================================================
 * L4: Complex Power Conservation (Steinmetz)
 * ========================================================================== */

int cp_verify_complex_power_balance(const double complex *s_complex,
                                     size_t n_branches, double tolerance)
{
    if (!s_complex) return -1;

    double complex sum_s = 0.0 + I * 0.0;
    for (size_t i = 0; i < n_branches; i++) {
        sum_s += s_complex[i];
    }

    /* Both real and imaginary parts must be zero within tolerance */
    double p_sum = creal(sum_s);
    double q_sum = cimag(sum_s);

    if (fabs(p_sum) > tolerance || fabs(q_sum) > tolerance) return 1;
    return 0;
}

/* ==========================================================================
 * L4: Series Complex Power Decomposition
 * ========================================================================== */

int cp_series_complex_power(const double complex *z_elements, size_t n_elements,
                            double complex i_phasor, double complex *s_total)
{
    if (!z_elements || !s_total || n_elements == 0) return -1;

    double i_mag_sq = creal(i_phasor * conj(i_phasor)); /* |I|² */

    *s_total = 0.0 + I * 0.0;
    for (size_t k = 0; k < n_elements; k++) {
        /* S_k = |I|² × Z_k  (since same I flows through all) */
        *s_total += i_mag_sq * z_elements[k];
    }
    return 0;
}

/* ==========================================================================
 * L4: Parallel Complex Power Decomposition
 * ========================================================================== */

int cp_parallel_complex_power(const double complex *y_elements, size_t n_elements,
                              double complex v_phasor, double complex *s_total)
{
    if (!y_elements || !s_total || n_elements == 0) return -1;

    double v_mag_sq = creal(v_phasor * conj(v_phasor));

    *s_total = 0.0 + I * 0.0;
    for (size_t k = 0; k < n_elements; k++) {
        /* S_k = |V|² × conj(Y_k)  (same V across all) */
        *s_total += v_mag_sq * conj(y_elements[k]);
    }
    return 0;
}

/* ==========================================================================
 * L5: Optimal Shunt Compensation Admittance
 * ========================================================================== */

double complex cp_compensation_admittance(double v_rms, double p_load,
                                           double pf_old, double pf_target,
                                           int is_lagging)
{
    if (v_rms < CP_EPS || p_load < 0.0) return 0.0 + I * 0.0;
    if (pf_old <= 0.0 || pf_target <= 0.0) return 0.0 + I * 0.0;
    if (pf_old >= 1.0 || pf_target > 1.0) return 0.0 + I * 0.0;

    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);

    /* B_c = (P/V²) × (tan φ_old - tan φ_new)
     *
     * For lagging load (inductive, φ > 0), we add capacitive susceptance B_c > 0.
     * Y_comp = j × B_c  (purely imaginary, no real part — ideal capacitor)
     */
    double tan_diff = tan(phi_old) - tan(phi_new);
    if (tan_diff < 0.0) tan_diff = 0.0; /* Don't over-correct */

    double b_c = (p_load / (v_rms * v_rms)) * tan_diff;

    /* If original is leading, we need inductive compensation (B_c < 0) */
    if (!is_lagging) b_c = -b_c;

    return 0.0 + I * b_c; /* pure susceptance */
}

/* ==========================================================================
 * L5: Capacitor Sizing for PF Correction (Single-Phase)
 * ========================================================================== */

double cp_capacitor_for_pf_correction(double v_rms, double p_load,
                                       double f_hz, double pf_old,
                                       double pf_target, int is_lagging)
{
    if (v_rms < CP_EPS || f_hz < CP_EPS || p_load <= 0.0) return -1.0;
    if (pf_old <= 0.0 || pf_target <= 0.0) return -1.0;
    if (pf_old > 1.0 || pf_target > 1.0) return -1.0;

    /* C = B_c / (2πf) = P(tan φ_old - tan φ_new) / (2πf V²) */
    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);
    double tan_diff = tan(phi_old) - tan(phi_new);
    if (tan_diff < 0.0) tan_diff = 0.0;

    double c = p_load * tan_diff / (2.0 * M_PI * f_hz * v_rms * v_rms);

    /* Leading loads need inductive compensation, not capacitors */
    if (!is_lagging) return -1.0;

    return c;
}

/* ==========================================================================
 * L7: Three-Phase Capacitor Bank Configuration
 * ========================================================================== */

int cp_three_phase_pf_capacitor(double v_ll_rms, double p_3phase,
                                 double f_hz, double pf_old, double pf_target,
                                 int is_lagging, double *c_delta_out,
                                 double *c_wye_out)
{
    if (!c_delta_out || !c_wye_out) return -1;
    if (v_ll_rms < CP_EPS || f_hz < CP_EPS || p_3phase <= 0.0) return -1;
    if (pf_old <= 0.0 || pf_target <= 0.0) return -1;
    if (pf_old > 1.0 || pf_target > 1.0) return -1;

    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);
    double tan_diff = tan(phi_old) - tan(phi_new);
    if (tan_diff < 0.0) tan_diff = 0.0;

    double q_total = p_3phase * tan_diff;
    double omega = 2.0 * M_PI * f_hz;

    /* Delta connection: each capacitor sees V_LL
     * Q_delta_per_phase = V_LL² × ω × C_delta
     * Q_3phase_total = 3 × Q_delta_per_phase
     * C_delta = Q_total / (3 × ω × V_LL²)
     */
    *c_delta_out = q_total / (3.0 * omega * v_ll_rms * v_ll_rms);

    /* Wye connection: each capacitor sees V_LL/√3
     * Q_wye_per_phase = (V_LL/√3)² × ω × C_wye = V_LL² × ω × C_wye / 3
     * Q_3phase_total = 3 × Q_wye_per_phase = V_LL² × ω × C_wye
     * C_wye = Q_total / (ω × V_LL²) = 3 × C_delta
     *
     * Delta gives 3× the reactive power per unit capacitance compared to wye.
     */
    *c_wye_out = q_total / (omega * v_ll_rms * v_ll_rms);

    if (!is_lagging) {
        /* For leading loads, capacitors not applicable — need reactors */
        *c_delta_out = -1.0;
        *c_wye_out   = -1.0;
    } else {
        if (tan_diff <= 0.0) {
            *c_delta_out = 0.0;
            *c_wye_out   = 0.0;
        }
    }

    return 0;
}

/* ==========================================================================
 * L5: kVAR Rating Required
 * ========================================================================== */

double cp_kvar_required(double p_load, double pf_old, double pf_target)
{
    if (p_load <= 0.0 || pf_old <= 0.0 || pf_target <= 0.0) return -1.0;
    if (pf_old > 1.0 || pf_target > 1.0 || pf_target < pf_old) return -1.0;

    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);
    double q_c = p_load * (tan(phi_old) - tan(phi_new));

    return q_c; /* In VAR */
}
