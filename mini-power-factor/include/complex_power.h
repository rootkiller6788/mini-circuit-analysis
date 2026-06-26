/**
 * @file complex_power.h
 * @brief Complex power theory and mathematical structures for AC power analysis
 *
 * Complex power S = P + jQ is the fundamental mathematical construct
 * that unifies real and reactive power in a single phasor quantity.
 * Introduced by C.P. Steinmetz in 1893, it revolutionized AC power
 * analysis by enabling algebraic (rather than graphical) treatment.
 *
 * S = V · I* = P + jQ  where I* is the complex conjugate of I
 *
 * References:
 *   - C.P. Steinmetz, "Complex Quantities and Their Use in Electrical
 *     Engineering" (1893) — AIEE Proceedings
 *   - W. Leonhard, "Control of Electrical Drives" (3rd ed, 2001)
 *   - IEEE Std 1459-2010 §3.1
 *   - MIT 6.061 "Introduction to Electric Power Systems"
 *   - Berkeley EE137A "Power Electronics"
 *
 * Knowledge coverage:
 *   L2 (Concepts):   Complex power, apparent S, power triangle, power angle
 *   L3 (Math):       Complex arithmetic, phasor rotation, quadrature decomposition
 *   L4 (Laws):       Steinmetz's formulation, conservation of complex power
 *   L5 (Algorithms): Numerical complex power computation from sampled data
 */

#ifndef COMPLEX_POWER_H
#define COMPLEX_POWER_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L3: Complex Phasor Operations
 * ========================================================================== */

/** Small tolerance for floating-point comparisons in phasor domain */
#define CP_EPS 1e-12

/**
 * @brief Create a complex number from polar form (magnitude, angle).
 *
 * z = r × e^{jθ} = r × (cos θ + j sin θ)
 *
 * @param magnitude  Magnitude (radius) r
 * @param angle_rad  Angle θ in radians
 * @return Complex number in rectangular form
 */
static inline double complex cp_polar(double magnitude, double angle_rad) {
    return magnitude * (cos(angle_rad) + I * sin(angle_rad));
}

/**
 * @brief Get magnitude of a complex number.
 *
 * |z| = sqrt(Re(z)² + Im(z)²)
 */
static inline double cp_abs(double complex z) {
    return cabs(z);
}

/**
 * @brief Get argument (angle) of a complex number in radians.
 *
 * arg(z) = atan2(Im(z), Re(z))
 */
static inline double cp_arg_rad(double complex z) {
    return carg(z);
}

/**
 * @brief Get argument in degrees.
 */
static inline double cp_arg_deg(double complex z) {
    return carg(z) * 180.0 / M_PI;
}

/**
 * @brief Complex conjugate.
 *
 * z* = Re(z) - j Im(z)
 */
static inline double complex cp_conj(double complex z) {
    return conj(z);
}

/**
 * @brief Phase rotation: multiply by e^{jθ}.
 *
 * Rotates a phasor by angle θ counterclockwise.
 * Used in dq0 transformation and symmetrical components.
 */
static inline double complex cp_rotate(double complex z, double angle_rad) {
    double complex rot = cos(angle_rad) + I * sin(angle_rad);
    return z * rot;
}

/**
 * @brief Compute real part squared for energy calculations.
 */
static inline double cp_real_sq(double complex z) {
    double r = creal(z);
    return r * r;
}

/**
 * @brief Compute the product z1 × conj(z2) used in complex power.
 */
static inline double complex cp_mul_conj(double complex z1, double complex z2) {
    return z1 * conj(z2);
}

/* ==========================================================================
 * L2: Complex Power Computation
 * ========================================================================== */

/**
 * @brief Complex power from voltage and current phasors.
 *
 * Fundamental equation of AC power theory:
 *   S = V_phasor × conj(I_phasor) = P + jQ
 *
 * where V_phasor = V ∠θv and I_phasor = I ∠θi
 *
 * Derivation:
 *   S = (V ∠θv) × (I ∠(-θi)) = V×I ∠(θv - θi)
 *     = V×I cos(θv-θi) + j V×I sin(θv-θi)
 *     = P + jQ
 *
 * @param v_phasor  Voltage phasor (complex, RMS)
 * @param i_phasor  Current phasor (complex, RMS)
 * @return Complex power S = P + jQ [VA]
 */
double complex cp_complex_power(double complex v_phasor, double complex i_phasor);

/**
 * @brief Extract real power P from complex power.
 */
static inline double cp_real_power(double complex s) {
    return creal(s);
}

/**
 * @brief Extract reactive power Q from complex power.
 */
static inline double cp_reactive_power(double complex s) {
    return cimag(s);
}

/**
 * @brief Compute apparent power |S| = sqrt(P² + Q²).
 *
 * The magnitude of complex power. Represents the total "vector" power
 * that the source must supply. For sinusoidal systems, |S| = V_rms × I_rms.
 */
static inline double cp_apparent_power(double complex s) {
    return cabs(s);
}

/**
 * @brief Compute power factor from complex power.
 *
 * PF = P / |S| = cos(φ)  where φ = atan2(Q, P)
 */
static inline double cp_power_factor(double complex s) {
    double abs_s = cabs(s);
    if (abs_s < CP_EPS) return 1.0;
    return creal(s) / abs_s;
}

/**
 * @brief Compute power angle φ from complex power.
 *
 * φ = atan2(Q, P) = power angle
 * φ > 0 → lagging PF (inductive)
 * φ < 0 → leading PF (capacitive)
 */
static inline double cp_power_angle_rad(double complex s) {
    return atan2(cimag(s), creal(s));
}

/**
 * @brief Compute the displacement power factor from complex power.
 *
 * Same as PF for sinusoidal waveforms, but distinct from true PF
 * when harmonics are present. dPF = P₁ / S₁.
 */
static inline double cp_displacement_pf(double complex s) {
    double abs_s = cabs(s);
    if (abs_s < CP_EPS) return 1.0;
    double pf = creal(s) / abs_s;
    if (pf > 1.0) return 1.0;
    if (pf < -1.0) return -1.0;
    return pf;
}

/* ==========================================================================
 * L3: Impedance-Based Power Analysis
 * ========================================================================== */

/**
 * @brief Compute complex power from voltage and impedance.
 *
 * S = V² / Z*
 * P = |V|² × Re(1/Z) = |V|² × G     where G = conductance
 * Q = |V|² × Im(1/Z) = -|V|² × B    where B = susceptance
 *
 * For a series RL load Z = R + jX_L:
 *   P = V² × R / (R² + X_L²)
 *   Q = V² × X_L / (R² + X_L²)
 *
 * @param v_magnitude  RMS voltage magnitude [V]
 * @param z_impedance  Load impedance [Ohm] (complex)
 * @return Complex power S = P + jQ
 */
double complex cp_power_from_v_z(double v_magnitude, double complex z_impedance);

/**
 * @brief Compute complex power from current and impedance.
 *
 * S = I² × Z
 * P = |I|² × R
 * Q = |I|² × X
 *
 * @param i_magnitude  RMS current magnitude [A]
 * @param z_impedance  Load impedance [Ohm] (complex)
 * @return Complex power S = P + jQ
 */
double complex cp_power_from_i_z(double i_magnitude, double complex z_impedance);

/**
 * @brief Compute load impedance from complex power and voltage.
 *
 * Z = V² / S* = V² / (P - jQ) = V²(P + jQ) / (P² + Q²)
 *
 * Useful for determining the equivalent circuit of an unknown load
 * from power measurements.
 *
 * @param v_rms   RMS voltage [V]
 * @param p       Real power [W]
 * @param q       Reactive power [VAR]
 * @return Load impedance Z = R + jX [Ohm]
 */
double complex cp_impedance_from_power(double v_rms, double p, double q);

/**
 * @brief Compute admittance Y = 1/Z = G + jB from power.
 *
 * Y = S* / V² = (P - jQ) / V²
 * G = P / V²   (conductance)
 * B = -Q / V²  (susceptance, negative for inductive)
 *
 * @param v_rms   RMS voltage [V]
 * @param p       Real power [W]
 * @param q       Reactive power [VAR]
 * @return Admittance Y = G + jB [Siemens]
 */
double complex cp_admittance_from_power(double v_rms, double p, double q);

/* ==========================================================================
 * L4: Conservation Laws in Complex Power Form
 * ========================================================================== */

/**
 * @brief Verify Steinmetz power balance: Σ S_k = 0 for all branches.
 *
 * In a network with sinusoidal sources at the same frequency, the
 * algebraic sum of complex powers over all branches is zero.
 *
 * This is a stronger statement than separate P and Q conservation
 * because it implies both simultaneously via complex arithmetic.
 *
 * Σ (P_k + jQ_k) = 0  →  Σ P_k = 0 AND Σ Q_k = 0
 *
 * @param s_complex    Array of complex powers [VA] for each branch
 * @param n_branches   Number of branches
 * @param tolerance    Numerical tolerance for zero check
 * @return 0 if balanced, 1 if imbalance detected
 */
int cp_verify_complex_power_balance(const double complex *s_complex,
                                     size_t n_branches, double tolerance);

/**
 * @brief Verify that series elements share the same complex power.
 *
 * In a series connection, all elements carry the same current I,
 * so S_k = I² × Z_k. The total S = Σ S_k = I² × Σ Z_k = I² × Z_eq.
 *
 * @param z_elements   Array of complex impedances [Ohm]
 * @param n_elements   Number of series elements
 * @param i_phasor     Series current phasor [A]
 * @param s_total      [out] Total complex power [VA]
 * @return 0 on success
 */
int cp_series_complex_power(const double complex *z_elements, size_t n_elements,
                            double complex i_phasor, double complex *s_total);

/**
 * @brief Verify that parallel elements share the same voltage.
 *
 * In a parallel connection, each element has the same voltage V across it,
 * so S_k = V² / Z_k*. The total S = Σ S_k = V² × Σ (1/Z_k*) = V² × Y*.
 *
 * @param y_elements   Array of complex admittances [Siemens]
 * @param n_elements   Number of parallel elements
 * @param v_phasor     Parallel voltage phasor [V]
 * @param s_total      [out] Total complex power [VA]
 * @return 0 on success
 */
int cp_parallel_complex_power(const double complex *y_elements, size_t n_elements,
                              double complex v_phasor, double complex *s_total);

/* ==========================================================================
 * L5: Complex Power Compensation Optimization
 * ========================================================================== */

/**
 * @brief Compute optimal shunt compensation admittance for PF correction.
 *
 * To raise PF from cos φ_old to cos φ_new by adding shunt susceptance B_c:
 *   B_c = (P/V²) × (tan φ_old - tan φ_new)
 *
 * Capacitive compensation (B_c > 0) corrects lagging PF.
 * Inductive compensation (B_c < 0) corrects leading PF.
 *
 * Result in complex admittance form: Y_comp = 0 + jB_c [Siemens]
 *
 * @param v_rms      RMS voltage [V]
 * @param p_load     Real power of load [W]
 * @param pf_old     Original power factor
 * @param pf_target  Target power factor
 * @param is_lagging Non-zero if original PF is lagging
 * @return Compensation admittance Y_c [Siemens] (purely imaginary)
 */
double complex cp_compensation_admittance(double v_rms, double p_load,
                                           double pf_old, double pf_target,
                                           int is_lagging);

/**
 * @brief Compute the capacitor value needed for PF correction.
 *
 * C = B_c / (2πf) = P(tan φ_old - tan φ_new) / (2πf V²)
 *
 * @param v_rms      RMS voltage [V]
 * @param p_load     Real power [W]
 * @param f_hz       System frequency [Hz]
 * @param pf_old     Original power factor
 * @param pf_target  Target power factor
 * @param is_lagging Non-zero if original PF is lagging
 * @return Capacitance [Farads], or -1 on error
 */
double cp_capacitor_for_pf_correction(double v_rms, double p_load,
                                       double f_hz, double pf_old,
                                       double pf_target, int is_lagging);

/**
 * @brief Compute the optimal capacitor bank configuration for three-phase.
 *
 * For three-phase PF correction, capacitors can be connected in delta
 * or wye. Delta connection provides 3× the per-phase capacitance
 * of wye for the same voltage rating.
 *
 * @param v_ll_rms    Line-to-line voltage [V]
 * @param p_3phase    Total three-phase real power [W]
 * @param f_hz        System frequency [Hz]
 * @param pf_old      Original PF
 * @param pf_target   Target PF
 * @param is_lagging  Non-zero if original PF is lagging
 * @param c_delta_out [out] Capacitance per phase (delta) [F]
 * @param c_wye_out   [out] Capacitance per phase (wye) [F]
 * @return 0 on success, -1 on error
 */
int cp_three_phase_pf_capacitor(double v_ll_rms, double p_3phase,
                                 double f_hz, double pf_old, double pf_target,
                                 int is_lagging, double *c_delta_out,
                                 double *c_wye_out);

/**
 * @brief Compute kVAR rating of capacitor bank for PF correction.
 *
 * Q_c = P × (tan(acos(pf_old)) - tan(acos(pf_target)))
 *
 * @param p_load     Real power [W]
 * @param pf_old     Original power factor
 * @param pf_target  Target power factor
 * @return Required kVAR rating
 */
double cp_kvar_required(double p_load, double pf_old, double pf_target);

#ifdef __cplusplus
}
#endif

#endif /* COMPLEX_POWER_H */
