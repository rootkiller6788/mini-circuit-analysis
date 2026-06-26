/**
 * @file pf_correction.c
 * @brief Power factor correction algorithms — capacitor sizing, step control,
 *        detuning reactor design, and active PFC control
 *
 * Power factor correction is the application of reactive power compensation
 * to improve the power factor of an electrical load. For industrial loads,
 * this typically means adding shunt capacitors to counteract inductive
 * reactive power from motors, transformers, and discharge lighting.
 *
 * Active PFC uses switched-mode power converters to shape the input current
 * to follow the sinusoidal voltage, achieving near-unity PF and low THD.
 *
 * Knowledge points:
 *   L1: Capacitor bank step configuration
 *   L1: Active PFC boost converter state model
 *   L5: Capacitor sizing (single-phase and three-phase)
 *   L5: Automatic step bank design
 *   L5: Detuning reactor resonance analysis
 *   L5: Active PFC average current mode control
 *   L5: Boost PFC THD estimation
 *   L6: Industrial PF correction end-to-end (sizing + resonance + payback)
 *   L6: Step switching optimization with hysteresis
 *   L6: Harmonic resonance risk assessment
 *
 * References:
 *   - IEEE Std 18-2012, "Shunt Power Capacitors"
 *   - IEEE Std 519-2014, "Harmonic Control"
 *   - Erickson & Maksimovic, "Power Electronics" (2001) Ch.18
 *   - Dixon, "Average Current Mode Control" (Unitrode AN)
 */

#include "pf_correction.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define PFC_EPS 1e-12

/* ==========================================================================
 * L5: Single-Phase Capacitor Sizing
 * ========================================================================== */

int pfc_size_capacitor(double p_load, double v_rms, double f_hz,
                       double pf_old, double pf_target, int is_lagging,
                       double *capacitance, double *kvar_rating)
{
    if (!capacitance || !kvar_rating) return -1;
    if (p_load <= 0.0 || v_rms < PFC_EPS || f_hz < PFC_EPS) return -1;
    if (pf_old <= 0.0 || pf_target <= 0.0 || pf_old > 1.0 || pf_target > 1.0)
        return -1;
    if (pf_target < pf_old) return -1;

    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);
    double tan_diff = tan(phi_old) - tan(phi_new);
    if (tan_diff < 0.0) tan_diff = 0.0;

    /* Q_c = P × (tan φ_old - tan φ_new) */
    double q_c = p_load * tan_diff;
    *kvar_rating = q_c;

    /* C = Q_c / (2πf V²) */
    double omega = 2.0 * M_PI * f_hz;
    *capacitance = q_c / (omega * v_rms * v_rms);

    if (!is_lagging) {
        /* Leading load: PF correction requires inductors, not capacitors */
        *capacitance = -1.0;
    }

    return 0;
}

/* ==========================================================================
 * L5: Three-Phase Capacitor Bank Sizing
 * ========================================================================== */

int pfc_size_capacitor_3phase(double p_3phase, double v_ll_rms, double f_hz,
                               double pf_old, double pf_target, int is_lagging,
                               int connection, double *c_per_phase,
                               double *kvar_total)
{
    if (!c_per_phase || !kvar_total) return -1;
    if (p_3phase <= 0.0 || v_ll_rms < PFC_EPS || f_hz < PFC_EPS) return -1;
    if (pf_old <= 0.0 || pf_target <= 0.0) return -1;

    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);
    double tan_diff = tan(phi_old) - tan(phi_new);
    if (tan_diff < 0.0) tan_diff = 0.0;

    double q_total = p_3phase * tan_diff;
    *kvar_total = q_total;

    double omega = 2.0 * M_PI * f_hz;

    if (connection == 0) {
        /* Delta: C_phase = Q_total / (3 × ω × V_LL²) */
        *c_per_phase = q_total / (3.0 * omega * v_ll_rms * v_ll_rms);
    } else {
        /* Wye: C_phase = Q_total / (ω × V_LL²) */
        *c_per_phase = q_total / (omega * v_ll_rms * v_ll_rms);
    }

    if (!is_lagging) *c_per_phase = -1.0;
    return 0;
}

/* ==========================================================================
 * L5: Automatic Capacitor Bank Step Design
 * ========================================================================== */

int pfc_design_step_bank(double kvar_required, uint32_t num_steps,
                          double resolution_kvar, pfc_capacitor_bank_t *bank)
{
    if (!bank || num_steps == 0 || num_steps > 16) return -1;
    if (kvar_required <= 0.0 || resolution_kvar <= 0.0) return -1;

    memset(bank, 0, sizeof(*bank));
    bank->num_steps   = num_steps;
    bank->total_kvar  = kvar_required;
    bank->v_rated     = 480.0; /* common default */
    bank->f_rated     = 60.0;
    bank->active_mask = 0;

    /* Binary-weighted step sizes for maximum range with minimum steps.
     * step_kvar[i] = min_step × 2^i, scaled to fit total. */
    double min_step = resolution_kvar;
    double sum = 0.0;
    for (uint32_t i = 0; i < num_steps; i++) {
        double step_size = min_step * (1U << i);
        if (step_size + sum > kvar_required) {
            step_size = kvar_required - sum;
        }
        bank->step_kvar[i] = step_size;
        sum += step_size;
        if (sum >= kvar_required - PFC_EPS) break;
    }

    /* If sum falls short, distribute remainder equally */
    if (sum < kvar_required - PFC_EPS) {
        double remainder = kvar_required - sum;
        double each = remainder / (double)num_steps;
        for (uint32_t i = 0; i < num_steps; i++) {
            bank->step_kvar[i] += each;
        }
    }

    bank->discharge_time_sec = 60.0; /* IEEE Std 18: discharge to <50V in 1 min */
    bank->step_delay_sec     = 15;   /* Typical contactor delay */
    return 0;
}

/* ==========================================================================
 * L5: Reactive Power Contribution of a Capacitor Step
 * ========================================================================== */

double pfc_step_kvar(double c_farads, double v_rms, double f_hz,
                     double detune_factor)
{
    if (c_farads < 0.0 || v_rms < PFC_EPS || f_hz < PFC_EPS) return 0.0;
    if (detune_factor < 0.0 || detune_factor >= 1.0) detune_factor = 0.0;

    /* Without detuning: Q = V² × ωC
     * With detuning reactor (p%): effective reactance changes.
     *
     * X_c = 1/(ωC)
     * X_L_detune = p × X_c (where p = detune_factor, e.g., 0.07)
     * X_net = X_c - X_L_detune = X_c × (1 - p)
     * Q_eff = V² / X_net = V² × ωC / (1 - p)
     *
     * The detuning reactor reduces the effective capacitance at fundamental
     * frequency, so the kVAR output is slightly higher than the capacitor rating
     * (because the reactor voltage boost adds to the capacitor voltage).
     */
    double omega = 2.0 * M_PI * f_hz;
    double x_c = 1.0 / (omega * c_farads);
    double x_net = x_c * (1.0 - detune_factor);
    if (x_net < PFC_EPS) return 0.0;

    return (v_rms * v_rms) / x_net;
}

/* ==========================================================================
 * L6: Resonant Frequency with Capacitor Bank
 * ========================================================================== */

double pfc_resonant_frequency(double s_sc_va, double q_c_var, double f_sys_hz)
{
    /* Parallel resonance between system inductance and capacitor bank:
     *
     * System short-circuit impedance: X_sys = V² / S_sc
     *   L_sys = X_sys / (2πf_sys) = V² / (2πf_sys × S_sc)
     *
     * Capacitor reactance: X_c = V² / Q_c
     *   C = Q_c / (2πf_sys × V²)
     *
     * Resonance: ω_r L_sys = 1/(ω_r C)
     *   f_r = 1/(2π√(L_sys C)) = f_sys × √(S_sc / Q_c)
     *
     * This simplified formula gives the resonant frequency in Hz.
     */
    if (s_sc_va < PFC_EPS || q_c_var < PFC_EPS || f_sys_hz < PFC_EPS)
        return 0.0;

    return f_sys_hz * sqrt(s_sc_va / q_c_var);
}

/* ==========================================================================
 * L6: Detuning Necessity Check
 * ========================================================================== */

int pfc_needs_detuning(double f_res_hz, double f_sys_hz,
                        int *harmonic_order)
{
    if (!harmonic_order) return -1;
    if (f_sys_hz < 1.0 || f_res_hz < 1.0) return -1;

    /* Characteristic harmonics for 6-pulse converters (most common):
     * h = 6k ± 1 = 5, 7, 11, 13, 17, 19, 23, 25, ...
     */
    static const int characteristic_h[] = {
        3, 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35, 37
    };
    int n_h = (int)(sizeof(characteristic_h) / sizeof(characteristic_h[0]));
    int closest_h = 0;
    double min_dist = 1e9;

    for (int i = 0; i < n_h; i++) {
        double f_h = f_sys_hz * (double)characteristic_h[i];
        double dist = fabs(f_res_hz - f_h);
        if (dist < min_dist) {
            min_dist = dist;
            closest_h = characteristic_h[i];
        }
    }

    *harmonic_order = closest_h;

    /* Detuning recommended if resonant frequency is within ±10% of a
     * characteristic harmonic. */
    double f_h = f_sys_hz * (double)closest_h;
    double rel_dist = min_dist / f_h;

    if (rel_dist < 0.10) return 1; /* Detuning recommended */
    return 0; /* Safe */
}

/* ==========================================================================
 * L5: Boost PFC Controller Initialization
 * ========================================================================== */

int pfc_boost_init(pfc_boost_state_t *state, double v_in_rms, double v_out,
                    double p_out, double f_sw, double f_line)
{
    if (!state) return -1;
    if (v_in_rms <= 0.0 || v_out <= 0.0 || p_out <= 0.0) return -1;
    if (f_sw <= 0.0 || f_line <= 0.0) return -1;

    memset(state, 0, sizeof(*state));

    state->v_in_rms_nom = v_in_rms;
    state->v_out_ref    = v_out;
    state->p_out_nom    = p_out;
    state->f_sw         = f_sw;
    state->f_line       = f_line;
    state->v_out        = v_out;  /* Start at target */

    /* Design boost inductor for CCM operation:
     * L = V_in_rms² × (V_out - V_in_peak) / (ΔI × f_sw × V_out)
     * where ΔI ≈ 20% of peak input current
     */
    double v_in_peak = v_in_rms * 1.41421356;
    double i_in_peak = p_out / v_in_rms * 1.41421356;
    double ripple_i = 0.20 * i_in_peak;
    /* d_min is the minimum duty cycle (at peak input voltage).
     * CCM requires L > (V_in_peak × d_min) / (ΔI × f_sw) */
    double d_min = 1.0 - v_in_peak / v_out;
    state->l_boost = v_in_peak * d_min / (ripple_i * f_sw);
    if (state->l_boost < 1e-6) state->l_boost = 500e-6;

    /* Output capacitor: for 5% voltage ripple at 2×f_line
     * C_out = P_out / (4π × f_line × V_out × ΔV)
     */
    double delta_v = 0.05 * v_out;
    state->c_out = p_out / (2.0 * M_PI * f_line * v_out * delta_v);
    if (state->c_out < 100e-6) state->c_out = 470e-6;

    /* Voltage loop: bandwidth ≈ 10 Hz (well below 2×f_line = 120 Hz)
     * Current loop: bandwidth ≈ 2 kHz (well below f_sw/2)
     *
     * Standard PI tuning for PFC (Erickson & Maksimovic §18.3):
     */
    double bwi_v = 2.0 * M_PI * 10.0;
    double bwi_c = 2.0 * M_PI * 2000.0;

    state->kp_voltage = bwi_v * state->c_out;
    state->ki_voltage = bwi_v * bwi_v * state->c_out / 4.0;
    state->kp_current = bwi_c * state->l_boost;
    state->ki_current = bwi_c * bwi_c * state->l_boost / 4.0;

    return 0;
}

/* ==========================================================================
 * L5: Boost PFC Control Step (Average Current Mode)
 * ========================================================================== */

int pfc_boost_control_step(pfc_boost_state_t *state, double v_in_inst,
                            double i_in_meas, double ts_sec)
{
    if (!state || ts_sec <= 0.0) return -1;

    double v_in_abs = fabs(v_in_inst);
    double v_err = state->v_out_ref - state->v_out;

    /* Integrate voltage error */
    state->v_err_integral += v_err * ts_sec;

    /* Anti-windup: clamp integrator */
    if (state->v_err_integral > 10.0)  state->v_err_integral = 10.0;
    if (state->v_err_integral < -10.0) state->v_err_integral = -10.0;

    /* Voltage loop PI: produces current reference amplitude */
    double i_amplitude = state->kp_voltage * v_err
                          + state->ki_voltage * state->v_err_integral;
    if (i_amplitude < 0.0) i_amplitude = 0.0;

    /* Current reference: shaped to follow |v_in| waveform */
    double v_in_peak = state->v_in_rms_nom * 1.41421356;
    double shape = 0.0;
    if (v_in_peak > PFC_EPS) shape = v_in_abs / v_in_peak;
    if (shape > 1.0) shape = 1.0;

    state->i_ref_peak = i_amplitude;
    double i_ref = i_amplitude * shape;

    /* Current loop PI */
    double i_err = i_ref - i_in_meas;
    state->i_err_integral += i_err * ts_sec;

    /* Anti-windup for current integrator */
    if (state->i_err_integral > 5.0)  state->i_err_integral = 5.0;
    if (state->i_err_integral < -5.0) state->i_err_integral = -5.0;

    double duty_raw = state->kp_current * i_err
                       + state->ki_current * state->i_err_integral;

    /* Feedforward: d_ff = 1 - v_in_abs / v_out (for boost) */
    double d_ff = 1.0 - (v_in_abs / (state->v_out + PFC_EPS));
    double duty = duty_raw + d_ff;

    /* Clamp duty cycle */
    if (duty > 0.95) duty = 0.95;
    if (duty < 0.0)  duty = 0.0;

    state->duty_cycle = duty;

    /* Update output voltage model (simplified):
     * dv_out/dt = (1-d)×I_L/C_out - P_out/(V_out×C_out) */
    double i_l = i_in_meas;
    double di_load = state->p_out_nom / (state->v_out + PFC_EPS);
    double dv_out = ((1.0 - duty) * i_l / state->c_out
                     - di_load / state->c_out) * ts_sec;
    state->v_out += dv_out;

    /* Estimate PF from duty cycle (simplified CCM model) */
    state->pf_achieved = 0.99; /* base CCM PF */

    return 0;
}

/* ==========================================================================
 * L5: Boost PFC PF Estimation
 * ========================================================================== */

double pfc_boost_estimate_pf(const pfc_boost_state_t *state)
{
    if (!state) return 0.0;

    /* For an ideal CCM boost PFC, PF ≈ 0.99.
     * Degradation factors:
     *   - Crossover distortion at zero-crossing (~0.5% reduction)
     *   - Limited current loop bandwidth (~0.2%)
     *   - DCM operation at light load (significant below 20% load)
     */

    double load_factor = state->p_out_nom
                          / (state->v_out_ref * state->v_out_ref
                             / (2.0 * M_PI * state->f_line * state->l_boost));
    if (load_factor > 1.0) load_factor = 1.0;

    double pf = 0.995 * load_factor + 0.85 * (1.0 - load_factor);
    if (pf > 0.999) pf = 0.999;
    return pf;
}

/* ==========================================================================
 * L5: Boost PFC THD Estimation
 * ========================================================================== */

double pfc_estimate_thd(const pfc_boost_state_t *state, double v_in_instant,
                         double delta_i_l_pp)
{
    if (!state || state->i_ref_peak < PFC_EPS) return 0.0;

    (void)v_in_instant;

    /* THD from inductor current ripple (high-frequency component):
     * THD_ripple ≈ ΔI_pp / (2√3 × I_fundamental_rms)
     *
     * Plus crossover distortion component (~1-2% at moderate load).
     */
    double i_fund_rms = state->p_out_nom / (state->v_in_rms_nom + PFC_EPS);
    double thd_ripple = delta_i_l_pp / (2.0 * 1.73205 * i_fund_rms);
    double thd_crossover = 0.015; /* 1.5% typical */

    double thd = sqrt(thd_ripple * thd_ripple + thd_crossover * thd_crossover);
    if (thd > 0.5) thd = 0.5;
    return thd;
}

/* ==========================================================================
 * L6: Industrial PF Correction Solution
 * ========================================================================== */

int pfc_solve_industrial(double p_load, double v_ll, double f_hz,
                          double pf_old, double pf_target, double s_sc_kva,
                          double *c_out, double *kvar_out,
                          double *payback_mon, double *f_res_out)
{
    if (!c_out || !kvar_out || !payback_mon || !f_res_out) return -1;
    if (p_load <= 0.0 || v_ll < 1.0 || f_hz < 1.0) return -1;
    if (pf_old <= 0.0 || pf_target <= 0.0 || pf_old > 1.0 || pf_target > 1.0)
        return -1;

    /* Step 1: Compute required kVAR */
    double phi_old = acos(pf_old);
    double phi_new = acos(pf_target);
    double q_c = p_load * (tan(phi_old) - tan(phi_new));
    *kvar_out = q_c;

    /* Step 2: Compute delta-connected capacitance */
    double omega = 2.0 * M_PI * f_hz;
    *c_out = q_c / (3.0 * omega * v_ll * v_ll);

    /* Step 3: Compute resonant frequency */
    *f_res_out = f_hz * sqrt((s_sc_kva * 1000.0) / q_c);
    if (q_c < PFC_EPS) *f_res_out = 0.0;

    /* Step 4: Compute payback period
     * Assume: capacitor cost = $30/kVAR installed
     *         PF penalty = $0.75/kVAR/month below 0.85
     *         Loss savings ≈ 2% of kVAR reduction
     */
    double kvar_billable = p_load * tan(acos(0.85));
    if (kvar_billable < 0.0) kvar_billable = 0.0;
    double penalty_saved_monthly = kvar_billable * 0.75;
    if (penalty_saved_monthly < 0.0) penalty_saved_monthly = 0.0;

    double cap_cost = q_c * 30.0 / 1000.0; /* $30/kVAR → $/VAR */
    if (penalty_saved_monthly < PFC_EPS) {
        *payback_mon = 60.0; /* no penalty → long payback */
    } else {
        *payback_mon = cap_cost / penalty_saved_monthly;
    }

    return 0;
}

/* ==========================================================================
 * L6: Capacitor Step Switching Optimization
 * ========================================================================== */

int pfc_optimize_steps(const pfc_capacitor_bank_t *bank,
                        double q_target_var, double q_measured_var,
                        double hysteresis_var, uint32_t *new_mask)
{
    if (!bank || !new_mask) return -1;

    double q_error = q_measured_var - q_target_var;
    uint32_t current = bank->active_mask;
    uint32_t next = current;
    int changes = 0;

    /* Only act if outside hysteresis band */
    if (fabs(q_error) < hysteresis_var) {
        *new_mask = current;
        return 0;
    }

    if (q_error > hysteresis_var) {
        /* Too much reactive power → remove capacitive steps */
        for (uint32_t i = 0; i < bank->num_steps; i++) {
            if (next & (1U << i)) {
                /* Check if removing this step helps */
                double q_after = q_measured_var - bank->step_kvar[i];
                if (q_after > q_target_var + hysteresis_var) {
                    next &= ~(1U << i);
                    changes++;
                } else {
                    break; /* further removal would undershoot */
                }
            }
        }
    } else if (q_error < -hysteresis_var) {
        /* Too little reactive power → add capacitive steps */
        for (uint32_t i = 0; i < bank->num_steps; i++) {
            if (!(next & (1U << i))) {
                double q_after = q_measured_var + bank->step_kvar[i];
                if (q_after < q_target_var - hysteresis_var) {
                    next |= (1U << i);
                    changes++;
                } else {
                    break;
                }
            }
        }
    }

    *new_mask = next;
    return changes;
}

/* ==========================================================================
 * L6: Harmonic Resonance Risk with Capacitor Bank
 * ========================================================================== */

int pfc_harmonic_risk(double s_sc_va, double q_c_var, double f_fund_hz,
                       int h_order, double *magnification)
{
    if (!magnification) return -1;
    if (s_sc_va <= 0.0 || q_c_var <= 0.0 || f_fund_hz <= 0.0) return -1;
    if (h_order < 2) return -1;

    /* Parallel resonance magnification factor:
     *
     * At harmonic h, the system impedance at the PCC with capacitor is:
     *   Z_h = (j h X_s) ∥ (-j X_c / h)
     *
     * where X_s = V²/S_sc, X_c = V²/Q_c
     *
     * Magnification at resonance: M ≈ Q_c / (h S_sc) at h = h_res
     * For h ≠ h_res: M = |Z_h| / |h X_s|
     *
     * The magnification factor tells us how much the harmonic voltage
     * at order h will be amplified by the capacitor bank.
     */
    double h = (double)h_order;
    double f_res = f_fund_hz * sqrt(s_sc_va / q_c_var);
    double h_res = f_res / f_fund_hz;

    /* Simplified magnification factor based on proximity to resonance */
    double h_ratio = h / h_res;
    double denom = fabs(1.0 - h_ratio * h_ratio);
    if (denom < 0.01) denom = 0.01; /* clamp near resonance */

    *magnification = 1.0 / denom;

    /* Risk if magnification > 2 (voltage harmonic amplified > 2×) */
    if (*magnification > 2.0) return 1;
    return 0;
}
