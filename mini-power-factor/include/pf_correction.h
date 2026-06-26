/**
 * @file pf_correction.h
 * @brief Power factor correction algorithms, capacitor sizing, and filter design
 *
 * Power factor correction (PFC) is the process of adjusting the reactive power
 * balance in an AC system to bring the power factor closer to unity. This
 * improves energy efficiency, reduces I²R losses, and avoids utility penalties.
 *
 * Passive PFC: capacitor banks, tuned filters, synchronous condensers
 * Active PFC:  switched-mode converters (boost PFC, bridgeless PFC)
 *
 * References:
 *   - Erickson & Maksimovic, "Fundamentals of Power Electronics" (2nd ed, 2001) Ch.18
 *   - IEEE Std 18-2012, "Shunt Power Capacitors"
 *   - IEC 61000-3-2, "Limits for Harmonic Current Emissions"
 *   - Dixon, "Average Current Mode Control of Switching Power Supplies" (Unitrode AN)
 *
 * Knowledge coverage:
 *   L1 (Definitions): PFC, kVAR, capacitor bank, detuning reactor
 *   L2 (Concepts):   Shunt compensation, series compensation, resonance avoidance
 *   L5 (Algorithms): Capacitor sizing, step control, active PFC control loop
 *   L6 (Problems):   Industrial PF correction, harmonic resonance avoidance
 */

#ifndef PF_CORRECTION_H
#define PF_CORRECTION_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: PFC System Types and Configurations
 * ========================================================================== */

/**
 * @brief Types of power factor correction systems.
 */
typedef enum {
    PFC_TYPE_NONE             = 0,  /**< No PFC */
    PFC_TYPE_FIXED_CAPACITOR  = 1,  /**< Fixed capacitor bank */
    PFC_TYPE_AUTO_CAPACITOR   = 2,  /**< Automatic (stepped) capacitor bank */
    PFC_TYPE_SYNCH_CONDENSER  = 3,  /**< Synchronous condenser */
    PFC_TYPE_STATIC_VAR       = 4,  /**< Static VAR compensator (SVC) */
    PFC_TYPE_STATCOM          = 5,  /**< Static synchronous compensator */
    PFC_TYPE_ACTIVE_PFC       = 6,  /**< Active PFC (boost converter) */
    PFC_TYPE_BRIDGELESS_PFC   = 7,  /**< Bridgeless totem-pole PFC */
    PFC_TYPE_HYBRID           = 8,  /**< Hybrid active/passive */
} pfc_type_t;

/**
 * @brief Capacitor bank step configuration.
 *
 * Automatic PF correction systems use multiple capacitor steps
 * switched by contactors or thyristors to achieve variable
 * reactive power compensation.
 */
typedef struct {
    uint32_t num_steps;              /**< Number of capacitor steps */
    double   step_kvar[16];          /**< kVAR rating of each step */
    double   total_kvar;             /**< Total installed kVAR */
    double   v_rated;                /**< Rated voltage [V] */
    double   f_rated;                /**< Rated frequency [Hz] */
    double   detune_percent;         /**< Detuning reactor percentage (e.g., 7% for 189 Hz) */
    uint8_t  is_detuned;             /**< Flag: anti-resonance detuning applied */
    uint8_t  step_delay_sec;         /**< Delay between step switching [seconds] */
    double   discharge_time_sec;     /**< Capacitor discharge time [seconds] */
    uint32_t active_mask;            /**< Bitmask of currently active steps */
} pfc_capacitor_bank_t;

/**
 * @brief Active PFC controller state for boost converter PFC.
 *
 * Implements average current mode control for a boost PFC converter.
 * The control objective: shape the input current to follow the
 * rectified sinusoidal voltage, achieving near-unity PF.
 *
 * Outer voltage loop: regulates DC output voltage V_out
 * Inner current loop: shapes inductor current to match V_in waveform
 */
typedef struct {
    /* Power stage parameters */
    double   l_boost;            /**< Boost inductor [Henries] */
    double   c_out;              /**< Output capacitor [Farads] */
    double   v_in_rms_nom;       /**< Nominal input RMS voltage [V] */
    double   v_out_ref;          /**< Reference output voltage [V] */
    double   p_out_nom;          /**< Nominal output power [W] */
    double   f_sw;               /**< Switching frequency [Hz] */
    double   f_line;             /**< Line frequency [Hz] */

    /* Controller gains */
    double   kp_voltage;         /**< Voltage loop proportional gain */
    double   ki_voltage;         /**< Voltage loop integral gain */
    double   kp_current;         /**< Current loop proportional gain */
    double   ki_current;         /**< Current loop integral gain */

    /* State variables */
    double   v_out;              /**< Actual output voltage [V] */
    double   i_l;                /**< Inductor current [A] */
    double   duty_cycle;         /**< Current duty cycle [0-1] */
    double   v_err_integral;     /**< Voltage error integral */
    double   i_err_integral;     /**< Current error integral */
    double   i_ref_peak;         /**< Peak current reference */

    /* Derived metrics */
    double   pf_achieved;        /**< Achieved power factor */
    double   thd_i_percent;      /**< Input current THD [%] */
    double   efficiency;         /**< Converter efficiency */
} pfc_boost_state_t;

/* ==========================================================================
 * L5: Capacitor Bank Sizing Algorithms
 * ========================================================================== */

/**
 * @brief Size capacitor bank for single-phase PF correction.
 *
 * Computes required capacitance to correct from pf_old to pf_target.
 *
 * Q_c = P × (tan φ_old - tan φ_new)
 * C   = Q_c / (2πf V²)
 *
 * For delta-connected three-phase capacitors:
 *   C_phase = Q_c / (6πf V_LL²)
 *
 * For wye-connected three-phase capacitors:
 *   C_phase = Q_c / (2πf V_LL²)
 *
 * @param p_load       Load real power [W]
 * @param v_rms        RMS voltage [V]
 * @param f_hz         System frequency [Hz]
 * @param pf_old       Original power factor [0-1]
 * @param pf_target    Target power factor [0-1, typically 0.95]
 * @param is_lagging   Non-zero for lagging (inductive) load
 * @param capacitance  [out] Required capacitance [Farads]
 * @param kvar_rating  [out] Capacitor reactive power rating [VAR]
 * @return 0 on success, -1 on invalid input
 */
int pfc_size_capacitor(double p_load, double v_rms, double f_hz,
                       double pf_old, double pf_target, int is_lagging,
                       double *capacitance, double *kvar_rating);

/**
 * @brief Size three-phase capacitor bank.
 *
 * Supports both delta and wye configurations. Delta is standard for
 * low-voltage (<600V) applications because capacitors see line-to-line
 * voltage and provide 3× the reactive power of wye.
 *
 * @param p_3phase      Total three-phase real power [W]
 * @param v_ll_rms      Line-to-line RMS voltage [V]
 * @param f_hz          System frequency [Hz]
 * @param pf_old        Original power factor
 * @param pf_target     Target power factor
 * @param is_lagging    Non-zero for inductive load
 * @param connection    0=delta, 1=wye
 * @param c_per_phase   [out] Capacitance per phase [F]
 * @param kvar_total    [out] Total three-phase kVAR
 * @return 0 on success, -1 on error
 */
int pfc_size_capacitor_3phase(double p_3phase, double v_ll_rms, double f_hz,
                               double pf_old, double pf_target, int is_lagging,
                               int connection, double *c_per_phase,
                               double *kvar_total);

/**
 * @brief Design an automatic capacitor bank step configuration.
 *
 * Creates a binary-weighted or equal-sized step configuration
 * to cover the required kVAR range with minimal steps.
 *
 * @param kvar_required   Total kVAR needed for full correction
 * @param num_steps       Desired number of steps
 * @param resolution_kvar Minimum kVAR per step (controls granularity)
 * @param bank            [out] Configured capacitor bank
 * @return 0 on success, -1 on error
 */
int pfc_design_step_bank(double kvar_required, uint32_t num_steps,
                          double resolution_kvar, pfc_capacitor_bank_t *bank);

/**
 * @brief Calculate the reactive power contribution of a capacitor bank step.
 *
 * Q_step = V² × 2πfC × (1 / (1 - k_detune))
 *
 * where k_detune is the detuning factor (e.g., 0.07 for 7% reactor).
 * Detuning shifts the LC resonance below the lowest characteristic harmonic.
 *
 * @param c_farads      Capacitance of this step [F]
 * @param v_rms         RMS voltage [V]
 * @param f_hz          System frequency [Hz]
 * @param detune_factor Detuning reactor fraction (0 = no detuning)
 * @return Reactive power [VAR] contributed by this step
 */
double pfc_step_kvar(double c_farads, double v_rms, double f_hz,
                     double detune_factor);

/**
 * @brief Predict resonant frequency of capacitor bank with system inductance.
 *
 * When a capacitor bank is connected to a system with short-circuit
 * inductance L_sys, parallel resonance occurs at:
 *
 *   f_res = f_sys × sqrt(S_sc / Q_c)
 *
 * where S_sc = system short-circuit power [VA]
 *       Q_c  = capacitor bank reactive power [VAR]
 *
 * Harmonics close to f_res will be amplified — detuning is required
 * if f_res is near an integer multiple of f_sys (3rd, 5th, 7th, etc.).
 *
 * @param s_sc_va       System short-circuit power [VA]
 * @param q_c_var       Capacitor bank reactive power [VAR]
 * @param f_sys_hz      System fundamental frequency [Hz]
 * @return Resonant frequency [Hz]
 */
double pfc_resonant_frequency(double s_sc_va, double q_c_var, double f_sys_hz);

/**
 * @brief Determine if detuning reactor is needed.
 *
 * A detuning reactor is recommended if the resonant frequency falls
 * near a characteristic harmonic. Common harmonic orders to avoid:
 * 3rd, 5th, 7th, 11th, 13th (for 6-pulse rectifiers).
 *
 * @param f_res_hz        Resonant frequency [Hz]
 * @param f_sys_hz        System frequency [Hz]
 * @param harmonic_order  [out] Nearest problematic harmonic order
 * @return 1 if detuning recommended, 0 if safe, -1 on error
 */
int pfc_needs_detuning(double f_res_hz, double f_sys_hz,
                        int *harmonic_order);

/* ==========================================================================
 * L5: Active PFC Control Algorithms
 * ========================================================================== */

/**
 * @brief Initialize a boost PFC controller state.
 *
 * Sets up the power stage parameters and initializes controller
 * gains based on desired bandwidth and phase margin.
 *
 * Voltage loop bandwidth: typically 10-20 Hz (well below 2×f_line)
 * Current loop bandwidth: typically 1-10 kHz (well below f_sw/2)
 *
 * @param state             [out] Initialized PFC state
 * @param v_in_rms          Nominal input voltage [V]
 * @param v_out             Target output voltage [V]
 * @param p_out             Nominal power [W]
 * @param f_sw              Switching frequency [Hz]
 * @param f_line            Line frequency [Hz]
 * @return 0 on success
 */
int pfc_boost_init(pfc_boost_state_t *state, double v_in_rms, double v_out,
                    double p_out, double f_sw, double f_line);

/**
 * @brief Run one control cycle of the boost PFC controller.
 *
 * Executes the average current mode control algorithm:
 * 1. Sample V_in (rectified sine)
 * 2. Voltage error → PI → current reference amplitude
 * 3. Current reference = amplitude × |sin(ωt)|
 * 4. Current error → PI → duty cycle
 * 5. Update state for next cycle
 *
 * @param state        PFC controller state (updated in place)
 * @param v_in_inst    Instantaneous rectified input voltage [V]
 * @param i_in_meas    Measured inductor current [A]
 * @param ts_sec       Sampling time [seconds]
 * @return 0 on success
 */
int pfc_boost_control_step(pfc_boost_state_t *state, double v_in_inst,
                            double i_in_meas, double ts_sec);

/**
 * @brief Estimate the power factor achieved by a boost PFC.
 *
 * For CCM boost PFC with average current mode control:
 * PF ≈ 1 / sqrt(1 + THD_i²)
 *
 * THD_i is affected by: crossover distortion, light load DCM operation,
 * input voltage distortion, and limited control bandwidth.
 *
 * @param state   PFC state
 * @return Estimated power factor [0-1]
 */
double pfc_boost_estimate_pf(const pfc_boost_state_t *state);

/**
 * @brief Compute the input current THD for a PFC converter.
 *
 * Simplified model based on inductor current ripple and
 * crossover distortion around zero-crossing.
 *
 * @param state            PFC state
 * @param v_in_instant     Instantaneous input voltage [V]
 * @param delta_i_l_pp     Peak-to-peak inductor current ripple [A]
 * @return Estimated THD [0-1] (not percentage)
 */
double pfc_estimate_thd(const pfc_boost_state_t *state, double v_in_instant,
                         double delta_i_l_pp);

/* ==========================================================================
 * L6: Canonical PF Correction Problems
 * ========================================================================== */

/**
 * @brief Solve industrial PF correction: size capacitor bank and verify savings.
 *
 * Typical industrial scenario:
 * - 480V, 60Hz, 500 kW motor load at PF=0.70 lagging
 * - Target PF = 0.95
 * - Compute: kVAR needed, capacitor value, payback period, resonance check
 *
 * @param p_load       Real power [W]
 * @param v_ll         Line-to-line voltage [V]
 * @param f_hz         Frequency [Hz]
 * @param pf_old       Current PF
 * @param pf_target    Target PF
 * @param s_sc_kva     Short-circuit kVA at PCC
 * @param c_out        [out] Required capacitance (delta) [F]
 * @param kvar_out     [out] Required kVAR
 * @param payback_mon  [out] Simple payback [months]
 * @param f_res_out    [out] Resonant frequency [Hz]
 * @return 0 on success
 */
int pfc_solve_industrial(double p_load, double v_ll, double f_hz,
                          double pf_old, double pf_target, double s_sc_kva,
                          double *c_out, double *kvar_out,
                          double *payback_mon, double *f_res_out);

/**
 * @brief Optimize capacitor step switching sequence.
 *
 * When multiple capacitor steps are available, this function determines
 * which steps to energize/de-energize to achieve a target reactive power
 * with minimal switching operations (hysteresis control).
 *
 * @param bank            Current capacitor bank state
 * @param q_target_var    Target reactive power (negative = absorb) [VAR]
 * @param q_measured_var  Currently measured reactive power [VAR]
 * @param hysteresis_var  Hysteresis band to prevent hunting [VAR]
 * @param new_mask        [out] New step activation bitmask
 * @return Number of steps changed, -1 on error
 */
int pfc_optimize_steps(const pfc_capacitor_bank_t *bank,
                        double q_target_var, double q_measured_var,
                        double hysteresis_var, uint32_t *new_mask);

/**
 * @brief Assess risk of harmonic resonance with capacitor bank.
 *
 * Evaluates the frequency response of the system impedance at the
 * point of common coupling (PCC) with the capacitor bank installed.
 * Computes the impedance magnification factor at each harmonic.
 *
 * @param s_sc_va        System short-circuit power [VA]
 * @param q_c_var        Capacitor reactive power [VAR]
 * @param f_fund_hz      Fundamental frequency [Hz]
 * @param h_order        Harmonic order to check (e.g., 5, 7, 11)
 * @param magnification  [out] Impedance magnification at harmonic h
 * @return 0 if safe, 1 if resonance risk detected
 */
int pfc_harmonic_risk(double s_sc_va, double q_c_var, double f_fund_hz,
                       int h_order, double *magnification);

#ifdef __cplusplus
}
#endif

#endif /* PF_CORRECTION_H */
