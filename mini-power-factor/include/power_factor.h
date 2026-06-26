/**
 * @file power_factor.h
 * @brief Core definitions and types for AC power factor analysis
 *
 * Power factor (PF) is the ratio of real power to apparent power in an AC circuit.
 * It quantifies how effectively electrical power is converted into useful work.
 *
 * Key formula:  PF = P / |S| = cos(φ)   (for sinusoidal waveforms)
 *              PF_true = PF_displacement × PF_distortion  (for non-sinusoidal)
 *
 * References:
 *   - IEEE Std 1459-2010, "Definitions for Power Quantities"
 *   - C.P. Steinmetz, "Theory and Calculation of Alternating Current Phenomena" (1897)
 *   - W. Shepherd, P. Zand, "Energy Flow and Power Factor" (1979)
 *   - MIT 6.061 / Berkeley EE105 / Tsinghua Power Electronics
 *
 * Knowledge coverage:
 *   L1 (Definitions): Real/Reactive/Apparent Power, PF, displacement PF, distortion PF
 *   L2 (Concepts):   Power triangle, leading/lagging, energy storage in reactance
 *   L3 (Math):       Complex power S = P + jQ, phasor arithmetic
 *   L4 (Laws):       Conservation of real power, Boucherot's theorem
 */

#ifndef POWER_FACTOR_H
#define POWER_FACTOR_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Core Power Type Definitions
 * ========================================================================== */

#define PF_MAX_HARMONICS    64
#define PF_MAX_NAME_LEN     64
#define PF_EPSILON          1e-12

/**
 * @brief Power factor type classification.
 *
 * PF can be leading (capacitive load, current leads voltage) or
 * lagging (inductive load, current lags voltage). Unity PF means
 * purely resistive load.
 */
typedef enum {
    PF_TYPE_UNITY     = 0,   /**< PF = 1.0, purely resistive */
    PF_TYPE_LAGGING   = 1,   /**< Inductive: current lags voltage (φ > 0) */
    PF_TYPE_LEADING   = -1,  /**< Capacitive: current leads voltage (φ < 0) */
    PF_TYPE_UNDEFINED = 99   /**< Cannot determine (zero power) */
} pf_type_t;

/**
 * @brief Core power quantities for single-phase sinusoidal AC.
 *
 * Real Power P:     Average rate of energy transfer [Watts]
 * Reactive Power Q:  Rate of energy oscillation between source and load [VAR]
 * Apparent Power S:  Product of RMS voltage and current [VA]
 *
 * Power Triangle: S² = P² + Q²
 *
 * Theta [theta_p] is the impedance angle: θ = arctan(Q/P)
 * For sinusoidal: P = Vrms * Irms * cos(θv - θi)
 *                 Q = Vrms * Irms * sin(θv - θi)
 *                 S = Vrms * Irms
 */
typedef struct {
    double   p_real;          /**< Real power P [Watts] — average (heating) power */
    double   q_reactive;      /**< Reactive power Q [VAR] — exchanged energy per cycle */
    double   s_apparent;      /**< Apparent power S [VA] — product Vrms × Irms */
    double   pf;              /**< Power factor = P / S  (range [0, 1]) */
    double   phi_rad;         /**< Phase angle φ between voltage and current [radians] */
    double   phi_deg;         /**< Phase angle φ between voltage and current [degrees] */
    pf_type_t type;           /**< Leading, lagging, or unity */
    double   v_rms;           /**< RMS voltage [V] */
    double   i_rms;           /**< RMS current [A] */
} pf_single_phase_t;

/**
 * @brief Three-phase power quantities.
 *
 * For balanced three-phase systems with line-to-line voltage V_LL
 * and line current I_L:
 *   P_3φ = √3 × V_LL × I_L × cos(φ)
 *   Q_3φ = √3 × V_LL × I_L × sin(φ)
 *   S_3φ = √3 × V_LL × I_L
 *
 * For unbalanced systems, total power is the sum of per-phase powers:
 *   P_total = P_a + P_b + P_c
 */
typedef struct {
    double   p_total;          /**< Total three-phase real power [W] */
    double   q_total;          /**< Total three-phase reactive power [VAR] */
    double   s_total;          /**< Total three-phase apparent power [VA] */
    double   p_a, p_b, p_c;    /**< Per-phase real power [W] */
    double   q_a, q_b, q_c;    /**< Per-phase reactive power [VAR] */
    double   v_ll_rms;         /**< Line-to-line RMS voltage [V] */
    double   i_l_rms;          /**< Line RMS current [A] */
    double   pf_total;         /**< Total power factor */
    double   pf_a, pf_b, pf_c; /**< Per-phase power factor */
    double   phi_deg;          /**< Phase angle [degrees] */
    pf_type_t type;            /**< Leading/lagging per total */
    uint8_t  is_balanced;      /**< Flag: phases are balanced */
    double   voltage_unbalance; /**< NEMA voltage unbalance % */
    double   current_unbalance; /**< NEMA current unbalance % */
} pf_three_phase_t;

/**
 * @brief Harmonic-level power data for distortion analysis.
 *
 * The spectral components of voltage and current beyond the fundamental
 * give rise to distortion power D (also called harmonic reactive power).
 *
 * From IEEE 1459-2010:
 *   V_rms² = V₁² + V_H²  where V_H² = Σ V_h² for h > 1
 *   I_rms² = I₁² + I_H²
 *   S² = P₁² + Q₁² + D_H² + D_V² + S_H²
 */
typedef struct {
    uint32_t harmonic_order;        /**< Harmonic number h (1 = fundamental) */
    double   v_magnitude;           /**< Voltage magnitude at harmonic h [V] */
    double   v_phase_deg;           /**< Voltage phase at harmonic h [degrees] */
    double   i_magnitude;           /**< Current magnitude at harmonic h [A] */
    double   i_phase_deg;           /**< Current phase at harmonic h [degrees] */
    double   p_h;                   /**< Real power at harmonic h [W] */
    double   q_h;                   /**< Reactive power at harmonic h [VAR] */
} pf_harmonic_component_t;

/**
 * @brief Complete non-sinusoidal power analysis per IEEE 1459-2010.
 *
 * Fundamental power quantities (subscript 1):
 *   P₁ = V₁ × I₁ × cos(φ₁)   — fundamental real power
 *   Q₁ = V₁ × I₁ × sin(φ₁)   — fundamental reactive power
 *   S₁ = V₁ × I₁             — fundamental apparent power
 *
 * Non-fundamental distortion:
 *   D_I = V₁ × I_H   — current distortion power
 *   D_V = V_H × I₁   — voltage distortion power
 *   S_H = V_H × I_H   — harmonic apparent power
 *
 * Total non-fundamental apparent power:
 *   S_N² = D_I² + D_V² + S_H²
 */
typedef struct {
    /* Fundamental (60/50 Hz) quantities */
    double   p1;                    /**< Fundamental real power [W] */
    double   q1;                    /**< Fundamental reactive power [VAR] */
    double   s1;                    /**< Fundamental apparent power [VA] */
    double   pf1;                   /**< Displacement power factor P₁/S₁ */
    double   phi1_rad;              /**< Fundamental phase angle [rad] */
    double   phi1_deg;              /**< Fundamental phase angle [deg] */
    pf_type_t dcp_type;             /**< Displacement PF type: leading/lagging */

    /* Total quantities (including all harmonics) */
    double   p_total;               /**< Total real power P = Σ P_h [W] */
    double   s_total;               /**< Total apparent power S = V_rms × I_rms [VA] */
    double   pf_true;               /**< True power factor = P_total / S_total */

    /* Distortion decomposition (IEEE 1459-2010) */
    double   v_rms;                 /**< Total RMS voltage */
    double   v1_rms;                /**< Fundamental RMS voltage */
    double   v_h_rms;               /**< Harmonic RMS voltage (h > 1) */
    double   i_rms;                 /**< Total RMS current */
    double   i1_rms;                /**< Fundamental RMS current */
    double   i_h_rms;               /**< Harmonic RMS current (h > 1) */
    double   d_i;                   /**< Current distortion power V₁ × I_H [var] */
    double   d_v;                   /**< Voltage distortion power V_H × I₁ [var] */
    double   s_h;                   /**< Harmonic apparent power V_H × I_H [VA] */
    double   s_n;                   /**< Non-fundamental apparent power [VA] */
    double   n_active;              /**< Non-active power N [var] (IEEE 1459) */

    /* Distortion metrics */
    double   thd_v_percent;         /**< Voltage THD [%]: (V_H/V₁) × 100 */
    double   thd_i_percent;         /**< Current THD [%]: (I_H/I₁) × 100 */
    double   pf_distortion;         /**< Distortion power factor = 1/√((1+THDv²)(1+THDi²)) */

    /* Harmonic components */
    uint32_t num_harmonics;         /**< Number of harmonic components stored */
    pf_harmonic_component_t harmonics[PF_MAX_HARMONICS];
} pf_ieee1459_t;

/**
 * @brief Power consumption time-series record for monitoring.
 *
 * Captures power measurements at regular intervals for energy auditing,
 * load profiling, and demand analysis. Used in smart grid applications.
 */
typedef struct {
    uint64_t timestamp_sec;         /**< Unix timestamp of measurement */
    double   p_instant;             /**< Instantaneous real power [W] */
    double   q_instant;             /**< Instantaneous reactive power [VAR] */
    double   s_instant;             /**< Instantaneous apparent power [VA] */
    double   pf_instant;            /**< Instantaneous power factor */
    double   v_rms;                 /**< RMS voltage at this sample */
    double   i_rms;                 /**< RMS current at this sample */
    double   frequency_hz;          /**< Measured line frequency */
    double   energy_wh;             /**< Cumulative energy [Watt-hours] */
    double   demand_w;              /**< Rolling demand [W] (typically 15-min) */
} pf_measurement_t;

/* ==========================================================================
 * L1: Core PF Calculation Functions
 * ========================================================================== */

/**
 * @brief Compute single-phase power quantities from RMS values and phase angle.
 *
 * @param v_rms    RMS voltage [V]
 * @param i_rms    RMS current [A]
 * @param phi_deg  Phase angle between V and I [degrees], positive=lagging
 * @param result   [out] Computed power quantities
 * @return 0 on success, -1 on invalid input
 *
 * Formula: P = V×I×cos(φ), Q = V×I×sin(φ), S = V×I, PF = cos(φ)
 */
int pf_compute_single_phase(double v_rms, double i_rms, double phi_deg,
                            pf_single_phase_t *result);

/**
 * @brief Compute PF from time-domain voltage and current samples.
 *
 * Uses sliding-window RMS and cross-correlation to determine
 * real power and phase angle from sampled waveforms.
 *
 * @param v_samples  Array of voltage samples [V]
 * @param i_samples  Array of current samples [A]
 * @param n_samples  Number of samples (must be ≥ 2)
 * @param dt_sec     Sampling interval [seconds]
 * @param result     [out] Computed power quantities
 * @return 0 on success, -1 on invalid input
 *
 * Algorithm:
 *   V_rms = sqrt( (1/N) Σ v[n]² )
 *   I_rms = sqrt( (1/N) Σ i[n]² )
 *   P = (1/N) Σ v[n] × i[n]         — average instantaneous power
 *   S = V_rms × I_rms
 *   Q = sqrt(S² - P²)                — sign from phase detection
 *   PF = P / S
 */
int pf_compute_from_samples(const double *v_samples, const double *i_samples,
                            size_t n_samples, double dt_sec,
                            pf_single_phase_t *result);

/**
 * @brief Determine PF type (leading/lagging) from cross-correlation analysis.
 *
 * Computes the phase relationship by finding the lag at which
 * cross-correlation between V and I is maximized.
 *
 * @param v_samples  Voltage waveform samples
 * @param i_samples  Current waveform samples
 * @param n_samples  Number of samples
 * @param phase_lag_deg  [out] Phase lag in degrees (positive = I lags V)
 * @return 0 on success, -1 on error
 */
int pf_detect_phase_lag(const double *v_samples, const double *i_samples,
                        size_t n_samples, double *phase_lag_deg);

/* ==========================================================================
 * L2: Power Factor Classification and Properties
 * ========================================================================== */

/**
 * @brief Classify power factor quality per utility standards.
 *
 * | PF Range         | Class    | Typical Load                |
 * |------------------|----------|-----------------------------|
 * | 0.95 - 1.00      | Good     | Resistive heaters, PFC PSUs |
 * | 0.85 - 0.95      | Fair     | Induction motors (loaded)   |
 * | 0.70 - 0.85      | Poor     | Lightly loaded motors       |
 * | 0.00 - 0.70      | Bad      | Arc furnaces, welders       |
 */
typedef enum {
    PF_CLASS_GOOD      = 0,   /**< PF ≥ 0.95 */
    PF_CLASS_FAIR      = 1,   /**< 0.85 ≤ PF < 0.95 */
    PF_CLASS_POOR      = 2,   /**< 0.70 ≤ PF < 0.85 */
    PF_CLASS_BAD       = 3,   /**< PF < 0.70 */
    PF_CLASS_LEADING   = 4,   /**< Leading PF (over-correction warning) */
} pf_class_t;

pf_class_t pf_classify(double pf, pf_type_t type);

/**
 * @brief Compute the reactive power needed to achieve a target PF.
 *
 * Given current real power P and current/desired power factors,
 * computes Q_c = P × (tan(acos(pf_old)) - tan(acos(pf_new)))
 *
 * @param p         Real power [W]
 * @param pf_old    Current power factor [0-1]
 * @param pf_target Target power factor [0-1]
 * @param lagging   Non-zero if original load is lagging (inductive)
 * @return Required reactive power compensation [VAR]. Positive = need capacitive.
 */
double pf_required_reactive_comp(double p, double pf_old, double pf_target,
                                 int lagging);

/**
 * @brief Compute the economic benefit of PF correction.
 *
 * Calculates annual savings from reduced I²R losses and/or
 * avoided utility penalty charges for low PF.
 *
 * @param p_load         Load real power [W]
 * @param pf_old         Original power factor
 * @param pf_new         Corrected power factor
 * @param line_r         Line resistance [Ohm]
 * @param hours_per_year Operating hours
 * @param cost_per_kwh   Electricity cost [$ per kWh]
 * @return Estimated annual savings [$]
 */
double pf_savings_estimate(double p_load, double pf_old, double pf_new,
                           double line_r, double hours_per_year,
                           double cost_per_kwh);

/* ==========================================================================
 * L3: Phasor-Based Power Computation
 * ========================================================================== */

/**
 * @brief Compute complex power S = V × I* from voltage and current phasors.
 *
 * S = P + jQ = V_rms ∠θv × (I_rms ∠θi)* = V_rms × I_rms ∠(θv - θi)
 *
 * @param v_mag_rms   RMS voltage magnitude [V]
 * @param v_ang_rad   Voltage phase angle [rad]
 * @param i_mag_rms   RMS current magnitude [A]
 * @param i_ang_rad   Current phase angle [rad]
 * @param result      [out] Computed single-phase power
 * @return 0 on success
 */
int pf_phasor_power(double v_mag_rms, double v_ang_rad,
                    double i_mag_rms, double i_ang_rad,
                    pf_single_phase_t *result);

/**
 * @brief Compute three-phase power from line quantities.
 *
 * For balanced three-phase:
 *   P_3φ = √3 × V_LL × I_L × cos(φ)
 *
 * For unbalanced, sums per-phase with symmetrical components.
 *
 * @param v_ll_rms   Line-to-line RMS voltage [V]
 * @param i_l_rms    Line RMS current [A]
 * @param phi_deg    Phase angle [degrees]
 * @param balanced   Non-zero if system is balanced
 * @param result     [out] Three-phase power quantities
 * @return 0 on success
 */
int pf_compute_three_phase(double v_ll_rms, double i_l_rms, double phi_deg,
                            int balanced, pf_three_phase_t *result);

/**
 * @brief Symmetrical components decomposition (Fortescue transform).
 *
 * Decomposes three unbalanced phasors into positive, negative, and
 * zero sequence components using the Fortescue transformation (1918):
 *
 *   [V0]    [1   1   1 ] [Va]
 *   [V1] =  [1   a   a²] [Vb]  where a = 1∠120° = e^{j2π/3}
 *   [V2]    [1   a²  a ] [Vc]
 *
 * @param va_mag,va_ang  Phase A magnitude/angle (voltage or current)
 * @param vb_mag,vb_ang  Phase B magnitude/angle
 * @param vc_mag,vc_ang  Phase C magnitude/angle
 * @param v0_mag,v0_ang  [out] Zero sequence
 * @param v1_mag,v1_ang  [out] Positive sequence
 * @param v2_mag,v2_ang  [out] Negative sequence
 * @return 0 on success
 */
int pf_symmetrical_components(double va_mag, double va_ang,
                               double vb_mag, double vb_ang,
                               double vc_mag, double vc_ang,
                               double *v0_mag, double *v0_ang,
                               double *v1_mag, double *v1_ang,
                               double *v2_mag, double *v2_ang);

/**
 * @brief Compute voltage unbalance per NEMA MG1 definition.
 *
 * NEMA unbalance = max deviation from average / average × 100%
 *
 * @param v_ab, v_bc, v_ca  Line-to-line voltages [V]
 * @return Voltage unbalance percentage
 */
double pf_voltage_unbalance_percent(double v_ab, double v_bc, double v_ca);

/* ==========================================================================
 * L4: Energy Conservation and Power Theorems
 * ========================================================================== */

/**
 * @brief Verify real power balance in a circuit node/network.
 *
 * Per conservation of energy: Σ P_in = Σ P_out + Σ P_losses
 *
 * @param p_sources    Array of source powers [W] (positive = supplying)
 * @param n_sources    Number of sources
 * @param p_loads      Array of load powers [W] (positive = consuming)
 * @param n_loads      Number of loads
 * @param tolerance    Allowed numerical tolerance
 * @return 0 if balanced, 1 if imbalance detected, -1 on error
 */
int pf_verify_power_balance(const double *p_sources, size_t n_sources,
                            const double *p_loads, size_t n_loads,
                            double tolerance);

/**
 * @brief Verify Boucherot's Theorem: total reactive power is conserved.
 *
 * Boucherot's theorem (1896): In any network with sinusoidal sources
 * at the same frequency, the sum of reactive powers over all branches
 * is zero: Σ Q_k = 0.
 *
 * This holds because reactive power is associated with energy storage
 * in L and C, not energy consumption.
 *
 * @param q_values   Array of reactive powers [VAR] for all branches
 * @param n_branches Number of branches
 * @param tolerance  Numerical tolerance
 * @return 0 if theorem holds, 1 if violation detected
 */
int pf_verify_boucherot(const double *q_values, size_t n_branches,
                         double tolerance);

/* ==========================================================================
 * L5: RMS and Statistical Methods
 * ========================================================================== */

/**
 * @brief Compute true RMS from sampled data using sliding window.
 *
 * Implements the standard RMS computation:
 *   V_rms[n] = sqrt( (1/W) Σ_{k=n-W+1}^{n} v[k]² )
 *
 * Uses a numerically stable Welford-style accumulator to avoid
 * catastrophic cancellation.
 *
 * @param samples      Input sample array
 * @param n_samples    Total number of samples
 * @param window_size  Number of samples per window
 * @param rms_out      [out] Array of RMS values (length = n_samples - window_size + 1)
 * @return 0 on success, -1 on error
 */
int pf_sliding_rms(const double *samples, size_t n_samples,
                   size_t window_size, double *rms_out);

/**
 * @brief Compute exponential moving average (EMA) RMS for real-time PF tracking.
 *
 * Uses infinite impulse response:  y[n] = α·x[n]² + (1-α)·y[n-1]
 * where α = 2/(tau_samples + 1)
 *
 * @param samples    Input sample array
 * @param n_samples  Number of samples
 * @param tau        Time constant in samples (e.g., 10 for fast, 100 for slow)
 * @param rms_out    [out] Computed EMA-RMS array
 * @return 0 on success
 */
int pf_ema_rms(const double *samples, size_t n_samples,
               size_t tau, double *rms_out);

/**
 * @brief Compute crest factor (ratio of peak to RMS).
 *
 * CF = V_peak / V_rms
 *
 * For a pure sine wave: CF = √2 ≈ 1.414
 * High CF indicates waveform distortion or impulsive loads.
 *
 * @param samples    Voltage/current samples
 * @param n_samples  Number of samples
 * @return Crest factor, or -1 on error
 */
double pf_crest_factor(const double *samples, size_t n_samples);

/**
 * @brief Compute form factor (ratio of RMS to average of absolute value).
 *
 * FF = V_rms / V_avg_rectified
 *
 * For a pure sine wave: FF = π/(2√2) ≈ 1.111
 * For a square wave: FF = 1.0
 * For a triangle wave: FF = 2/√3 ≈ 1.155
 *
 * Deviations from 1.111 indicate harmonic content.
 *
 * @param samples    Voltage/current samples
 * @param n_samples  Number of samples
 * @return Form factor, or -1 on error
 */
double pf_form_factor(const double *samples, size_t n_samples);

/**
 * @brief Compute demand interval from measurement history.
 *
 * Rolling demand over a specified interval, per utility metering practice.
 * Typically 15-minute or 30-minute sliding window.
 *
 * @param measurements Array of power measurements
 * @param n_meas       Number of measurements
 * @param interval_sec Demand interval in seconds
 * @param demand_out   [out] Computed demand array
 * @return 0 on success, -1 on error
 */
int pf_demand_interval(const pf_measurement_t *measurements, size_t n_meas,
                       double interval_sec, double *demand_out);

#ifdef __cplusplus
}
#endif

#endif /* POWER_FACTOR_H */
