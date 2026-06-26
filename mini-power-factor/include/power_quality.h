/**
 * @file power_quality.h
 * @brief Power quality metrics, standards, and monitoring for AC systems
 *
 * Power quality encompasses voltage stability, frequency regulation,
 * harmonic distortion, transient overvoltages, sags, swells, and
 * interruptions — all of which affect the power factor and
 * operational reliability of electrical equipment.
 *
 * References:
 *   - IEEE Std 1159-2019, "Monitoring Electric Power Quality"
 *   - IEEE Std 1459-2010, "Definitions for Power Quantities"
 *   - IEC 61000-4-30, "Power Quality Measurement Methods"
 *   - Bollen, "Understanding Power Quality Problems" (2000)
 *   - Dugan et al., "Electrical Power Systems Quality" (3rd ed, 2012)
 *   - MIT 6.691 "Seminar in Electric Power Systems"
 *
 * Knowledge coverage:
 *   L1 (Definitions): Sag, swell, flicker, unbalance, notch, transient
 *   L2 (Concepts):   Voltage tolerance (ITIC/CBEMA), power acceptability
 *   L4 (Laws):       IEEE 1159 event classification, ITIC curve
 *   L7 (Applications): Smart grid PQ monitoring, data center PQ, EV charging PQ
 */

#ifndef POWER_QUALITY_H
#define POWER_QUALITY_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PQ_MAX_EVENTS    256
#define PQ_MAX_HARMONICS 64

/* ==========================================================================
 * L1: Power Quality Event Types (IEEE 1159-2019)
 * ========================================================================== */

/**
 * @brief IEEE 1159-2019 power quality event categories.
 *
 * | Category      | Duration              | Voltage Magnitude          |
 * |---------------|-----------------------|----------------------------|
 * | Transient     | < 0.5 cycle           | 0 - 800%                   |
 * | Short-duration| 0.5 cycle - 1 min     | 0.1 - 1.8 pu               |
 * | Long-duration | > 1 min               | 0.0 - 1.2 pu               |
 * | Imbalance     | Steady-state          | Deviation from symmetry    |
 * | Waveform dist.| Steady-state          | Harmonic/interharmonic     |
 * | Flicker       | Intermittent          | 0.1 - 7% ΔV               |
 */
typedef enum {
    PQ_EVENT_NONE               = 0,
    /* Transients */
    PQ_EVENT_IMPULSE_TRANSIENT  = 1,   /**< Lightning, switching surge (<50ns-1ms) */
    PQ_EVENT_OSCILLATORY_TRANSIENT = 2, /**< Capacitor switching (0.3-50ms, 0-8pu) */
    /* Short-duration RMS variations */
    PQ_EVENT_SAG_INSTANTANEOUS  = 3,   /**< 0.5-30 cycles, 0.1-0.9 pu */
    PQ_EVENT_SWELL_INSTANTANEOUS = 4,  /**< 0.5-30 cycles, 1.1-1.8 pu */
    PQ_EVENT_INTERRUPTION_MOMENTARY = 5, /**< 0.5 cycles - 3 sec, < 0.1 pu */
    PQ_EVENT_SAG_MOMENTARY      = 6,   /**< 30 cycles - 3 sec, 0.1-0.9 pu */
    PQ_EVENT_SWELL_MOMENTARY    = 7,   /**< 30 cycles - 3 sec, 1.1-1.4 pu */
    PQ_EVENT_SAG_TEMPORARY      = 8,   /**< 3 sec - 1 min, 0.1-0.9 pu */
    PQ_EVENT_SWELL_TEMPORARY    = 9,   /**< 3 sec - 1 min, 1.1-1.2 pu */
    /* Long-duration */
    PQ_EVENT_INTERRUPTION_SUSTAINED = 10, /**< > 1 min, 0.0 pu */
    PQ_EVENT_UNDERVOLTAGE       = 11,  /**< > 1 min, 0.8-0.9 pu */
    PQ_EVENT_OVERVOLTAGE        = 12,  /**< > 1 min, 1.1-1.2 pu */
    /* Steady-state */
    PQ_EVENT_VOLTAGE_UNBALANCE  = 13,  /**< Negative sequence ≠ 0 */
    PQ_EVENT_HARMONIC_DISTORTION = 14, /**< THD exceeds limit */
    PQ_EVENT_FLICKER            = 15,  /**< Pst > 1.0 or Plt > 0.8 */
    PQ_EVENT_FREQUENCY_DEV      = 16,  /**< Δf beyond ±0.5 Hz (60Hz) or ±0.2 Hz (50Hz) */
    PQ_EVENT_DC_OFFSET          = 17,  /**< DC component in AC system */
    PQ_EVENT_NOTCHING           = 18,  /**< Commutation notches from converters */
    PQ_EVENT_NOISE              = 19,  /**< Broadband conducted EMI */
} pq_event_type_t;

/**
 * @brief Single power quality event record.
 *
 * Compliant with IEEE 1159-2019 event classification and
 * IEC 61000-4-30 measurement aggregation intervals.
 */
typedef struct {
    uint64_t timestamp_epoch;           /**< Unix epoch timestamp of event start */
    pq_event_type_t event_type;         /**< IEEE 1159 event category */
    double   duration_sec;              /**< Event duration [seconds] */
    double   v_min_pu;                  /**< Minimum voltage during event [pu of nominal] */
    double   v_max_pu;                  /**< Maximum voltage during event [pu] */
    double   v_avg_pu;                  /**< Average voltage during event [pu] */
    double   i_max_pu;                  /**< Maximum current [pu] */
    double   frequency_hz;              /**< Frequency at event onset */
    double   thd_v_pre;                 /**< Voltage THD before event [%] */
    double   pf_pre;                    /**< Power factor before event */
    uint8_t  phase_involved;            /**< Bitmask: b0=PhaseA, b1=PhaseB, b2=PhaseC */
    char     location[32];              /**< Measurement point identifier */
    char     likely_cause[64];          /**< Probable cause analysis */
} pq_event_t;

/**
 * @brief Power quality monitoring aggregator (IEC 61000-4-30 Class A).
 *
 * Maintains 10-min, 2-hour, and daily PQ statistics per
 * IEC 61000-4-30 measurement aggregation recommendations.
 */
typedef struct {
    double   v_rms_10min_avg;           /**< 10-minute RMS voltage average [V] */
    double   v_rms_10min_min;           /**< 10-minute minimum */
    double   v_rms_10min_max;           /**< 10-minute maximum */
    double   thd_v_10min_avg;           /**< 10-minute THD average [%] */
    double   thd_v_10min_95pct;         /**< 10-minute THD 95th percentile */
    double   freq_10min_avg;            /**< 10-minute frequency average [Hz] */
    double   pf_10min_avg;              /**< 10-minute average PF */
    double   unbalance_10min_avg;       /**< 10-minute voltage unbalance [%] */
    double   flicker_pst;               /**< Short-term flicker severity Pst (10 min) */
    double   flicker_plt;               /**< Long-term flicker severity Plt (2 h) */
    uint32_t sag_count_24h;             /**< Number of sags in last 24 hours */
    uint32_t swell_count_24h;           /**< Number of swells in last 24 hours */
    uint32_t interruption_count_24h;    /**< Number of interruptions in last 24 hours */
    pq_event_t event_log[PQ_MAX_EVENTS];
    uint32_t event_count;
} pq_monitor_t;

/* ==========================================================================
 * L2: ITIC/CBEMA Voltage Tolerance Curve
 * ========================================================================== */

/**
 * @brief ITIC (CBEMA) curve evaluation point.
 *
 * The ITIC curve defines the voltage tolerance envelope for
 * electronic equipment. Events inside the envelope are acceptable;
 * events outside may cause malfunction or damage.
 *
 * | Region        | Duration         | Voltage  |
 * |---------------|------------------|----------|
 * | Acceptable    | —                | Envelope |
 * | Prohibited    | Any              | > envelope upper |
 * | Prohibited    | > 0.5 sec        | 0%       |
 */
typedef enum {
    PQ_ITIC_ACCEPTABLE      = 0,  /**< Event within ITIC envelope */
    PQ_ITIC_ABOVE_UPPER     = 1,  /**< Voltage above upper envelope (damage risk) */
    PQ_ITIC_BELOW_LOWER     = 2,  /**< Voltage below lower envelope (malfunction) */
    PQ_ITIC_BOUNDARY        = 3,  /**< Exactly on boundary */
} pq_itic_result_t;

/**
 * @brief Evaluate an event against the ITIC (CBEMA) curve.
 *
 * The ITIC curve (revised 2000) defines acceptable voltage
 * at 60Hz equipment terminals:
 *
 * Upper envelope:
 *   500% for 0.001s, 200% for 0.01s, 120% for 0.5s,
 *   110% for steady-state (> 10s)
 *
 * Lower envelope:
 *   0% for 0.02s (one cycle), 70% for 0.5s,
 *   80% for 10s, 90% for steady-state
 *
 * @param duration_sec  Event duration
 * @param v_percent     Voltage as percent of nominal
 * @return ITIC acceptability classification
 */
pq_itic_result_t pq_check_itic(double duration_sec, double v_percent);

/**
 * @brief Get the ITIC upper envelope voltage for a given duration.
 *
 * @param duration_sec  Duration [seconds]
 * @return Upper envelope voltage as percent of nominal
 */
double pq_itic_upper_envelope(double duration_sec);

/**
 * @brief Get the ITIC lower envelope voltage for a given duration.
 *
 * @param duration_sec  Duration [seconds]
 * @return Lower envelope voltage as percent of nominal
 */
double pq_itic_lower_envelope(double duration_sec);

/* ==========================================================================
 * L4: Power Quality Quantification
 * ========================================================================== */

/**
 * @brief Compute short-term flicker severity Pst per IEC 61000-4-15.
 *
 * Flicker is the impression of unsteadiness of visual sensation
 * induced by a light stimulus whose luminance fluctuates with time.
 * Caused by voltage fluctuations from arc furnaces, welders, etc.
 *
 * Pst is computed over 10 minutes using a lamp-eye-brain model.
 * Pst > 1.0 is considered irritating.
 *
 * @param v_rms_samples  Array of RMS voltage samples (≥ 1 sample/sec for 10 min)
 * @param n_samples      Number of samples (should cover 10 minutes)
 * @param nominal_v      Nominal RMS voltage [V]
 * @return Pst value, or -1 on error
 */
double pq_compute_pst(const double *v_rms_samples, size_t n_samples,
                       double nominal_v);

/**
 * @brief Compute long-term flicker severity Plt from Pst values.
 *
 * Plt = ³√( (1/12) Σ_{i=1}^{12} Pst_i³ )
 *
 * Computed over 2 hours from twelve 10-minute Pst values.
 * Plt > 0.8 requires mitigation per IEC 61000-3-7.
 *
 * @param pst_values     Array of 12 Pst values
 * @return Plt value, or -1 on error
 */
double pq_compute_plt(const double *pst_values);

/**
 * @brief Compute the System Average Interruption Frequency Index (SAIFI).
 *
 * SAIFI = Σ (customers interrupted) / (total customers served)
 *
 * Standard IEEE 1366-2012 reliability index for distribution systems.
 *
 * @param customers_per_event  Array of customer counts per interruption
 * @param n_events             Number of interruption events
 * @param total_customers      Total customers served
 * @return SAIFI [interruptions per customer per year]
 */
double pq_compute_saifi(const uint32_t *customers_per_event,
                         size_t n_events, uint32_t total_customers);

/**
 * @brief Compute the System Average Interruption Duration Index (SAIDI).
 *
 * SAIDI = Σ (customer-minutes interrupted) / (total customers)
 *
 * @param customer_min_per_event  Customer-minutes per event
 * @param n_events                Number of events
 * @param total_customers         Total customers served
 * @return SAIDI [minutes per customer per year]
 */
double pq_compute_saidi(const double *customer_min_per_event,
                         size_t n_events, uint32_t total_customers);

/**
 * @brief Compute the Momentary Average Interruption Frequency Index (MAIFI).
 *
 * MAIFI accounts for momentary interruptions (< 5 minutes).
 *
 * @param momentary_events  Number of momentary events per customer group
 * @param n_groups          Number of customer groups
 * @param total_customers   Total customers
 * @return MAIFI
 */
double pq_compute_maifi(const uint32_t *momentary_events,
                         size_t n_groups, uint32_t total_customers);

/* ==========================================================================
 * L7: Application-Level Power Quality Assessment
 * ========================================================================== */

/**
 * @brief Assess power quality for data center operation (80 PLUS / Energy Star).
 *
 * Data centers require:
 * - PF ≥ 0.95 at IT load
 * - THD_i ≤ 5% (IEEE 519 at PCC)
 * - Voltage within ±10% of nominal
 * - 80 PLUS Titanium: ≥96% efficiency at 50% load
 *
 * @param pf_measured       Measured PF at PDU level
 * @param thd_i_percent     Current THD [%]
 * @param v_deviation_pct   Voltage deviation from nominal [%]
 * @param efficiency_pct    Power supply efficiency [%]
 * @param tier              [out] 80 PLUS tier (0=bronze..4=titanium, >=6=no tier)
 * @param score             [out] PQ score 0-100
 * @return 0 if meets requirements, bitmask of issues found
 */
int pq_assess_datacenter(double pf_measured, double thd_i_percent,
                          double v_deviation_pct, double efficiency_pct,
                          int *tier, double *score);

/**
 * @brief Assess power quality for EV charging station (SAE J2894 / IEC 61851).
 *
 * EV chargers must meet:
 * - PF ≥ 0.95 at rated power
 * - THD_i ≤ 5%
 * - Voltage sag ride-through per SAE J2894
 *
 * @param pf_rated         PF at rated charging power
 * @param thd_i_percent    Current THD [%]
 * @param v_sag_remaining  Voltage during sag [% of nominal]
 * @param sag_duration_ms  Sag duration [ms]
 * @param compliance       [out] 1 if SAE J2894 compliant, 0 if not
 * @return 0 on success
 */
int pq_assess_ev_charger(double pf_rated, double thd_i_percent,
                          double v_sag_remaining, double sag_duration_ms,
                          int *compliance);

/**
 * @brief Compute the cost of poor power quality.
 *
 * Estimates annual losses from:
 * - Increased I²R losses due to low PF
 * - Harmonic losses in transformers and cables
 * - Production downtime from voltage sags
 * - Utility PF penalty charges
 *
 * Based on EPRI cost-of-PQ study methodology.
 *
 * @param avg_pf           Average power factor
 * @param annual_kwh       Annual energy consumption [kWh]
 * @param cost_per_kwh     Electricity rate [$/kWh]
 * @param pf_penalty_rate  Utility PF penalty rate [$/kVAR/month]
 * @param sag_events_yr    Number of disruptive sag events per year
 * @param cost_per_sag     Average cost per sag event [$]
 * @return Total annual PQ cost [$/year]
 */
double pq_cost_of_poor_quality(double avg_pf, double annual_kwh,
                                double cost_per_kwh, double pf_penalty_rate,
                                int sag_events_yr, double cost_per_sag);

/**
 * @brief Generate a power quality report card.
 *
 * Summarizes PQ metrics against industry benchmarks.
 *
 * @param monitor     PQ monitor data
 * @param report_buf  [out] Report text buffer
 * @param buf_len     Buffer length
 * @return 0 on success
 */
int pq_generate_report(const pq_monitor_t *monitor, char *report_buf,
                        size_t buf_len);

#ifdef __cplusplus
}
#endif

#endif /* POWER_QUALITY_H */
