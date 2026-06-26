#include "transient_defs.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*==================================================================
 * L1: Component Model Initialization
 * Each function initializes a circuit component model with
 * parasitic parameters critical for accurate transient simulation.
 *==================================================================*/

int capacitor_model_init(CapacitorModel_t *c, int id, const char *name,
                         double C, double ESR, double ESL)
{
    if (!c) return -1;
    if (C < 0.0 || ESR < 0.0 || ESL < 0.0) return -1;
    memset(c, 0, sizeof(*c));
    c->id = id;
    c->C = C;
    c->ESR = ESR;
    c->ESL = ESL;
    c->R_leakage = 1e6;
    c->V_rated = 50.0;
    c->V_initial = 0.0;
    c->temp_coeff_ppm = 200.0;
    c->dissipation_factor = 0.02;
    if (name) {
        strncpy(c->name, name, MAX_NAME_TA - 1);
        c->name[MAX_NAME_TA - 1] = '\0';
    } else {
        snprintf(c->name, MAX_NAME_TA, "C%d", id);
    }
    return 0;
}

int inductor_model_init(InductorModel_t *l, int id, const char *name,
                        double L, double DCR)
{
    if (!l) return -1;
    if (L < 0.0 || DCR < 0.0) return -1;
    memset(l, 0, sizeof(*l));
    l->id = id;
    l->L = L;
    l->DCR = DCR;
    l->EPR = 1e6;
    l->EPC = 1e-12;
    l->I_rated = 1.0;
    l->I_saturation = 2.0;
    l->I_initial = 0.0;
    if (name) {
        strncpy(l->name, name, MAX_NAME_TA - 1);
        l->name[MAX_NAME_TA - 1] = '\0';
    } else {
        snprintf(l->name, MAX_NAME_TA, "L%d", id);
    }
    return 0;
}

int switch_model_init(SwitchModel_t *sw, int id, const char *name,
                      double R_on, double R_off)
{
    if (!sw) return -1;
    if (R_on < 0.0 || R_off <= R_on) return -1;
    memset(sw, 0, sizeof(*sw));
    sw->id = id;
    sw->R_on = R_on;
    sw->R_off = R_off;
    sw->V_threshold = 0.7;
    sw->t_rise = 1e-9;
    sw->t_fall = 1e-9;
    sw->t_delay_on = 5e-9;
    sw->t_delay_off = 10e-9;
    sw->state = 0;
    if (name) {
        strncpy(sw->name, name, MAX_NAME_TA - 1);
        sw->name[MAX_NAME_TA - 1] = '\0';
    } else {
        snprintf(sw->name, MAX_NAME_TA, "SW%d", id);
    }
    return 0;
}

int diode_model_init(DiodeModel_t *d, int id, const char *name,
                     double I_s, double V_t, double C_j0)
{
    if (!d) return -1;
    if (I_s <= 0.0 || V_t <= 0.0 || C_j0 < 0.0) return -1;
    memset(d, 0, sizeof(*d));
    d->id = id;
    d->I_s = I_s;
    d->V_t = V_t;
    d->n = 1.0;
    d->C_j0 = C_j0;
    d->V_bi = 0.7;
    d->m = 0.5;
    d->t_rr = 4e-9;
    d->Q_rr = 50e-12;
    if (name) {
        strncpy(d->name, name, MAX_NAME_TA - 1);
        d->name[MAX_NAME_TA - 1] = '\0';
    } else {
        snprintf(d->name, MAX_NAME_TA, "D%d", id);
    }
    return 0;
}

void transient_sim_init(TransientSim_t *sim, int n_states)
{
    if (!sim || n_states <= 0) return;
    memset(sim, 0, sizeof(*sim));
    sim->n_states = n_states;
    sim->state = (double*)calloc((size_t)n_states, sizeof(double));
    sim->state_deriv = (double*)calloc((size_t)n_states, sizeof(double));
    sim->solver = SOLVER_RK4;
    sim->tolerance = 1e-6;
    sim->dt_min = 1e-15;
    sim->dt_max = 1.0;
    sim->stiffness = STIFFNESS_NONSTIFF;
    if (!sim->state || !sim->state_deriv) {
        free(sim->state);
        free(sim->state_deriv);
        sim->n_states = 0;
    }
}

void transient_sim_free(TransientSim_t *sim)
{
    if (!sim) return;
    free(sim->state);
    free(sim->state_deriv);
    sim->state = NULL;
    sim->state_deriv = NULL;
    sim->n_states = 0;
}

/*==================================================================
 * L1: Fundamental Parameter Computations
 * Each function computes a basic characteristic parameter
 * from circuit element values.
 *==================================================================*/

/* L2: Time constant tau = R*C (seconds)
 * This is the time for the response to reach 63.2% of final value.
 * After 5*tau, response is within 0.67% of final value. */
double compute_time_constant_rc(double R, double C)
{
    if (R < 0.0 || C < 0.0) return -1.0;
    return R * C;
}

/* L2: Time constant tau = L/R (seconds)
 * Energy stored in inductor decays with this time constant. */
double compute_time_constant_rl(double R, double L)
{
    if (R <= 0.0 || L < 0.0) return -1.0;
    return L / R;
}

/* L3-L4: Compute complete second-order characteristic parameters
 * from R, L, C values for series or parallel RLC circuits.
 *
 * Series RLC:
 *   alpha = R/(2L), omega_n = 1/sqrt(LC)
 *   zeta = alpha/omega_n = (R/2)*sqrt(C/L)
 *
 * Parallel RLC:
 *   alpha = 1/(2RC), omega_n = 1/sqrt(LC)
 *   zeta = alpha/omega_n = (1/(2R))*sqrt(L/C) */
SecondOrderParams_t compute_second_order_params(double R, double L, double C,
                                                 int is_parallel)
{
    SecondOrderParams_t p;
    memset(&p, 0, sizeof(p));

    if (L <= 0.0 || C <= 0.0) {
        p.omega_n = 0.0;
        p.damping = DAMPING_OVERDAMPED;
        return p;
    }

    p.omega_n = 1.0 / sqrt(L * C);

    if (is_parallel) {
        if (R <= 0.0) {
            p.alpha = INFINITY;
            p.zeta = INFINITY;
        } else {
            p.alpha = 1.0 / (2.0 * R * C);
            p.zeta = p.alpha / p.omega_n;
        }
    } else {
        if (L <= 0.0) {
            p.alpha = INFINITY;
            p.zeta = INFINITY;
        } else {
            p.alpha = R / (2.0 * L);
            p.zeta = p.alpha / p.omega_n;
        }
    }

    p.sigma = -p.alpha;

    if (p.zeta > 1.0) {
        p.damping = DAMPING_OVERDAMPED;
        double disc = p.omega_n * sqrt(p.zeta * p.zeta - 1.0);
        p.pole_real_1 = -p.alpha + disc;
        p.pole_real_2 = -p.alpha - disc;
        p.pole_imag = 0.0;
        p.omega_d = 0.0;
    } else if (fabs(p.zeta - 1.0) < 1e-10) {
        p.damping = DAMPING_CRITICALLY_DAMPED;
        p.pole_real_1 = -p.alpha;
        p.pole_real_2 = -p.alpha;
        p.pole_imag = 0.0;
        p.omega_d = 0.0;
    } else if (p.zeta > 0.0) {
        p.damping = DAMPING_UNDERDAMPED;
        p.omega_d = p.omega_n * sqrt(1.0 - p.zeta * p.zeta);
        p.pole_real_1 = -p.alpha;
        p.pole_real_2 = -p.alpha;
        p.pole_imag = p.omega_d;
    } else if (fabs(p.zeta) < 1e-10) {
        p.damping = DAMPING_UNDAMPED;
        p.omega_d = p.omega_n;
        p.pole_real_1 = 0.0;
        p.pole_real_2 = 0.0;
        p.pole_imag = p.omega_n;
    } else {
        p.damping = DAMPING_NEGATIVE;
        p.pole_real_1 = -p.alpha;
        p.pole_real_2 = -p.alpha;
        p.pole_imag = 0.0;
        p.omega_d = 0.0;
    }

    return p;
}

/* L2: Classify damping based on zeta value */
DampingClass_t classify_damping(double zeta)
{
    if (zeta > 1.0 + 1e-9) return DAMPING_OVERDAMPED;
    if (fabs(zeta - 1.0) < 1e-9) return DAMPING_CRITICALLY_DAMPED;
    if (zeta > 1e-9) return DAMPING_UNDERDAMPED;
    if (fabs(zeta) < 1e-9) return DAMPING_UNDAMPED;
    return DAMPING_NEGATIVE;
}

/* L2: Quality factor Q for series RLC.
 * Q = (1/R)*sqrt(L/C) = omega_n * L / R
 * Q measures the sharpness of resonance. Higher Q = sharper peak,
 * longer ringing in transient response.
 * For R->0, Q->infinity. For large R, Q->0. */
double compute_quality_factor_series(double R, double L, double C)
{
    if (R <= 0.0) return INFINITY;
    if (L <= 0.0 || C <= 0.0) return 0.0;
    return sqrt(L / C) / R;
}

/* L2: Quality factor Q for parallel RLC.
 * Q = R*sqrt(C/L) = R/(omega_n*L)
 * Large R => high Q. Small R => low Q. */
double compute_quality_factor_parallel(double R, double L, double C)
{
    if (R <= 0.0) return 0.0;
    if (L <= 0.0 || C <= 0.0) return 0.0;
    return R * sqrt(C / L);
}

/* L4: Resonant frequency omega_0 = 1/sqrt(LC)
 * Reference: Thomson (Lord Kelvin), 1853.
 * At resonance, inductive and capacitive reactances cancel.
 * This is the undamped natural frequency. */
double compute_resonant_frequency(double L, double C)
{
    if (L <= 0.0 || C <= 0.0) return 0.0;
    return 1.0 / sqrt(L * C);
}

/* L4: Damped natural frequency omega_d = omega_n * sqrt(1 - zeta^2)
 * Only defined for underdamped systems (0 <= zeta < 1).
 * This is the frequency of oscillation observed in the transient response. */
double compute_damped_frequency(double omega_n, double zeta)
{
    if (omega_n <= 0.0) return 0.0;
    if (zeta >= 1.0) return 0.0;
    return omega_n * sqrt(1.0 - zeta * zeta);
}

/* L4: Neper frequency (attenuation constant) alpha
 * Series: alpha = R/(2L)
 * Parallel: alpha = 1/(2RC) */
double compute_neper_frequency(double R, double L, double C, int is_parallel)
{
    if (is_parallel) {
        if (R <= 0.0 || C <= 0.0) return INFINITY;
        return 1.0 / (2.0 * R * C);
    } else {
        if (L <= 0.0) return INFINITY;
        return R / (2.0 * L);
    }
}

/*==================================================================
 * L1: Energy Storage Computations
 * Reference: Maxwell 1861 (capacitor energy), Faraday 1831 (inductor energy)
 *
 * Capacitor stores energy in electric field: E = 1/2 * C * V^2
 * Inductor stores energy in magnetic field: E = 1/2 * L * I^2
 *
 * Energy continuity: v_c(0+) = v_c(0-), i_l(0+) = i_l(0-)
 * These continuity conditions are fundamental to transient analysis.
 *==================================================================*/

EnergyState_t compute_energy_state(double C, double L, double vc, double il)
{
    EnergyState_t e;
    e.v_capacitor = vc;
    e.i_inductor = il;
    e.energy_c = 0.5 * C * vc * vc;
    e.energy_l = 0.5 * L * il * il;
    e.total_energy = e.energy_c + e.energy_l;
    return e;
}

/*==================================================================
 * L3: Pole-Zero Analysis
 *
 * The characteristic equation determines the form of the natural response.
 * First-order: s + 1/tau = 0  =>  s = -1/tau (single real pole)
 * Second-order: s^2 + 2*zeta*omega_n*s + omega_n^2 = 0
 *==================================================================*/

/* L3: First-order pole location.
 * Single pole at s = -1/tau on negative real axis.
 * More negative => faster decay. */
int compute_poles_first_order(double tau)
{
    if (tau <= 0.0) return 0; /* unstable or degenerate */
    return 1; /* one pole at s = -1/tau */
}

/* L3: Second-order pole locations.
 * Overdamped: two distinct real poles
 * Critically damped: repeated real pole
 * Underdamped: complex conjugate pair
 * Undamped: purely imaginary pair
 *
 * Returns a ComplexFrequency_t representing the pole location.
 * For complex poles, sigma = real part, omega = imag part. */
ComplexFrequency_t compute_poles_second_order(double zeta, double omega_n)
{
    ComplexFrequency_t cf;
    double alpha = zeta * omega_n;

    if (zeta >= 1.0) {
        /* Real poles */
        cf.sigma = -alpha - omega_n * sqrt(zeta * zeta - 1.0);
        cf.omega = 0.0;
    } else if (zeta > 0.0) {
        cf.sigma = -alpha;
        cf.omega = omega_n * sqrt(1.0 - zeta * zeta);
    } else if (fabs(zeta) < 1e-10) {
        cf.sigma = 0.0;
        cf.omega = omega_n;
    } else {
        cf.sigma = -alpha;
        cf.omega = 0.0;
    }
    return cf;
}

/*==================================================================
 * L6: Transient Performance Metrics
 *
 * These metrics quantify transient response quality from measured
 * or simulated waveforms. They are fundamental to control system
 * design and circuit characterization.
 *==================================================================*/

/* L6: Compute performance metrics from a time-sampled waveform.
 * v_initial and v_final are the initial and final steady-state values.
 * Returns overshoot, settling time, rise time, etc. */
TransientMetrics_t compute_metrics_from_waveform(const double *t, const double *v,
                                                  int n, double v_initial, double v_final)
{
    TransientMetrics_t m;
    memset(&m, 0, sizeof(m));

    if (!t || !v || n < 2) {
        m.final_value = v_final;
        return m;
    }

    double delta_v = v_final - v_initial;
    if (fabs(delta_v) < 1e-15) delta_v = 1.0;

    m.final_value = v_final;
    m.peak_value = v_initial;
    m.steady_state_error = 0.0;

    int i_peak = 0;
    int settling_2_idx = -1;
    int settling_5_idx = -1;
    int rise_10_idx = -1;
    int rise_90_idx = -1;
    int delay_50_idx = -1;

    double v_10 = v_initial + 0.10 * delta_v;
    double v_90 = v_initial + 0.90 * delta_v;
    double v_50 = v_initial + 0.50 * delta_v;

    /* Find peak */
    for (int i = 0; i < n; i++) {
        if (v[i] > m.peak_value) {
            m.peak_value = v[i];
            i_peak = i;
        }
        if (v[i] < m.peak_value * 0.0) {
            /* track undershoot */
        }
    }

    /* Find 10%, 50%, 90% crossing times */
    for (int i = 0; i < n; i++) {
        if (rise_10_idx < 0 && v[i] >= v_10) rise_10_idx = i;
        if (delay_50_idx < 0 && v[i] >= v_50) delay_50_idx = i;
        if (rise_90_idx < 0 && v[i] >= v_90) rise_90_idx = i;
    }

    if (rise_10_idx >= 0 && rise_90_idx >= 0) {
        m.rise_time = t[rise_90_idx] - t[rise_10_idx];
    }
    if (delay_50_idx >= 0) {
        m.delay_time = t[delay_50_idx] - t[0];
    }
    if (i_peak > 0) {
        m.peak_time = t[i_peak];
    }

    /* Overshoot percentage */
    if (fabs(v_final) > 1e-15) {
        m.overshoot_pct = ((m.peak_value - v_final) / fabs(v_final)) * 100.0;
        if (m.overshoot_pct < 0.0) m.overshoot_pct = 0.0;
    }

    /* Settling time: find when signal stays within 2% band */
    double band_2 = 0.02 * fabs(v_final);
    double band_5 = 0.05 * fabs(v_final);
    for (int i = n - 1; i >= 0; i--) {
        if (settling_2_idx < 0 && fabs(v[i] - v_final) > band_2) {
            settling_2_idx = (i + 1 < n) ? i + 1 : i;
        }
        if (settling_5_idx < 0 && fabs(v[i] - v_final) > band_5) {
            settling_5_idx = (i + 1 < n) ? i + 1 : i;
        }
    }
    if (settling_2_idx >= 0) m.settling_time_2pct = t[settling_2_idx];
    if (settling_5_idx >= 0) m.settling_time_5pct = t[settling_5_idx];

    /* Steady-state error */
    m.steady_state_error = v[n-1] - v_final;

    return m;
}

/* L6: Overshoot as a function of damping ratio.
 * PO = 100 * exp(-pi*zeta / sqrt(1 - zeta^2))
 * Valid for 0 <= zeta < 1. For zeta >= 1, overshoot is 0.
 * This formula is derived from the peak of the underdamped step response. */
double overshoot_from_zeta(double zeta)
{
    if (zeta >= 1.0) return 0.0;
    if (zeta <= 0.0) return 100.0;
    double exponent = -M_PI * zeta / sqrt(1.0 - zeta * zeta);
    return 100.0 * exp(exponent);
}

/* L6: Inverse relation: zeta from measured overshoot.
 * zeta = -ln(PO/100) / sqrt(pi^2 + ln(PO/100)^2)
 * Used to identify damping from experimental step response data. */
double zeta_from_overshoot(double overshoot_pct)
{
    if (overshoot_pct <= 0.0) return 1.0;
    if (overshoot_pct >= 100.0) return 0.0;
    double ln_os = log(overshoot_pct / 100.0);
    return -ln_os / sqrt(M_PI * M_PI + ln_os * ln_os);
}

/* L6: 2% settling time: t_s = 4 / (zeta * omega_n)
 * Conservative estimate: response stays within 2% after 4 time constants.
 * More precisely: t_s = -ln(0.02) / (zeta * omega_n) ≈ 3.91 / (zeta*omega_n) */
double settling_time_2pct(double zeta, double omega_n)
{
    if (zeta <= 0.0 || omega_n <= 0.0) return INFINITY;
    return 4.0 / (zeta * omega_n);
}

/* L6: 5% settling time: t_s = 3 / (zeta * omega_n) */
double settling_time_5pct(double zeta, double omega_n)
{
    if (zeta <= 0.0 || omega_n <= 0.0) return INFINITY;
    return 3.0 / (zeta * omega_n);
}

/* L6: Peak time: t_p = pi / (omega_n * sqrt(1 - zeta^2)) = pi / omega_d
 * Time at which the first peak occurs for underdamped step response. */
double peak_time_from_params(double omega_n, double zeta)
{
    if (omega_n <= 0.0 || zeta >= 1.0 || zeta < 0.0) return INFINITY;
    double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
    if (omega_d <= 0.0) return INFINITY;
    return M_PI / omega_d;
}

/* L6: 10-90% rise time approximation for second-order systems.
 * t_r ≈ (1.8) / omega_n for zeta ≈ 0.5
 * More generally: t_r ≈ (1 + 1.1*zeta + 1.4*zeta^2) / omega_n
 * (Empirical fit from Ogata, Modern Control Engineering) */
double rise_time_10_90_from_params(double omega_n, double zeta)
{
    if (omega_n <= 0.0) return INFINITY;
    if (zeta < 0.0) zeta = 0.0;
    double numerator = 1.0 + 1.1 * zeta + 1.4 * zeta * zeta;
    return numerator / omega_n;
}

/*==================================================================
 * L4: Energy and Power in Transient Circuits
 *==================================================================*/

/* L4: Capacitor energy: E_c = 0.5 * C * v^2 (Joules) */
double capacitor_energy(double C, double v)
{
    return 0.5 * C * v * v;
}

/* L4: Inductor energy: E_l = 0.5 * L * i^2 (Joules) */
double inductor_energy(double L, double i)
{
    return 0.5 * L * i * i;
}

/* L4: Power dissipated in resistor: P_R = i^2 * R (Watts) */
double resistor_power(double i, double R)
{
    return i * i * R;
}

/* L4: Energy dissipated over time interval dt:
 * E_R = integral(i^2 * R * dt) ≈ i^2 * R * dt */
double resistor_energy_approx(double i, double R, double dt)
{
    return i * i * R * dt;
}
