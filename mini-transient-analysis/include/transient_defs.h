#ifndef TRANSIENT_DEFS_H
#define TRANSIENT_DEFS_H
#include <stddef.h>
#include <math.h>
#include <complex.h>

/*==================================================================
 * L1: Core Definitions - Transient Analysis
 * Reference: Hayt, Kemmerly & Durbin (2019), Nilsson & Riedel (2019),
 *           Dorf & Svoboda (2014), Sedra & Smith (2020)
 *
 * Transient analysis studies circuit behaviour during the transition
 * from one steady state to another, governed by differential equations
 * arising from energy storage elements (L, C).
 *
 * Key equation: d/dt(state) = f(state, input, t)
 *==================================================================*/

typedef enum {
    ORDER_FIRST = 1,
    ORDER_SECOND = 2,
    ORDER_THIRD = 3,
    ORDER_FOURTH = 4,
    ORDER_HIGHER = 5
} SystemOrder_t;

typedef enum {
    DAMPING_OVERDAMPED = 0,
    DAMPING_CRITICALLY_DAMPED = 1,
    DAMPING_UNDERDAMPED = 2,
    DAMPING_UNDAMPED = 3,
    DAMPING_NEGATIVE = 4
} DampingClass_t;

typedef enum {
    RESPONSE_NATURAL = 0,
    RESPONSE_FORCED = 1,
    RESPONSE_COMPLETE = 2,
    RESPONSE_STEP = 3,
    RESPONSE_IMPULSE = 4,
    RESPONSE_RAMP = 5,
    RESPONSE_SINUSOIDAL = 6,
    RESPONSE_PULSE = 7,
    RESPONSE_EXPONENTIAL = 8
} ResponseType_t;

typedef struct {
    double rise_time;
    double fall_time;
    double delay_time;
    double settling_time_2pct;
    double settling_time_5pct;
    double peak_time;
    double overshoot_pct;
    double undershoot_pct;
    double steady_state_error;
    double peak_value;
    double final_value;
} TransientMetrics_t;

typedef struct {
    double tau;
    double tau_effective;
    double thermal_time_constant;
    double mechanical_tau;
} TimeConstants_t;

typedef struct {
    double zeta;
    double omega_n;
    double omega_d;
    double alpha;
    double sigma;
    double pole_real_1;
    double pole_real_2;
    double pole_imag;
    DampingClass_t damping;
} SecondOrderParams_t;

typedef struct {
    double v_capacitor;
    double i_inductor;
    double energy_c;
    double energy_l;
    double total_energy;
} EnergyState_t;

typedef struct {
    double vc0;
    double il0;
    double dvdt0;
    double didt0;
    int has_initial_charge;
} InitialConditions_t;

typedef enum {
    SOLVER_ANALYTIC = 0,
    SOLVER_EULER_FWD = 1,
    SOLVER_EULER_BWD = 2,
    SOLVER_TRAPEZOIDAL = 3,
    SOLVER_RK4 = 4,
    SOLVER_RK45 = 5,
    SOLVER_GEAR2 = 6,
    SOLVER_AB2 = 7,
    SOLVER_AM2 = 8
} SolverType_t;

typedef enum {
    STIFFNESS_NONSTIFF = 0,
    STIFFNESS_MILDLY_STIFF = 1,
    STIFFNESS_STIFF = 2,
    STIFFNESS_VERY_STIFF = 3
} StiffnessClass_t;

typedef enum {
    SWITCH_ON = 0,
    SWITCH_OFF = 1,
    SWITCH_COMMUTATION = 2,
    SWITCH_FAULT = 3,
    SWITCH_ZVS = 4,
    SWITCH_ZCS = 5
} SwitchEvent_t;

typedef enum {
    BJT_CUTOFF = 0,
    BJT_FORWARD_ACTIVE = 1,
    BJT_SATURATION = 2,
    BJT_REVERSE_ACTIVE = 3
} BJTRegion_t;

typedef enum {
    MOSFET_CUTOFF = 0,
    MOSFET_LINEAR = 1,
    MOSFET_SATURATION = 2,
    MOSFET_SUBTHRESHOLD = 3
} MOSFETRegion_t;

#define MAX_COMPONENTS 64
#define MAX_NODES_TA 128
#define MAX_NAME_TA 48

typedef struct {
    int id;
    char name[MAX_NAME_TA];
    double C;
    double ESR;
    double ESL;
    double R_leakage;
    double V_rated;
    double V_initial;
    double temp_coeff_ppm;
    double dissipation_factor;
} CapacitorModel_t;

typedef struct {
    int id;
    char name[MAX_NAME_TA];
    double L;
    double DCR;
    double EPR;
    double EPC;
    double I_rated;
    double I_saturation;
    double I_initial;
} InductorModel_t;

typedef struct {
    int id;
    char name[MAX_NAME_TA];
    double R_on;
    double R_off;
    double V_threshold;
    double t_rise;
    double t_fall;
    double t_delay_on;
    double t_delay_off;
    int state;
} SwitchModel_t;

typedef struct {
    int id;
    char name[MAX_NAME_TA];
    double I_s;
    double V_t;
    double n;
    double C_j0;
    double V_bi;
    double m;
    double t_rr;
    double Q_rr;
} DiodeModel_t;

typedef struct {
    double t;
    double t_start;
    double t_end;
    double dt;
    double dt_min;
    double dt_max;
    int n_states;
    double *state;
    double *state_deriv;
    SolverType_t solver;
    StiffnessClass_t stiffness;
    double tolerance;
    int step_count;
    int rejected_steps;
    double cpu_time_ms;
} TransientSim_t;

typedef struct {
    double t_event;
    int event_type;
    double value_before;
    double value_after;
    void (*callback)(void*);
    void *user_data;
} TransientEvent_t;

typedef struct {
    double sigma;
    double omega;
} ComplexFrequency_t;

typedef struct {
    int n_poles;
    int n_zeros;
    ComplexFrequency_t poles[32];
    ComplexFrequency_t zeros[32];
    double dc_gain;
} PoleZeroMap_t;

typedef void (*ODERHSFunc)(double t, const double *y, double *dydt, void *params);

typedef struct {
    int id;
    char name[MAX_NAME_TA];
    double (*i_v_func)(double v, void *params);
    double (*v_i_func)(double i, void *params);
    double (*di_dv_func)(double v, void *params);
    void *params;
    double v_op;
    double i_op;
    double g_op;
} NonlinearElement_t;

typedef struct {
    double Vin;
    double Vout;
    double L;
    double C;
    double R_load;
    double f_sw;
    double duty;
    double I_load;
    double ripple_pct;
} BuckConverterParams_t;

typedef struct {
    double Vin;
    double Vout;
    double L;
    double C;
    double R_load;
    double f_sw;
    double duty;
} BoostConverterParams_t;

int capacitor_model_init(CapacitorModel_t *c, int id, const char *name,
                         double C, double ESR, double ESL);
int inductor_model_init(InductorModel_t *l, int id, const char *name,
                        double L, double DCR);
int switch_model_init(SwitchModel_t *sw, int id, const char *name,
                      double R_on, double R_off);
int diode_model_init(DiodeModel_t *d, int id, const char *name,
                     double I_s, double V_t, double C_j0);
void transient_sim_init(TransientSim_t *sim, int n_states);
void transient_sim_free(TransientSim_t *sim);

double compute_time_constant_rc(double R, double C);
double compute_time_constant_rl(double R, double L);
SecondOrderParams_t compute_second_order_params(double R, double L, double C, int is_parallel);
DampingClass_t classify_damping(double zeta);
double compute_quality_factor_series(double R, double L, double C);
double compute_quality_factor_parallel(double R, double L, double C);
double compute_resonant_frequency(double L, double C);
double compute_damped_frequency(double omega_n, double zeta);
double compute_neper_frequency(double R, double L, double C, int is_parallel);

EnergyState_t compute_energy_state(double C, double L, double vc, double il);

int compute_poles_first_order(double tau);
ComplexFrequency_t compute_poles_second_order(double zeta, double omega_n);

TransientMetrics_t compute_metrics_from_waveform(const double *t, const double *v,
                                                  int n, double v_initial, double v_final);
double overshoot_from_zeta(double zeta);
double zeta_from_overshoot(double overshoot_pct);
double settling_time_2pct(double zeta, double omega_n);
double settling_time_5pct(double zeta, double omega_n);
double peak_time_from_params(double omega_n, double zeta);
double rise_time_10_90_from_params(double omega_n, double zeta);

#endif /* TRANSIENT_DEFS_H */
