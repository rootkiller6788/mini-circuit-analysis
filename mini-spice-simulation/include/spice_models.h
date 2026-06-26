/**
 * @file spice_models.h
 * @brief Device model evaluation and MNA stamping functions
 *
 * Knowledge coverage:
 *   L1 (Definitions): Shockley diode equation, MOSFET Shichman-Hodges, BJT Ebers-Moll
 *   L2 (Core Concepts): Companion models for numerical integration
 *   L4 (Fundamental Laws): Ohm's Law, capacitor/inductor I-V relationships
 *   L5 (Algorithms): Nonlinear device linearization (Jacobian), companion model stamping
 *
 * Each device type provides:
 *   1. A "load" function that evaluates I-V and partial derivatives
 *   2. A "stamp" function that contributes to the MNA matrix and RHS
 *
 * Course alignment:
 *   Berkeley EE105 (Analog IC): Device physics and SPICE models
 *   MIT 6.630 (EM): Lumped-element modeling
 *   TU Munich (High-Frequency): Nonlinear device characterization
 */

#ifndef SPICE_MODELS_H
#define SPICE_MODELS_H

#include <stdint.h>
#include "spice_netlist.h"
#include "spice_matrix.h"

/* ── Physical Constants ────────────────────────────────────────────── */

/** Boltzmann constant k_B (J/K) */
#define SPICE_BOLTZMANN      1.380649e-23
/** Electron charge q (C) */
#define SPICE_E_CHARGE       1.602176634e-19
/** Absolute zero offset (K = °C + 273.15) */
#define SPICE_KELVIN_OFFSET  273.15

/**
 * @brief Compute thermal voltage V_T = k*T/q
 *
 * @param temp_celsius Temperature in degrees Celsius
 * @return V_T in volts
 *
 * At 27°C: V_T ≈ 25.85 mV
 *
 * Reference: Sedra & Smith §3.5
 */
double spice_thermal_voltage(double temp_celsius);

/* ── L1/L2: Resistor Model ─────────────────────────────────────────── */

/**
 * @brief Stamp a linear resistor into the MNA system
 *
 * Conductance G = 1/R between nodes n+ and n-:
 *   [ +G  -G ] [v_n+]   [ 0 ]
 *   [ -G  +G ] [v_n-] = [ 0 ]
 *
 * @param G  Dense conductance matrix (MNA)
 * @param rhs Right-hand side vector
 * @param nplus  Positive node index
 * @param nminus Negative node index
 * @param resistance R in ohms
 * @param mna_size Total MNA system size
 *
 * Reference: Kielkowski (1998) §2.2
 * Complexity: O(1)
 */
void spice_stamp_resistor(spice_dense_matrix_t *G, spice_vector_t *rhs,
                           int32_t nplus, int32_t nminus,
                           double resistance, int32_t mna_size);

/* ── L1/L2: Capacitor Model ────────────────────────────────────────── */

/**
 * @brief Stamp a capacitor companion model for transient analysis
 *
 * Trapezoidal rule discretizes I = C * dV/dt:
 *   I_{n+1} = G_eq * V_{n+1} + I_eq
 * where G_eq = 2*C/dt, I_eq = -(G_eq * V_n + I_n)
 *
 * This is a Norton equivalent: conductance G_eq in parallel with
 * current source I_eq connected between n+ and n-.
 *
 * @param G  Conductance matrix
 * @param rhs RHS vector
 * @param nplus  Positive node
 * @param nminus Negative node
 * @param capacitance C in farads
 * @param dt Time step (seconds)
 * @param prev_voltage V across cap at previous time step
 * @param prev_current I through cap at previous time step
 * @param mna_size MNA system size
 *
 * Reference: Chua & Lin (1975) §7.2
 */
void spice_stamp_capacitor_tran(spice_dense_matrix_t *G, spice_vector_t *rhs,
                                 int32_t nplus, int32_t nminus,
                                 double capacitance, double dt,
                                 double prev_voltage, double prev_current,
                                 int32_t mna_size);

/**
 * @brief Stamp capacitor for DC analysis (open circuit)
 *
 * In DC, a capacitor is an open circuit → no contribution to MNA.
 * This is a no-op for documentation completeness.
 */
void spice_stamp_capacitor_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                               int32_t nplus, int32_t nminus,
                               int32_t mna_size);

/**
 * @brief Stamp capacitor for AC analysis (admittance jωC)
 *
 * Contributes s*C = jωC to the complex admittance matrix.
 *
 * @param G Complex conductance matrix
 * @param rhs Complex RHS
 * @param nplus  Positive node
 * @param nminus Negative node
 * @param capacitance C in farads
 * @param omega Angular frequency (rad/s), ω = 2πf
 * @param mna_size MNA size
 *
 * Complex stamp: Y = jωC between n+ and n-
 *   [ +jωC  -jωC ]
 *   [ -jωC  +jωC ]
 */
void spice_stamp_capacitor_ac(spice_complex_t *G, spice_complex_t *rhs,
                               int32_t nplus, int32_t nminus,
                               double capacitance, double omega,
                               int32_t mna_size);

/* ── L1/L2: Inductor Model ─────────────────────────────────────────── */

/**
 * @brief Stamp an inductor for DC analysis (short circuit)
 *
 * In DC, an inductor is a short circuit. In MNA, this adds
 * a branch current variable I_L and a KVL constraint:
 *   V(n+) - V(n-) = 0 (DC steady-state: dI/dt = 0)
 *
 * This effectively ties n+ and n- together.
 *
 * @param G  Conductance matrix
 * @param rhs RHS vector
 * @param nplus  Positive node
 * @param nminus Negative node
 * @param branch_idx Index of inductor branch current variable
 * @param mna_size MNA size
 */
void spice_stamp_inductor_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                              int32_t nplus, int32_t nminus,
                              int32_t branch_idx, int32_t mna_size);

/**
 * @brief Stamp inductor companion model for transient analysis
 *
 * Trapezoidal rule: V = L * dI/dt  →
 *   V_{n+1} = R_eq * I_{n+1} + V_eq
 * where R_eq = 2*L/dt, V_eq = -(R_eq * I_n + V_n)
 *
 * For MNA with branch current I_L:
 *   KVL: V(n+) - V(n-) - R_eq * I_L = V_eq
 *
 * @param G  Conductance matrix
 * @param rhs RHS vector
 * @param nplus  Positive node
 * @param nminus Negative node
 * @param inductance L in henries
 * @param dt Time step
 * @param prev_current I_L at previous step
 * @param prev_voltage V_L at previous step
 * @param branch_idx Index of I_L in MNA system
 * @param mna_size MNA size
 */
void spice_stamp_inductor_tran(spice_dense_matrix_t *G, spice_vector_t *rhs,
                                int32_t nplus, int32_t nminus,
                                double inductance, double dt,
                                double prev_current, double prev_voltage,
                                int32_t branch_idx, int32_t mna_size);

/**
 * @brief Stamp inductor for AC analysis (impedance jωL)
 *
 * Adds the branch equation: V(n+) - V(n-) - jωL * I_L = 0
 *
 * @param G Complex conductance matrix
 * @param rhs Complex RHS
 * @param nplus/nminus Terminal nodes
 * @param inductance L in henries
 * @param omega Angular frequency
 * @param branch_idx Branch current index
 * @param mna_size MNA size
 */
void spice_stamp_inductor_ac(spice_complex_t *G, spice_complex_t *rhs,
                              int32_t nplus, int32_t nminus,
                              double inductance, double omega,
                              int32_t branch_idx, int32_t mna_size);

/* ── L1/L2: Independent Sources ────────────────────────────────────── */

/**
 * @brief Stamp an independent voltage source
 *
 * V-source adds a branch current variable and KVL constraint:
 *   V(n+) - V(n-) = V_dc
 *
 * @param G  Conductance matrix (row/col added for branch current)
 * @param rhs RHS vector
 * @param nplus/nminus Terminal nodes
 * @param dc_value DC voltage value
 * @param branch_idx Index of branch current in MNA
 * @param mna_size MNA size
 */
void spice_stamp_vsource_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                             int32_t nplus, int32_t nminus,
                             double dc_value, int32_t branch_idx,
                             int32_t mna_size);

/**
 * @brief Stamp independent current source
 *
 * I-source contributes to RHS directly:
 *   rhs[n+] -= I_dc, rhs[n-] += I_dc
 *
 * No extra variables needed.
 */
void spice_stamp_isource_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                             int32_t nplus, int32_t nminus,
                             double dc_value, int32_t mna_size);

/* ── L2/L5: Diode Model (Nonlinear) ────────────────────────────────── */

/**
 * @brief Evaluate diode current and conductance
 *
 * Shockley equation: I_D = I_S * (exp(V_D / (n*V_T)) - 1)
 * Small-signal conductance: g_D = dI_D/dV_D = I_S/(n*V_T) * exp(V_D/(n*V_T))
 *
 * @param vd Voltage across diode (V_D = V_anode - V_cathode)
 * @param model Diode model parameters
 * @param temp_celsius Temperature
 * @param out_current Output: diode current I_D
 * @param out_conduct Output: small-signal conductance g_D
 *
 * Includes series resistance effects if model->series_r > 0.
 *
 * Reference: Sedra & Smith §4.3 (Shockley equation)
 */
void spice_diode_evaluate(double vd, const spice_diode_model_t *model,
                           double temp_celsius,
                           double *out_current, double *out_conduct);

/**
 * @brief Stamp a diode for DC operating point (Newton-Raphson)
 *
 * Nonlinear stamp: contributes I_D(V_D) to RHS
 *                and g_D to conductance matrix (companion model)
 *
 * I_D flows from anode (nplus) to cathode (nminus).
 * The companion model is: I = I_D(V_D^{k}) + g_D^{k} * (V_D - V_D^{k})
 *
 * @param G  Conductance matrix (stamped with g_D)
 * @param rhs RHS vector (stamped with -I_D + g_D*V_D)
 * @param nplus  Anode node
 * @param nminus Cathode node
 * @param vd Current estimate of V_D
 * @param model Diode model
 * @param temp_celsius Temperature
 * @param mna_size MNA size
 *
 * Reference: Kielkowski (1998) §4.1
 */
void spice_stamp_diode_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                           int32_t nplus, int32_t nminus,
                           double vd, const spice_diode_model_t *model,
                           double temp_celsius, int32_t mna_size);

/* ── L2/L5: MOSFET Level 1 Model (Nonlinear) ───────────────────────── */

/**
 * @brief Evaluate MOSFET Level 1 drain current and transconductances
 *
 * Shichman-Hodges model:
 *   Cutoff:   V_GS < V_TH  →  I_D = 0
 *   Triode:   V_DS < V_GS - V_TH  →
 *     I_D = K * [(V_GS - V_TH)*V_DS - V_DS^2/2] * (1 + λ*V_DS)
 *   Saturation: V_DS >= V_GS - V_TH  →
 *     I_D = (K/2) * (V_GS - V_TH)^2 * (1 + λ*V_DS)
 *
 * K = KP * (W/L)_eff
 *
 * @param vgs Gate-source voltage
 * @param vds Drain-source voltage
 * @param vbs Bulk-source voltage (for V_TH body effect)
 * @param model MOSFET model parameters
 * @param width Gate width (m)
 * @param length Gate length (m)
 * @param out_id Output: drain current I_D
 * @param out_gm Output: transconductance gm = dI_D/dV_GS
 * @param out_gds Output: output conductance gds = dI_D/dV_DS
 * @param out_gmb Output: body transconductance gmb = dI_D/dV_BS
 *
 * Reference: Shichman & Hodges (1968), "Modeling and Simulation of
 *            Insulated-Gate Field-Effect Transistor Switching Circuits"
 *            IEEE JSSC SC-3(3)
 */
void spice_mosfet_evaluate(double vgs, double vds, double vbs,
                            const spice_mos_model_t *model,
                            double width, double length,
                            double *out_id, double *out_gm,
                            double *out_gds, double *out_gmb);

/**
 * @brief Stamp NMOS Level 1 for DC analysis
 *
 * MOSFET has 4 terminals: D(rain), G(ate), S(ource), B(ulk).
 * NMOS: current flows from drain to source, controlled by V_GS.
 *
 * Companion model is a VCCS: I_DS = f(V_GS, V_DS) with
 * linearization around current operating point.
 *
 * @param G  Conductance matrix
 * @param rhs RHS vector
 * @param nd Drain node
 * @param ng Gate node
 * @param ns Source node
 * @param nb Bulk node
 * @param model MOSFET model
 * @param width Gate width
 * @param length Gate length
 * @param vgs/vds/vbs Current voltage estimates
 * @param temp_celsius Temperature
 * @param mna_size MNA size
 */
void spice_stamp_mosfet_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                            int32_t nd, int32_t ng, int32_t ns, int32_t nb,
                            const spice_mos_model_t *model,
                            double width, double length,
                            double vgs, double vds, double vbs,
                            double temp_celsius, int32_t mna_size);

/* ── L2/L5: BJT Model (Nonlinear) ──────────────────────────────────── */

/**
 * @brief Evaluate BJT terminal currents (Ebers-Moll simplified)
 *
 * Forward-active region:
 *   I_C = I_S * exp(V_BE / (NF * V_T))
 *   I_B = I_C / BF
 *   I_E = -(I_C + I_B)
 *
 * With Early effect: I_C *= (1 + V_CE / VAF)
 *
 * @param vbe Base-emitter voltage
 * @param vbc Base-collector voltage
 * @param model BJT model parameters
 * @param temp_celsius Temperature
 * @param out_ic Collector current
 * @param out_ib Base current
 * @param out_gm Transconductance dI_C/dV_BE
 * @param out_go Output conductance dI_C/dV_CE
 * @param out_gpi Input conductance dI_B/dV_BE
 *
 * Reference: Massobrio & Antognetti (1998) §2
 */
void spice_bjt_evaluate(double vbe, double vbc,
                         const spice_bjt_model_t *model,
                         double temp_celsius,
                         double *out_ic, double *out_ib,
                         double *out_gm, double *out_go, double *out_gpi);

/**
 * @brief Stamp BJT for DC analysis (simplified π-model)
 *
 * BJT has 3 terminals: C(ollector), B(ase), E(mitter).
 * NPN: I_C flows into collector, I_B into base, I_E out of emitter.
 *
 * Companion model uses the hybrid-π small-signal equivalent:
 *   g_m between C and E (controlled by V_BE)
 *   g_π between B and E
 *   g_o  between C and E
 *
 * @param G  Conductance matrix
 * @param rhs RHS vector
 * @param nc Collector node
 * @param nb Base node
 * @param ne Emitter node
 * @param model BJT model
 * @param vbe/vbc Current voltage estimates
 * @param temp_celsius Temperature
 * @param mna_size MNA size
 */
void spice_stamp_bjt_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                         int32_t nc, int32_t nb, int32_t ne,
                         const spice_bjt_model_t *model,
                         double vbe, double vbc,
                         double temp_celsius, int32_t mna_size);

#endif /* SPICE_MODELS_H */
