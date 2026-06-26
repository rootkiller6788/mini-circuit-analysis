/**
 * @file spice_models.c
 * @brief Device model evaluation and MNA stamp implementations
 *
 * Implements the nonlinear device equations and their linearized
 * companion models for MNA stamping.
 *
 * Core device physics:
 *   Diode:   Shockley equation (1949)
 *   MOSFET:  Shichman-Hodges Level 1 model (1968)
 *   BJT:     Simplified Ebers-Moll model (1954)
 *
 * Each nonlinear device contributes:
 *   - I_k (current vector) to the RHS (function evaluation)
 *   - G_k (conductance Jacobian) to the MNA matrix (linearization)
 *
 * Reference texts:
 *   Sedra & Smith, "Microelectronic Circuits" 8th ed. (2020)
 *   Kielkowski, "Inside SPICE" 2nd ed. (1998)
 *   Massobrio & Antognetti, "Semiconductor Device Modeling with SPICE" 2nd ed. (1998)
 */

#include "spice_models.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ═══════════════════════════════════════════════════════════════════════
 * Physical Constants and Utilities
 * ═══════════════════════════════════════════════════════════════════════ */

double spice_thermal_voltage(double temp_celsius) {
    double T_kelvin = temp_celsius + SPICE_KELVIN_OFFSET;
    return SPICE_BOLTZMANN * T_kelvin / SPICE_E_CHARGE;
}

/**
 * @brief Exponential with limiting to prevent overflow in Newton iteration
 *
 * SPICE uses pn-junction limiting: for forward bias > 0, the exponential
 * is limited to prevent divergence during Newton iteration.
 *
 * If V_D > V_crit, use linear extrapolation from the critical point.
 * This is the standard SPICE "pn junction limiting" technique.
 *
 * Reference: Nagel (1975) §6.2.2, "Exponential Limiting"
 */
static double safe_exp(double x, double *out_deriv, double v_crit) {
    if (x < -100.0) {
        /* Deep reverse bias: exp(x) ≈ 0, derivative ≈ 0 */
        if (out_deriv) *out_deriv = 0.0;
        return 0.0;
    }
    if (x > v_crit) {
        /* Forward bias limiting: linear extrapolation */
        double exp_crit = exp(v_crit);
        double val = exp_crit * (1.0 + x - v_crit);
        if (out_deriv) *out_deriv = exp_crit;
        return val;
    }
    double val = exp(x);
    if (out_deriv) *out_deriv = val;
    return val;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Linear Device Stamps
 * ═══════════════════════════════════════════════════════════════════════ */

void spice_stamp_resistor(spice_dense_matrix_t *G, spice_vector_t *rhs,
                           int32_t nplus, int32_t nminus,
                           double resistance, int32_t mna_size) {
    (void)rhs;
    (void)mna_size;
    if (resistance <= 0.0) return;

    double g = 1.0 / resistance;

    /* Stamp conductance into MNA matrix:
     *   G[n+][n+] += g, G[n+][n-] -= g
     *   G[n-][n+] -= g, G[n-][n-] += g
     */
    if (nplus > 0) {
        spice_dense_add(G, nplus - 1, nplus - 1, g);
        if (nminus > 0) spice_dense_add(G, nplus - 1, nminus - 1, -g);
    }
    if (nminus > 0) {
        if (nplus > 0) spice_dense_add(G, nminus - 1, nplus - 1, -g);
        spice_dense_add(G, nminus - 1, nminus - 1, g);
    }
}

void spice_stamp_capacitor_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                               int32_t nplus, int32_t nminus,
                               int32_t mna_size) {
    /* DC: capacitor is an open circuit — no contribution to MNA */
    (void)G; (void)rhs; (void)nplus; (void)nminus; (void)mna_size;
}

void spice_stamp_capacitor_tran(spice_dense_matrix_t *G, spice_vector_t *rhs,
                                 int32_t nplus, int32_t nminus,
                                 double capacitance, double dt,
                                 double prev_voltage, double prev_current,
                                 int32_t mna_size) {
    (void)mna_size;
    if (capacitance <= 0.0 || dt <= 0.0) return;

    /* Trapezoidal companion model:
     *   I_{n+1} = (2C/dt) * V_{n+1} + (-2C/dt * V_n - I_n)
     *            = G_eq * V_{n+1} + I_eq
     *
     * G_eq = 2*C/dt
     * I_eq = -(2*C/dt) * V_n - I_n
     */
    double g_eq = 2.0 * capacitance / dt;
    double i_eq = -g_eq * prev_voltage - prev_current;

    /* Stamp G_eq */
    if (nplus > 0) {
        spice_dense_add(G, nplus - 1, nplus - 1, g_eq);
        if (nminus > 0) spice_dense_add(G, nplus - 1, nminus - 1, -g_eq);
        /* RHS contribution: -I_eq flows out of n+, into n- */
        if (rhs) rhs->data[nplus - 1] += i_eq;
    }
    if (nminus > 0) {
        if (nplus > 0) spice_dense_add(G, nminus - 1, nplus - 1, -g_eq);
        spice_dense_add(G, nminus - 1, nminus - 1, g_eq);
        if (rhs) rhs->data[nminus - 1] -= i_eq;
    }
}

void spice_stamp_capacitor_ac(spice_complex_t *G, spice_complex_t *rhs,
                               int32_t nplus, int32_t nminus,
                               double capacitance, double omega,
                               int32_t mna_size) {
    (void)rhs;
    (void)mna_size;
    if (capacitance <= 0.0) return;

    /* AC admittance: Y = jωC between n+ and n- */
    spice_complex_t Y = capacitance * omega * I;

    if (nplus > 0) {
        G[(nplus - 1) + (nplus - 1) * (size_t)mna_size] += Y;
        if (nminus > 0)
            G[(nplus - 1) + (nminus - 1) * (size_t)mna_size] -= Y;
    }
    if (nminus > 0) {
        if (nplus > 0)
            G[(nminus - 1) + (nplus - 1) * (size_t)mna_size] -= Y;
        G[(nminus - 1) + (nminus - 1) * (size_t)mna_size] += Y;
    }
}

void spice_stamp_inductor_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                              int32_t nplus, int32_t nminus,
                              int32_t branch_idx, int32_t mna_size) {
    (void)mna_size;
    /* DC: inductor acts as short circuit.
     * In MNA, this means V(n+) - V(n-) = 0, enforced by branch current I_L.
     *
     * KVL row (branch_idx): V(n+) - V(n-) = 0
     * KCL columns: I_L contributes +1 at n+, -1 at n-.
     *
     * Stamp:
     *   Row branch_idx:  +1 at col n+, -1 at col n-  → RHS = 0
     *   Col branch_idx:  +1 at row n+, -1 at row n-  → 0 (for DC, R=0)
     */
    int32_t br = branch_idx; /* Row/col index in MNA (0-based) */
    if (nplus > 0) {
        spice_dense_set(G, br, nplus - 1, 1.0);
        spice_dense_set(G, nplus - 1, br, 1.0);
    }
    if (nminus > 0) {
        spice_dense_set(G, br, nminus - 1, -1.0);
        spice_dense_set(G, nminus - 1, br, -1.0);
    }
    /* RHS: 0 for short circuit */
    if (rhs) rhs->data[br] = 0.0;
}

void spice_stamp_inductor_tran(spice_dense_matrix_t *G, spice_vector_t *rhs,
                                int32_t nplus, int32_t nminus,
                                double inductance, double dt,
                                double prev_current, double prev_voltage,
                                int32_t branch_idx, int32_t mna_size) {
    (void)mna_size;
    if (inductance <= 0.0 || dt <= 0.0) return;

    /* Trapezoidal companion:
     *   V_L = L * dI/dt
     *   V_{n+1} = R_eq * I_{n+1} + V_eq
     *   R_eq = 2*L/dt,  V_eq = -(R_eq * I_n + V_n)
     *
     * KVL constraint (row branch_idx):
     *   V(n+) - V(n-) - R_eq * I_L = V_eq
     *
     * KCL (column branch_idx): I_L adds to node equations
     */
    double r_eq = 2.0 * inductance / dt;
    double v_eq = -(r_eq * prev_current + prev_voltage);

    int32_t br = branch_idx;

    /* Row br: KVL equation */
    if (nplus > 0)  spice_dense_set(G, br, nplus - 1,  1.0);
    if (nminus > 0) spice_dense_set(G, br, nminus - 1, -1.0);
    spice_dense_set(G, br, br, -r_eq);
    if (rhs) rhs->data[br] = v_eq;

    /* Column br: KCL contribution */
    if (nplus > 0)  spice_dense_set(G, nplus - 1, br,  1.0);
    if (nminus > 0) spice_dense_set(G, nminus - 1, br, -1.0);
}

void spice_stamp_inductor_ac(spice_complex_t *G, spice_complex_t *rhs,
                              int32_t nplus, int32_t nminus,
                              double inductance, double omega,
                              int32_t branch_idx, int32_t mna_size) {
    if (inductance <= 0.0) return;

    /* AC: impedance Z_L = jωL
     * KVL: V(n+) - V(n-) - jωL * I_L = 0
     */
    spice_complex_t Z = inductance * omega * I;
    int32_t br = branch_idx;
    (void)rhs;

    if (nplus > 0)  G[br + (nplus - 1) * (size_t)mna_size]  =  1.0;
    if (nminus > 0) G[br + (nminus - 1) * (size_t)mna_size] = -1.0;
    G[br + br * (size_t)mna_size] = -Z;

    if (nplus > 0)  G[(nplus - 1) + br * (size_t)mna_size]  =  1.0;
    if (nminus > 0) G[(nminus - 1) + br * (size_t)mna_size] = -1.0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Independent Source Stamps
 * ═══════════════════════════════════════════════════════════════════════ */

void spice_stamp_vsource_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                             int32_t nplus, int32_t nminus,
                             double dc_value, int32_t branch_idx,
                             int32_t mna_size) {
    int32_t br = branch_idx;
    (void)mna_size;

    /* KVL: V(n+) - V(n-) = V_dc
     * Row br: +1 at n+, -1 at n-
     * RHS: V_dc */
    if (nplus > 0)  spice_dense_set(G, br, nplus - 1,  1.0);
    if (nminus > 0) spice_dense_set(G, br, nminus - 1, -1.0);
    if (rhs) rhs->data[br] = dc_value;

    /* KCL: I_V flows into n+, out of n-
     * Col br: +1 at n+, -1 at n- */
    if (nplus > 0)  spice_dense_set(G, nplus - 1, br,  1.0);
    if (nminus > 0) spice_dense_set(G, nminus - 1, br, -1.0);
}

void spice_stamp_isource_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                             int32_t nplus, int32_t nminus,
                             double dc_value, int32_t mna_size) {
    (void)G;
    (void)mna_size;

    /* Current source: injects I_dc into n+, extracts from n-.
     * KCL at n+:   ... + I_source = 0  → RHS contribution
     * KCL at n-:   ... - I_source = 0
     *
     * In MNA, the convention is that RHS = -I_source for nodes
     * where current is injected. So:
     *   RHS[n+] += I_dc  (current source pushes current into node)
     *   RHS[n-] -= I_dc
     *
     * Actually: The MNA equation is G * V = J, where J contains
     * independent current sources flowing INTO each node.
     * If I_dc flows from n+ to n-, then J[n+] -= I_dc, J[n-] += I_dc.
     *
     * Wait - standard SPICE convention:
     * Current flowing INTO a node is positive for RHS.
     * A current source I from n+ to n- means:
     *   Into n+: 0 (current flows out)
     *   Into n-: I (current flows in)
     * So RHS[n+] stays unchanged, RHS[n-] += I_dc.
     * But we must account for sign: the I-source has positive terminal at n+.
     * Current leaves n+ and enters n-. So:
     *   RHS[n+] -= I_dc (current leaving n+)
     *   RHS[n-] += I_dc (current entering n-)
     */
    if (nplus > 0 && rhs)  rhs->data[nplus - 1]  -= dc_value;
    if (nminus > 0 && rhs) rhs->data[nminus - 1] += dc_value;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Diode Model: Shockley Equation
 * ═══════════════════════════════════════════════════════════════════════ */

void spice_diode_evaluate(double vd, const spice_diode_model_t *model,
                           double temp_celsius,
                           double *out_current, double *out_conduct) {
    if (!model) {
        if (out_current) *out_current = 0.0;
        if (out_conduct) *out_conduct = 0.0;
        return;
    }

    double vt = spice_thermal_voltage(temp_celsius);
    double nvt = model->n_ideality * vt;

    if (nvt <= 0.0) {
        if (out_current) *out_current = 0.0;
        if (out_conduct) *out_conduct = 0.0;
        return;
    }

    /* If series resistance R_S > 0, the internal junction voltage
     * V_J = V_D - I_D * R_S. We solve iteratively for simplicity:
     * use the external V_D as first estimate. */
    double vj = vd;
    double id, gd;

    /* Simple iteration for series resistance effect */
    for (int iter = 0; iter < 5; iter++) {
        double exp_deriv;
        /* Exponential with limiting at 40/nVT to prevent overflow */
        double v_crit = 40.0 / (1.0 / nvt); /* ≈ 40 * nVT */
        if (v_crit > 80.0) v_crit = 80.0;

        double exp_v = safe_exp(vj / nvt, &exp_deriv, v_crit / nvt);

        id = model->is_saturation * (exp_v - 1.0);
        gd = model->is_saturation * exp_deriv / nvt;

        if (model->series_r <= 0.0) break;

        /* Update V_J = V_D - I_D * R_S */
        vj = vd - id * model->series_r;
    }

    /* Breakdown region: for V_D < -BV, current increases sharply */
    if (vd < -model->breakdown_v && model->breakdown_v < 1e29) {
        double v_bd = -vd - model->breakdown_v;
        if (v_bd > 0) {
            double ibv_contrib = model->breakdown_i *
                (exp(v_bd / (model->n_ideality * vt)) - 1.0);
            id -= ibv_contrib;
            gd += model->breakdown_i *
                exp(v_bd / (model->n_ideality * vt)) /
                (model->n_ideality * vt);
        }
    }

    if (out_current) *out_current = id;
    if (out_conduct) *out_conduct = gd;
}

void spice_stamp_diode_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                           int32_t nplus, int32_t nminus,
                           double vd, const spice_diode_model_t *model,
                           double temp_celsius, int32_t mna_size) {
    (void)mna_size;
    if (!model) return;

    double id, gd;
    spice_diode_evaluate(vd, model, temp_celsius, &id, &gd);

    /* Companion model (Newton linearization at V_D^k):
     *   I = I_D(V_D^k) + g_D^k * (V_D - V_D^k)
     *     = (I_D^k - g_D^k * V_D^k) + g_D^k * V_D
     *
     * MNA stamp: conductance gd between n+ and n-
     *            current source I_eq = I_D^k - g_D^k * V_D^k
     *
     * Current flows from anode (n+) to cathode (n-).
     * Into n+: -(I_D)  → RHS contribution
     * Into n-: +(I_D)  → RHS contribution
     *
     * Newton form: RHS gets -(I_D - g_D*V_D) which = -I_D + g_D*V_D
     *
     * Actually in SPICE: the stamp adds to the RHS:
     *   RHS[n+] += -(I_D - g_D * V_D)
     *   RHS[n-] += +(I_D - g_D * V_D)
     *
     * Because the KCL equation at n+ is:
     *   Σ G*V + g_D*(V_n+ - V_n-) = -(I_D - g_D*V_D)
     * → g_D*(V_n+ - V_n-) appears on LHS (stamped to G matrix)
     * → -(I_D - g_D*V_D) appears on RHS
     */
    double i_eq = id - gd * vd;  /* Companion current source */

    /* Stamp conductance g_D */
    if (nplus > 0) {
        spice_dense_add(G, nplus - 1, nplus - 1, gd);
        if (nminus > 0) spice_dense_add(G, nplus - 1, nminus - 1, -gd);
        if (rhs) rhs->data[nplus - 1] -= i_eq;
    }
    if (nminus > 0) {
        if (nplus > 0) spice_dense_add(G, nminus - 1, nplus - 1, -gd);
        spice_dense_add(G, nminus - 1, nminus - 1, gd);
        if (rhs) rhs->data[nminus - 1] += i_eq;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * MOSFET Level 1 Model: Shichman-Hodges
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Compute effective threshold voltage with body effect
 *
 * V_TH = V_TO + GAMMA * (sqrt(PHI - V_BS) - sqrt(PHI))
 *
 * Only valid for V_BS <= 0 (standard body bias).
 * For PMOS, all signs are reversed in the calling code.
 *
 * Reference: Shichman & Hodges (1968), Eq. (4)
 */
static double mos_vth(const spice_mos_model_t *model, double vbs) {
    if (!model) return 0.0;
    double phi = model->phi_surface;
    if (phi <= 0.0) phi = 0.7;
    /* Prevent negative argument to sqrt */
    double arg = phi - vbs;
    if (arg <= 0.0) arg = 0.0;
    return model->vth0 + model->gamma_body * (sqrt(arg) - sqrt(phi));
}

void spice_mosfet_evaluate(double vgs, double vds, double vbs,
                            const spice_mos_model_t *model,
                            double width, double length,
                            double *out_id, double *out_gm,
                            double *out_gds, double *out_gmb) {
    /* Default outputs */
    double id = 0.0, gm = 0.0, gds = 0.0, gmb = 0.0;

    if (!model || length <= 0.0 || width <= 0.0) goto done;

    double vth = mos_vth(model, vbs);
    double k = model->kp * (width / length);

    /* For PMOS, we use absolute values internally.
     * The model stores VTO as negative for PMOS; we use fabs and swap signs. */
    int is_pmos = (model->type == 1);
    double vgs_eff = is_pmos ? -vgs : vgs;
    double vds_eff = is_pmos ? -vds : vds;
    double vth_eff = is_pmos ? -vth : vth;

    if (vgs_eff <= vth_eff) {
        /* Cutoff region */
        id = 0.0;
        gm = 0.0;
        gds = 0.0;
        gmb = 0.0;
        goto done;
    }

    double vdsat = vgs_eff - vth_eff;

    if (vds_eff < vdsat) {
        /* Triode (linear) region:
         * I_D = K * [(V_GS - V_TH)*V_DS - V_DS^2/2] * (1 + λ*V_DS)
         */
        double triode_term = (vgs_eff - vth_eff) * vds_eff - 0.5 * vds_eff * vds_eff;
        double chan_mod = 1.0 + model->lambda_channel * vds_eff;
        id = k * triode_term * chan_mod;

        /* gm = dI_D/dV_GS = K * V_DS * (1 + λ*V_DS) */
        gm = k * vds_eff * chan_mod;

        /* gds = dI_D/dV_DS = K * [(V_GS - V_TH - V_DS)*(1+λ*V_DS) + λ*triode_term] */
        gds = k * ((vgs_eff - vth_eff - vds_eff) * chan_mod +
                   model->lambda_channel * triode_term);

    } else {
        /* Saturation region:
         * I_D = (K/2) * (V_GS - V_TH)^2 * (1 + λ*V_DS)
         */
        double vov = vgs_eff - vth_eff;  /* Overdrive voltage */
        double chan_mod = 1.0 + model->lambda_channel * vds_eff;
        id = 0.5 * k * vov * vov * chan_mod;

        /* gm = K * (V_GS - V_TH) * (1 + λ*V_DS) */
        gm = k * vov * chan_mod;

        /* gds = (K/2) * (V_GS - V_TH)^2 * λ */
        gds = 0.5 * k * vov * vov * model->lambda_channel;
    }

    /* gmb = dI_D/dV_BS = dI_D/dV_TH * dV_TH/dV_BS = -gm * dV_TH/dV_BS
     *
     * dV_TH/dV_BS = -GAMMA / (2 * sqrt(PHI - V_BS))
     * So gmb = gm * GAMMA / (2 * sqrt(PHI - V_BS))
     */
    {
        double phi = model->phi_surface;
        if (phi <= 0.0) phi = 0.7;
        double arg = phi - vbs;
        if (arg <= 0.0) arg = 1e-3;  /* Avoid div by zero */
        gmb = gm * model->gamma_body / (2.0 * sqrt(arg));
    }

    /* For PMOS, flip current direction and transconductance signs */
    if (is_pmos) {
        id = -id;
        /* For PMOS: gm and gmb are defined with respect to V_SG and V_SB, keep positive */
    }

done:
    if (out_id)  *out_id  = id;
    if (out_gm)  *out_gm  = gm;
    if (out_gds) *out_gds = gds;
    if (out_gmb) *out_gmb = gmb;
}

void spice_stamp_mosfet_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                            int32_t nd, int32_t ng, int32_t ns, int32_t nb,
                            const spice_mos_model_t *model,
                            double width, double length,
                            double vgs, double vds, double vbs,
                            double temp_celsius, int32_t mna_size) {
    (void)temp_celsius;
    (void)mna_size;

    if (!model) return;

    double id, gm, gds, gmb;
    spice_mosfet_evaluate(vgs, vds, vbs, model, width, length,
                           &id, &gm, &gds, &gmb);

    /* MOSFET stamp uses a VCCS companion model:
     *
     * I_DS flows from drain (nd) to source (ns).
     * Linearized: I_DS = I_D^0 + gm*(v_GS - v_GS^0) + gds*(v_DS - v_DS^0) + gmb*(v_BS - v_BS^0)
     *
     * = I_eq + gm*v_G + (gds - gm)*v_S - gds*v_D + gmb*(v_B - v_S)
     *   ... where I_eq = I_D^0 - gm*V_GS^0 - gds*V_DS^0 - gmb*V_BS^0
     *
     * Stamp to G matrix (conductance contributions):
     *   G[d][d] += gds
     *   G[d][s] -= gds
     *   G[s][d] -= gds
     *   G[s][s] += gds
     *
     * gm between d-s controlled by g: +gm at (d,g), -gm at (s,g)
     *   G[d][g] += gm
     *   G[s][g] -= gm
     *
     * gmb between d-s controlled by b: +gmb at (d,b), -gmb at (s,b)
     *   G[d][b] += gmb
     *   G[s][b] -= gmb
     *
     * RHS: companion current source I_eq
     *   Into drain: +I_eq
     *   Into source: -I_eq
     *
     * But in SPICE MNA convention:
     *   The RHS represents -I (net current into node is 0 → Σ = 0)
     *   So if I_DS flows d→s:
     *     RHS[d] -= I_DS (current leaves drain)
     *     RHS[s] += I_DS (current enters source)
     *
     * For Newton: I_DS is replaced by linearized form,
     * the constant part goes to RHS, the voltage-dependent part to G.
     */
    double i_eq = id - gm * vgs - gds * vds - gmb * vbs;

    /* Drain node stamp */
    if (nd > 0) {
        /* G matrix */
        spice_dense_add(G, nd - 1, nd - 1, gds);
        if (ns > 0) spice_dense_add(G, nd - 1, ns - 1, -gds);
        if (ng > 0) spice_dense_add(G, nd - 1, ng - 1, gm);
        if (nb > 0) spice_dense_add(G, nd - 1, nb - 1, gmb);
        /* RHS: current leaving drain */
        if (rhs) rhs->data[nd - 1] -= i_eq;
    }

    /* Source node stamp */
    if (ns > 0) {
        if (nd > 0) spice_dense_add(G, ns - 1, nd - 1, -gds);
        spice_dense_add(G, ns - 1, ns - 1, gds);
        if (ng > 0) spice_dense_add(G, ns - 1, ng - 1, -gm);
        if (nb > 0) spice_dense_add(G, ns - 1, nb - 1, -gmb);
        /* RHS: current entering source */
        if (rhs) rhs->data[ns - 1] += i_eq;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * BJT Model: Simplified Ebers-Moll (Forward Active)
 * ═══════════════════════════════════════════════════════════════════════ */

void spice_bjt_evaluate(double vbe, double vbc,
                         const spice_bjt_model_t *model,
                         double temp_celsius,
                         double *out_ic, double *out_ib,
                         double *out_gm, double *out_go, double *out_gpi) {
    double ic = 0.0, ib = 0.0, gm = 0.0, gpi = 0.0, go = 0.0;

    if (!model) goto done;

    double vt = spice_thermal_voltage(temp_celsius);

    /* Forward-active region:
     *   I_C = I_S * exp(V_BE / (NF * V_T)) * (1 + V_CE / VAF)
     *   I_B = I_C / BF
     *
     * where V_CE = V_BE - V_BC (approximately)
     */
    double nf_vt = model->nf_coeff * vt;
    if (nf_vt <= 0.0) goto done;

    double exp_deriv;
    double vbe_crit = 80.0;
    double exp_be = safe_exp(vbe / nf_vt, &exp_deriv, vbe_crit / nf_vt);

    double vce = vbe - vbc;
    double early_factor = 1.0;
    if (model->vaf_early > 0.0 && model->vaf_early < 1e29) {
        early_factor = 1.0 + vce / model->vaf_early;
        if (early_factor < 0.1) early_factor = 0.1;  /* Prevent negative */
    }

    ic = model->is_saturation * exp_be * early_factor;
    ib = ic / model->bf_forward;

    /* Transconductance gm = dI_C/dV_BE ≈ I_C / (NF * V_T) */
    gm = (model->is_saturation / nf_vt) * exp_deriv * early_factor;

    /* Output conductance go = dI_C/dV_CE
     *   I_C = I_S * exp(V_BE/NF·VT) * (1 + V_CE/VAF)
     *   go = I_S * exp(...) / VAF ≈ I_C / VAF for large VAF */
    if (model->vaf_early > 0.0 && model->vaf_early < 1e29) {
        go = ic / (model->vaf_early + vce);
    }

    /* Input conductance gπ = dI_B/dV_BE = gm / BF */
    gpi = gm / model->bf_forward;

done:
    if (out_ic)  *out_ic  = ic;
    if (out_ib)  *out_ib  = ib;
    if (out_gm)  *out_gm  = gm;
    if (out_go)  *out_go  = go;
    if (out_gpi) *out_gpi = gpi;
}

void spice_stamp_bjt_dc(spice_dense_matrix_t *G, spice_vector_t *rhs,
                         int32_t nc, int32_t nb, int32_t ne,
                         const spice_bjt_model_t *model,
                         double vbe, double vbc,
                         double temp_celsius, int32_t mna_size) {
    (void)mna_size;
    if (!model) return;

    double ic, ib, gm, go, gpi;
    spice_bjt_evaluate(vbe, vbc, model, temp_celsius,
                        &ic, &ib, &gm, &go, &gpi);

    /* Hybrid-π companion model (simplified for NPN):
     *
     * Between B and E: gπ (input conductance), with companion current
     * Between C and E: go (output conductance)
     * VCCS: gm * V_BE from C to E (controlled by V_BE)
     *
     * I_B flows into base: I_B = gπ*(V_B - V_E) + I_B_eq
     * I_C flows into collector: I_C = go*(V_C - V_E) + gm*(V_B - V_E) + I_C_eq
     *
     * Companion currents (Newton linearization):
     *   I_B_eq = I_B^0 - gπ * V_BE^0
     *   I_C_eq = I_C^0 - gm * V_BE^0 - go * V_CE^0
     *
     * KCL at each node:
     *   Node B: I_B + ... = 0 → gπ*(V_B-V_E) + ... = -I_B_eq
     *   Node C: I_C + ... = 0 → go*(V_C-V_E) + gm*(V_B-V_E) + ... = -I_C_eq
     *   Node E: -(I_B+I_C) + ... = 0
     */
    double ib_eq = ib - gpi * vbe;
    double vce = vbe - vbc;
    double ic_eq = ic - gm * vbe - go * vce;

    /* ── Base node (nb) ── */
    if (nb > 0) {
        spice_dense_add(G, nb - 1, nb - 1, gpi);
        if (ne > 0) spice_dense_add(G, nb - 1, ne - 1, -gpi);
        if (rhs) rhs->data[nb - 1] -= ib_eq;
    }

    /* ── Collector node (nc) ── */
    if (nc > 0) {
        /* go between C and E */
        spice_dense_add(G, nc - 1, nc - 1, go);
        if (ne > 0) spice_dense_add(G, nc - 1, ne - 1, -go);
        /* gm from base to collector (VCCS): I_C += gm * (V_B - V_E) */
        if (nb > 0) spice_dense_add(G, nc - 1, nb - 1, gm);
        if (ne > 0) spice_dense_add(G, nc - 1, ne - 1, -gm);
        /* RHS */
        if (rhs) rhs->data[nc - 1] -= ic_eq;
    }

    /* ── Emitter node (ne) ── */
    if (ne > 0) {
        /* go and gm contributions (reverse sign) */
        if (nc > 0) spice_dense_add(G, ne - 1, nc - 1, -go);
        spice_dense_add(G, ne - 1, ne - 1, go);
        if (nb > 0) {
            spice_dense_add(G, ne - 1, nb - 1, -gpi - gm);
        }
        spice_dense_add(G, ne - 1, ne - 1, gpi + go + gm);
        /* RHS: -(I_B+I_C) flows out of emitter */
        if (rhs) rhs->data[ne - 1] += (ib_eq + ic_eq);
    }
}
