#include "circuit_elements.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*==================================================================
 * Element Initialization Functions (L1: Core Definitions)
 * Each function initializes a circuit element with the defining
 * physical parameters and terminal connections.
 *==================================================================*/

int resistor_init(Resistor_t *r, int id, const char *name,
                  int n1, int n2, double R_ohms)
{
    if (!r) return -1;
    if (R_ohms < 0.0) return -1;
    if (n1 < 0 || n2 < 0) return -1;

    memset(r, 0, sizeof(*r));
    r->id = id;
    r->type = ELEM_RESISTOR;
    r->R = R_ohms;
    r->tolerance = 0.05;   /* default 5% */
    r->temp_coeff = 100.0; /* default 100 ppm/C for carbon film */
    r->power_rating = 0.25;/* default 0.25W */

    if (name) {
        strncpy(r->name, name, MAX_NAME_LEN - 1);
        r->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(r->name, MAX_NAME_LEN, "R%d", id);
    }

    r->t[0].node_id = n1;
    r->t[1].node_id = n2;
    snprintf(r->t[0].label, MAX_NAME_LEN, "R%d_t1", id);
    snprintf(r->t[1].label, MAX_NAME_LEN, "R%d_t2", id);

    return 0;
}

int capacitor_init(Capacitor_t *c, int id, const char *name,
                   int n1, int n2, double C_farads)
{
    if (!c) return -1;
    if (C_farads < 0.0) return -1;
    if (n1 < 0 || n2 < 0) return -1;

    memset(c, 0, sizeof(*c));
    c->id = id;
    c->type = ELEM_CAPACITOR;
    c->C = C_farads;
    c->tolerance = 0.10;    /* default 10% for electrolytic */
    c->temp_coeff = 200.0;
    c->voltage_rating = 50.0;
    c->esr = 0.1;           /* typical 0.1 Ohm ESR */
    c->esl = 1e-9;          /* typical 1 nH ESL */
    c->leakage_r = 1e6;     /* typical 1 MOhm leakage */

    if (name) {
        strncpy(c->name, name, MAX_NAME_LEN - 1);
        c->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(c->name, MAX_NAME_LEN, "C%d", id);
    }

    c->t[0].node_id = n1;
    c->t[1].node_id = n2;
    snprintf(c->t[0].label, MAX_NAME_LEN, "C%d_t1", id);
    snprintf(c->t[1].label, MAX_NAME_LEN, "C%d_t2", id);

    return 0;
}

int inductor_init(Inductor_t *l, int id, const char *name,
                  int n1, int n2, double L_henries)
{
    if (!l) return -1;
    if (L_henries < 0.0) return -1;
    if (n1 < 0 || n2 < 0) return -1;

    memset(l, 0, sizeof(*l));
    l->id = id;
    l->type = ELEM_INDUCTOR;
    l->L = L_henries;
    l->tolerance = 0.10;
    l->temp_coeff = 500.0;
    l->current_rating = 1.0;
    l->dcr = 0.5;           /* typical 0.5 Ohm DCR */
    l->core_material = 1.0; /* ferrite */

    if (name) {
        strncpy(l->name, name, MAX_NAME_LEN - 1);
        l->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(l->name, MAX_NAME_LEN, "L%d", id);
    }

    l->t[0].node_id = n1;
    l->t[1].node_id = n2;
    snprintf(l->t[0].label, MAX_NAME_LEN, "L%d_t1", id);
    snprintf(l->t[1].label, MAX_NAME_LEN, "L%d_t2", id);

    return 0;
}

int dc_voltage_source_init(DCVoltageSource_t *s, int id, const char *name,
                           int n_plus, int n_minus, double V)
{
    if (!s) return -1;
    if (n_plus < 0 || n_minus < 0) return -1;

    memset(s, 0, sizeof(*s));
    s->id = id;
    s->type = ELEM_DC_VOLTAGE_SRC;
    s->V_dc = V;
    s->R_internal = 0.0; /* ideal */

    if (name) {
        strncpy(s->name, name, MAX_NAME_LEN - 1);
        s->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(s->name, MAX_NAME_LEN, "Vs_dc%d", id);
    }

    s->t[0].node_id = n_plus;   /* positive terminal */
    s->t[1].node_id = n_minus;  /* negative terminal */
    snprintf(s->t[0].label, MAX_NAME_LEN, "Vs%d_p", id);
    snprintf(s->t[1].label, MAX_NAME_LEN, "Vs%d_n", id);

    return 0;
}

int dc_current_source_init(DCCurrentSource_t *s, int id, const char *name,
                           int n_from, int n_to, double I_val)
{
    if (!s) return -1;
    if (n_from < 0 || n_to < 0) return -1;

    memset(s, 0, sizeof(*s));
    s->id = id;
    s->type = ELEM_DC_CURRENT_SRC;
    s->I_dc = I_val;
    s->R_internal = 1e12; /* near-infinite internal resistance */

    if (name) {
        strncpy(s->name, name, MAX_NAME_LEN - 1);
        s->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(s->name, MAX_NAME_LEN, "Is_dc%d", id);
    }

    /* Current flows from node n_from into node n_to through the source */
    s->t[0].node_id = n_from;
    s->t[1].node_id = n_to;
    snprintf(s->t[0].label, MAX_NAME_LEN, "Is%d_f", id);
    snprintf(s->t[1].label, MAX_NAME_LEN, "Is%d_t", id);

    return 0;
}

int ac_voltage_source_init(ACVoltageSource_t *s, int id, const char *name,
                           int n_plus, int n_minus,
                           double V_amp, double freq, double phase_deg,
                           double V_offset)
{
    if (!s) return -1;
    if (V_amp < 0.0) return -1;
    if (freq < 0.0) return -1;
    if (n_plus < 0 || n_minus < 0) return -1;

    memset(s, 0, sizeof(*s));
    s->id = id;
    s->type = ELEM_AC_VOLTAGE_SRC;
    s->V_amplitude = V_amp;
    s->V_rms = V_amp / sqrt(2.0);
    s->freq_hz = freq;
    s->phase_deg = phase_deg;
    s->V_offset = V_offset;
    s->R_internal = 0.0;

    if (name) {
        strncpy(s->name, name, MAX_NAME_LEN - 1);
        s->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(s->name, MAX_NAME_LEN, "Vs_ac%d", id);
    }

    s->t[0].node_id = n_plus;
    s->t[1].node_id = n_minus;

    return 0;
}

int ac_current_source_init(ACCurrentSource_t *s, int id, const char *name,
                           int n_from, int n_to,
                           double I_amp, double freq, double phase_deg,
                           double I_offset)
{
    if (!s) return -1;
    if (I_amp < 0.0) return -1;
    if (n_from < 0 || n_to < 0) return -1;

    memset(s, 0, sizeof(*s));
    s->id = id;
    s->type = ELEM_AC_CURRENT_SRC;
    s->I_amplitude = I_amp;
    s->I_rms = I_amp / sqrt(2.0);
    s->freq_hz = freq;
    s->phase_deg = phase_deg;
    s->I_offset = I_offset;
    s->R_internal = 1e12;

    if (name) {
        strncpy(s->name, name, MAX_NAME_LEN - 1);
        s->name[MAX_NAME_LEN - 1] = '\0';
    }

    s->t[0].node_id = n_from;
    s->t[1].node_id = n_to;

    return 0;
}

int opamp_ideal_init(OpAmpIdeal_t *oa, int id, const char *name,
                     int n_plus, int n_minus, int n_out, double A_ol)
{
    if (!oa) return -1;
    if (A_ol <= 0.0) return -1;

    memset(oa, 0, sizeof(*oa));
    oa->id = id;
    oa->type = ELEM_OPAMP_IDEAL;
    oa->A_ol = A_ol;
    oa->gbw = 1e6; /* default 1 MHz GBW */

    if (name) {
        strncpy(oa->name, name, MAX_NAME_LEN - 1);
        oa->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(oa->name, MAX_NAME_LEN, "OA%d", id);
    }

    oa->t[0].node_id = n_plus;
    oa->t[1].node_id = n_minus;
    oa->t[2].node_id = n_out;

    return 0;
}

int transformer_ideal_init(TransformerIdeal_t *tf, int id, const char *name,
                           int n_p1, int n_p2, int n_s1, int n_s2, double n_turns)
{
    if (!tf) return -1;
    if (n_turns <= 0.0) return -1;

    memset(tf, 0, sizeof(*tf));
    tf->id = id;
    tf->type = ELEM_TRANSFORMER_IDEAL;
    tf->n = n_turns;
    tf->Lm = 1e-3;   /* default 1 mH magnetizing inductance */
    tf->L_leak_pri = 1e-6;
    tf->L_leak_sec = 1e-6;

    if (name) {
        strncpy(tf->name, name, MAX_NAME_LEN - 1);
        tf->name[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(tf->name, MAX_NAME_LEN, "TF%d", id);
    }

    tf->t[0].node_id = n_p1;
    tf->t[1].node_id = n_p2;
    tf->t[2].node_id = n_s1;
    tf->t[3].node_id = n_s2;

    return 0;
}

/*==================================================================
 * Circuit Netlist Management (L1: graph theory for circuits)
 *==================================================================*/

void circuit_init(Circuit_t *ckt, const char *name)
{
    if (!ckt) return;

    memset(ckt, 0, sizeof(*ckt));

    if (name) {
        strncpy(ckt->name, name, MAX_NAME_LEN - 1);
        ckt->name[MAX_NAME_LEN - 1] = '\0';
    }

    /* Create ground node (node 0) by default */
    circuit_add_node(ckt, "GND", 1);
    ckt->ground_node = 0;
}

int circuit_add_node(Circuit_t *ckt, const char *label, int is_ground)
{
    if (!ckt) return -1;
    if (ckt->n_nodes >= MAX_NODES) return -1;

    int idx = ckt->n_nodes;
    Node_t *n = &ckt->nodes[idx];

    n->id = idx;
    if (label) {
        strncpy(n->label, label, MAX_NAME_LEN - 1);
        n->label[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(n->label, MAX_NAME_LEN, "N%d", idx);
    }
    n->is_ground = is_ground;
    n->branch_count = 0;

    if (is_ground) {
        ckt->ground_node = idx;
    }

    ckt->n_nodes++;
    return idx;
}

int circuit_find_node(const Circuit_t *ckt, const char *label)
{
    if (!ckt || !label) return -1;

    for (int i = 0; i < ckt->n_nodes; i++) {
        if (strcmp(ckt->nodes[i].label, label) == 0) {
            return i;
        }
    }
    return -1;
}

int circuit_add_branch(Circuit_t *ckt, int from, int to,
                       ElementType_t type, int elem_idx, double value,
                       const char *label)
{
    if (!ckt) return -1;
    if (ckt->n_branches >= MAX_BRANCHES) return -1;
    if (from < 0 || from >= ckt->n_nodes) return -1;
    if (to < 0 || to >= ckt->n_nodes) return -1;

    int idx = ckt->n_branches;
    Branch_t *b = &ckt->branches[idx];

    b->id = idx;
    b->node_from = from;
    b->node_to = to;
    b->elem_type = type;
    b->elem_index = elem_idx;
    b->value = value;

    if (label) {
        strncpy(b->label, label, MAX_NAME_LEN - 1);
        b->label[MAX_NAME_LEN - 1] = '\0';
    } else {
        snprintf(b->label, MAX_NAME_LEN, "B%d", idx);
    }

    /* Update node incidence */
    Node_t *nf = &ckt->nodes[from];
    Node_t *nt = &ckt->nodes[to];
    if (nf->branch_count < MAX_BRANCHES) {
        nf->branches[nf->branch_count++] = idx;
    }
    if (nt->branch_count < MAX_BRANCHES) {
        nt->branches[nt->branch_count++] = idx;
    }

    ckt->n_branches++;
    return idx;
}

int circuit_validate(const Circuit_t *ckt)
{
    if (!ckt) return -1;

    /* Check that ground node exists */
    int has_ground = 0;
    for (int i = 0; i < ckt->n_nodes; i++) {
        if (ckt->nodes[i].is_ground) {
            has_ground = 1;
            break;
        }
    }
    if (!has_ground) return -1;

    /* Check that all nodes have at least one connection (no floating nodes) */
    for (int i = 0; i < ckt->n_nodes; i++) {
        if (ckt->nodes[i].is_ground) continue;
        if (ckt->nodes[i].branch_count == 0) return -2;
    }

    /* Check that at least one branch exists */
    if (ckt->n_branches == 0) return -3;

    /* Check no duplicate branch IDs */
    for (int i = 0; i < ckt->n_branches; i++) {
        for (int j = i + 1; j < ckt->n_branches; j++) {
            if (ckt->branches[i].id == ckt->branches[j].id) return -4;
        }
    }

    return 0;
}

/*==================================================================
 * L3 ? Complex Impedance Computation Functions
 * Mathematical Structure: Z = R + jX, Y = G + jB
 *==================================================================*/

Impedance_t impedance_resistor(double R_ohms)
{
    Impedance_t Z;
    Z.re = R_ohms;
    Z.im = 0.0;
    return Z;
}

Impedance_t impedance_capacitor(double C_farads, double omega_rad_per_s)
{
    Impedance_t Z;
    if (C_farads <= 0.0 || omega_rad_per_s <= 0.0) {
        /* DC: open circuit; infinite frequency: short circuit */
        Z.re = (omega_rad_per_s <= 0.0) ? INFINITY : 0.0;
        Z.im = 0.0;
        return Z;
    }
    /* Z_C = 1/(j*omega*C) = 0 - j/(omega*C) */
    Z.re = 0.0;
    Z.im = -1.0 / (omega_rad_per_s * C_farads);
    return Z;
}

Impedance_t impedance_inductor(double L_henries, double omega_rad_per_s)
{
    Impedance_t Z;
    if (L_henries <= 0.0 || omega_rad_per_s <= 0.0) {
        /* DC: short circuit; infinite frequency: open circuit */
        Z.re = 0.0;
        Z.im = (omega_rad_per_s <= 0.0) ? 0.0 : INFINITY;
        return Z;
    }
    /* Z_L = j*omega*L */
    Z.re = 0.0;
    Z.im = omega_rad_per_s * L_henries;
    return Z;
}

Impedance_t impedance_rlc_series(double R, double L, double C, double omega)
{
    /* Z = R + j(omega*L - 1/(omega*C)) */
    Impedance_t Z;
    Z.re = R;

    double X_L = omega * L;
    double X_C_div = (omega * C);
    double X_C;

    if (X_C_div == 0.0) {
        X_C = -INFINITY;
    } else {
        X_C = -1.0 / X_C_div;
    }

    Z.im = X_L + X_C;

    /* Resonance condition: X_L + X_C = 0 => omega*L = 1/(omega*C)
     * => omega_0 = 1/sqrt(L*C) */
    return Z;
}

Impedance_t impedance_rlc_parallel(double R, double L, double C, double omega)
{
    /* Admittance: Y = 1/R + j(omega*C - 1/(omega*L)) */
    /* Impedance: Z = 1/Y */
    Impedance_t Z;

    double G = (R > 0.0) ? (1.0 / R) : INFINITY;
    double B_C = omega * C;
    double B_L_div = omega * L;

    double B_L;
    if (B_L_div == 0.0) {
        B_L = -INFINITY;
    } else {
        B_L = -1.0 / B_L_div;
    }

    double B = B_C + B_L;

    /* Y = G + jB, Z = 1/Y = (G - jB)/(G^2 + B^2) */
    double denom = G * G + B * B;
    if (denom == 0.0) {
        Z.re = INFINITY;
        Z.im = 0.0;
    } else {
        Z.re = G / denom;
        Z.im = -B / denom;
    }
    return Z;
}

Admittance_t impedance_to_admittance(Impedance_t Z)
{
    /* Y = 1/Z = 1/(R+jX) = (R-jX)/(R^2+X^2) */
    Admittance_t Y;
    double denom = Z.re * Z.re + Z.im * Z.im;

    if (denom == 0.0) {
        Y.re = INFINITY;
        Y.im = 0.0;
    } else {
        Y.re = Z.re / denom;
        Y.im = -Z.im / denom;
    }
    return Y;
}

Impedance_t admittance_to_impedance(Admittance_t Y)
{
    /* Z = 1/Y = 1/(G+jB) = (G-jB)/(G^2+B^2) */
    Impedance_t Z;
    double denom = Y.re * Y.re + Y.im * Y.im;

    if (denom == 0.0) {
        Z.re = INFINITY;
        Z.im = 0.0;
    } else {
        Z.re = Y.re / denom;
        Z.im = -Y.im / denom;
    }
    return Z;
}

Impedance_t impedance_series(Impedance_t Z1, Impedance_t Z2)
{
    /* Z_eq = Z1 + Z2 */
    Impedance_t Z;
    Z.re = Z1.re + Z2.re;
    Z.im = Z1.im + Z2.im;
    return Z;
}

Impedance_t impedance_parallel(Impedance_t Z1, Impedance_t Z2)
{
    /* Z_eq = (Z1 * Z2) / (Z1 + Z2) */
    /* Product: (R1+jX1)*(R2+jX2) = (R1*R2-X1*X2) + j(R1*X2+R2*X1) */
    double prod_re = Z1.re * Z2.re - Z1.im * Z2.im;
    double prod_im = Z1.re * Z2.im + Z2.re * Z1.im;

    /* Sum: (R1+R2) + j(X1+X2) */
    double sum_re = Z1.re + Z2.re;
    double sum_im = Z1.im + Z2.im;

    /* Division: (prod_re+j*prod_im)/(sum_re+j*sum_im) */
    double denom = sum_re * sum_re + sum_im * sum_im;

    Impedance_t Z;
    if (denom == 0.0) {
        Z.re = INFINITY;
        Z.im = 0.0;
    } else {
        Z.re = (prod_re * sum_re + prod_im * sum_im) / denom;
        Z.im = (prod_im * sum_re - prod_re * sum_im) / denom;
    }
    return Z;
}

Admittance_t admittance_parallel(Admittance_t Y1, Admittance_t Y2)
{
    /* Y_eq = Y1 + Y2 (admittances add in parallel) */
    Admittance_t Y;
    Y.re = Y1.re + Y2.re;
    Y.im = Y1.im + Y2.im;
    return Y;
}

Impedance_t phasor_to_impedance(const Phasor_t *p)
{
    /* A angle(phi) -> A*cos(phi) + j*A*sin(phi) */
    Impedance_t Z;
    if (!p) {
        Z.re = 0.0;
        Z.im = 0.0;
        return Z;
    }

    double phi_rad = p->phase_deg * M_PI / 180.0;
    Z.re = p->magnitude * cos(phi_rad);
    Z.im = p->magnitude * sin(phi_rad);
    return Z;
}

Phasor_t impedance_to_phasor(Impedance_t Z, double freq_hz)
{
    /* R+jX -> |Z| angle(atan2(X,R)) */
    Phasor_t p;
    p.magnitude = sqrt(Z.re * Z.re + Z.im * Z.im);
    p.phase_deg = atan2(Z.im, Z.re) * 180.0 / M_PI;
    p.freq_hz = freq_hz;
    return p;
}

/*==================================================================
 * L2/L4 ? Quality Factor, Resonance, Bandwidth, Damping
 *==================================================================*/

double q_factor_series(double R, double L, double C)
{
    /* Q = (1/R) * sqrt(L/C) = omega_0 * L / R
     * For R=0 (ideal LC): Q -> infinity */
    if (R <= 0.0) return INFINITY;
    if (C <= 0.0 || L <= 0.0) return 0.0;
    return sqrt(L / C) / R;
}

double q_factor_parallel(double R, double L, double C)
{
    /* Q = R * sqrt(C/L) = R/(omega_0 * L)
     * For R=infinity (ideal LC): Q -> infinity */
    if (R <= 0.0) return 0.0;
    if (L <= 0.0 || C <= 0.0) return 0.0;
    return R * sqrt(C / L);
}

double resonance_omega(double L, double C)
{
    /* omega_0 = 1/sqrt(L*C) (rad/s)
     * This is the undamped natural frequency.
     * Reference: Thomson (1853), formula for LC oscillation. */
    if (L <= 0.0 || C <= 0.0) return 0.0;
    return 1.0 / sqrt(L * C);
}

double resonance_freq_hz(double L, double C)
{
    /* f_0 = omega_0/(2*pi) = 1/(2*pi*sqrt(L*C)) (Hz) */
    return resonance_omega(L, C) / (2.0 * M_PI);
}

double bandwidth_series_radps(double R, double L)
{
    /* BW = R/L (rad/s) for series RLC.
     * BW = omega_0 / Q */
    if (L <= 0.0) return INFINITY;
    return R / L;
}

double bandwidth_parallel_radps(double R, double C)
{
    /* BW = 1/(R*C) (rad/s) for parallel RLC */
    if (R <= 0.0 || C <= 0.0) return INFINITY;
    return 1.0 / (R * C);
}

double damping_factor_series(double R, double L, double C)
{
    /* zeta = (R/2) * sqrt(C/L) = alpha/omega_0
     * alpha = R/(2L), omega_0 = 1/sqrt(LC) */
    if (L <= 0.0 || C <= 0.0) return 0.0;
    double alpha = R / (2.0 * L);
    double omega_0 = 1.0 / sqrt(L * C);
    if (omega_0 <= 0.0) return INFINITY;
    return alpha / omega_0;
}

double damping_factor_parallel(double R, double L, double C)
{
    /* zeta = (1/(2*R)) * sqrt(L/C) for parallel RLC */
    if (R <= 0.0 || L <= 0.0 || C <= 0.0) return 0.0;
    double alpha = 1.0 / (2.0 * R * C);
    double omega_0 = 1.0 / sqrt(L * C);
    if (omega_0 <= 0.0) return INFINITY;
    return alpha / omega_0;
}
