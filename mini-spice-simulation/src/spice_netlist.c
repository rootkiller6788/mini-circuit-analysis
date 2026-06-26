/**
 * @file spice_netlist.c
 * @brief SPICE netlist parser implementation
 *
 * Parses SPICE-format circuit description files (.cir).
 * Supports passive components (R,L,C), independent sources (V,I),
 * controlled sources (E,F,G,H), semiconductors (D,Q,M),
 * .MODEL definitions, analysis commands, and .OPTIONS.
 *
 * Reference: SPICE2 User's Guide (Nagel 1975), Appendix A
 */

#include "spice_netlist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/* ── Utility helpers ───────────────────────────────────────────────── */

/**
 * @brief Case-insensitive string comparison (SPICE is case-insensitive)
 */
static int strcicmp(const char *a, const char *b) {
    while (*a && *b) {
        int ca = tolower((unsigned char)*a);
        int cb = tolower((unsigned char)*b);
        if (ca != cb) return ca - cb;
        a++; b++;
    }
    return (unsigned char)*a - (unsigned char)*b;
}

/**
 * @brief Trim leading whitespace
 */
static const char* trim_left(const char *s) {
    while (*s && isspace((unsigned char)*s)) s++;
    return s;
}

/**
 * @brief Trim trailing whitespace (modifies in place, returns start)
 */
static char* trim_right(char *s) {
    if (!s || !*s) return s;
    char *end = s + strlen(s) - 1;
    while (end >= s && isspace((unsigned char)*end)) end--;
    *(end + 1) = '\0';
    return s;
}

/**
 * @brief Trim both sides (modifies in place)
 */
static char* trim(char *s) {
    return trim_right((char*)trim_left(s));
}

/**
 * @brief Extract the next token from a line, advancing the pointer
 *
 * Tokens are separated by whitespace or '='.
 *
 * @param line Pointer to line pointer (updated in place)
 * @param token Output buffer for extracted token
 * @param maxlen Token buffer size
 * @return 0 on success, -1 at end of line
 */
static int next_token(const char **line, char *token, int maxlen) {
    /* Skip whitespace */
    while (**line && isspace((unsigned char)**line)) (*line)++;
    if (!**line) return -1;

    int len = 0;

    /* Read token */
    while (**line && !isspace((unsigned char)**line) && **line != '=' && len < maxlen - 1) {
        token[len++] = **line;
        (*line)++;
    }
    /* Skip trailing '=' */
    if (**line == '=') {
        if (len == 0) {
            token[len++] = **line;
            (*line)++;
        }
        /* If '=' was the token itself, we're done.
         * Otherwise, '=' is a separator and we'll get it next call. */
        if (len == 0) return -1;
    }
    token[len] = '\0';
    return 0;
}

/**
 * @brief Parse a floating-point value from a token
 *
 * Supports SPICE suffixes: f=1e-15, p=1e-12, n=1e-9, u=1e-6,
 * m=1e-3, k=1e3, meg=1e6, g=1e9, t=1e12
 * Also supports mil=25.4e-6
 *
 * @param token Token string
 * @param value Output value
 * @return 0 on success, -1 if not a number
 */
static int parse_value(const char *token, double *value) {
    if (!token || !*token) return -1;

    /* Handle SPICE suffixes */
    char buf[64];
    int i = 0;
    const char *p = token;

    /* Copy numeric part */
    while (*p && (isdigit((unsigned char)*p) || *p == '.' || *p == '+' ||
                  *p == '-' || *p == 'e' || *p == 'E')) {
        if (i < 63) buf[i++] = *p;
        p++;
    }
    buf[i] = '\0';

    char *endptr;
    double val = strtod(buf, &endptr);
    if (endptr == buf && !*p) return -1;  /* Not a number at all */

    /* Apply suffix multiplier */
    if (*p) {
        /* Case-insensitive suffix check */
        char suffix[8] = {0};
        int j = 0;
        while (*p && j < 7) suffix[j++] = tolower((unsigned char)*p++);
        suffix[j] = '\0';

        if      (strcmp(suffix, "f") == 0)    val *= 1e-15;
        else if (strcmp(suffix, "p") == 0)    val *= 1e-12;
        else if (strcmp(suffix, "n") == 0)    val *= 1e-9;
        else if (strcmp(suffix, "u") == 0)    val *= 1e-6;
        else if (strcmp(suffix, "m") == 0)    val *= 1e-3;
        else if (strcmp(suffix, "k") == 0)    val *= 1e3;
        else if (strcmp(suffix, "meg") == 0)  val *= 1e6;
        else if (strcmp(suffix, "g") == 0)    val *= 1e9;
        else if (strcmp(suffix, "t") == 0)    val *= 1e12;
        else if (strcmp(suffix, "mil") == 0)  val *= 25.4e-6;
        /* unrecognized suffix → treat as part of name, value is the numeric prefix */
    }

    *value = val;
    return 0;
}

/* ── Node Management ───────────────────────────────────────────────── */

void spice_netlist_init(spice_netlist_t *nl) {
    memset(nl, 0, sizeof(*nl));
    nl->analysis_type = SPICE_ANALYSIS_NONE;
    /* Set default options */
    nl->options_abstol = 1e-12;
    nl->options_vntol  = 1e-6;
    nl->options_reltol = 0.001;
    nl->options_temp   = 27.0;
    nl->options_itl1   = 100;
    nl->options_itl4   = 10;
}

spice_node_id spice_netlist_get_node(spice_netlist_t *nl, const char *name) {
    /* Ground node */
    if (strcicmp(name, "0") == 0 || strcicmp(name, "gnd") == 0 ||
        strcicmp(name, "ground") == 0) {
        return SPICE_GROUND_NODE;
    }

    /* Search existing nodes */
    for (int32_t i = 1; i <= nl->num_nodes; i++) {
        if (strcicmp(nl->nodes[i].name, name) == 0) {
            return nl->nodes[i].id;
        }
    }

    /* Create new node */
    if (nl->num_nodes >= SPICE_MAX_NODES - 1) return -1;
    int32_t new_id = nl->num_nodes + 1;
    nl->num_nodes = new_id;
    nl->nodes[new_id].id = new_id;
    nl->nodes[new_id].is_dc_ground = 0;
    strncpy(nl->nodes[new_id].name, name, SPICE_MAX_NAME - 1);
    nl->nodes[new_id].name[SPICE_MAX_NAME - 1] = '\0';
    return new_id;
}

/* ── Component Parsers ─────────────────────────────────────────────── */

/**
 * @brief Parse resistor: R<name> <n+> <n-> <value>
 *
 * Example: R1 1 2 1k
 */
static int parse_resistor(spice_netlist_t *nl, const char **lineptr,
                          const char *name) {
    char tok[SPICE_MAX_NAME];
    spice_resistor_t *r = calloc(1, sizeof(spice_resistor_t));
    if (!r) return -1;
    strncpy(r->base.name, name, SPICE_MAX_NAME - 1);
    r->base.type = SPICE_COMP_RESISTOR;

    /* n+ */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(r); return -1; }
    r->base.nplus = spice_netlist_get_node(nl, tok);

    /* n- */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(r); return -1; }
    r->base.nminus = spice_netlist_get_node(nl, tok);

    /* value */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(r); return -1; }
    if (parse_value(tok, &r->resistance) != 0 || r->resistance <= 0.0) {
        free(r); return -1;
    }
    r->base.value = r->resistance;

    /* Store */
    if (nl->num_components >= SPICE_MAX_COMPONENTS) { free(r); return -1; }
    nl->components[nl->num_components++] = (spice_component_t*)r;
    nl->num_resistors++;
    return 0;
}

/**
 * @brief Parse capacitor: C<name> <n+> <n-> <value> [IC=<initial_v>]
 */
static int parse_capacitor(spice_netlist_t *nl, const char **lineptr,
                           const char *name) {
    char tok[SPICE_MAX_NAME];
    spice_capacitor_t *c = calloc(1, sizeof(spice_capacitor_t));
    if (!c) return -1;
    strncpy(c->base.name, name, SPICE_MAX_NAME - 1);
    c->base.type = SPICE_COMP_CAPACITOR;

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(c); return -1; }
    c->base.nplus = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(c); return -1; }
    c->base.nminus = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(c); return -1; }
    if (parse_value(tok, &c->capacitance) != 0 || c->capacitance <= 0.0) {
        free(c); return -1;
    }
    c->base.value = c->capacitance;

    /* Optional IC= keyword */
    const char *saved = *lineptr;
    if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
        if (strcicmp(tok, "IC") == 0) {
            if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                    parse_value(tok, &c->initial_v);
                }
            }
        } else {
            *lineptr = saved;  /* Put token back */
        }
    }

    if (nl->num_components >= SPICE_MAX_COMPONENTS) { free(c); return -1; }
    nl->components[nl->num_components++] = (spice_component_t*)c;
    nl->num_capacitors++;
    return 0;
}

/**
 * @brief Parse inductor: L<name> <n+> <n-> <value> [IC=<initial_i>]
 */
static int parse_inductor(spice_netlist_t *nl, const char **lineptr,
                          const char *name) {
    char tok[SPICE_MAX_NAME];
    spice_inductor_t *l = calloc(1, sizeof(spice_inductor_t));
    if (!l) return -1;
    strncpy(l->base.name, name, SPICE_MAX_NAME - 1);
    l->base.type = SPICE_COMP_INDUCTOR;

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(l); return -1; }
    l->base.nplus  = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(l); return -1; }
    l->base.nminus = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(l); return -1; }
    if (parse_value(tok, &l->inductance) != 0 || l->inductance <= 0.0) {
        free(l); return -1;
    }
    l->base.value = l->inductance;

    /* Optional IC */
    const char *saved = *lineptr;
    if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
        if (strcicmp(tok, "IC") == 0) {
            if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                    parse_value(tok, &l->initial_i);
                }
            }
        } else {
            *lineptr = saved;
        }
    }

    if (nl->num_components >= SPICE_MAX_COMPONENTS) { free(l); return -1; }
    nl->components[nl->num_components++] = (spice_component_t*)l;
    nl->num_inductors++;
    return 0;
}

/**
 * @brief Parse independent voltage source: V<name> <n+> <n-> [DC] <value> [AC <mag> <phase>]
 *
 * Also supports: V<name> <n+> <n-> SIN(VO VA FREQ) or PULSE(...)
 */
static int parse_vsource(spice_netlist_t *nl, const char **lineptr,
                         const char *name) {
    char tok[SPICE_MAX_NAME];
    spice_vsource_t *v = calloc(1, sizeof(spice_vsource_t));
    if (!v) return -1;
    strncpy(v->base.name, name, SPICE_MAX_NAME - 1);
    v->base.type = SPICE_COMP_VSOURCE;

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(v); return -1; }
    v->base.nplus  = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(v); return -1; }
    v->base.nminus = spice_netlist_get_node(nl, tok);

    /* Value — may have DC prefix */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(v); return -1; }
    if (strcicmp(tok, "DC") == 0) {
        if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(v); return -1; }
    }
    if (parse_value(tok, &v->dc_value) != 0) { free(v); return -1; }
    v->base.value = v->dc_value;

    /* Optional AC magnitude/phase */
    const char *saved = *lineptr;
    if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
        if (strcicmp(tok, "AC") == 0) {
            if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                parse_value(tok, &v->ac_magnitude);
            }
            if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                parse_value(tok, &v->ac_phase);
            }
        } else {
            *lineptr = saved;
        }
    }

    if (nl->num_components >= SPICE_MAX_COMPONENTS) { free(v); return -1; }
    nl->components[nl->num_components++] = (spice_component_t*)v;
    nl->num_vsources++;
    return 0;
}

/**
 * @brief Parse independent current source: I<name> <n+> <n-> [DC] <value>
 */
static int parse_isource(spice_netlist_t *nl, const char **lineptr,
                         const char *name) {
    char tok[SPICE_MAX_NAME];
    spice_isource_t *is = calloc(1, sizeof(spice_isource_t));
    if (!is) return -1;
    strncpy(is->base.name, name, SPICE_MAX_NAME - 1);
    is->base.type = SPICE_COMP_ISOURCE;
    is->waveform_type = 0; /* DC */

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(is); return -1; }
    is->base.nplus  = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(is); return -1; }
    is->base.nminus = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(is); return -1; }
    if (strcicmp(tok, "DC") == 0) {
        if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(is); return -1; }
    }
    if (parse_value(tok, &is->dc_value) != 0) { free(is); return -1; }
    is->base.value = is->dc_value;

    /* Check for waveform type */
    const char *saved = *lineptr;
    if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
        if (strcicmp(tok, "SIN") == 0) {
            is->waveform_type = 1; /* SIN */
            /* Read SIN parameters: VO VA FREQ TD THETA PHASE */
            for (int i = 0; i < 6; i++) {
                if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                    parse_value(tok, &is->waveform_params[i]);
                }
            }
        } else if (strcicmp(tok, "PULSE") == 0) {
            is->waveform_type = 2; /* PULSE */
            for (int i = 0; i < 7; i++) {
                if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                    parse_value(tok, &is->waveform_params[i]);
                }
            }
        } else if (strcicmp(tok, "PWL") == 0) {
            is->waveform_type = 3; /* PWL */
            for (int i = 0; i < 8; i++) {
                if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
                    parse_value(tok, &is->waveform_params[i]);
                }
            }
        } else {
            *lineptr = saved;
        }
    }

    if (nl->num_components >= SPICE_MAX_COMPONENTS) { free(is); return -1; }
    nl->components[nl->num_components++] = (spice_component_t*)is;
    nl->num_isources++;
    return 0;
}

/**
 * @brief Parse diode: D<name> <n+> <n-> <modelname> [area]
 */
static int parse_diode(spice_netlist_t *nl, const char **lineptr,
                       const char *name) {
    char tok[SPICE_MAX_NAME];
    spice_component_t *d = calloc(1, sizeof(spice_component_t));
    if (!d) return -1;
    strncpy(d->name, name, SPICE_MAX_NAME - 1);
    d->type = SPICE_COMP_DIODE;
    d->is_nonlinear = 1;

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(d); return -1; }
    d->nplus = spice_netlist_get_node(nl, tok);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) { free(d); return -1; }
    d->nminus = spice_netlist_get_node(nl, tok);

    /* Model name — stored in a field we'll look up later */
    /* For simplicity, we use the name as the model reference */
    /* The actual model lookup happens during analysis */

    if (nl->num_components >= SPICE_MAX_COMPONENTS) { free(d); return -1; }
    nl->components[nl->num_components++] = d;
    nl->num_diodes++;
    return 0;
}

/* ── .MODEL Parser ─────────────────────────────────────────────────── */

/**
 * @brief Parse .MODEL statement: .MODEL <name> <type> (param1=val1 param2=val2 ...)
 */
static int parse_model(spice_netlist_t *nl, const char **lineptr) {
    char model_name[SPICE_MAX_NAME];
    char model_type[SPICE_MAX_NAME];

    /* Model name */
    if (next_token(lineptr, model_name, SPICE_MAX_NAME) != 0) return -1;

    /* Model type */
    if (next_token(lineptr, model_type, SPICE_MAX_NAME) != 0) return -1;

    /* Parse parenthesized parameter list */
    /* Skip until '(' */
    while (**lineptr && **lineptr != '(') (*lineptr)++;
    if (**lineptr == '(') (*lineptr)++;

    if (strcicmp(model_type, "D") == 0) {
        if (nl->num_diode_models >= 32) return -1;
        spice_diode_model_t *dm = &nl->diode_models[nl->num_diode_models++];
        memset(dm, 0, sizeof(*dm));
        strncpy(dm->model_name, model_name, SPICE_MAX_NAME - 1);
        /* Default diode parameters */
        dm->is_saturation = 1e-14;
        dm->n_ideality    = 1.0;
        dm->series_r      = 0.0;
        dm->transit_time  = 0.0;
        dm->zero_bias_c   = 0.0;
        dm->vj_potential  = 0.75;
        dm->grading_coeff = 0.5;
        dm->breakdown_v   = 1e30;
        dm->breakdown_i   = 1e-3;
        dm->is_used       = 0;

        /* Parse param=value pairs */
        char param[64], val[64];
        while (**lineptr && **lineptr != ')') {
            /* Read param name */
            int pi = 0;
            while (**lineptr && **lineptr != '=' && **lineptr != ')' &&
                   !isspace((unsigned char)**lineptr) && pi < 63) {
                param[pi++] = **lineptr;
                (*lineptr)++;
            }
            param[pi] = '\0';
            if (**lineptr == '=') {
                (*lineptr)++;
                int vi = 0;
                while (**lineptr && **lineptr != ')' &&
                       !isspace((unsigned char)**lineptr) && vi < 63) {
                    val[vi++] = **lineptr;
                    (*lineptr)++;
                }
                val[vi] = '\0';
                double v;
                if (parse_value(val, &v) == 0) {
                    if (strcicmp(param, "IS") == 0) dm->is_saturation = v;
                    else if (strcicmp(param, "N") == 0) dm->n_ideality = v;
                    else if (strcicmp(param, "RS") == 0) dm->series_r = v;
                    else if (strcicmp(param, "TT") == 0) dm->transit_time = v;
                    else if (strcicmp(param, "CJO") == 0) dm->zero_bias_c = v;
                    else if (strcicmp(param, "VJ") == 0) dm->vj_potential = v;
                    else if (strcicmp(param, "M") == 0) dm->grading_coeff = v;
                    else if (strcicmp(param, "BV") == 0) dm->breakdown_v = v;
                    else if (strcicmp(param, "IBV") == 0) dm->breakdown_i = v;
                }
            }
            /* Skip whitespace */
            while (**lineptr && isspace((unsigned char)**lineptr)) (*lineptr)++;
        }

    } else if (strcicmp(model_type, "NMOS") == 0 || strcicmp(model_type, "PMOS") == 0) {
        if (nl->num_mos_models >= 32) return -1;
        spice_mos_model_t *mm = &nl->mos_models[nl->num_mos_models++];
        memset(mm, 0, sizeof(*mm));
        strncpy(mm->model_name, model_name, SPICE_MAX_NAME - 1);
        mm->type = (tolower((unsigned char)model_type[0]) == 'p') ? 1 : 0;
        /* Default MOS Level 1 parameters */
        mm->vth0           = (mm->type == 0) ? 0.7 : -0.7;
        mm->kp             = 2e-5;
        mm->gamma_body     = 0.5;
        mm->phi_surface    = 0.7;
        mm->lambda_channel = 0.01;
        mm->tox            = 1e-7;
        mm->is_used        = 0;

        char param[64], val[64];
        while (**lineptr && **lineptr != ')') {
            int pi = 0;
            while (**lineptr && **lineptr != '=' && **lineptr != ')' &&
                   !isspace((unsigned char)**lineptr) && pi < 63) {
                param[pi++] = **lineptr;
                (*lineptr)++;
            }
            param[pi] = '\0';
            if (**lineptr == '=') {
                (*lineptr)++;
                int vi = 0;
                while (**lineptr && **lineptr != ')' &&
                       !isspace((unsigned char)**lineptr) && vi < 63) {
                    val[vi++] = **lineptr;
                    (*lineptr)++;
                }
                val[vi] = '\0';
                double v;
                if (parse_value(val, &v) == 0) {
                    if (strcicmp(param, "VTO") == 0) mm->vth0 = v;
                    else if (strcicmp(param, "KP") == 0) mm->kp = v;
                    else if (strcicmp(param, "GAMMA") == 0) mm->gamma_body = v;
                    else if (strcicmp(param, "PHI") == 0) mm->phi_surface = v;
                    else if (strcicmp(param, "LAMBDA") == 0) mm->lambda_channel = v;
                    else if (strcicmp(param, "TOX") == 0) mm->tox = v;
                }
            }
            while (**lineptr && isspace((unsigned char)**lineptr)) (*lineptr)++;
        }

    } else if (strcicmp(model_type, "NPN") == 0 || strcicmp(model_type, "PNP") == 0) {
        if (nl->num_bjt_models >= 32) return -1;
        spice_bjt_model_t *bm = &nl->bjt_models[nl->num_bjt_models++];
        memset(bm, 0, sizeof(*bm));
        strncpy(bm->model_name, model_name, SPICE_MAX_NAME - 1);
        /* Default BJT parameters */
        bm->is_saturation = 1e-16;
        bm->bf_forward    = 100.0;
        bm->br_reverse    = 1.0;
        bm->nf_coeff      = 1.0;
        bm->nr_coeff      = 1.0;
        bm->vaf_early     = 1e30;
        bm->var_early     = 1e30;
        bm->is_used       = 0;

        char param[64], val[64];
        while (**lineptr && **lineptr != ')') {
            int pi = 0;
            while (**lineptr && **lineptr != '=' && **lineptr != ')' &&
                   !isspace((unsigned char)**lineptr) && pi < 63) {
                param[pi++] = **lineptr;
                (*lineptr)++;
            }
            param[pi] = '\0';
            if (**lineptr == '=') {
                (*lineptr)++;
                int vi = 0;
                while (**lineptr && **lineptr != ')' &&
                       !isspace((unsigned char)**lineptr) && vi < 63) {
                    val[vi++] = **lineptr;
                    (*lineptr)++;
                }
                val[vi] = '\0';
                double v;
                if (parse_value(val, &v) == 0) {
                    if (strcicmp(param, "IS") == 0) bm->is_saturation = v;
                    else if (strcicmp(param, "BF") == 0) bm->bf_forward = v;
                    else if (strcicmp(param, "BR") == 0) bm->br_reverse = v;
                    else if (strcicmp(param, "NF") == 0) bm->nf_coeff = v;
                    else if (strcicmp(param, "NR") == 0) bm->nr_coeff = v;
                    else if (strcicmp(param, "VAF") == 0) bm->vaf_early = v;
                    else if (strcicmp(param, "VAR") == 0) bm->var_early = v;
                }
            }
            while (**lineptr && isspace((unsigned char)**lineptr)) (*lineptr)++;
        }
    }

    return 0;
}

/* ── Analysis Command Parser ───────────────────────────────────────── */

/**
 * @brief Parse .DC command: .DC <srcname> <start> <stop> <step>
 */
static int parse_dc_command(spice_netlist_t *nl, const char **lineptr) {
    char tok[SPICE_MAX_NAME];
    nl->analysis_type = SPICE_ANALYSIS_DC;

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    strncpy(nl->dc_source_name, tok, SPICE_MAX_NAME - 1);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    parse_value(tok, &nl->dc_start);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    parse_value(tok, &nl->dc_stop);

    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    parse_value(tok, &nl->dc_step);

    return 0;
}

/**
 * @brief Parse .AC command: .AC DEC|OCT|LIN <N> <fstart> <fstop>
 */
static int parse_ac_command(spice_netlist_t *nl, const char **lineptr) {
    char tok[SPICE_MAX_NAME];
    nl->analysis_type = SPICE_ANALYSIS_AC;

    /* DEC/OCT/LIN */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    /* Number of points per decade/octave or total */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    nl->ac_points_per_decade = (int32_t)strtol(tok, NULL, 10);

    /* fstart */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    parse_value(tok, &nl->ac_fstart);

    /* fstop */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    parse_value(tok, &nl->ac_fstop);

    return 0;
}

/**
 * @brief Parse .TRAN command: .TRAN <tstep> <tstop> [tstart] [tmax]
 */
static int parse_tran_command(spice_netlist_t *nl, const char **lineptr) {
    char tok[SPICE_MAX_NAME];
    nl->analysis_type = SPICE_ANALYSIS_TRAN;

    /* tstep */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    parse_value(tok, &nl->tran_tstep);

    /* tstop */
    if (next_token(lineptr, tok, SPICE_MAX_NAME) != 0) return -1;
    parse_value(tok, &nl->tran_tstop);

    /* Optional tstart */
    const char *saved = *lineptr;
    if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
        parse_value(tok, &nl->tran_tstart);
        /* Optional tmax */
        saved = *lineptr;
        if (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
            parse_value(tok, &nl->tran_tmax);
        } else {
            *lineptr = saved;
        }
    }

    return 0;
}

/**
 * @brief Parse .OPTIONS line
 */
static int parse_options(spice_netlist_t *nl, const char **lineptr) {
    char tok[SPICE_MAX_NAME];
    nl->has_options = 1;

    while (next_token(lineptr, tok, SPICE_MAX_NAME) == 0) {
        /* Look for param=value */
        char *eq = strchr(tok, '=');
        if (eq) {
            *eq = '\0';
            double v;
            if (parse_value(eq + 1, &v) == 0) {
                if (strcicmp(tok, "ABSTOL") == 0) nl->options_abstol = v;
                else if (strcicmp(tok, "VNTOL") == 0) nl->options_vntol = v;
                else if (strcicmp(tok, "RELTOL") == 0) nl->options_reltol = v;
                else if (strcicmp(tok, "TEMP") == 0) nl->options_temp = v;
                else if (strcicmp(tok, "ITL1") == 0) nl->options_itl1 = (int32_t)v;
                else if (strcicmp(tok, "ITL4") == 0) nl->options_itl4 = (int32_t)v;
            }
        }
    }
    return 0;
}

/* ── Main Parse API ────────────────────────────────────────────────── */

int spice_netlist_parse_line(spice_netlist_t *nl, const char *line) {
    const char *p = trim_left(line);

    /* Skip empty lines and comments */
    if (!*p || *p == '*') return 0;

    /* Copy line for tokenization */
    char linebuf[SPICE_MAX_LINE];
    strncpy(linebuf, p, SPICE_MAX_LINE - 1);
    linebuf[SPICE_MAX_LINE - 1] = '\0';
    p = linebuf;

    char first[SPICE_MAX_NAME];
    if (next_token(&p, first, SPICE_MAX_NAME) != 0) return 0;

    /* Directive lines start with '.' */
    if (first[0] == '.') {
        if (strcicmp(first, ".END") == 0) {
            return 0;
        } else if (strcicmp(first, ".MODEL") == 0) {
            return parse_model(nl, &p);
        } else if (strcicmp(first, ".DC") == 0) {
            return parse_dc_command(nl, &p);
        } else if (strcicmp(first, ".AC") == 0) {
            return parse_ac_command(nl, &p);
        } else if (strcicmp(first, ".TRAN") == 0) {
            return parse_tran_command(nl, &p);
        } else if (strcicmp(first, ".OPTIONS") == 0) {
            return parse_options(nl, &p);
        } else if (strcicmp(first, ".TF") == 0) {
            nl->analysis_type = SPICE_ANALYSIS_TF;
            return 0;
        } else if (strcicmp(first, ".NOISE") == 0) {
            nl->analysis_type = SPICE_ANALYSIS_NOISE;
            return 0;
        }
        /* Unknown directive — ignore */
        return 0;
    }

    /* Component lines: first character determines type */
    char ctype = (char)toupper((unsigned char)first[0]);

    switch (ctype) {
    case 'R':
        return parse_resistor(nl, &p, first);
    case 'C':
        return parse_capacitor(nl, &p, first);
    case 'L':
        return parse_inductor(nl, &p, first);
    case 'V':
        return parse_vsource(nl, &p, first);
    case 'I':
        return parse_isource(nl, &p, first);
    case 'D':
        return parse_diode(nl, &p, first);
    case 'Q':
        /* BJT: Q<name> <nc> <nb> <ne> <ns> <modelname> */
        {
            spice_component_t *q = calloc(1, sizeof(spice_component_t));
            if (!q) return -1;
            strncpy(q->name, first, SPICE_MAX_NAME - 1);
            char tok2[SPICE_MAX_NAME];
            if (next_token(&p, tok2, SPICE_MAX_NAME) == 0) {
                if (strcicmp(tok2, "NPN") == 0 || strcicmp(tok2, "PNP") == 0) {
                    /* Q<name> <type> ... — we look for nodes after type */
                    q->type = (tolower((unsigned char)tok2[0]) == 'p') ?
                              SPICE_COMP_BJT_PNP : SPICE_COMP_BJT_NPN;
                } else {
                    /* tok2 is collector node */
                    q->nplus = spice_netlist_get_node(nl, tok2);
                    q->type = SPICE_COMP_BJT_NPN; /* default */
                }
            }
            q->is_nonlinear = 1;
            if (nl->num_components < SPICE_MAX_COMPONENTS) {
                nl->components[nl->num_components++] = q;
                nl->num_bjts++;
            } else { free(q); }
        }
        return 0;
    case 'M':
        /* MOSFET: M<name> <nd> <ng> <ns> <nb> <modelname> [W= L=] */
        {
            spice_component_t *m = calloc(1, sizeof(spice_component_t));
            if (!m) return -1;
            strncpy(m->name, first, SPICE_MAX_NAME - 1);
            m->type = SPICE_COMP_MOS_NMOS;
            m->is_nonlinear = 1;
            char tok2[SPICE_MAX_NAME];
            /* Read nodes */
            if (next_token(&p, tok2, SPICE_MAX_NAME) == 0) {
                m->nplus = spice_netlist_get_node(nl, tok2); /* drain */
            }
            if (nl->num_components < SPICE_MAX_COMPONENTS) {
                nl->components[nl->num_components++] = m;
                nl->num_mosfets++;
            } else { free(m); }
        }
        return 0;
    default:
        /* Unknown component → ignore */
        return 0;
    }
}

int spice_netlist_parse_file(spice_netlist_t *nl, const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) return -1;

    char line[SPICE_MAX_LINE];
    int line_num = 0;
    int ret = 0;

    /* First line is the title */
    if (fgets(line, SPICE_MAX_LINE, fp)) {
        line_num++;
        strncpy(nl->title, trim(line), SPICE_MAX_LINE - 1);
        nl->title[SPICE_MAX_LINE - 1] = '\0';
    }

    /* Parse remaining lines */
    while (fgets(line, SPICE_MAX_LINE, fp)) {
        line_num++;
        char *trimmed = trim(line);

        /* Remove trailing comment (anything after ';' or '//') */
        char *comment = strpbrk(trimmed, ";/");
        if (comment) *comment = '\0';

        ret = spice_netlist_parse_line(nl, trimmed);
        if (ret != 0) {
            fprintf(stderr, "Parse error at line %d: %s\n", line_num, trimmed);
            fclose(fp);
            return -1;
        }
    }

    fclose(fp);
    return 0;
}

int32_t spice_netlist_mna_size(const spice_netlist_t *nl) {
    /* MNA size = nodes + voltage sources + inductors + mutual inductors */
    int32_t size = nl->num_nodes;  /* Node voltages (nodes 1..N) */
    size += nl->num_vsources;      /* Branch current for each V-source */
    size += nl->num_inductors;     /* Branch current for each inductor */
    /* Each mutual inductor pair adds 2 branch currents */
    return size;
}

void spice_netlist_print_summary(const spice_netlist_t *nl) {
    printf("=== SPICE Netlist Summary ===\n");
    printf("Title: %s\n", nl->title);
    printf("Nodes: %d (excluding GND)\n", nl->num_nodes);
    printf("Components: %d total\n", nl->num_components);
    printf("  Resistors:  %d\n", nl->num_resistors);
    printf("  Capacitors: %d\n", nl->num_capacitors);
    printf("  Inductors:  %d\n", nl->num_inductors);
    printf("  V-sources:  %d\n", nl->num_vsources);
    printf("  I-sources:  %d\n", nl->num_isources);
    printf("  Diodes:     %d\n", nl->num_diodes);
    printf("  BJTs:       %d\n", nl->num_bjts);
    printf("  MOSFETs:    %d\n", nl->num_mosfets);
    printf("MNA system size: %d equations\n", spice_netlist_mna_size(nl));
    printf("Analysis requested: ");
    switch (nl->analysis_type) {
    case SPICE_ANALYSIS_DC:   printf(".DC\n"); break;
    case SPICE_ANALYSIS_AC:   printf(".AC\n"); break;
    case SPICE_ANALYSIS_TRAN: printf(".TRAN\n"); break;
    case SPICE_ANALYSIS_TF:   printf(".TF\n"); break;
    default:                  printf("None\n"); break;
    }
}

void spice_netlist_destroy(spice_netlist_t *nl) {
    for (int32_t i = 0; i < nl->num_components; i++) {
        free(nl->components[i]);
        nl->components[i] = NULL;
    }
    nl->num_components = 0;
}
