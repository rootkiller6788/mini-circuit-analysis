/**
 * @file spice_netlist.h
 * @brief SPICE netlist parser — data structures and API
 *
 * Knowledge coverage:
 *   L1 (Definitions): Component types enum, Node/Component/Netlist structs
 *   L2 (Core Concepts): Circuit topology representation, node numbering
 *   Reference: SPICE2 (Nagel 1975), "Inside SPICE" (Kielkowski 1998)
 *
 * Course alignment:
 *   Berkeley EE105 (Analog IC): netlist-based circuit description
 *   MIT 6.630 (EM Waves): lumped element circuit topology
 *   Stanford EE247 (Optical): circuit modeling interface
 */

#ifndef SPICE_NETLIST_H
#define SPICE_NETLIST_H

#include <stddef.h>
#include <stdint.h>

/* ── L1: Core Definitions ─────────────────────────────────────────── */

/** Maximum number of nodes in a circuit */
#define SPICE_MAX_NODES      1024
/** Maximum number of components */
#define SPICE_MAX_COMPONENTS 2048
/** Maximum identifier string length */
#define SPICE_MAX_NAME       64
/** Maximum characters per netlist line */
#define SPICE_MAX_LINE       256
/** Ground (reference) node number — always node 0 */
#define SPICE_GROUND_NODE    0

/**
 * @brief Component type enumeration
 *
 * Each type corresponds to a unique two-terminal or multi-terminal
 * device with its own stamp (contribution to the MNA matrix).
 */
typedef enum {
    SPICE_COMP_RESISTOR     = 0,  /**< Linear resistor        R */
    SPICE_COMP_CAPACITOR    = 1,  /**< Linear capacitor       C */
    SPICE_COMP_INDUCTOR     = 2,  /**< Linear inductor        L */
    SPICE_COMP_VSOURCE      = 3,  /**< Independent voltage source  V */
    SPICE_COMP_ISOURCE      = 4,  /**< Independent current source  I */
    SPICE_COMP_VCVS         = 5,  /**< Voltage-controlled voltage source E */
    SPICE_COMP_VCCS         = 6,  /**< Voltage-controlled current source G */
    SPICE_COMP_CCCS         = 7,  /**< Current-controlled current source F */
    SPICE_COMP_CCVS         = 8,  /**< Current-controlled voltage source H */
    SPICE_COMP_DIODE        = 9,  /**< PN junction diode       D */
    SPICE_COMP_BJT_NPN      = 10, /**< NPN bipolar transistor  Q */
    SPICE_COMP_BJT_PNP      = 11, /**< PNP bipolar transistor  Q */
    SPICE_COMP_MOS_NMOS     = 12, /**< N-channel MOSFET        M */
    SPICE_COMP_MOS_PMOS     = 13, /**< P-channel MOSFET        M */
    SPICE_COMP_MUTUAL_IND   = 14, /**< Mutual inductance pair  K */
    SPICE_COMP_UNKNOWN      = 99  /**< Unrecognized component         */
} spice_component_type_t;

/**
 * @brief Analysis type requested in the netlist (.DC, .AC, .TRAN)
 */
typedef enum {
    SPICE_ANALYSIS_NONE   = 0,  /**< No analysis requested          */
    SPICE_ANALYSIS_DC     = 1,  /**< DC sweep / operating point      */
    SPICE_ANALYSIS_AC     = 2,  /**< AC small-signal frequency sweep */
    SPICE_ANALYSIS_TRAN   = 3,  /**< Transient time-domain           */
    SPICE_ANALYSIS_TF     = 4,  /**< Transfer function (DC small-sig) */
    SPICE_ANALYSIS_NOISE  = 5   /**< Noise analysis                  */
} spice_analysis_cmd_t;

/** @brief Single node identifier (index into the node list) */
typedef int32_t spice_node_id;

/**
 * @brief Structure representing one named circuit node
 *
 * Node 0 is always ground (GND). Named nodes (e.g., "n1", "vout")
 * are mapped to positive integer IDs during parsing.
 */
typedef struct {
    char      name[SPICE_MAX_NAME]; /**< Node name string (e.g. "vout") */
    int32_t   id;                   /**< Unique integer ID (0 = GND)    */
    int32_t   is_dc_ground;         /**< 1 if this node is at DC ground */
} spice_node_t;

/**
 * @brief Generic two-terminal component
 *
 * All components inherit these fields. Specific parameters
 * (resistance, capacitance, model name) are stored in the
 * union within each component type.
 */
typedef struct {
    char                  name[SPICE_MAX_NAME]; /**< Component label (e.g. "R1")    */
    spice_component_type_t type;                /**< Component type enum             */
    spice_node_id          nplus;               /**< Positive terminal node ID       */
    spice_node_id          nminus;              /**< Negative terminal node ID       */
    double                 value;               /**< Primary parameter (R, C, L, V, I) */
    int32_t                is_nonlinear;        /**< 1 if component requires Newton iteration */
} spice_component_t;

/**
 * @brief Resistor: V = I * R
 *
 * The resistor stamp contributes to the MNA conductance matrix G.
 * For a resistor between nodes n+ and n- with resistance R:
 *   G[n+][n+] += 1/R,  G[n+][n-] -= 1/R
 *   G[n-][n+] -= 1/R,  G[n-][n-] += 1/R
 */
typedef struct {
    spice_component_t base;
    double            resistance;    /**< Resistance in ohms                */
    double            temp_coeff;    /**< Temperature coefficient (1/degC)  */
} spice_resistor_t;

/**
 * @brief Capacitor: I = C * dV/dt
 *
 * In transient analysis, the capacitor is discretized via companion models.
 * Trapezoidal rule: I_{n+1} = (2C/h)*(V_{n+1} - V_n) - I_n
 * This gives a conductance 2C/h in parallel with a current source.
 */
typedef struct {
    spice_component_t base;
    double            capacitance;   /**< Capacitance in farads             */
    double            initial_v;     /**< Initial voltage across capacitor  */
    double            prev_voltage;  /**< Voltage at previous time step     */
    double            prev_current;  /**< Current at previous time step     */
} spice_capacitor_t;

/**
 * @brief Inductor: V = L * dI/dt
 *
 * In MNA, the inductor adds a branch current variable and a KVL equation.
 * Trapezoidal companion: V_{n+1} = (2L/h)*I_{n+1} - [(2L/h)*I_n + V_n]
 */
typedef struct {
    spice_component_t base;
    double            inductance;    /**< Inductance in henries             */
    double            initial_i;     /**< Initial current through inductor  */
    double            prev_current;  /**< Current at previous time step     */
    double            prev_voltage;  /**< Voltage at previous time step     */
} spice_inductor_t;

/**
 * @brief Independent voltage source
 *
 * In MNA, voltage sources add a branch current variable and a constraint
 * equation: V(n+) - V(n-) = V_source. This increases the matrix size.
 */
typedef struct {
    spice_component_t base;
    double            dc_value;      /**< DC voltage in volts               */
    double            ac_magnitude;  /**< AC magnitude for AC analysis      */
    double            ac_phase;      /**< AC phase in degrees               */
    int32_t           branch_index;  /**< Index of branch current in MNA    */
} spice_vsource_t;

/**
 * @brief Independent current source
 *
 * Current sources contribute directly to the RHS vector in MNA.
 * They do not add extra variables.
 */
typedef struct {
    spice_component_t base;
    double            dc_value;      /**< DC current in amperes             */
    double            ac_magnitude;  /**< AC magnitude for AC analysis      */
    double            ac_phase;      /**< AC phase in degrees               */
    int32_t           waveform_type; /**< 0=DC, 1=SIN, 2=PULSE, 3=PWL      */
    double            waveform_params[8]; /**< Parameters for time-varying sources */
} spice_isource_t;

/**
 * @brief Diode model parameters (Shockley equation)
 *
 * I_D = I_S * (exp(V_D / (n*V_T)) - 1)
 *
 * V_T = k*T/q ≈ 25.85 mV at 300 K (thermal voltage)
 *
 * Reference: Sedra & Smith §4.3
 */
typedef struct {
    char      model_name[SPICE_MAX_NAME]; /**< Model reference name        */
    double    is_saturation;  /**< Saturation current I_S (A)              */
    double    n_ideality;     /**< Ideality factor n (typically 1-2)       */
    double    series_r;       /**< Series resistance R_S (ohms)            */
    double    transit_time;   /**< Transit time TT (seconds)               */
    double    zero_bias_c;    /**< Zero-bias junction capacitance CJO (F)  */
    double    vj_potential;   /**< Junction built-in potential VJ (V)      */
    double    grading_coeff;  /**< Grading coefficient M                   */
    double    breakdown_v;    /**< Reverse breakdown voltage BV (V)        */
    double    breakdown_i;    /**< Current at breakdown IBV (A)            */
    int32_t   is_used;        /**< Model reference count                   */
} spice_diode_model_t;

/**
 * @brief MOSFET Level 1 model (Shichman-Hodges)
 *
 * Triode region:  I_D = K * (W/L) * [(V_GS - V_TH)*V_DS - V_DS^2/2] * (1 + λ*V_DS)
 * Saturation:     I_D = (K/2) * (W/L) * (V_GS - V_TH)^2 * (1 + λ*V_DS)
 *
 * Reference: "MOSFET Models for VLSI Circuit Simulation" (Arora 1993)
 */
typedef struct {
    char      model_name[SPICE_MAX_NAME];
    double    vth0;           /**< Zero-bias threshold voltage VTO (V)     */
    double    kp;             /**< Transconductance parameter KP (A/V^2)   */
    double    gamma_body;     /**< Body effect parameter GAMMA (V^0.5)    */
    double    phi_surface;    /**< Surface potential PHI (V)               */
    double    lambda_channel; /**< Channel-length modulation LAMBDA (1/V)  */
    double    tox;            /**< Gate oxide thickness TOX (m)            */
    double    cgso;           /**< Gate-source overlap capacitance (F/m)   */
    double    cgdo;           /**< Gate-drain overlap capacitance (F/m)    */
    double    cj;             /**< Bottom junction capacitance (F/m^2)     */
    double    cjsw;           /**< Sidewall junction capacitance (F/m)     */
    double    mj;             /**< Bottom junction grading coefficient     */
    double    mjsw;           /**< Sidewall junction grading coefficient   */
    int32_t   type;           /**< 0=NMOS, 1=PMOS                          */
    int32_t   is_used;
} spice_mos_model_t;

/**
 * @brief BJT model parameters (simplified Ebers-Moll)
 *
 * Forward active: I_C = I_S * exp(V_BE / V_T)
 *                I_B = I_C / β_F
 *
 * Reference: "Semiconductor Device Modeling with SPICE" (Massobrio & Antognetti 1998)
 */
typedef struct {
    char      model_name[SPICE_MAX_NAME];
    double    is_saturation;  /**< Transport saturation current IS (A)     */
    double    bf_forward;     /**< Forward current gain β_F / BF           */
    double    br_reverse;     /**< Reverse current gain β_R / BR           */
    double    nf_coeff;       /**< Forward emission coefficient NF         */
    double    nr_coeff;       /**< Reverse emission coefficient NR         */
    double    vaf_early;      /**< Forward Early voltage VAF (V)           */
    double    var_early;      /**< Reverse Early voltage VAR (V)           */
    double    cje;            /**< B-E zero-bias depletion capacitance (F) */
    double    cjc;            /**< B-C zero-bias depletion capacitance (F) */
    double    cjs;            /**< C-S zero-bias depletion capacitance (F) */
    double    tf_transit;     /**< Forward transit time TF (s)             */
    double    tr_transit;     /**< Reverse transit time TR (s)             */
    int32_t   is_used;
} spice_bjt_model_t;

/**
 * @brief Parsed SPICE netlist
 *
 * Contains all nodes, components, models, and analysis commands
 * extracted from the input file.
 */
typedef struct {
    char                    title[SPICE_MAX_LINE]; /**< First line of netlist */
    spice_node_t            nodes[SPICE_MAX_NODES];      /**< Node table       */
    int32_t                 num_nodes;                    /**< Node count (excl. GND) */
    spice_component_t*      components[SPICE_MAX_COMPONENTS]; /**< Component ptr array */
    int32_t                 num_components;               /**< Total component count  */
    int32_t                 num_resistors;
    int32_t                 num_capacitors;
    int32_t                 num_inductors;
    int32_t                 num_vsources;
    int32_t                 num_isources;
    int32_t                 num_diodes;
    int32_t                 num_bjts;
    int32_t                 num_mosfets;
    spice_diode_model_t     diode_models[32];    /**< .MODEL D cards            */
    int32_t                 num_diode_models;
    spice_mos_model_t       mos_models[32];      /**< .MODEL NMOS/PMOS cards    */
    int32_t                 num_mos_models;
    spice_bjt_model_t       bjt_models[32];      /**< .MODEL NPN/PNP cards      */
    int32_t                 num_bjt_models;
    spice_analysis_cmd_t    analysis_type;       /**< Requested analysis         */
    double                  dc_start, dc_stop, dc_step;   /**< .DC parameters   */
    char                    dc_source_name[SPICE_MAX_NAME];
    double                  ac_fstart, ac_fstop;   /**< .AC parameters (Hz)     */
    int32_t                 ac_points_per_decade;
    double                  tran_tstep, tran_tstop; /**< .TRAN parameters (s)   */
    double                  tran_tstart;
    double                  tran_tmax;              /**< Maximum time step       */
    int32_t                 has_options;
    double                  options_abstol;    /**< Absolute current tolerance   */
    double                  options_vntol;     /**< Absolute voltage tolerance   */
    double                  options_reltol;    /**< Relative tolerance            */
    double                  options_temp;      /**< Circuit temperature (degC)   */
    int32_t                 options_itl1;      /**< DC iteration limit           */
    int32_t                 options_itl4;      /**< Transient iteration limit    */
} spice_netlist_t;

/* ── Parser API ────────────────────────────────────────────────────── */

/**
 * @brief Initialize an empty netlist structure
 * @param nl Pointer to netlist to initialize
 *
 * Complexity: O(1)
 */
void spice_netlist_init(spice_netlist_t *nl);

/**
 * @brief Parse a SPICE netlist from file
 *
 * @param nl       Pointer to netlist structure to fill
 * @param filename Path to SPICE netlist file (.cir)
 * @return 0 on success, negative on parse error
 *
 * Supports: R, L, C, V, I, D, Q (bipolar), M (MOSFET),
 *           .MODEL, .DC, .AC, .TRAN, .OPTIONS, .END
 *
 * Complexity: O(N_lines * N_components)
 */
int spice_netlist_parse_file(spice_netlist_t *nl, const char *filename);

/**
 * @brief Parse a single netlist line
 *
 * @param nl   Netlist to add component to
 * @param line One line of SPICE input
 * @return 0 on success, negative on error
 */
int spice_netlist_parse_line(spice_netlist_t *nl, const char *line);

/**
 * @brief Find or create a node by name
 *
 * @param nl   Netlist to search
 * @param name Node name (case-insensitive match)
 * @return Node ID (>=0), or -1 if node table is full
 *
 * GND is always node 0. All other nodes get sequential IDs starting at 1.
 */
spice_node_id spice_netlist_get_node(spice_netlist_t *nl, const char *name);

/**
 * @brief Free all memory associated with a netlist
 *
 * @param nl Netlist to destroy
 *
 * Frees all dynamically allocated component structures.
 */
void spice_netlist_destroy(spice_netlist_t *nl);

/**
 * @brief Count total number of MNA variables
 *
 * MNA_size = num_nodes + num_vsources + num_inductors + num_mutual_inductors
 *
 * Each voltage source, inductor, and mutual inductor adds one
 * branch current variable to the system.
 *
 * @param nl Netlist
 * @return Total number of equations in MNA system
 */
int32_t spice_netlist_mna_size(const spice_netlist_t *nl);

/**
 * @brief Print netlist summary to stdout
 *
 * @param nl Parsed netlist
 */
void spice_netlist_print_summary(const spice_netlist_t *nl);

#endif /* SPICE_NETLIST_H */
