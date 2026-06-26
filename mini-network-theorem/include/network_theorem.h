/**
 * @file network_theorem.h
 * @brief Network Theorems — Core Definitions and Data Structures
 *
 * This header defines all fundamental data types for circuit network theorem
 * analysis, including Thevenin, Norton, Superposition, Maximum Power Transfer,
 * Reciprocity, Millman, Tellegen, and Two-Port network parameter representations.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: typedef struct for all core circuit entities
 *   L2 - Core Concepts: Equivalent circuits, linearity, source transformation
 *   L3 - Mathematical Structures: Complex impedance, matrices for nodal/mesh
 *   L4 - Fundamental Laws: Theorems validated via implemented algorithms
 *
 * Reference: Sedra & Smith "Microelectronic Circuits" (2020) Ch 1, App D
 *            Desoer & Kuh "Basic Circuit Theory" (1969)
 */

#ifndef NETWORK_THEOREM_H
#define NETWORK_THEOREM_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Core Circuit Element Definitions
 * ========================================================================== */

/** Complex impedance: Z = R + jX (ohms); admittance: Y = G + jB (siemens) */
typedef struct {
    double real;    /**< Resistance R (Ω) */
    double imag;    /**< Reactance X (Ω), positive = inductive, negative = capacitive */
} ComplexImpedance;

/** Complex admittance: Y = G + jB */
typedef struct {
    double real;    /**< Conductance G (S) */
    double imag;    /**< Susceptance B (S) */
} ComplexAdmittance;

/** Independent voltage source: v(t) = V_dc + V_ac * cos(ωt + φ) */
typedef struct {
    double dc_offset;       /**< DC component (V) */
    double ac_amplitude;    /**< AC peak amplitude (V) */
    double angular_freq;    /**< Angular frequency ω (rad/s) */
    double phase;           /**< Phase angle φ (radians) */
} VoltageSource;

/** Independent current source: i(t) = I_dc + I_ac * cos(ωt + φ) */
typedef struct {
    double dc_offset;       /**< DC component (A) */
    double ac_amplitude;    /**< AC peak amplitude (A) */
    double angular_freq;    /**< Angular frequency ω (rad/s) */
    double phase;           /**< Phase angle φ (radians) */
} CurrentSource;

/** Circuit element types for unified element representation */
typedef enum {
    ELEM_RESISTOR = 0,
    ELEM_CAPACITOR,
    ELEM_INDUCTOR,
    ELEM_VOLTAGE_SOURCE,
    ELEM_CURRENT_SOURCE,
    ELEM_VCVS,
    ELEM_VCCS,
    ELEM_CCVS,
    ELEM_CCCS,
    ELEM_OPEN,
    ELEM_SHORT
} ElementType;

/** A single circuit branch element */
typedef struct {
    ElementType type;
    uint32_t node_from;
    uint32_t node_to;
    double value;
    uint32_t ctrl_elem;
} CircuitElement;

/*
 * L1: Thevenin Equivalent Circuit
 * Thevenin: Any linear two-terminal network can be replaced by
 * a single voltage source V_th in series with an impedance Z_th.
 * V_th = open-circuit voltage at the terminals
 * Z_th = equivalent impedance with all independent sources deactivated
 */
typedef struct {
    ComplexImpedance Z_th;
    VoltageSource    V_th;
    double           frequency;
} TheveninEquivalent;

/*
 * L1: Norton Equivalent Circuit
 * Norton: Any linear two-terminal network can be replaced by
 * a single current source I_n in parallel with admittance Y_n.
 * I_n = short-circuit current, Y_n = 1/Z_th
 */
typedef struct {
    ComplexAdmittance Y_n;
    CurrentSource     I_n;
    double            frequency;
} NortonEquivalent;

/*
 * L1: Superposition Principle
 * In a linear circuit with N independent sources, the total response
 * is the sum of responses due to each source acting alone.
 */
typedef struct {
    double    total_response;
    double   *individual_contributions;
    uint32_t  num_sources;
    uint8_t   is_voltage;
} SuperpositionResult;

/*
 * L1: Maximum Power Transfer
 * For a source Z_s = R_s + jX_s delivering to load Z_L = R_L + jX_L:
 *   Conjugate match (AC): Z_L = Z_s* -> P_max = |V_s|^2/(8R_s)
 *   Resistive match (DC): R_L = R_s -> P_max = V_s^2/(4R_s)
 */
typedef struct {
    double source_impedance_real;
    double source_impedance_imag;
    double optimal_load_real;
    double optimal_load_imag;
    double max_power;
    double source_voltage_rms;
    uint8_t is_ac;
} MaxPowerTransferResult;

/*
 * L1: Reciprocity Theorem
 * In a linear bilateral network, excitation/response ratio is unchanged
 * when positions are interchanged. For two-ports: Z_12 = Z_21.
 */
typedef struct {
    double z11, z12, z21, z22;
    double y11, y12, y21, y22;
    double h11, h12, h21, h22;
    uint8_t is_reciprocal;
    double reciprocity_error;
} ReciprocityParams;

/* L1: Z-parameters: [V]=[z][I], measured with output open */
typedef struct {
    ComplexImpedance z11, z12, z21, z22;
} ZParameters;

/* L1: Y-parameters: [I]=[y][V], measured with output short */
typedef struct {
    ComplexAdmittance y11, y12, y21, y22;
} YParameters;

/* L1: H-parameters: hybrid, [V1,I2]' = H*[I1,V2]' */
typedef struct {
    ComplexImpedance  h11;
    double            h12_real, h12_imag;
    double            h21_real, h21_imag;
    ComplexAdmittance h22;
} HParameters;

/* L1: ABCD transmission parameters: [V1,I1]' = ABCD*[V2,-I2]' */
typedef struct {
    double A_real, A_imag;
    double B_real, B_imag;
    double C_real, C_imag;
    double D_real, D_imag;
} ABCDParameters;

/* L1: S-parameters for RF/microwave, normalized to Z0 */
typedef struct {
    ComplexImpedance s11, s12, s21, s22;
    double Z0;
} SParameters;

/* L1: Network Topology */
typedef struct {
    uint32_t id;
    double   voltage;
    uint8_t  is_ground;
} CircuitNode;

typedef struct {
    uint32_t id;
    uint32_t node_from;
    uint32_t node_to;
    uint32_t element_idx;
    double   current;
    double   voltage;
} CircuitBranch;

typedef struct {
    uint32_t  id;
    uint32_t *branch_ids;
    uint32_t  num_branches;
    double    mesh_current;
} CircuitMesh;

typedef struct {
    CircuitNode    *nodes;
    uint32_t        num_nodes;
    CircuitBranch  *branches;
    uint32_t        num_branches;
    CircuitElement *elements;
    uint32_t        num_elements;
    CircuitMesh    *meshes;
    uint32_t        num_meshes;
    uint32_t        ground_node;
} CircuitTopology;

/*
 * L1: Millman's Theorem
 * For N parallel branches: V_common = (Σ V_k/Z_k) / (Σ 1/Z_k)
 */
typedef struct {
    double   *branch_voltages;
    double   *branch_impedances;
    uint32_t  num_branches;
    double    common_node_voltage;
} MillmanResult;

/*
 * L1: Tellegen's Theorem
 * For any lumped network: Σ(v_k * i'_k) = 0 across any two states
 */
typedef struct {
    double   *branch_voltages_state1;
    double   *branch_currents_state2;
    uint32_t  num_branches;
    double    power_sum;
    double    tolerance;
} TellegenResult;

/* L1: Wye-Delta Transformation */
typedef struct {
    double R1, R2, R3;
} WyeNetwork;

typedef struct {
    double R12, R23, R31;
} DeltaNetwork;

/* L1: Compensation Theorem */
typedef struct {
    uint32_t  modified_branch;
    double    original_impedance;
    double    new_impedance;
    double    delta_z;
    double    original_current;
    double    compensation_voltage;
} CompensationResult;

/* L1: Substitution Theorem */
typedef struct {
    uint32_t  branch_id;
    double    known_voltage;
    double    known_current;
    uint8_t   replace_with_voltage;
} SubstitutionResult;

/* L3: Modified Nodal Analysis (MNA) Matrix — SPICE-style */
typedef struct {
    double  *A;
    double  *b;
    double  *x;
    uint32_t dim;
    uint32_t n_nodes;
    uint32_t n_vsources;
} MNASystem;

/* L5: Network Reduction State */
typedef struct {
    CircuitTopology *topology;
    uint32_t          step;
    uint32_t          max_steps;
    double            tolerance;
} NetworkReduction;

/* ==========================================================================
 * API Function Declarations
 * ========================================================================== */

/* L2: Core Concept Implementations */
void thevenin_to_norton(const TheveninEquivalent *thev, NortonEquivalent *nort);
void norton_to_thevenin(const NortonEquivalent *nort, TheveninEquivalent *thev);
ComplexImpedance impedance_resistor(double R);
ComplexImpedance impedance_capacitor(double C, double frequency);
ComplexImpedance impedance_inductor(double L, double frequency);
ComplexImpedance impedance_series(ComplexImpedance z1, ComplexImpedance z2);
ComplexImpedance impedance_parallel(ComplexImpedance z1, ComplexImpedance z2);
ComplexAdmittance impedance_to_admittance(ComplexImpedance z);
ComplexImpedance admittance_to_impedance(ComplexAdmittance y);

/* L4: Fundamental Theorem Implementations */
int compute_thevenin(const CircuitTopology *ckt, uint32_t a, uint32_t b, TheveninEquivalent *result);
int compute_norton(const CircuitTopology *ckt, uint32_t a, uint32_t b, NortonEquivalent *result);
int superposition_solve(const CircuitTopology *ckt, uint32_t target, SuperpositionResult *result);
int max_power_transfer(const TheveninEquivalent *source, MaxPowerTransferResult *result);
int verify_reciprocity(const ZParameters *zp, ReciprocityParams *result);
int millman_solve(const double *v, const double *z, uint32_t n, MillmanResult *result);
int tellegen_verify(const double *v1, const double *i2, uint32_t n, double tol, TellegenResult *result);

/* L3: Matrix-Based Solvers */
double *build_nodal_admittance_matrix(const CircuitTopology *ckt, uint32_t *n_out);
int solve_nodal_voltages(double *Y, double *I, uint32_t n, double *V);
int build_current_vector(const CircuitTopology *ckt, double *I, uint32_t n);
double *build_mesh_impedance_matrix(const CircuitTopology *ckt, uint32_t *m_out);
int solve_mesh_currents(double *Z, double *V, uint32_t n, double *I);
int lu_decompose(double *A, uint32_t n);
void lu_solve(const double *LU, const double *b, uint32_t n, double *x);
int gaussian_elimination(double *A, double *b, uint32_t n, double *x);

/* L5: Transformation Methods */
void wye_to_delta(const WyeNetwork *wye, DeltaNetwork *delta);
void delta_to_wye(const DeltaNetwork *delta, WyeNetwork *wye);
void source_transform_thevenin_to_norton(double V, double Z, double *I_out, double *Z_out);
void source_transform_norton_to_thevenin(double I, double Z, double *V_out, double *Z_out);
int network_reduce(CircuitTopology *ckt, NetworkReduction *state);

/* L2: Two-Port Parameter Conversions */
int z_to_y(const ZParameters *z, YParameters *y);
int y_to_z(const YParameters *y, ZParameters *z);
int z_to_h(const ZParameters *z, HParameters *h);
int h_to_z(const HParameters *h, ZParameters *z);
int z_to_abcd(const ZParameters *z, ABCDParameters *abcd);
void abcd_cascade(const ABCDParameters *a, const ABCDParameters *b, ABCDParameters *result);
int z_to_s(const ZParameters *z, double Z0, SParameters *s);
int s_to_z(const SParameters *s, double Z0, ZParameters *z);

/* L5: Advanced Solvers */
int mna_build_system(const CircuitTopology *ckt, MNASystem *sys);
int mna_solve(MNASystem *sys);
void mna_free(MNASystem *sys);

/* L2: Utility Functions */
double impedance_magnitude(ComplexImpedance z);
double impedance_phase(ComplexImpedance z);
double power_dissipated(double voltage, double current, double resistance);
double reactive_power(double current, double reactance);
double apparent_power(double v_rms, double i_rms);
double power_factor(double v_phase, double i_phase);
CircuitTopology *topology_create(uint32_t nn, uint32_t nb, uint32_t ne, uint32_t nm);
void topology_free(CircuitTopology *ckt);

/* --- L3: Additional Matrix Operations --- */
double compute_residual_norm(const double *A, const double *x, const double *b, uint32_t n);
double matrix_condition_estimate(const double *A, uint32_t n);
double matrix_determinant(const double *A, uint32_t n);
int matrix_inverse(const double *A, uint32_t n, double *inv_A);

/* --- L5: Series/Parallel Combination Functions --- */
double combine_series_resistors(double R1, double R2);
double combine_parallel_resistors(double R1, double R2);
double combine_series_capacitors(double C1, double C2);
double combine_parallel_capacitors(double C1, double C2);
double combine_series_inductors(double L1, double L2);
double combine_parallel_inductors(double L1, double L2);

/* --- L5: Source Deactivation --- */
int deactivate_independent_sources(CircuitElement *elements, uint32_t count);

/* --- L5: Superposition Verification --- */
double superposition_verify(const CircuitTopology *ckt, uint32_t target_node);

/* --- L4: Compensation Theorem --- */
int compensation_compute(uint32_t branch_idx, double Z_orig, double Z_new,
                          double I_orig, CompensationResult *result);

/* --- L6: Wheatstone Bridge --- */
typedef struct {
    double R1, R2, R3, R4;
    double R_g;
    double V_supply;
    double V_output;
    double I_g;
    int    is_balanced;
    double balance_ratio;
} WheatstoneBridge;

int wheatstone_analyze(double R1, double R2, double R3, double R4,
                        double R_g, double V_supply, WheatstoneBridge *bridge);

/* --- L6: Voltage Divider --- */
typedef struct {
    double R1, R2, R_load, V_in, V_out;
    double I_total, I_load, output_impedance;
} VoltageDivider;

int voltage_divider_analyze(double R1, double R2, double R_load,
                             double V_in, VoltageDivider *div);

/* --- L6: Current Divider --- */
typedef struct {
    double R1, R2, I_in, I1, I2, V_parallel;
} CurrentDivider;

int current_divider_analyze(double R1, double R2, double I_in,
                             CurrentDivider *div);

/* --- L2: Two-Port Analysis --- */
ComplexImpedance two_port_input_impedance(const ZParameters *zp, ComplexImpedance ZL);
ComplexImpedance two_port_output_impedance(const ZParameters *zp, ComplexImpedance ZS);
double two_port_voltage_gain(const ZParameters *zp, double ZL);

/* --- L7: Audio Amplifier Matching --- */
typedef struct {
    double amplifier_Zout;
    double speaker_Znom;
    double speaker_Zmin;
    double amplifier_Vrms;
    double matched_power;
    double actual_power;
    double efficiency;
    double damping_factor;
} AudioMatchingResult;

int audio_amplifier_match(double amp_Zout, double speaker_Z,
                           double amp_Vrms, AudioMatchingResult *result);

/* --- L8: Sparse Matrix CSR --- */
typedef struct {
    double   *values;
    uint32_t *col_idx;
    uint32_t *row_ptr;
    uint32_t  nnz;
    uint32_t  n_rows;
    uint32_t  n_cols;
} SparseMatrixCSR;

SparseMatrixCSR *dense_to_csr(const double *A, uint32_t n);
void csr_mat_vec_mul(const SparseMatrixCSR *csr, const double *x, double *y);
void csr_free(SparseMatrixCSR *csr);

/* --- L8: Newton-Raphson --- */
int newton_raphson_step(double *x, uint32_t n,
                        int (*eval_F)(const double *x, double *F, void *ctx),
                        int (*eval_J)(const double *x, double *J, void *ctx),
                        void *ctx, double *dx_norm);

#ifdef __cplusplus
}
#endif

#endif /* NETWORK_THEOREM_H */
