#ifndef CIRCUIT_ELEMENTS_H
#define CIRCUIT_ELEMENTS_H
#include <stddef.h>
#include <stdint.h>
#include <complex.h>
#ifdef I
#undef I
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* L1 ? Core Definitions: Circuit Elements, Units, Component Models
 * Reference: Ulaby & Ravaioli, Sedra & Smith, Hayt et al.
 * Each element maps to a physical component with governing equation:
 *   Resistor: v=R*i (Ohm 1827), Capacitor: i=C*dv/dt (Maxwell 1861)
 *   Inductor: v=L*di/dt (Faraday 1831) */

typedef enum { UNIT_VOLT=0,UNIT_AMPERE=1,UNIT_OHM=2,UNIT_FARAD=3,UNIT_HENRY=4,UNIT_WATT=5,UNIT_VAR=6,UNIT_VA=7,UNIT_RAD_PER_S=8,UNIT_HERTZ=9,UNIT_SECOND=10,UNIT_DEGREE=11,UNIT_DB=12,UNIT_NEPER=13,UNIT_DIMENSIONLESS=14 } Unit_t;

typedef enum { ELEM_RESISTOR=0,ELEM_CAPACITOR=1,ELEM_INDUCTOR=2,ELEM_DC_VOLTAGE_SRC=3,ELEM_DC_CURRENT_SRC=4,ELEM_AC_VOLTAGE_SRC=5,ELEM_AC_CURRENT_SRC=6,ELEM_VCVS=7,ELEM_VCCS=8,ELEM_CCVS=9,ELEM_CCCS=10,ELEM_SHORT=11,ELEM_OPEN=12,ELEM_DIODE_IDEAL=13,ELEM_TRANSFORMER_IDEAL=14,ELEM_OPAMP_IDEAL=15,ELEM_GYRATOR=16 } ElementType_t;

#define MAX_TERMINALS 4
#define MAX_NAME_LEN 32

typedef struct { int node_id; char label[MAX_NAME_LEN]; } Terminal_t;

/* Resistor: L1 v=R*i, L4 Ohm's Law 1827 */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; double R; double tolerance; double temp_coeff; double power_rating; } Resistor_t;

/* Capacitor: L1 i=C*dv/dt, Maxwell 1861 */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; double C; double tolerance; double temp_coeff; double voltage_rating; double esr; double esl; double leakage_r; } Capacitor_t;

/* Inductor: L1 v=L*di/dt, Faraday 1831 */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; double L; double tolerance; double temp_coeff; double current_rating; double dcr; double core_material; } Inductor_t;

/* DC Voltage Source */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; double V_dc; double R_internal; } DCVoltageSource_t;

/* DC Current Source */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; double I_dc; double R_internal; } DCCurrentSource_t;

/* AC Voltage Source: v(t)=V_offset+V_amplitude*cos(2*pi*f*t+phase) */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; double V_amplitude; double V_rms; double freq_hz; double phase_deg; double V_offset; double R_internal; } ACVoltageSource_t;

/* AC Current Source */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; double I_amplitude; double I_rms; double freq_hz; double phase_deg; double I_offset; double R_internal; } ACCurrentSource_t;

typedef enum { DEP_GAIN_VOLTAGE=0,DEP_GAIN_TRANSCONDUCT=1,DEP_GAIN_TRANSRESIST=2,DEP_GAIN_CURRENT=3 } DependentGainType_t;

typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[2]; int ctrl_node_plus; int ctrl_node_minus; int ctrl_branch; double gain; DependentGainType_t gain_type; } DependentSource_t;

/* Ideal Op-Amp (L2): v+=v-, i+=i-=0 */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[3]; double A_ol; double gbw; } OpAmpIdeal_t;

/* Ideal Transformer (L4): v2=n*v1, i1=-n*i2, n=N2/N1 */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[4]; double n; double Lm; double L_leak_pri; double L_leak_sec; } TransformerIdeal_t;

/* Gyrator (L8): Tellegen 1948, Z1=r^2/Z2 */
typedef struct { int id; char name[MAX_NAME_LEN]; ElementType_t type; Terminal_t t[4]; double r; } Gyrator_t;

/* L3 Math: Complex impedance Z=R+jX */
typedef struct { double re; double im; } Impedance_t;

/* Complex admittance Y=G+jB */
typedef struct { double re; double im; } Admittance_t;

/* Phasor (L2: Steinmetz 1893) */
typedef struct { double magnitude; double phase_deg; double freq_hz; } Phasor_t;

#define MAX_NODES 256
#define MAX_BRANCHES 512
#define MAX_LOOPS 128

typedef struct { int id; int node_from; int node_to; ElementType_t elem_type; int elem_index; double value; char label[MAX_NAME_LEN]; } Branch_t;

typedef struct { int id; char label[MAX_NAME_LEN]; int is_ground; int branch_count; int branches[MAX_BRANCHES]; } Node_t;

typedef struct { int id; char label[MAX_NAME_LEN]; int branch_count; int branches[MAX_BRANCHES]; int directions[MAX_BRANCHES]; } Loop_t;

/* Complete circuit netlist */
typedef struct {
    char name[MAX_NAME_LEN];
    int n_resistors; Resistor_t resistors[MAX_BRANCHES];
    int n_capacitors; Capacitor_t capacitors[MAX_BRANCHES];
    int n_inductors; Inductor_t inductors[MAX_BRANCHES];
    int n_dc_vsrc; DCVoltageSource_t dc_vsrc[MAX_BRANCHES];
    int n_dc_isrc; DCCurrentSource_t dc_isrc[MAX_BRANCHES];
    int n_ac_vsrc; ACVoltageSource_t ac_vsrc[MAX_BRANCHES];
    int n_ac_isrc; ACCurrentSource_t ac_isrc[MAX_BRANCHES];
    int n_dep_src; DependentSource_t dep_src[MAX_BRANCHES];
    int n_opamps; OpAmpIdeal_t opamps[MAX_BRANCHES];
    int n_transformers; TransformerIdeal_t transformers[MAX_BRANCHES];
    int n_nodes; Node_t nodes[MAX_NODES];
    int n_branches; Branch_t branches[MAX_BRANCHES];
    int n_loops; Loop_t loops[MAX_LOOPS];
    int ground_node;
} Circuit_t;

/* Element creation API */
int resistor_init(Resistor_t *r,int id,const char *name,int n1,int n2,double R_ohms);
int capacitor_init(Capacitor_t *c,int id,const char *name,int n1,int n2,double C_farads);
int inductor_init(Inductor_t *l,int id,const char *name,int n1,int n2,double L_henries);
int dc_voltage_source_init(DCVoltageSource_t *s,int id,const char *name,int np,int nm,double V);
int dc_current_source_init(DCCurrentSource_t *s,int id,const char *name,int nf,int nt,double I_val);
int ac_voltage_source_init(ACVoltageSource_t *s,int id,const char *name,int np,int nm,double Va,double f,double ph,double Vo);
int ac_current_source_init(ACCurrentSource_t *s,int id,const char *name,int nf,int nt,double Ia,double f,double ph,double Io);
int opamp_ideal_init(OpAmpIdeal_t *oa,int id,const char *name,int np,int nm,int no,double A);
int transformer_ideal_init(TransformerIdeal_t *tf,int id,const char *name,int p1,int p2,int s1,int s2,double n);
void circuit_init(Circuit_t *ckt,const char *name);
int circuit_add_branch(Circuit_t *ckt,int from,int to,ElementType_t type,int elem_idx,double value,const char *label);
int circuit_add_node(Circuit_t *ckt,const char *label,int is_ground);
int circuit_find_node(const Circuit_t *ckt,const char *label);
int circuit_validate(const Circuit_t *ckt);

/* L3 Complex impedance functions */
Impedance_t impedance_resistor(double R);
Impedance_t impedance_capacitor(double C,double omega);
Impedance_t impedance_inductor(double L,double omega);
Impedance_t impedance_rlc_series(double R,double L,double C,double omega);
Impedance_t impedance_rlc_parallel(double R,double L,double C,double omega);
Admittance_t impedance_to_admittance(Impedance_t Z);
Impedance_t admittance_to_impedance(Admittance_t Y);
Impedance_t impedance_series(Impedance_t Z1,Impedance_t Z2);
Impedance_t impedance_parallel(Impedance_t Z1,Impedance_t Z2);
Admittance_t admittance_parallel(Admittance_t Y1,Admittance_t Y2);
Impedance_t phasor_to_impedance(const Phasor_t *p);
Phasor_t impedance_to_phasor(Impedance_t Z,double freq_hz);
double q_factor_series(double R,double L,double C);
double q_factor_parallel(double R,double L,double C);
double resonance_omega(double L,double C);
double resonance_freq_hz(double L,double C);
double bandwidth_series_radps(double R,double L);
double bandwidth_parallel_radps(double R,double C);
double damping_factor_series(double R,double L,double C);
double damping_factor_parallel(double R,double L,double C);
#endif
