#ifndef CIRCUIT_UTILS_H
#define CIRCUIT_UTILS_H
#include "circuit_elements.h"
#include "ac_analysis.h"
#include <stdio.h>

/* Utility functions for circuit analysis.
 * Unit conversions, formatting, validation, and helper functions.
 * Reference: IEEE Std 315, IEC 60027 */

/* Unit conversion helpers (L1) */
double deg_to_rad(double deg);
double rad_to_deg(double rad);
double hz_to_radps(double hz);
double radps_to_hz(double radps);
double celsius_to_kelvin(double celsius);
double db_to_linear(double db);
double linear_to_db(double linear);
double dbm_to_watts(double dbm);
double watts_to_dbm(double watts);
double dbm_to_volts_rms(double dbm, double Z0);
double volts_rms_to_dbm(double V_rms, double Z0);

/* SI prefix scaling */
double pico_to_base(double p);   /* 1e-12 */
double nano_to_base(double n);   /* 1e-9 */
double micro_to_base(double u);  /* 1e-6 */
double milli_to_base(double m);  /* 1e-3 */
double kilo_to_base(double k);   /* 1e3 */
double mega_to_base(double M);   /* 1e6 */
double giga_to_base(double G);   /* 1e9 */

double base_to_milli(double v);
double base_to_micro(double v);
double base_to_nano(double v);
double base_to_pico(double v);

/* Standard resistor values (E12, E24 series) */
double nearest_e12_value(double R);
double nearest_e24_value(double R);
int is_standard_e12(double R, double tol);
int is_standard_e24(double R, double tol);

/* Round to significant digits */
double round_to_sigdig(double val, int n_digits);
double round_to_decimals(double val, int n_decimals);

/* String formatting for display */
const char* unit_to_string(Unit_t unit);
const char* element_type_to_string(ElementType_t type);
void print_impedance(Impedance_t Z, FILE *fp);
void print_complex_power(ComplexPower_t cp, FILE *fp);
void print_resonance(Resonance_t res, FILE *fp);

/* Error handling */
typedef enum {
    ERR_NONE = 0,
    ERR_NULL_POINTER,
    ERR_DIVISION_BY_ZERO,
    ERR_SINGULAR_MATRIX,
    ERR_NO_CONVERGENCE,
    ERR_INVALID_ARGUMENT,
    ERR_OUT_OF_MEMORY,
    ERR_CIRCUIT_NOT_CONNECTED,
    ERR_FLOATING_NODE,
    ERR_VOLTAGE_LOOP,
    ERR_CURRENT_CUTSET,
    ERR_MAX_ELEMENTS,
    ERR_MAX_NODES,
    ERR_UNKNOWN_ELEMENT,
} ErrorCode_t;

const char* error_string(ErrorCode_t err);

/* Tolerance analysis (L8: Monte Carlo support) */
double random_normal(double mean, double stddev);
double random_uniform(double a, double b);

/* Resistor value with tolerance applied */
double resistor_with_tolerance(double nominal, double tol_pct);

/* Comparator for double near-equality (within relative tolerance) */
int double_near(double a, double b, double rel_tol);

/* Comparator for double near-equality with absolute tolerance */
int double_near_abs(double a, double b, double abs_tol);

/* Check if a floating-point value is valid (not NaN, not Inf) */
int is_valid_double(double x);

/* CRC-16 for netlist verification (detect data corruption) */
unsigned short crc16_netlist(const Circuit_t *ckt);

/* Memory-safe string copy (ensures null termination) */
size_t safe_strcpy(char *dst, const char *src, size_t dst_size);

/* Compute parallel equivalent: returns (a*b)/(a+b) safely */
double parallel_equiv(double a, double b);

/* Three-element parallel: (a*b*c)/(a*b + b*c + c*a) */
double parallel_equiv_3(double a, double b, double c);

/* AWG wire gauge functions */
double awg_to_diameter_mm(int awg);
double awg_to_resistance_per_meter(int awg);
int nearest_awg_from_diameter_mm(double dia_mm);

#endif /* CIRCUIT_UTILS_H */
