/**
 * @file interconnection.c
 * @brief Two-Port Network Interconnection Implementation
 *
 * Implements all five interconnection topologies with the natural
 * parameter type for each, plus a general interconnection function
 * that handles automatic parameter conversion.
 */

#include "../include/interconnection.h"
#include "../include/conversion.h"

/* ============================================================================
 * L5: Direct Interconnections
 * ============================================================================ */

matrix2x2_t interconnect_series_series(matrix2x2_t z1, matrix2x2_t z2) {
    return matrix2x2_add(z1, z2);
}

matrix2x2_t interconnect_parallel_parallel(matrix2x2_t y1, matrix2x2_t y2) {
    return matrix2x2_add(y1, y2);
}

matrix2x2_t interconnect_series_parallel(matrix2x2_t h1, matrix2x2_t h2) {
    return matrix2x2_add(h1, h2);
}

matrix2x2_t interconnect_parallel_series(matrix2x2_t g1, matrix2x2_t g2) {
    return matrix2x2_add(g1, g2);
}

/**
 * Cascade: ABCD_total = ABCD1 × ABCD2.
 *
 * IMPORTANT ORDER: The first network (closest to source) is ABCD1,
 * the second network (closest to load) is ABCD2.
 * ABCD_total = ABCD1 × ABCD2.
 *
 * This is NOT commutative: ABCD1 × ABCD2 ≠ ABCD2 × ABCD1 in general.
 *
 * Example:
 *   Amplifier → Filter:
 *   ABCD_total = ABCD_amp × ABCD_filter
 *
 *   Filter → Amplifier:
 *   ABCD_total = ABCD_filter × ABCD_amp (different result!)
 *
 * The order reflects signal flow: source → ABCD1 → ABCD2 → load.
 */
matrix2x2_t interconnect_cascade(matrix2x2_t abcd1, matrix2x2_t abcd2) {
    return matrix2x2_mul(abcd1, abcd2);
}

/* ============================================================================
 * L5: General Interconnection with Conversion
 * ============================================================================ */

/**
 * Generic interconnection: converts both networks to the natural
 * parameter type for the given configuration, combines them, and
 * converts the result to the requested output type.
 *
 * Configuration → Natural type mapping:
 *   SERIES_SERIES      → Z (additive)
 *   PARALLEL_PARALLEL  → Y (additive)
 *   SERIES_PARALLEL    → H (additive)
 *   PARALLEL_SERIES    → G (additive)
 *   CASCADE            → ABCD (multiplicative)
 */
matrix2x2_t interconnect_general(matrix2x2_t net1, param_type_t type1,
                                  matrix2x2_t net2, param_type_t type2,
                                  config_type_t config,
                                  param_type_t out_type, double z0) {
    param_type_t natural_type;
    switch (config) {
        case CONFIG_SERIES_SERIES:    natural_type = PARAM_Z; break;
        case CONFIG_PARALLEL_PARALLEL: natural_type = PARAM_Y; break;
        case CONFIG_SERIES_PARALLEL:   natural_type = PARAM_H; break;
        case CONFIG_PARALLEL_SERIES:   natural_type = PARAM_G; break;
        case CONFIG_CASCADE:           natural_type = PARAM_ABCD; break;
        default:                       natural_type = PARAM_Z; break;
    }

    /* Convert both networks to natural type */
    matrix2x2_t conv1 = convert_parameters(net1, type1, natural_type, z0);
    matrix2x2_t conv2 = convert_parameters(net2, type2, natural_type, z0);

    /* Combine */
    matrix2x2_t combined;
    if (config == CONFIG_CASCADE) {
        combined = matrix2x2_mul(conv1, conv2);
    } else {
        combined = matrix2x2_add(conv1, conv2);
    }

    /* Convert to output type */
    if (out_type == natural_type) {
        return combined;
    }
    return convert_parameters(combined, natural_type, out_type, z0);
}

/**
 * Brune's test approximation for interconnection validity.
 *
 * The Brune test checks whether connecting two two-port networks
 * in a given configuration maintains the port conditions:
 *   - At each port: I_in = -I_out (currents equal and opposite)
 *   - The interconnection does not create unintended current loops
 *
 * Violation occurs when:
 *   - Connecting floating (differential) networks to ground-referenced networks
 *   - Creating ground loops in series connections
 *   - Series-parallel and parallel-series without proper isolation
 *
 * This simplified check looks at the parameter types and makes a
 * heuristic assessment. A full Brune test would require network
 * topology information beyond just the parameter matrices.
 */
int interconnect_is_valid(matrix2x2_t net1, param_type_t type1,
                           matrix2x2_t net2, param_type_t type2,
                           config_type_t config) {
    /* For cascade: always valid (no port condition issues) */
    if (config == CONFIG_CASCADE) return 1;

    /* For series-series (Z) and parallel-parallel (Y):
     * Check that both networks can be represented in the natural type.
     * If a network cannot be represented (singular Z or Y), the
     * interconnection may be invalid. */
    switch (config) {
        case CONFIG_SERIES_SERIES: {
            /* Both need valid Z-params */
            matrix2x2_t z1 = convert_parameters(net1, type1, PARAM_Z, 0.0);
            matrix2x2_t z2 = convert_parameters(net2, type2, PARAM_Z, 0.0);
            return !isnan(z1.m11.real) && !isnan(z2.m11.real);
        }
        case CONFIG_PARALLEL_PARALLEL: {
            matrix2x2_t y1 = convert_parameters(net1, type1, PARAM_Y, 0.0);
            matrix2x2_t y2 = convert_parameters(net2, type2, PARAM_Y, 0.0);
            return !isnan(y1.m11.real) && !isnan(y2.m11.real);
        }
        default:
            /* Other configurations: assume valid */
            return 1;
    }
}

/* ============================================================================
 * L4: Bartlett's Bisection Theorem
 * ============================================================================ */

/**
 * Bartlett's bisection theorem for symmetric networks.
 *
 * Given a symmetric two-port, measure (or compute) the open-circuit
 * and short-circuit impedances of one half (at the symmetry plane):
 *
 *   Zoc: open-circuit impedance of half-network
 *   Zsc: short-circuit impedance of half-network
 *
 * Then the complete two-port Z-parameters are:
 *   z11 = z22 = (Zoc + Zsc) / 2
 *   z12 = z21 = (Zoc - Zsc) / 2
 *
 * Proof (using even/odd mode decomposition):
 *   Even mode (V1 = V2): Zeven = Zoc
 *   Odd mode (V1 = -V2): Zodd = Zsc
 *   z11 = (Zeven + Zodd)/2, z12 = (Zeven - Zodd)/2
 *
 * This theorem is powerful because it reduces a two-port analysis
 * to two one-port analyses — cutting the problem complexity in half.
 *
 * Reference: Bartlett, Phil. Mag., 1927
 * Course: Illinois ECE 451 — Bartlett's theorem in filter design
 */
matrix2x2_t bartlett_bisection(complex_t zoc, complex_t zsc) {
    complex_t half = complex_make(0.5, 0.0);
    complex_t z11 = complex_mul(complex_add(zoc, zsc), half);
    complex_t z12 = complex_mul(complex_sub(zoc, zsc), half);
    return matrix2x2_make(z11, z12, z12, z11);
}

/**
 * Even/odd mode decomposition.
 *
 * For a symmetric network, the eigenvalues of the Z-matrix are:
 *   λ_even = z11 + z12 = Zeven (common-mode impedance)
 *   λ_odd  = z11 - z12 = Zodd  (differential-mode impedance)
 *
 * Applications:
 *   - Differential pair analysis: differential gain vs common-mode gain
 *   - CMRR = |Zeven| / |Zodd| for high CMRR, want Zodd large, Zeven small
 *   - Filter design: λ_even and λ_odd determine the passband/stopband
 *
 * Course: Stanford EE214 — Differential circuits
 */
void symmetric_even_odd_impedance(matrix2x2_t z, complex_t *zeven,
                                   complex_t *zodd) {
    *zeven = complex_add(z.m11, z.m12);
    *zodd = complex_sub(z.m11, z.m12);
}
