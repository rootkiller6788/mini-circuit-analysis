/**
 * @file phasor_operations.h
 * @brief Phasor arithmetic, dq0 transforms, and three-phase system analysis
 *
 * Phasors are complex-number representations of sinusoidal quantities
 * that simplify AC circuit analysis by converting differential equations
 * to algebraic ones.
 *
 * A sinusoidal quantity v(t) = V_m cos(ωt + φ) is represented by the phasor
 *   V = (V_m/√2) ∠φ = V_rms × e^{jφ}
 *
 * References:
 *   - C.P. Steinmetz, "Theory and Calculation of Alternating Current
 *     Phenomena" (1897)
 *   - R.H. Park, "Two-Reaction Theory of Synchronous Machines" (1929)
 *   - P.C. Krause, "Analysis of Electric Machinery" (3rd ed, 2013)
 *   - Fortescue, "Method of Symmetrical Coordinates" (1918, AIEE Trans.)
 *   - MIT 6.685 "Electric Machines" / Berkeley EE137A
 *
 * Knowledge coverage:
 *   L1 (Definitions): Phasor, RMS phasor, sequence components, dq0 frame
 *   L2 (Concepts):   Phasor diagrams, leading/lagging, balanced/unbalanced
 *   L3 (Math):       Fortescue transform, Park/Clarke transform, matrix ops
 *   L5 (Algorithms): Symmetrical components computation, dq0 reference frame
 */

#ifndef PHASOR_OPERATIONS_H
#define PHASOR_OPERATIONS_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ==========================================================================
 * L1: Phasor Representation Types
 * ========================================================================== */

/**
 * @brief Single phasor in polar and rectangular forms.
 *
 * V = V_mag ∠φ = V_mag × e^{jφ} = V_re + j·V_im
 *
 * For RMS phasors, V_mag = V_rms (already divided by √2).
 * For peak phasors, V_mag = V_peak.
 */
typedef struct {
    double   magnitude;       /**< Phasor magnitude (RMS or peak) */
    double   angle_rad;       /**< Phase angle [radians] */
    double   angle_deg;       /**< Phase angle [degrees] */
    double complex rect;      /**< Rectangular form: Re + j Im */
} phasor_t;

/**
 * @brief Three-phase phasor set (ABC frame).
 */
typedef struct {
    phasor_t phase_a;         /**< Phase A phasor */
    phasor_t phase_b;         /**< Phase B phasor */
    phasor_t phase_c;         /**< Phase C phasor */
    uint8_t  is_balanced;     /**< Flag: magnitudes equal & 120° separation */
    uint8_t  is_positive_seq; /**< Flag: ABC rotation (A-B-C = positive) */
} phasor_abc_t;

/**
 * @brief Symmetrical components (012 frame).
 *
 * Fortescue transform (1918):
 *   [V0]    [1   1   1 ] [Va]
 *   [V1] =  [1   a   a²] [Vb]   where a = 1∠120° = e^{j2π/3}
 *   [V2]    [1   a²  a ] [Vc]
 *
 * V0 = zero sequence (in-phase, same magnitude)
 * V1 = positive sequence (balanced ABC rotation)
 * V2 = negative sequence (balanced ACB rotation)
 */
typedef struct {
    phasor_t zero_seq;        /**< Zero-sequence component V0 */
    phasor_t pos_seq;         /**< Positive-sequence component V1 */
    phasor_t neg_seq;         /**< Negative-sequence component V2 */
    double   unbalance_percent; /**< V2/V1 × 100% — voltage unbalance factor */
    double   zero_seq_percent;  /**< V0/V1 × 100% */
} phasor_012_t;

/**
 * @brief dq0 reference frame quantities (Park transform).
 *
 * The Park transform (1929) converts ABC stationary-frame quantities
 * to rotating dq0 reference frame, making AC quantities appear DC
 * in steady state — essential for machine control and grid-sync PLL.
 *
 *   [Vd]   [cos(θ)  cos(θ-2π/3)  cos(θ+2π/3)] [Va]
 *   [Vq] = [-sin(θ) -sin(θ-2π/3) -sin(θ+2π/3)] [Vb]
 *   [V0]   [1/2      1/2          1/2        ] [Vc]
 *
 * where θ = ωt + θ₀ is the reference frame angle.
 */
typedef struct {
    double   d;               /**< Direct-axis component */
    double   q;               /**< Quadrature-axis component */
    double   zero;            /**< Zero-sequence component */
    double   magnitude;       /**< sqrt(d² + q²) */
    double   angle_rad;       /**< atan2(q, d) — frame angle offset */
} dq0_t;

/* ==========================================================================
 * L1: Phasor Construction and Conversion
 * ========================================================================== */

/**
 * @brief Create a phasor from polar coordinates.
 *
 * @param magnitude  Magnitude (RMS or peak)
 * @param angle_rad  Angle in radians
 * @return Populated phasor
 */
phasor_t phasor_from_polar(double magnitude, double angle_rad);

/**
 * @brief Create a phasor from rectangular coordinates.
 *
 * @param re  Real part
 * @param im  Imaginary part
 * @return Populated phasor
 */
phasor_t phasor_from_rect(double re, double im);

/**
 * @brief Create a phasor from peak magnitude and degrees.
 *
 * Converts peak to RMS by dividing by √2.
 *
 * @param v_peak     Peak voltage/current
 * @param angle_deg  Phase angle in degrees
 * @return RMS phasor
 */
phasor_t phasor_from_peak(double v_peak, double angle_deg);

/**
 * @brief Convert a phasor to its complex conjugate.
 *
 * Useful for complex power computation: S = V × I*
 */
phasor_t phasor_conjugate(const phasor_t *p);

/* ==========================================================================
 * L2: Phasor Arithmetic Operations
 * ========================================================================== */

/**
 * @brief Add two phasors (KVL: series voltages add).
 *
 * V_total = V1 + V2  (rectangular addition)
 */
phasor_t phasor_add(const phasor_t *a, const phasor_t *b);

/**
 * @brief Subtract two phasors.
 */
phasor_t phasor_sub(const phasor_t *a, const phasor_t *b);

/**
 * @brief Multiply two phasors (magnitudes multiply, angles add).
 */
phasor_t phasor_mul(const phasor_t *a, const phasor_t *b);

/**
 * @brief Divide two phasors (magnitudes divide, angles subtract).
 *
 * Used for impedance computation: Z = V / I
 *
 * @param a  Numerator phasor
 * @param b  Denominator phasor (must have non-zero magnitude)
 * @return Quotient phasor, magnitude=0 if division by zero
 */
phasor_t phasor_div(const phasor_t *a, const phasor_t *b);

/**
 * @brief Scale a phasor by a real scalar.
 */
phasor_t phasor_scale(const phasor_t *p, double scalar);

/* ==========================================================================
 * L3: Symmetrical Components (Fortescue Transform)
 * ========================================================================== */

/**
 * @brief Compute symmetrical components from ABC phasors.
 *
 * Implements the Fortescue transformation:
 *   [V0]   [1   1   1 ] [Va]
 *   [V1] = [1   a   a²] [Vb] × 1/3
 *   [V2]   [1   a²  a ] [Vc]
 *
 * @param abc  Three-phase phasors in ABC frame
 * @param zpn  [out] Symmetrical components (zero/pos/neg sequence)
 * @return 0 on success
 */
int phasor_abc_to_012(const phasor_abc_t *abc, phasor_012_t *zpn);

/**
 * @brief Reconstruct ABC phasors from symmetrical components.
 *
 * Inverse Fortescue transform:
 *   [Va]   [1   1   1 ] [V0]
 *   [Vb] = [1   a²  a ] [V1]
 *   [Vc]   [1   a   a²] [V2]
 *
 * @param zpn  Symmetrical components
 * @param abc  [out] Reconstructed ABC phasors
 * @return 0 on success
 */
int phasor_012_to_abc(const phasor_012_t *zpn, phasor_abc_t *abc);

/**
 * @brief Compute three-phase power from symmetrical components.
 *
 * S_3phase = 3 × (V1 × I1* + V2 × I2* + V0 × I0*)
 *
 * For balanced systems (V2=V0=I2=I0=0): S = 3 × V1 × I1*
 *
 * @param v_zpn  Voltage symmetrical components
 * @param i_zpn  Current symmetrical components
 * @param s_total [out] Total complex power (P + jQ)
 * @return 0 on success
 */
int phasor_power_from_components(const phasor_012_t *v_zpn,
                                  const phasor_012_t *i_zpn,
                                  double complex *s_total);

/* ==========================================================================
 * L3: Clarke and Park Transforms (αβ0 and dq0)
 * ========================================================================== */

/**
 * @brief Clarke transform: ABC → αβ0 (stationary reference frame).
 *
 *   [Vα]   [ 1  -1/2  -1/2 ] [Va]
 *   [Vβ] = [ 0  √3/2 -√3/2 ] [Vb] × 2/3
 *   [V0]   [1/2  1/2   1/2 ] [Vc]
 *
 * αβ components are orthogonal and stationary (not rotating).
 * Widely used in vector control of AC machines.
 *
 * @param abc       ABC phasors (instantaneous values)
 * @param alpha     [out] α-axis component
 * @param beta      [out] β-axis component
 * @param zero      [out] Zero-sequence component
 * @return 0 on success
 */
int phasor_clarke_transform(double a, double b, double c,
                             double *alpha, double *beta, double *zero);

/**
 * @brief Inverse Clarke transform: αβ0 → ABC.
 *
 * @param alpha  α-axis component
 * @param beta   β-axis component
 * @param zero   Zero-sequence component
 * @param a,b,c  [out] ABC components
 * @return 0 on success
 */
int phasor_inverse_clarke(double alpha, double beta, double zero,
                           double *a, double *b, double *c);

/**
 * @brief Park transform: ABC → dq0 (rotating reference frame).
 *
 * Converts stationary ABC quantities to a reference frame rotating
 * at angular frequency ω (angle θ = ωt + θ₀).
 *
 *   [Vd]   [ cos(θ)   cos(θ-2π/3)   cos(θ+2π/3) ] [Va]
 *   [Vq] = [-sin(θ)  -sin(θ-2π/3)  -sin(θ+2π/3) ] [Vb] × 2/3
 *   [V0]   [ 1/2       1/2           1/2         ] [Vc]
 *
 * In steady state balanced conditions, d and q are DC constants.
 *
 * @param a,b,c     Instantaneous ABC values
 * @param theta_rad Rotor/reference frame angle [radians]
 * @param dq        [out] dq0 quantities
 * @return 0 on success
 */
int phasor_park_transform(double a, double b, double c, double theta_rad,
                           dq0_t *dq);

/**
 * @brief Inverse Park transform: dq0 → ABC.
 *
 * @param dq         dq0 quantities
 * @param theta_rad  Reference frame angle [radians]
 * @param a,b,c      [out] ABC components
 * @return 0 on success
 */
int phasor_inverse_park(const dq0_t *dq, double theta_rad,
                         double *a, double *b, double *c);

/* ==========================================================================
 * L5: Phase-Locked Loop for Grid Synchronization
 * ========================================================================== */

/**
 * @brief Synchronous Reference Frame Phase-Locked Loop (SRF-PLL) state.
 *
 * The SRF-PLL is the standard method for grid synchronization in
 * three-phase converters. It uses the Park transform to align the
 * dq reference frame with the grid voltage vector.
 *
 * Control law: regulate Vq → 0 by adjusting ω, then integrate ω → θ.
 *   ω = ω_ff + (kp + ki/s) × Vq
 *   θ = ∫ ω dt
 *
 * Bandwidth tradeoff: higher BW → faster response but more noise
 *                     lower BW → better filtering but slower response
 *
 * Ref: Kaura & Blasko, "Operation of a Phase Locked Loop System
 *      Under Distorted Utility Conditions" (IEEE Trans. IA, 1997)
 */
typedef struct {
    double   kp;               /**< Proportional gain (typically 2·ζ·ωn / Vm) */
    double   ki;               /**< Integral gain (typically ωn² / Vm) */
    double   omega_ff;         /**< Nominal angular frequency feedforward [rad/s] */
    double   omega_est;        /**< Estimated frequency [rad/s] */
    double   theta;            /**< Estimated phase angle [rad] */
    double   vq_integral;      /**< PI integrator state */
    double   vd_filtered;      /**< Filtered d-axis voltage (magnitude estimate) */
    double   ts;               /**< Sampling time [seconds] */
} srf_pll_t;

/**
 * @brief Initialize an SRF-PLL.
 *
 * @param pll           [out] PLL state
 * @param omega_nom     Nominal grid frequency [rad/s] (e.g., 2π×60)
 * @param bandwidth_hz  Desired PLL bandwidth [Hz]
 * @param ts            Sampling time [seconds]
 * @return 0 on success
 */
int srf_pll_init(srf_pll_t *pll, double omega_nom, double bandwidth_hz,
                  double ts);

/**
 * @brief Run one iteration of the SRF-PLL.
 *
 * @param pll    PLL state (updated in place)
 * @param v_a    Phase A instantaneous voltage [V]
 * @param v_b    Phase B instantaneous voltage [V]
 * @param v_c    Phase C instantaneous voltage [V]
 * @return 0 on success
 */
int srf_pll_step(srf_pll_t *pll, double v_a, double v_b, double v_c);

/**
 * @brief Get the estimated frequency from the PLL.
 *
 * @param pll  PLL state
 * @return Estimated frequency [rad/s]
 */
double srf_pll_get_frequency(const srf_pll_t *pll);

/**
 * @brief Get the estimated phase angle from the PLL.
 *
 * @param pll  PLL state
 * @return Estimated angle [rad] (wrapped to [0, 2π])
 */
double srf_pll_get_angle(const srf_pll_t *pll);

/**
 * @brief Get the estimated voltage magnitude from the PLL.
 *
 * Under balanced conditions, Vd converges to the phase voltage magnitude
 * and Vq converges to zero.
 *
 * @param pll  PLL state
 * @return Estimated voltage magnitude [V]
 */
double srf_pll_get_magnitude(const srf_pll_t *pll);

#ifdef __cplusplus
}
#endif

#endif /* PHASOR_OPERATIONS_H */
