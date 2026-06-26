/**
 * @file phasor_operations.c
 * @brief Phasor arithmetic, Fortescue symmetrical components, Clarke/Park
 *        transforms, and SRF-PLL grid synchronization
 *
 * Phasors reduce sinusoidal steady-state AC analysis to complex algebra.
 * The Fortescue transform decomposes unbalanced three-phase systems into
 * balanced sequence components. The Park transform enables vector control
 * of AC machines by converting stationary AC to rotating DC quantities.
 *
 * Knowledge points:
 *   L1: Phasor construction (polar, rectangular, peak-to-RMS)
 *   L2: Phasor arithmetic (add, sub, mul, div, scale)
 *   L3: Fortescue symmetrical components transform (1918)
 *   L3: Inverse Fortescue transform
 *   L3: Clarke transform (αβ0 stationary frame)
 *   L3: Park transform (dq0 rotating frame)
 *   L5: SRF-PLL for three-phase grid synchronization
 *   L5: Complex power from sequence components
 *
 * References:
 *   - C.L. Fortescue, "Method of Symmetrical Co-ordinates" (AIEE, 1918)
 *   - R.H. Park, "Two-Reaction Theory of Synchronous Machines" (AIEE, 1929)
 *   - E. Clarke, "Circuit Analysis of AC Power Systems" (1943)
 *   - Kaura & Blasko, "Operation of a PLL Under Distorted Utility
 *     Conditions" (IEEE Trans. IA, 1997)
 */

#include "phasor_operations.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define PH_EPS 1e-12
#define SQRT3  1.7320508075688772

/* ==========================================================================
 * L1: Phasor Construction
 * ========================================================================== */

phasor_t phasor_from_polar(double magnitude, double angle_rad)
{
    phasor_t p;
    p.magnitude  = magnitude;
    p.angle_rad  = angle_rad;
    p.angle_deg  = angle_rad * 180.0 / M_PI;
    p.rect       = magnitude * (cos(angle_rad) + I * sin(angle_rad));
    return p;
}

phasor_t phasor_from_rect(double re, double im)
{
    phasor_t p;
    p.rect       = re + I * im;
    p.magnitude  = sqrt(re * re + im * im);
    p.angle_rad  = atan2(im, re);
    p.angle_deg  = p.angle_rad * 180.0 / M_PI;
    return p;
}

phasor_t phasor_from_peak(double v_peak, double angle_deg)
{
    double v_rms = v_peak / 1.4142135623730951;
    double ang_rad = angle_deg * M_PI / 180.0;
    return phasor_from_polar(v_rms, ang_rad);
}

phasor_t phasor_conjugate(const phasor_t *p)
{
    if (!p) { phasor_t zero = {0}; return zero; }
    phasor_t result;
    result.magnitude = p->magnitude;
    result.angle_rad = -p->angle_rad;
    result.angle_deg = -p->angle_deg;
    result.rect = conj(p->rect);
    return result;
}

/* ==========================================================================
 * L2: Phasor Arithmetic
 * ========================================================================== */

phasor_t phasor_add(const phasor_t *a, const phasor_t *b)
{
    if (!a || !b) { phasor_t zero = {0}; return zero; }
    double complex sum = a->rect + b->rect;
    return phasor_from_rect(creal(sum), cimag(sum));
}

phasor_t phasor_sub(const phasor_t *a, const phasor_t *b)
{
    if (!a || !b) { phasor_t zero = {0}; return zero; }
    double complex diff = a->rect - b->rect;
    return phasor_from_rect(creal(diff), cimag(diff));
}

phasor_t phasor_mul(const phasor_t *a, const phasor_t *b)
{
    if (!a || !b) { phasor_t zero = {0}; return zero; }
    /* Multiply: magnitudes multiply, angles add */
    return phasor_from_polar(a->magnitude * b->magnitude,
                             a->angle_rad + b->angle_rad);
}

phasor_t phasor_div(const phasor_t *a, const phasor_t *b)
{
    if (!a || !b) { phasor_t zero = {0}; return zero; }
    if (b->magnitude < PH_EPS) {
        phasor_t zero = {0};
        return zero;
    }
    /* Divide: magnitudes divide, angles subtract */
    return phasor_from_polar(a->magnitude / b->magnitude,
                             a->angle_rad - b->angle_rad);
}

phasor_t phasor_scale(const phasor_t *p, double scalar)
{
    if (!p) { phasor_t zero = {0}; return zero; }
    return phasor_from_polar(p->magnitude * scalar, p->angle_rad);
}

/* ==========================================================================
 * L3: Fortescue Transform — ABC to 012
 * ========================================================================== */

int phasor_abc_to_012(const phasor_abc_t *abc, phasor_012_t *zpn)
{
    if (!abc || !zpn) return -1;

    double complex va = abc->phase_a.rect;
    double complex vb = abc->phase_b.rect;
    double complex vc = abc->phase_c.rect;

    /* Fortescue operator a = e^{j·2π/3} = -0.5 + j·0.866 */
    double complex a  = -0.5 + I * 0.8660254037844386;
    double complex a2 = -0.5 - I * 0.8660254037844386; /* a² = a* */

    double complex v0 = (va + vb + vc) / 3.0;
    double complex v1 = (va + a * vb + a2 * vc) / 3.0;
    double complex v2 = (va + a2 * vb + a * vc) / 3.0;

    zpn->zero_seq = phasor_from_rect(creal(v0), cimag(v0));
    zpn->pos_seq  = phasor_from_rect(creal(v1), cimag(v1));
    zpn->neg_seq  = phasor_from_rect(creal(v2), cimag(v2));

    /* Unbalance: |V2|/|V1| × 100% */
    if (zpn->pos_seq.magnitude > PH_EPS) {
        zpn->unbalance_percent = zpn->neg_seq.magnitude / zpn->pos_seq.magnitude
                                  * 100.0;
    } else {
        zpn->unbalance_percent = 0.0;
    }

    /* Zero-sequence ratio */
    if (zpn->pos_seq.magnitude > PH_EPS) {
        zpn->zero_seq_percent = zpn->zero_seq.magnitude / zpn->pos_seq.magnitude
                                 * 100.0;
    } else {
        zpn->zero_seq_percent = 0.0;
    }

    return 0;
}

/* ==========================================================================
 * L3: Inverse Fortescue Transform — 012 to ABC
 * ========================================================================== */

int phasor_012_to_abc(const phasor_012_t *zpn, phasor_abc_t *abc)
{
    if (!zpn || !abc) return -1;

    double complex v0 = zpn->zero_seq.rect;
    double complex v1 = zpn->pos_seq.rect;
    double complex v2 = zpn->neg_seq.rect;

    double complex a  = -0.5 + I * 0.8660254037844386;
    double complex a2 = -0.5 - I * 0.8660254037844386;

    double complex va = v0 + v1 + v2;
    double complex vb = v0 + a2 * v1 + a * v2;
    double complex vc = v0 + a * v1 + a2 * v2;

    abc->phase_a = phasor_from_rect(creal(va), cimag(va));
    abc->phase_b = phasor_from_rect(creal(vb), cimag(vb));
    abc->phase_c = phasor_from_rect(creal(vc), cimag(vc));

    /* Check balance */
    double tol = 0.02; /* 2% tolerance for magnitude match */
    double avg_mag = (abc->phase_a.magnitude + abc->phase_b.magnitude
                       + abc->phase_c.magnitude) / 3.0;
    double dev_a = fabs(abc->phase_a.magnitude - avg_mag) / avg_mag;
    double dev_b = fabs(abc->phase_b.magnitude - avg_mag) / avg_mag;
    double dev_c = fabs(abc->phase_c.magnitude - avg_mag) / avg_mag;

    abc->is_balanced = (dev_a < tol && dev_b < tol && dev_c < tol) ? 1 : 0;
    abc->is_positive_seq = (zpn->neg_seq.magnitude < zpn->pos_seq.magnitude
                             * 0.05) ? 1 : 0;

    return 0;
}

/* ==========================================================================
 * L3: Complex Power from Symmetrical Components
 * ========================================================================== */

int phasor_power_from_components(const phasor_012_t *v_zpn,
                                  const phasor_012_t *i_zpn,
                                  double complex *s_total)
{
    if (!v_zpn || !i_zpn || !s_total) return -1;

    /* Total complex power from sequence components:
     * S_3phase = 3 × (V1·I1* + V2·I2* + V0·I0*)
     *
     * The factor 3 accounts for all three phases. Each sequence
     * network contributes independently (they are orthogonal).
     */

    double complex s1 = v_zpn->pos_seq.rect * conj(i_zpn->pos_seq.rect);
    double complex s2 = v_zpn->neg_seq.rect * conj(i_zpn->neg_seq.rect);
    double complex s0 = v_zpn->zero_seq.rect * conj(i_zpn->zero_seq.rect);

    *s_total = 3.0 * (s1 + s2 + s0);
    return 0;
}

/* ==========================================================================
 * L3: Clarke Transform — ABC → αβ0
 * ========================================================================== */

int phasor_clarke_transform(double a, double b, double c,
                             double *alpha, double *beta, double *zero)
{
    if (!alpha || !beta || !zero) return -1;

    /* Power-invariant Clarke transform:
     *   [α]   [ 1   -1/2   -1/2 ] [a]
     *   [β] = [ 0   √3/2  -√3/2 ] [b] × 2/3
     *   [0]   [1/2   1/2    1/2 ] [c]
     *
     * α is aligned with phase A.
     * β is orthogonal (90° lagging α).
     */

    *alpha = (2.0/3.0) * (a - 0.5 * b - 0.5 * c);
    *beta  = (2.0/3.0) * (0.8660254037844386 * b - 0.8660254037844386 * c);
    *zero  = (1.0/3.0) * (a + b + c);

    return 0;
}

/* ==========================================================================
 * L3: Inverse Clarke Transform — αβ0 → ABC
 * ========================================================================== */

int phasor_inverse_clarke(double alpha, double beta, double zero,
                           double *a, double *b, double *c)
{
    if (!a || !b || !c) return -1;

    /* Inverse power-invariant Clarke:
     *   [a]   [ 1       0      1 ] [α]
     *   [b] = [-1/2   √3/2    1 ] [β]
     *   [c]   [-1/2  -√3/2    1 ] [0]
     */
    *a = alpha + zero;
    *b = -0.5 * alpha + 0.8660254037844386 * beta + zero;
    *c = -0.5 * alpha - 0.8660254037844386 * beta + zero;

    return 0;
}

/* ==========================================================================
 * L3: Park Transform — ABC → dq0
 * ========================================================================== */

int phasor_park_transform(double a, double b, double c, double theta_rad,
                           dq0_t *dq)
{
    if (!dq) return -1;

    double cos_th = cos(theta_rad);
    double sin_th = sin(theta_rad);
    double cos_120 = cos(theta_rad - 2.0 * M_PI / 3.0);
    double sin_120 = -sin(theta_rad - 2.0 * M_PI / 3.0);
    double cos_240 = cos(theta_rad + 2.0 * M_PI / 3.0);
    double sin_240 = -sin(theta_rad + 2.0 * M_PI / 3.0);

    /* Standard Park transform:
     *   [Vd]   [ cos(θ)   cos(θ-120°)   cos(θ+120°) ] [Va]
     *   [Vq] = [-sin(θ)  -sin(θ-120°)  -sin(θ+120°) ] [Vb] × 2/3
     *   [V0]   [ 1/2       1/2           1/2         ] [Vc]
     */
    dq->d = (2.0/3.0) * (cos_th * a + cos_120 * b + cos_240 * c);
    dq->q = (2.0/3.0) * (sin_th * a + sin_120 * b + sin_240 * c);
    dq->zero = (a + b + c) / 3.0;

    dq->magnitude = sqrt(dq->d * dq->d + dq->q * dq->q);
    dq->angle_rad = atan2(dq->q, dq->d);

    return 0;
}

/* ==========================================================================
 * L3: Inverse Park Transform — dq0 → ABC
 * ========================================================================== */

int phasor_inverse_park(const dq0_t *dq, double theta_rad,
                         double *a, double *b, double *c)
{
    if (!dq || !a || !b || !c) return -1;

    double cos_th = cos(theta_rad);
    double sin_th = sin(theta_rad);

    /* Inverse Park:
     *   [Va]   [ cos(θ)  -sin(θ)  1 ] [Vd]
     *   [Vb] = [cos(θ-120°)  -sin(θ-120°)  1 ] [Vq]
     *   [Vc]   [cos(θ+120°)  -sin(θ+120°)  1 ] [V0]
     */
    *a = cos_th * dq->d - sin_th * dq->q + dq->zero;
    *b = cos(theta_rad - 2.0 * M_PI / 3.0) * dq->d
          - sin(theta_rad - 2.0 * M_PI / 3.0) * dq->q + dq->zero;
    *c = cos(theta_rad + 2.0 * M_PI / 3.0) * dq->d
          - sin(theta_rad + 2.0 * M_PI / 3.0) * dq->q + dq->zero;

    return 0;
}

/* ==========================================================================
 * L5: SRF-PLL Initialization
 * ========================================================================== */

int srf_pll_init(srf_pll_t *pll, double omega_nom, double bandwidth_hz,
                  double ts)
{
    if (!pll || omega_nom <= 0.0 || bandwidth_hz <= 0.0 || ts <= 0.0)
        return -1;

    memset(pll, 0, sizeof(*pll));

    /* SRF-PLL design (Kaura & Blasko, 1997):
     *
     * The closed-loop transfer function from phase error to phase estimate
     * is a second-order system:
     *   H(s) = (Kp + Ki/s) × Vm / (s + (Kp + Ki/s) × Vm)
     *        = (Vm·Kp·s + Vm·Ki) / (s² + Vm·Kp·s + Vm·Ki)
     *
     * With desired natural frequency ωn = 2π·BW and damping ζ = 1/√2:
     *   Kp = 2·ζ·ωn / Vm  ≈ 2·0.707·ωn / Vm
     *   Ki = ωn² / Vm
     *
     * where Vm is the expected phase voltage magnitude (we estimate 1.0 pu,
     * actual Vm is measured online from Vd).
     */

    double v_m_est = 1.0; /* normalized: actual Vm measured by PLL */
    double omega_n = 2.0 * M_PI * bandwidth_hz;
    double zeta = 0.70710678; /* 1/√2 = optimal damping */

    pll->kp = 2.0 * zeta * omega_n / v_m_est;
    pll->ki = omega_n * omega_n / v_m_est;
    pll->omega_ff = omega_nom;
    pll->omega_est = omega_nom;
    pll->theta = 0.0;
    pll->vq_integral = 0.0;
    pll->vd_filtered = v_m_est;
    pll->ts = ts;

    return 0;
}

/* ==========================================================================
 * L5: SRF-PLL Iteration Step
 * ========================================================================== */

int srf_pll_step(srf_pll_t *pll, double v_a, double v_b, double v_c)
{
    if (!pll) return -1;

    /* 1. Clarke transform */
    double alpha, beta, zero;
    phasor_clarke_transform(v_a, v_b, v_c, &alpha, &beta, &zero);

    /* 2. Park transform using estimated angle */
    double cos_th = cos(pll->theta);
    double sin_th = sin(pll->theta);
    double vd =  cos_th * alpha + sin_th * beta;
    double vq = -sin_th * alpha + cos_th * beta;

    /* 3. Filter Vd (voltage magnitude estimate) */
    double alpha_vd = 0.01; /* low-pass filter gain */
    pll->vd_filtered = (1.0 - alpha_vd) * pll->vd_filtered + alpha_vd * vd;

    /* 4. PI controller on Vq
     *    Regulate Vq → 0 to align d-axis with voltage vector.
     *    If Vq > 0: frequency is too low → increase ω
     *    If Vq < 0: frequency is too high → decrease ω
     */
    double vq_norm = (pll->vd_filtered > PH_EPS)
                      ? vq / pll->vd_filtered : vq;

    /* PI */
    double omega_correction = pll->kp * vq_norm + pll->ki * pll->vq_integral;
    pll->vq_integral += vq_norm * pll->ts;

    /* Anti-windup */
    double max_correction = 2.0 * M_PI * 5.0; /* ±5 Hz */
    if (pll->vq_integral > max_correction / pll->ki)
        pll->vq_integral = max_correction / pll->ki;
    if (pll->vq_integral < -max_correction / pll->ki)
        pll->vq_integral = -max_correction / pll->ki;

    /* 5. Update frequency */
    pll->omega_est = pll->omega_ff + omega_correction;

    /* 6. Integrate frequency to get angle */
    pll->theta += pll->omega_est * pll->ts;

    /* Wrap to [0, 2π] */
    while (pll->theta >= 2.0 * M_PI) pll->theta -= 2.0 * M_PI;
    while (pll->theta < 0.0) pll->theta += 2.0 * M_PI;

    return 0;
}

/* ==========================================================================
 * L5: SRF-PLL Accessors
 * ========================================================================== */

double srf_pll_get_frequency(const srf_pll_t *pll)
{
    if (!pll) return -1.0;
    return pll->omega_est;
}

double srf_pll_get_angle(const srf_pll_t *pll)
{
    if (!pll) return -1.0;
    return pll->theta;
}

double srf_pll_get_magnitude(const srf_pll_t *pll)
{
    if (!pll) return -1.0;
    return pll->vd_filtered;
}
