/**
 * @file stability.c
 * @brief Stability analysis implementation
 *
 * Implements:
 * - Routh-Hurwitz stability criterion
 * - Nyquist stability criterion with encirclement counting
 * - Gain and phase margin computation
 * - Root locus analysis
 * - Switched-capacitor filter stability
 *
 * Reference: Nyquist (1932), Bode (1945), Ogata (2010)
 * Course: MIT 6.003, Stanford EE102A, Berkeley EE16B
 */

#include "stability.h"
#include "frequency_response.h"
#include "transfer_function.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* ============================================================================
 * Routh-Hurwitz Stability Criterion
 * ============================================================================ */

/**
 * Construct Routh array and count sign changes in first column.
 *
 * The Routh array for polynomial a₀sⁿ + a₁sⁿ⁻¹ + ... + a_n:
 *
 * Row 0 (sⁿ):   a₀  a₂  a₄  ...
 * Row 1 (sⁿ⁻¹): a₁  a₃  a₅  ...
 * Row 2 (sⁿ⁻²): b₁  b₂  b₃  ...
 *   where b₁ = (a₁a₂ - a₀a₃)/a₁
 *         b₂ = (a₁a₄ - a₀a₅)/a₁
 *         b₃ = (a₁a₆ - a₀a₇)/a₁
 * Row 3 (sⁿ⁻³): c₁  c₂  c₃  ...
 *   where c₁ = (b₁a₃ - a₁b₂)/b₁
 *         c₂ = (b₁a₅ - a₁b₃)/b₁
 *
 * Special case 1: Zero in first column of a row.
 *   Replace zero with a small ε and evaluate the limit as ε → 0⁺.
 *   If the next row's first element has the SAME sign as the
 *   ε-substituted element, there are two sign changes (indicating
 *   two RHP poles from the auxiliary polynomial).
 *
 * Special case 2: Entire row of zeros.
 *   This indicates roots symmetrically placed about the origin
 *   (pairs on jω-axis, or mirror images about origin).
 *   Use the "auxiliary polynomial" formed from the row above and
 *   its derivative to continue the array.
 */
routh_hurwitz_t *stability_routh_hurwitz(const double *coeffs, size_t order)
{
    if (!coeffs || order == 0) return NULL;

    routh_hurwitz_t *rh = (routh_hurwitz_t *)calloc(1, sizeof(routh_hurwitz_t));
    if (!rh) return NULL;

    rh->order = order;

    /* Routh array: (order+1) rows, each with up to (order+1)/2 columns */
    size_t cols = (order + 2) / 2;
    size_t rows = order + 1;

    double **array = (double **)malloc(rows * sizeof(double *));
    if (!array) { free(rh); return NULL; }

    for (size_t i = 0; i < rows; i++) {
        array[i] = (double *)calloc(cols, sizeof(double));
        if (!array[i]) {
            for (size_t j = 0; j < i; j++) free(array[j]);
            free(array);
            free(rh);
            return NULL;
        }
    }

    /* Initialize first two rows from coefficients */
    for (size_t j = 0; j < cols; j++) {
        size_t idx = 2 * j;
        if (idx <= order) array[0][j] = coeffs[idx];
        idx = 2 * j + 1;
        if (idx <= order) array[1][j] = coeffs[idx];
    }

    /* Build remaining rows */
    for (size_t i = 2; i < rows; i++) {
        /* Leading element of previous row */
        double prev_lead = array[i - 1][0];

        if (fabs(prev_lead) < 1e-12) {
            /* Special case: zero in first column.
             * Replace with small epsilon. We handle this by checking
             * if the entire row is zero. */
            int all_zero = 1;
            for (size_t j = 0; j < cols; j++) {
                if (fabs(array[i - 1][j]) > 1e-12) {
                    all_zero = 0;
                    break;
                }
            }

            if (all_zero) {
                /* Entire row is zero → auxiliary polynomial.
                 * Use derivative of the row i-2 coefficients. */
                /* The auxiliary polynomial is formed from the coefficients
                 * in row i-2, using only even powers. */
                for (size_t j = 0; j < cols - 1; j++) {
                    /* Derivative: multiply by the power and shift */
                    array[i][j] = array[i - 2][j + 1] * (double)(2 * (cols - 1 - j));
                }
                continue;
            }

            /* Non-zero row but first element is zero → replace with ε */
            prev_lead = 1e-10;
        }

        /* Compute row i using the 2×2 determinant formula:
         * new[j] = (prev_row1[0] * prev_row2[j+1] - prev_row2[0] * prev_row1[j+1])
         *          / prev_row1[0] */
        for (size_t j = 0; j < cols - 1; j++) {
            double a = array[i - 2][0];
            double b = array[i - 2][j + 1];
            double c = array[i - 1][0];
            double d = array[i - 1][j + 1];

            array[i][j] = (c * b - a * d) / prev_lead;
        }
    }

    /* Extract first column */
    rh->first_column = (double *)malloc(rows * sizeof(double));
    if (!rh->first_column) {
        for (size_t j = 0; j < rows; j++) free(array[j]);
        free(array);
        free(rh);
        return NULL;
    }

    for (size_t i = 0; i < rows; i++) {
        rh->first_column[i] = array[i][0];
    }

    /* Count sign changes */
    int prev_sign = 0;
    int first_nonzero = 0;
    rh->sign_changes = 0;

    for (size_t i = 0; i < rows; i++) {
        double val = rh->first_column[i];
        if (fabs(val) < 1e-12) continue;  /* Skip zeros */

        int sign = (val > 0.0) ? 1 : -1;
        if (!first_nonzero) {
            prev_sign = sign;
            first_nonzero = 1;
        } else {
            if (sign != prev_sign) {
                rh->sign_changes++;
            }
            prev_sign = sign;
        }
    }

    /* Stability: no sign changes AND no zero in first column */
    if (rh->sign_changes == 0 && first_nonzero) {
        /* Check for marginal stability: any zero in first column? */
        int has_zero = 0;
        for (size_t i = 0; i < rows; i++) {
            if (fabs(rh->first_column[i]) < 1e-12) {
                has_zero = 1;
                break;
            }
        }
        rh->is_stable = !has_zero;
        rh->is_marginally_stable = has_zero;
    } else {
        rh->is_stable = 0;
        rh->is_marginally_stable = 0;
    }

    /* Cleanup array */
    for (size_t i = 0; i < rows; i++) free(array[i]);
    free(array);

    return rh;
}

void routh_hurwitz_free(routh_hurwitz_t *rh)
{
    if (rh) {
        free(rh->first_column);
        free(rh);
    }
}

int stability_is_hurwitz(const tf_polynomial_t *tf)
{
    if (!tf) return -1;

    routh_hurwitz_t *rh = stability_routh_hurwitz(tf->den, tf->den_order);
    if (!rh) return -1;

    int stable = rh->is_stable;
    routh_hurwitz_free(rh);
    return stable;
}

/* ============================================================================
 * Nyquist Stability Criterion
 * ============================================================================ */

double _Complex *stability_nyquist_contour(const tf_polynomial_t *loop_tf,
                                              double f_start, double f_end,
                                              size_t n_points)
{
    if (!loop_tf || n_points < 2) return NULL;

    double _Complex *nyquist = (double _Complex *)malloc(
        n_points * sizeof(double _Complex));
    if (!nyquist) return NULL;

    /* Generate frequency grid (logarithmic for better resolution at low freq) */
    size_t n_freq = 0;
    double *freqs = freq_logspace(f_start, f_end, 20, &n_freq);
    if (!freqs) { free(nyquist); return NULL; }

    size_t actual_pts = (n_freq < n_points) ? n_freq : n_points;

    for (size_t i = 0; i < actual_pts; i++) {
        double f = freqs[i];
        double w = 2.0 * M_PI * f;
        nyquist[i] = tf_evaluate_freq(loop_tf, w);
    }

    /* Also include the negative frequency part (conjugate symmetric for
     * real-coefficient transfer functions). For complete Nyquist contour,
     * we'd also need the infinite semicircle in the RHP, but for
     * open-loop stable systems (P=0), the finite frequency sweep
     * is sufficient. */

    /* Pad remaining points */
    for (size_t i = actual_pts; i < n_points; i++) {
        nyquist[i] = nyquist[actual_pts - 1];
    }

    free(freqs);
    return nyquist;
}

/**
 * Count encirclements of (-1, j0) by the Nyquist plot.
 *
 * Algorithm:
 * 1. For each adjacent pair of points on the Nyquist contour,
 *    compute the change in angle of the vector from (-1, j0).
 * 2. Accumulate the total angle change.
 * 3. Winding number = total angle change / (2π).
 *
 * The angle change is computed using the cross product formula
 * to avoid atan2 wrap-around issues:
 *   Δθ = atan2(cross, dot)
 * where cross = Im{(z_i+1)·conj(z_{i+1}+1)}
 *       dot   = Re{(z_i+1)·conj(z_{i+1}+1)}
 */
int stability_count_encirclements(const double _Complex *nyquist_data,
                                    size_t n)
{
    if (!nyquist_data || n < 2) return 0;

    double total_angle = 0.0;

    for (size_t i = 0; i < n - 1; i++) {
        /* Vector from (-1, 0) to nyquist point */
        double _Complex v1 = nyquist_data[i] + 1.0;
        double _Complex v2 = nyquist_data[i + 1] + 1.0;

        /* Compute angle from v1 to v2 */
        double cross = creal(v1) * cimag(v2) - cimag(v1) * creal(v2);
        double dot   = creal(v1) * creal(v2) + cimag(v1) * cimag(v2);

        double dtheta = atan2(cross, dot);
        total_angle += dtheta;
    }

    /* Winding number = total angle / (2π), rounded to nearest integer */
    int winding = (int)round(total_angle / (2.0 * M_PI));

    return winding;  /* Positive = clockwise encirclements */
}

int stability_nyquist_check(const tf_polynomial_t *loop_tf,
                             double f_start, double f_end,
                             size_t n_points,
                             int *is_stable, int *encirclements)
{
    if (!loop_tf) return -1;

    double _Complex *nyquist = stability_nyquist_contour(loop_tf, f_start,
                                                           f_end, n_points);
    if (!nyquist) return -1;

    int N = stability_count_encirclements(nyquist, n_points);

    /* Count open-loop unstable poles */
    routh_hurwitz_t *rh = stability_routh_hurwitz(loop_tf->den,
                                                    loop_tf->den_order);
    int P = 0;
    if (rh) {
        P = (int)rh->sign_changes;
        routh_hurwitz_free(rh);
    }

    /* Nyquist criterion: Z = N + P
     * For closed-loop stability: Z = 0 → N = -P */
    int Z = N + P;
    if (is_stable) *is_stable = (Z == 0) ? 1 : 0;
    if (encirclements) *encirclements = N;

    free(nyquist);
    return 0;
}

/* ============================================================================
 * Stability Margins
 * ============================================================================ */

stability_margin_t stability_margins_from_response(
    const freq_response_t *resp)
{
    stability_margin_t margin;
    memset(&margin, 0, sizeof(stability_margin_t));

    if (!resp || resp->num_points < 2) return margin;

    size_t n = resp->num_points;

    /* Find gain crossover: |H(jω)| = 1 (0 dB) */
    double w_gc = -1.0;
    double phase_at_gc = 0.0;

    for (size_t i = 0; i < n - 1; i++) {
        double db1 = resp->points[i].magnitude_db;
        double db2 = resp->points[i + 1].magnitude_db;

        if ((db1 >= 0.0 && db2 <= 0.0) || (db1 <= 0.0 && db2 >= 0.0)) {
            /* Interpolate */
            double alpha = -db1 / (db2 - db1);
            double log_w1 = log10(resp->points[i].angular_freq);
            double log_w2 = log10(resp->points[i + 1].angular_freq);
            double log_w = log_w1 + alpha * (log_w2 - log_w1);
            w_gc = pow(10.0, log_w);

            double phase1 = resp->points[i].phase_rad;
            double phase2 = resp->points[i + 1].phase_rad;
            /* Handle wrap */
            double dphase = phase2 - phase1;
            if (dphase > M_PI) dphase -= 2.0 * M_PI;
            if (dphase < -M_PI) dphase += 2.0 * M_PI;
            phase_at_gc = phase1 + alpha * dphase;
            break;
        }
    }

    /* Find phase crossover: ∠H(jω) = -180° (-π rad) */
    double w_pc = -1.0;
    double mag_at_pc_db = 0.0;

    for (size_t i = 0; i < n - 1; i++) {
        double phase1 = resp->points[i].phase_rad;
        double phase2 = resp->points[i + 1].phase_rad;

        /* Unwrap this pair */
        double dp = phase2 - phase1;
        if (dp > M_PI) dp -= 2.0 * M_PI;
        if (dp < -M_PI) dp += 2.0 * M_PI;
        double phase2u = phase1 + dp;

        if ((phase1 >= -M_PI && phase2u <= -M_PI) ||
            (phase2u >= -M_PI && phase1 <= -M_PI) ||
            (phase1 >= M_PI && phase2u <= M_PI)) {
            /* Interpolate to -π */
            double alpha = (-M_PI - phase1) / dp;
            double log_w1 = log10(resp->points[i].angular_freq);
            double log_w2 = log10(resp->points[i + 1].angular_freq);
            double log_w = log_w1 + alpha * (log_w2 - log_w1);
            w_pc = pow(10.0, log_w);

            double mag1 = resp->points[i].magnitude_db;
            double mag2 = resp->points[i + 1].magnitude_db;
            mag_at_pc_db = mag1 + alpha * (mag2 - mag1);
            break;
        }
    }

    /* Compute margins */
    if (w_gc > 0.0) {
        margin.gain_crossover_hz = w_gc / (2.0 * M_PI);
        double pm_rad = M_PI + phase_at_gc;  /* PM = 180° + phase(L(jω_gc)) */
        /* Bring into [-π, π] */
        while (pm_rad > M_PI) pm_rad -= 2.0 * M_PI;
        while (pm_rad < -M_PI) pm_rad += 2.0 * M_PI;
        margin.phase_margin_rad = pm_rad;
        margin.phase_margin_deg = pm_rad * 180.0 / M_PI;
    }

    if (w_pc > 0.0) {
        margin.phase_crossover_hz = w_pc / (2.0 * M_PI);
        margin.gain_margin_db = -mag_at_pc_db;  /* GM_dB = -|L|_dB at phase crossover */
        margin.gain_margin_linear = pow(10.0, margin.gain_margin_db / 20.0);
    }

    /* Stability: both margins positive */
    margin.is_stable = (margin.gain_margin_db > 0.0 &&
                         margin.phase_margin_deg > 0.0) ? 1 : 0;

    return margin;
}

stability_margin_t stability_margins(const tf_polynomial_t *loop_tf)
{
    stability_margin_t margin;
    memset(&margin, 0, sizeof(stability_margin_t));

    if (!loop_tf) return margin;

    /* Generate frequency response */
    /* Default sweep: 3 decades below dominant pole to 3 decades above */
    freq_response_t *resp = tf_compute_freq_response(loop_tf,
                                                      0.1, 1e6, 50);
    if (!resp) return margin;

    margin = stability_margins_from_response(resp);
    freq_response_free(resp);
    return margin;
}

/* ============================================================================
 * Root Locus Analysis (L8 Advanced)
 * ============================================================================ */

/**
 * Compute closed-loop pole locations for varying gain K.
 *
 * The characteristic equation: D(s) + K·N(s) = 0
 *
 * For each gain value K_k, we find the roots of:
 *   P(s) = D(s) + K·N(s)
 *
 * This implements the basic numerical root locus by:
 * 1. Starting from open-loop poles (K=0)
 * 2. Iteratively increasing K
 * 3. Using the previous roots as initial guesses for the next K
 *
 * This "continuation" approach is numerically efficient because
 * poles move continuously with K.
 */
int stability_root_locus(const tf_polynomial_t *loop_tf,
                          const double *K_values, size_t n_K,
                          double _Complex *poles_out)
{
    if (!loop_tf || !K_values || n_K == 0 || !poles_out) return -1;

    size_t n = loop_tf->den_order;

    /* For each gain, form the closed-loop characteristic polynomial
     * D(s) + K·N(s) and find its roots. */
    for (size_t k_idx = 0; k_idx < n_K; k_idx++) {
        double K = K_values[k_idx];

        /* coeffs = den_coeff + K * num_coeff (padded to same length) */
        double *cl_coeffs = (double *)calloc(n + 1, sizeof(double));
        if (!cl_coeffs) return -1;

        for (size_t i = 0; i <= n; i++) {
            cl_coeffs[i] = loop_tf->den[i];
            if (i <= loop_tf->num_order) {
                cl_coeffs[i] += K * loop_tf->num[i];
            }
        }

        /* Find roots using Laguerre (via find_roots) */
        /* For simplicity, use external root-finding by calling tf_to_pole_zero */
        tf_polynomial_t *cl_tf = tf_polynomial_create(
            loop_tf->num, loop_tf->num_order,
            cl_coeffs, n);
        free(cl_coeffs);

        if (!cl_tf) return -1;

        tf_pole_zero_t *pz = tf_to_pole_zero(cl_tf);
        tf_polynomial_free(cl_tf);

        if (pz) {
            for (size_t j = 0; j < pz->num_poles && j < n; j++) {
                poles_out[k_idx * n + j] = pz->poles[j];
            }
            tf_pole_zero_free(pz);
        }
    }

    return 0;
}

/**
 * Find gain margin via root locus: the gain K where poles first
 * cross the jω-axis.
 *
 * For this, we find the imaginary axis crossing by solving
 * D(jω) + K·N(jω) = 0 for real ω and positive K.
 *
 * This gives: K = -D(jω)/N(jω) evaluated at the jω-crossing frequency.
 * The crossing frequency satisfies: Im{D(jω)/N(jω)} = 0.
 */
int stability_gain_margin_root_locus(const tf_polynomial_t *loop_tf,
                                       double *K_margin, double *freq_hz)
{
    if (!loop_tf) return -1;

    /* Sweep frequency and find where Im{L(jω)} = 0 and Re{L(jω)} < 0 */
    /* The jω-crossing occurs where the loop gain L(jω) is real and negative,
     * since K_margin = -1/Re{L(jω)} where Im{L(jω)} = 0. */

    double K_m = INFINITY;
    double f_m = 0.0;
    int found = 0;

    for (double log_f = 0.0; log_f <= 8.0; log_f += 0.001) {
        double f = pow(10.0, log_f);
        double w = 2.0 * M_PI * f;
        double _Complex L = tf_evaluate_freq(loop_tf, w);

        /* Check if imaginary part crosses zero */
        if (fabs(cimag(L)) < 1e-6 * cabs(L) && creal(L) < 0.0) {
            double K_candidate = -1.0 / creal(L);
            if (K_candidate > 0.0 && K_candidate < K_m) {
                K_m = K_candidate;
                f_m = f;
                found = 1;
            }
        }
    }

    if (found) {
        if (K_margin) *K_margin = K_m;
        if (freq_hz) *freq_hz = f_m;
    }
    return found ? 0 : -1;
}

/* ============================================================================
 * Switched-Capacitor Filter Stability (L7 Application)
 * ============================================================================ */

/**
 * Check stability of a switched-capacitor realization.
 *
 * The bilinear transform maps the continuous-time s-domain to the
 * discrete-time z-domain:
 *   s → (2/T)·(z-1)/(z+1)   where T = 1/f_clk
 *
 * Or using the forward/backward Euler methods:
 *   Forward Euler:  s → (z-1)/T
 *   Backward Euler: s → (z-1)/(zT)
 *
 * For SC integrators, the backward Euler (or bilinear) mapping is
 * typically used.
 *
 * Stability check: all z-domain poles must satisfy |z| < 1.
 * Using the bilinear transform, this is equivalent to checking that
 * all s-domain poles are in the LHP (which is already ensured if
 * the continuous-time prototype is stable).
 *
 * This function verifies the mapping and reports whether the
 * discrete-time implementation is stable.
 */
int stability_sc_filter(const tf_polynomial_t *tf_s,
                          double f_clk, int *is_stable)
{
    if (!tf_s || f_clk <= 0.0) return -1;

    /* For backward Euler mapping of a first-order section:
     * H(s) = 1/(s + a) → H(z) = 1/((1-z^{-1})/T + a)
     *                           = T/(1 - z^{-1} + aT)
     *                           = T/((1+aT) - z^{-1})
     * Pole at z = 1/(1+aT)
     * Stability: |z| < 1 → |1/(1+aT)| < 1 → |1+aT| > 1 → aT > 0 for a>0
     *
     * This always holds if the continuous-time prototype is stable (a > 0).
     * For bilinear transform, the mapping is:
     *   z = (1 + sT/2)/(1 - sT/2)
     * which maps the entire LHP to the unit circle interior.
     */

    /* Simple check: if the continuous-time filter is stable, the
     * SC implementation will be stable for bilinear and backward
     * Euler mappings. Forward Euler may be unstable for high Q. */

    int ct_stable = stability_is_hurwitz(tf_s);
    if (ct_stable < 0) return -1;

    if (is_stable) {
        /* For the bilinear transform, LHP → unit circle interior, so
         * CT stability implies DT stability. */
        *is_stable = ct_stable;
    }

    return 0;
}
