/**
 * @file example_stability_analysis.c
 * @brief Example: Stability Analysis of Feedback Amplifier
 *
 * This example demonstrates:
 * 1. Open-loop transfer function stability assessment
 * 2. Routh-Hurwitz criterion on closed-loop poles
 * 3. Nyquist stability criterion
 * 4. Gain margin and phase margin computation
 * 5. Effect of feedback on bandwidth (gain-bandwidth trade-off)
 * 6. Root locus: how poles migrate with increasing gain
 *
 * L6 Canonical Problem: Stability analysis of a feedback amplifier is
 * essential to ensure the circuit won't oscillate. The Bode criterion
 * (checking phase when gain crosses 0 dB) is the most widely used
 * method in practice.
 *
 * Reference: Bode (1945), Ogata (2010) Ch. 8, Sedra & Smith (2020) Ch. 10
 * Course: MIT 6.003, Stanford EE102A, Berkeley EE105
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "frequency_response.h"
#include "transfer_function.h"
#include "bode_plot.h"
#include "stability.h"

int main(void)
{
    printf("╔══════════════════════════════════════════╗\n");
    printf("║  Feedback Amplifier Stability Analysis  ║\n");
    printf("╚══════════════════════════════════════════╝\n\n");

    /* ================================================================
     * Part 1: Third-order amplifier open-loop gain
     *
     * A(s) = A₀ / ((1 + s/ω₁)(1 + s/ω₂)(1 + s/ω₃))
     *
     * A₀ = 1000 (60 dB), ω₁ = 10⁴ rad/s, ω₂ = 10⁵ rad/s, ω₃ = 10⁶ rad/s
     * This represents a three-stage amplifier with dominant-pole
     * compensation at ω₁.
     *
     * Loop gain with feedback β = 0.1 (closed-loop gain ≈ 10):
     * L(s) = β·A(s) = 0.1 × 1000 / ((1+s/10⁴)(1+s/10⁵)(1+s/10⁶))
     *       = 100 / ((1+s/10⁴)(1+s/10⁵)(1+s/10⁶))
     * ================================================================ */

    printf("Amplifier Model:\n");
    printf("  DC gain A₀ = 1000 (60 dB)\n");
    printf("  Poles: f₁ = %.1f kHz, f₂ = %.1f kHz, f₃ = %.1f MHz\n",
           10e3 / (2.0 * M_PI), 100e3 / (2.0 * M_PI),
           1e6 / (2.0 * M_PI));
    printf("  Feedback β = 0.1 (closed-loop gain ≈ 10 = 20 dB)\n\n");

    /* Poles in rad/s */
    double w1 = 1e4, w2 = 1e5, w3 = 1e6;

    /* Denominator expansion: (1+s/w1)(1+s/w2)(1+s/w3)
     * = 1 + s(1/w1+1/w2+1/w3) + s²(1/(w1·w2)+1/(w2·w3)+1/(w3·w1))
     *   + s³/(w1·w2·w3) */
    double a0 = 1.0;
    double a1 = 1.0/w1 + 1.0/w2 + 1.0/w3;
    double a2 = 1.0/(w1*w2) + 1.0/(w2*w3) + 1.0/(w3*w1);
    double a3 = 1.0/(w1*w2*w3);

    /* Loop gain L(s) = K / (a0 + a1·s + a2·s² + a3·s³)
     * where K = β·A₀ = 100 */
    double K = 100.0;
    double num[1] = { K };
    double den[4] = { a0, a1, a2, a3 };

    tf_polynomial_t *L = tf_polynomial_create(num, 0, den, 3);
    if (!L) {
        printf("Error: Failed to create loop transfer function!\n");
        return 1;
    }

    printf("Loop gain L(s) = %.0f / (1 + %.2e·s + %.2e·s² + %.2e·s³)\n\n",
           K, a1, a2, a3);

    /* ================================================================
     * Part 2: Routh-Hurwitz on closed-loop characteristic equation
     *
     * Closed-loop: H(s) = A(s)/(1 + L(s))
     * Characteristic eq: 1 + L(s) = 0
     * → a₃s³ + a₂s² + a₁s + (a₀ + K) = 0
     * ================================================================ */

    double cl_coeffs[4] = { a0 + K, a1, a2, a3 };
    routh_hurwitz_t *rh = stability_routh_hurwitz(cl_coeffs, 3);

    printf("Closed-Loop Characteristic Polynomial:\n");
    printf("  s³ term: %.4e\n", cl_coeffs[3]);
    printf("  s² term: %.4e\n", cl_coeffs[2]);
    printf("  s¹ term: %.4e\n", cl_coeffs[1]);
    printf("  s⁰ term: %.4e\n", cl_coeffs[0]);

    printf("\nRouth-Hurwitz Analysis:\n");
    if (rh) {
        printf("  Sign changes in 1st column: %zu\n", rh->sign_changes);
        printf("  Is stable: %s\n", rh->is_stable ? "YES" : "NO");
        printf("  Is marginally stable: %s\n",
               rh->is_marginally_stable ? "YES" : "NO");

        if (!rh->is_stable) {
            printf("  → CLOSED-LOOP SYSTEM IS UNSTABLE!\n");
            printf("  → Number of RHP poles: %zu\n", rh->sign_changes);
        }

        routh_hurwitz_free(rh);
    }

    /* ================================================================
     * Part 3: Stability margins from Bode plot
     * ================================================================ */

    printf("\n═══════ Stability Margins from Bode Plot ═══════\n");

    stability_margin_t margins = stability_margins(L);

    printf("  Gain crossover frequency:  %.1f Hz\n",
           margins.gain_crossover_hz);
    printf("  Phase crossover frequency: %.1f Hz\n",
           margins.phase_crossover_hz);
    printf("  Phase margin:  PM = %+.1f°\n", margins.phase_margin_deg);
    printf("  Gain margin:   GM = %+.1f dB  (linear: %.2f)\n",
           margins.gain_margin_db, margins.gain_margin_linear);
    printf("  Is stable:     %s\n", margins.is_stable ? "YES" : "NO");

    printf("\nInterpretation:\n");
    if (margins.phase_margin_deg > 45.0) {
        printf("  PM > 45° → Well-damped response, good stability.\n");
    } else if (margins.phase_margin_deg > 0.0) {
        printf("  PM > 0° but < 45° → Stable but may have excessive ringing.\n");
    } else {
        printf("  PM < 0° → UNSTABLE! Will oscillate.\n");
    }

    if (margins.gain_margin_db > 6.0) {
        printf("  GM > 6 dB → Good gain margin.\n");
    } else if (margins.gain_margin_db > 0.0) {
        printf("  GM > 0 dB but < 6 dB → Marginal gain margin.\n");
    } else {
        printf("  GM < 0 dB → UNSTABLE!\n");
    }

    /* ================================================================
     * Part 4: Nyquist stability check
     * ================================================================ */

    printf("\n══════════ Nyquist Stability Criterion ══════════\n");

    int nyq_stable = 0, encirclements = 0;
    int ret = stability_nyquist_check(L, 1e-2, 1e7, 1000,
                                        &nyq_stable, &encirclements);
    if (ret == 0) {
        printf("  Encirclements of (-1, j0): N = %d\n", encirclements);
        printf("  Open-loop RHP poles:       P = ? (check Routh-Hurwitz)\n");
        printf("  Closed-loop RHP poles:     Z = N + P\n");
        printf("  Nyquist says:              %s\n",
               nyq_stable ? "STABLE" : "UNSTABLE");
    }

    /* ================================================================
     * Part 5: Root Locus — pole migration with gain
     * ================================================================ */

    printf("\n══════════════ Root Locus ══════════════\n");

    double K_values[] = { 1.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0 };
    size_t n_K = 7;
    size_t n_poles = L->den_order;

    double _Complex *poles = (double _Complex *)calloc(n_K * n_poles,
                                                         sizeof(double _Complex));
    if (poles) {
        ret = stability_root_locus(L, K_values, n_K, poles);
        if (ret == 0) {
            printf("\n  Closed-loop pole migration:\n");
            printf("  %8s  ", "Gain K");
            for (size_t j = 0; j < n_poles; j++)
                printf("  Pole %zu          ", j + 1);
            printf("\n");

            for (size_t k = 0; k < n_K; k++) {
                printf("  %8.1f  ", K_values[k]);
                for (size_t j = 0; j < n_poles; j++) {
                    double _Complex p = poles[k * n_poles + j];
                    if (cimag(p) > 1e-6) {
                        printf("  %.1f + j%.1f  ",
                               creal(p) / (2.0 * M_PI),
                               cimag(p) / (2.0 * M_PI));
                    } else {
                        printf("  %.1f Hz        ",
                               creal(p) / (2.0 * M_PI));
                    }
                }
                printf("\n");
            }
        }
        free(poles);
    }

    /* Gain margin from root locus */
    double K_margin = 0.0, f_margin = 0.0;
    ret = stability_gain_margin_root_locus(L, &K_margin, &f_margin);
    if (ret == 0) {
        printf("\n  Root locus gain margin: K_margin = %.1f\n", K_margin);
        printf("  Oscillation frequency at margin: %.1f Hz\n", f_margin);
    }

    /* ================================================================
     * Part 6: Effect of reducing feedback (increasing closed-loop gain)
     * ================================================================ */

    printf("\n══════ Effect of Feedback Factor ══════\n");
    printf("\n  %8s  %12s  %12s  %12s\n",
           "β", "CL Gain", "PM (°)", "GM (dB)");
    printf("  %8s  %12s  %12s  %12s\n",
           "------", "----------", "----------", "----------");

    double beta_values[] = { 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01 };
    for (int i = 0; i < 7; i++) {
        double beta = beta_values[i];
        double loop_K = beta * 1000.0;

        double num_b[1] = { loop_K };
        tf_polynomial_t *Lb = tf_polynomial_create(num_b, 0, den, 3);
        stability_margin_t m = stability_margins(Lb);

        printf("  %8.3f  %12.1f  %+12.1f  %+12.1f\n",
               beta, 1.0 / beta, m.phase_margin_deg, m.gain_margin_db);

        tf_polynomial_free(Lb);
    }

    printf("\n  Note: Reducing feedback β (higher CL gain) improves stability\n");
    printf("  but reduces bandwidth (gain-bandwidth trade-off).\n");

    /* Cleanup */
    tf_polynomial_free(L);

    printf("\nStability analysis complete.\n");
    return 0;
}
