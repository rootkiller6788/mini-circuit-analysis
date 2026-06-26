/**
 * @file example_three_phase_pq.c
 * @brief Three-phase power quality analysis with symmetrical components
 *
 * Demonstrates three-phase unbalanced power analysis using Fortescue
 * symmetrical components, Clarke/Park transforms, and SRF-PLL grid
 * synchronization for a motor drive system.
 *
 * This example (L7) simulates a real-world application: detecting
 * voltage unbalance in an industrial motor bus and assessing its
 * impact on power factor and motor performance.
 */

#include <stdio.h>
#include <math.h>
#include "../include/power_factor.h"
#include "../include/phasor_operations.h"
#include "../include/power_quality.h"
#include "../include/complex_power.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("========================================\n");
    printf("  Three-Phase Power Quality Analysis\n");
    printf("  Motor Bus: 480V, 60Hz, 200HP Induction Motor\n");
    printf("========================================\n\n");

    /* ==================================================================
     * Part 1: Balanced Three-Phase System (baseline)
     * ================================================================== */
    printf("--- Part 1: Balanced System Baseline ---\n");

    pf_three_phase_t balanced;
    pf_compute_three_phase(480.0, 200.0, 25.8, 1, &balanced);

    printf("  P_total  = %.1f kW\n", balanced.p_total / 1000.0);
    printf("  Q_total  = %.1f kVAR\n", balanced.q_total / 1000.0);
    printf("  S_total  = %.1f kVA\n", balanced.s_total / 1000.0);
    printf("  PF       = %.3f %s\n",
           balanced.pf_total,
           balanced.type == PF_TYPE_LAGGING ? "lagging" : "leading");
    printf("  Per-phase powers: Pa=%.1f, Pb=%.1f, Pc=%.1f kW\n\n",
           balanced.p_a / 1000.0, balanced.p_b / 1000.0, balanced.p_c / 1000.0);

    /* ==================================================================
     * Part 2: Unbalanced Condition (Phase B sag to 90%)
     * ================================================================== */
    printf("--- Part 2: Unbalanced System (Phase B at 90%%) ---\n");

    double v_unbalance = pf_voltage_unbalance_percent(480.0, 432.0, 480.0);
    printf("  Line voltages: Vab=480, Vbc=432, Vca=480 V\n");
    printf("  NEMA unbalance = %.1f %%\n\n", v_unbalance);

    /* Fortescue symmetrical components */
    /* Phase voltages: Va=277∠0°, Vb=249∠-120°, Vc=277∠120° */
    double va_ln = 277.0, vb_ln = 249.0, vc_ln = 277.0;
    double va_ang = 0.0;
    double vb_ang = -120.0 * M_PI / 180.0;
    double vc_ang =  120.0 * M_PI / 180.0;

    double v0_m, v0_a, v1_m, v1_a, v2_m, v2_a;
    pf_symmetrical_components(va_ln, va_ang,
                               vb_ln, vb_ang,
                               vc_ln, vc_ang,
                               &v0_m, &v0_a, &v1_m, &v1_a, &v2_m, &v2_a);

    printf("  Symmetrical Components (Fortescue 1918):\n");
    printf("    V0 = %.2f ∠ %.1f° V (zero seq)\n", v0_m, v0_a * 180.0 / M_PI);
    printf("    V1 = %.2f ∠ %.1f° V (positive seq)\n", v1_m, v1_a * 180.0 / M_PI);
    printf("    V2 = %.2f ∠ %.1f° V (negative seq)\n", v2_m, v2_a * 180.0 / M_PI);
    printf("    Unbalance (V2/V1) = %.2f %%\n\n",
           (v1_m > 0.0) ? (v2_m / v1_m * 100.0) : 0.0);

    /* Negative sequence → reverse rotating field → motor heating */
    printf("  Impact Assessment:\n");
    double neg_seq_pct = (v1_m > 0.0) ? (v2_m / v1_m * 100.0) : 0.0;
    if (neg_seq_pct > 5.0) {
        printf("    *** CRITICAL: >5%% negative sequence → motor derating needed\n");
    } else if (neg_seq_pct > 2.0) {
        printf("    *** WARNING: >2%% negative sequence → increased heating\n");
    } else {
        printf("    Acceptable: <2%% negative sequence\n");
    }
    printf("    NEMA MG1: Derate motor by %.1f%% for %.1f%% unbalance\n\n",
           neg_seq_pct * neg_seq_pct * 2.0, v_unbalance);

    /* ==================================================================
     * Part 3: Clarke and Park Transforms
     * ================================================================== */
    printf("--- Part 3: Clarke & Park Transforms ---\n");

    /* At t=0 with θ=0, balanced three-phase: */
    double a = 1.0, b = -0.5, c = -0.5;
    double alpha, beta, zero;
    phasor_clarke_transform(a, b, c, &alpha, &beta, &zero);

    printf("  Clarke (αβ0) at t=0:\n");
    printf("    α = %.4f, β = %.4f, 0 = %.4f\n\n", alpha, beta, zero);

    dq0_t dq;
    phasor_park_transform(a, b, c, 0.0, &dq);
    printf("  Park (dq0) at θ=0:\n");
    printf("    d = %.4f, q = %.4f, 0 = %.4f\n", dq.d, dq.q, dq.zero);
    printf("    |dq| = %.4f, angle = %.4f rad\n\n",
           dq.magnitude, dq.angle_rad);

    /* ==================================================================
     * Part 4: SRF-PLL Grid Synchronization
     * ================================================================== */
    printf("--- Part 4: SRF-PLL Grid Synchronization ---\n");

    srf_pll_t pll;
    double omega_grid = 2.0 * M_PI * 60.0;
    srf_pll_init(&pll, omega_grid, 10.0, 50e-6);  /* 50μs, 20kHz control */

    /* Simulate 100ms of PLL tracking */
    double omega_error_max = 0.0;
    for (int i = 0; i < 2000; i++) {
        double t = i * 50e-6;
        double va = 277.0 * cos(omega_grid * t);
        double vb = 277.0 * cos(omega_grid * t - 2.0 * M_PI / 3.0);
        double vc = 277.0 * cos(omega_grid * t + 2.0 * M_PI / 3.0);
        srf_pll_step(&pll, va, vb, vc);

        double omega_err = fabs(srf_pll_get_frequency(&pll) - omega_grid);
        if (omega_err > omega_error_max) omega_error_max = omega_err;
    }

    printf("  After 100ms lock:\n");
    printf("    Frequency = %.3f rad/s (%.2f Hz)\n",
           srf_pll_get_frequency(&pll),
           srf_pll_get_frequency(&pll) / (2.0 * M_PI));
    printf("    Angle = %.4f rad\n", srf_pll_get_angle(&pll));
    printf("    V_mag = %.2f V (true=277.0)\n", srf_pll_get_magnitude(&pll));
    printf("    Max frequency error = %.3f rad/s\n\n", omega_error_max);

    /* ==================================================================
     * Part 5: PF Correction for Motor Load
     * ================================================================== */
    printf("--- Part 5: PF Correction for 200HP Motor ---\n");

    double p_motor = 150000.0;  /* 150 kW (≈200HP) */
    double pf_motor = 0.82;
    double pf_target = 0.95;

    /* Current reactive power */
    double q_motor = p_motor * tan(acos(pf_motor));
    printf("  Before correction: P=%.1f kW, Q=%.1f kVAR, PF=%.2f\n",
           p_motor / 1000.0, q_motor / 1000.0, pf_motor);

    double c_phase, kvar_total;
    cp_three_phase_pf_capacitor(480.0, p_motor, 60.0,
                                 pf_motor, pf_target, 1,
                                 &c_phase, &kvar_total);

    printf("  Delta capacitor bank: %.1f μF per phase\n", c_phase * 1e6);
    printf("  Total kVAR: %.1f kVAR\n", kvar_total / 1000.0);

    double q_after = q_motor - kvar_total;
    double pf_after = cos(atan2(q_after, p_motor));
    printf("  After correction:  Q=%.1f kVAR, PF=%.3f\n",
           q_after / 1000.0, pf_after);

    /* ==================================================================
     * Part 6: Data Center PQ Assessment
     * ================================================================== */
    printf("\n--- Part 6: Data Center PQ Assessment ---\n");

    int tier; double score;
    int issues = pq_assess_datacenter(0.97, 3.0, 4.0, 96.0, &tier, &score);
    printf("  PF=0.97, THDi=3%%, ΔV=4%%, η=96%%\n");
    printf("  80 PLUS Tier: %d (Titanium)\n", tier);
    printf("  PQ Score: %.0f/100\n", score);
    printf("  Issues: %s\n\n", (issues == 0) ? "None" : "See bitmask");

    /* ==================================================================
     * Part 7: Cost of Poor PQ for this facility
     * ================================================================== */
    printf("--- Part 7: Cost of Power Quality (EPRI Method) ---\n");

    double annual_kwh = p_motor * 8000.0 / 1000.0;  /* 8000 hours/year */
    double pq_cost = pq_cost_of_poor_quality(0.82, annual_kwh,
                                              0.10, 0.75, 5, 5000.0);
    printf("  Annual energy: %.0f MWh\n", annual_kwh / 1000.0);
    printf("  Annual PQ cost: $%.0f\n\n", pq_cost);

    printf("========================================\n");
    return 0;
}
