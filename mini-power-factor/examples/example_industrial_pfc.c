/**
 * @file example_industrial_pfc.c
 * @brief Complete industrial power factor correction analysis
 *
 * Demonstrates end-to-end PF correction for a typical industrial facility:
 * - 500 kW motor load at PF=0.70 lagging (480V, 60Hz)
 * - Target PF = 0.95
 * - Computes capacitor sizing, resonance check, and payback period
 *
 * This is a canonical problem (L6) covering all aspects of
 * industrial PF correction: sizing, economics, and harmonic risk.
 */

#include <stdio.h>
#include <math.h>
#include "../include/power_factor.h"
#include "../include/complex_power.h"
#include "../include/pf_correction.h"
#include "../include/harmonic_power.h"

int main(void)
{
    printf("========================================\n");
    printf("  Industrial PF Correction — Complete Analysis\n");
    printf("  Facility: 500 kW, 480V, 60Hz, PF 0.70→0.95\n");
    printf("========================================\n\n");

    /* ==================================================================
     * Step 1: Analyze current power consumption
     * ================================================================== */
    double p_load  = 500000.0;   /* 500 kW */
    double v_ll    = 480.0;      /* 480V line-to-line */
    double f_hz    = 60.0;
    double pf_old  = 0.70;
    double pf_target = 0.95;
    double s_sc_kva = 5000.0;    /* 5 MVA short-circuit at PCC */

    pf_single_phase_t current_power;
    pf_compute_single_phase(v_ll, p_load / (v_ll * pf_old * 1.732), 0.0,
                            &current_power);
    /* The above uses per-phase approximation. Better: use three-phase directly. */

    double phi_old = acos(pf_old);
    double q_old = p_load * tan(phi_old);

    printf("Current Load Profile:\n");
    printf("  Real Power:       %8.1f kW\n", p_load / 1000.0);
    printf("  Reactive Power:   %8.1f kVAR\n", q_old / 1000.0);
    printf("  Apparent Power:   %8.1f kVA\n", p_load / (pf_old * 1000.0));
    printf("  Power Factor:     %.2f %s\n",
           pf_old, "lagging (inductive)");
    printf("  Line Current:     %8.1f A\n\n",
           p_load / (1.732 * v_ll * pf_old));

    /* ==================================================================
     * Step 2: Compute required compensation
     * ================================================================== */
    double c_delta, kvar_needed, payback_mon, f_res;
    pfc_solve_industrial(p_load, v_ll, f_hz,
                          pf_old, pf_target, s_sc_kva,
                          &c_delta, &kvar_needed,
                          &payback_mon, &f_res);

    printf("PF Correction Solution:\n");
    printf("  Required kVAR:     %8.1f kVAR\n", kvar_needed / 1000.0);
    printf("  Capacitance (Δ):   %8.2f μF per phase\n",
           c_delta * 1e6);
    printf("  Resonant freq:     %8.1f Hz\n", f_res);
    printf("  Payback period:    %8.1f months\n\n", payback_mon);

    /* ==================================================================
     * Step 3: Resonance and detuning assessment
     * ================================================================== */
    int harmonic_order;
    int needs_detune = pfc_needs_detuning(f_res, f_hz, &harmonic_order);
    printf("Harmonic Resonance Assessment:\n");
    if (needs_detune) {
        printf("  *** WARNING: Resonance near %dth harmonic! ***\n",
               harmonic_order);
        printf("  Recommend: 7%% detuning reactor (tuned to 189 Hz)\n");
    } else {
        printf("  Resonance at %.0f Hz — no detuning required.\n", f_res);
    }
    printf("  Nearest characteristic harmonic: %dth\n\n", harmonic_order);

    /* ==================================================================
     * Step 4: Design capacitor bank steps
     * ================================================================== */
    pfc_capacitor_bank_t bank;
    pfc_design_step_bank(kvar_needed / 1000.0, 5, 50.0, &bank);
    printf("Capacitor Bank Configuration (%u steps):\n", bank.num_steps);
    for (uint32_t i = 0; i < bank.num_steps; i++) {
        printf("  Step %u: %8.1f kVAR\n", i + 1, bank.step_kvar[i]);
    }
    printf("  Total installed: %.1f kVAR\n\n", bank.total_kvar);

    /* ==================================================================
     * Step 5: Verify power after correction
     * ================================================================== */
    double q_new = q_old - kvar_needed;
    double pf_new = cos(atan2(q_new, p_load));
    double s_new = p_load / pf_new;
    double i_new = s_new / (1.732 * v_ll);
    double i_old = p_load / (1.732 * v_ll * pf_old);

    printf("Post-Correction Profile:\n");
    printf("  New Reactive Power: %8.1f kVAR\n", q_new / 1000.0);
    printf("  New Power Factor:   %.3f\n", pf_new);
    printf("  New Line Current:   %8.1f A\n", i_new);
    printf("  Current Reduction:  %8.1f A (%.1f%%)\n\n",
           i_old - i_new, (i_old - i_new) / i_old * 100.0);

    /* ==================================================================
     * Step 6: Economic analysis
     * ================================================================== */
    double annual_savings = pf_savings_estimate(p_load, pf_old, pf_new,
                                                 0.02, 8000.0, 0.10);
    double cap_cost = kvar_needed * 30.0 / 1000.0; /* $30/kVAR installed */

    printf("Economic Analysis:\n");
    printf("  Capacitor bank cost:   $%8.0f\n", cap_cost);
    printf("  Annual energy savings: $%8.0f\n", annual_savings);
    printf("  Simple payback:        %8.1f months\n", payback_mon);
    printf("  == INVESTMENT RECOMMENDED ==\n\n");

    printf("========================================\n");
    return 0;
}
