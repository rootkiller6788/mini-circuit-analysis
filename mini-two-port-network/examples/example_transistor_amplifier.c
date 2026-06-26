/**
 * @file example_transistor_amplifier.c
 * @brief Example: BJT Common-Emitter Amplifier Two-Port Analysis
 *
 * Demonstrates how to use Y-parameters and H-parameters for analyzing
 * a common-emitter transistor amplifier stage.
 *
 * This example covers the complete analysis flow taught in Berkeley EE105
 * and Stanford EE214:
 *   1. Model the BJT using hybrid-π parameters
 *   2. Compute H-parameters for the CE configuration
 *   3. Calculate the loaded voltage gain
 *   4. Compute input and output impedances
 *   5. Demonstrate the Miller effect on input capacitance
 *   6. Convert to other parameter types for cascade analysis
 *
 * Course: Berkeley EE105 — BJT Amplifier Design
 * Ref: Sedra & Smith §7.3, §9.4
 */

#include <stdio.h>
#include <math.h>
#include "../include/two_port.h"
#include "../include/h_params.h"
#include "../include/y_params.h"
#include "../include/z_params.h"
#include "../include/conversion.h"
#include "../include/abcd_params.h"

int main(void) {
    printf("========================================\n");
    printf("BJT Common-Emitter Amplifier Analysis\n");
    printf("========================================\n\n");

    /* ============================================================
     * Step 1: Define BJT parameters (2N3904-like small-signal)
     * ============================================================ */
    double rbb = 50.0;       /* Base spreading resistance (Ω) */
    double rpi = 2500.0;     /* Base-emitter resistance (Ω) */
    double beta = 150.0;     /* DC current gain hFE */
    double ro = 50000.0;     /* Output resistance (Ω) */
    double gm = beta / rpi;  /* Transconductance = 150/2500 = 0.06 S */
    double ic = gm * 0.02585; /* Ic ≈ 1.55 mA (at VT=25.85mV) */

    printf("BJT Operating Point:\n");
    printf("  Ic = %.2f mA\n", ic * 1000.0);
    printf("  gm = %.1f mS\n", gm * 1000.0);
    printf("  rπ = %.1f kΩ\n", rpi / 1000.0);
    printf("  ro = %.1f kΩ\n", ro / 1000.0);
    printf("  β  = %.0f\n\n", beta);

    /* ============================================================
     * Step 2: Compute H-parameters (the classic BJT model)
     * ============================================================ */
    matrix2x2_t h_ce = hparams_bjt_ce(rbb, rpi, beta, ro);
    printf("CE H-Parameters (traditional BJT model):\n");
    printf("  h_ie = %.1f Ω     (input impedance)\n", h_ce.m11.real);
    printf("  h_re = %.6f       (reverse voltage, ≈ 0)\n", h_ce.m12.real);
    printf("  h_fe = %.1f        (forward current gain)\n", h_ce.m21.real);
    printf("  h_oe = %.6f S     (output admittance)\n\n", h_ce.m22.real);

    /* ============================================================
     * Step 3: Amplifier design — bias and load
     * ============================================================ */
    double rc = 4700.0;     /* Collector resistor (Ω) */
    double rl = 10000.0;    /* Load resistor (Ω) */
    double rload = (rc * rl) / (rc + rl);  /* AC load = RC || RL */
    complex_t yl = complex_make(1.0 / rload, 0.0);  /* Load admittance */
    complex_t zl = complex_make(rload, 0.0);

    printf("Amplifier Design:\n");
    printf("  RC = %.1f kΩ\n", rc / 1000.0);
    printf("  RL = %.1f kΩ\n", rl / 1000.0);
    printf("  AC Load = %.1f kΩ\n\n", rload / 1000.0);

    /* ============================================================
     * Step 4: Compute voltage gain using H-parameters
     * ============================================================ */
    complex_t av_h = hparams_voltage_gain(h_ce, yl);
    double av_h_mag = complex_mag(av_h);
    double av_h_db = 20.0 * log10(av_h_mag);

    printf("Voltage Gain Analysis:\n");
    printf("  Av = %.1f (%.1f dB)\n", av_h_mag, av_h_db);
    printf("  Phase = %.1f°\n", complex_arg(av_h) * 180.0 / M_PI);

    /* Compare with simplified formula: Av ≈ -gm * RL */
    double av_simple = -gm * rload;
    printf("  Simplified (-gm*RL) = %.1f\n\n", av_simple);

    /* ============================================================
     * Step 5: Input and output impedance
     * ============================================================ */
    complex_t zin_h = hparams_input_impedance(h_ce, yl);
    complex_t zout_h = hparams_output_admittance(h_ce,
        complex_make(50.0, 0.0));  /* Source impedance = 50Ω */

    printf("Impedance Analysis:\n");
    printf("  Zin  = %.1f Ω\n", zin_h.real);
    printf("  Zout = %.1f Ω\n\n", 1.0 / zout_h.real);

    /* ============================================================
     * Step 6: Miller effect — input capacitance
     * ============================================================ */
    double cf = 2e-12;  /* Cμ = Cbc ≈ 2 pF (typical) */
    double av_mag = fabs(av_simple);
    double cmiller = miller_input_capacitance(cf, av_mag);
    double cpi = 8e-12;  /* Cπ = Cbe ≈ 8 pF */
    double ctotal = cpi + cmiller;

    printf("Miller Effect:\n");
    printf("  Cμ (feedback)   = %.1f pF\n", cf * 1e12);
    printf("  Cπ (base-emitter) = %.1f pF\n", cpi * 1e12);
    printf("  C_miller (input)  = %.1f pF\n", cmiller * 1e12);
    printf("  C_total (input)   = %.1f pF\n", ctotal * 1e12);

    /* Compute the dominant pole frequency */
    double rth = zin_h.real;  /* Thevenin resistance at base */
    double f_h = 1.0 / (2.0 * M_PI * rth * ctotal);
    printf("  Dominant pole fh ≈ %.1f MHz\n\n", f_h / 1e6);

    /* ============================================================
     * Step 7: Convert to ABCD for cascade analysis
     * ============================================================ */
    matrix2x2_t abcd_amp = convert_parameters(h_ce, PARAM_H, PARAM_ABCD, 0.0);
    printf("ABCD Parameters (for cascade):\n");
    printf("  A = %.1f, B = %.1f Ω\n", abcd_amp.m11.real, abcd_amp.m12.real);
    printf("  C = %.6f S, D = %.1f\n", abcd_amp.m21.real, abcd_amp.m22.real);

    /* Cascade: Source(ABCD) → Amplifier(ABCD) → Load(ABCD) */
    matrix2x2_t abcd_load = abcd_series_z(complex_make(rload, 0.0));
    matrix2x2_t abcd_total = matrix2x2_mul(abcd_amp, abcd_load);

    complex_t av_cascade = abcd_voltage_transfer(abcd_total, complex_make(50, 0));
    printf("\nCascade analysis verification:\n");
    printf("  Av (from ABCD cascade) = %.1f\n", complex_mag(av_cascade));
    printf("\n========================================\n");
    printf("Analysis Complete (L4: Amplifier Two-Port)\n");
    printf("========================================\n");

    return 0;
}
