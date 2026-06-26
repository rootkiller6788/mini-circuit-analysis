#include <stdio.h>
#include <math.h>
#include "../include/transient_defs.h"
#include "../include/second_order.h"

int main(void) {
    printf("=== RLC Series Step Response Example ===\n\n");

    double L = 10e-3, C = 10e-6;
    printf("L=%.1f mH, C=%.1f uF\n\n", L*1e3, C*1e6);

    /* Demonstrate three damping regimes */
    double R_vals[] = {2.0, 63.245, 200.0};
    const char *regimes[] = {"Underdamped", "Critically Damped", "Overdamped"};

    for (int k = 0; k < 3; k++) {
        double R = R_vals[k];
        double zeta = rlc_series_zeta(R, L, C);
        double omega_n = rlc_series_omega_n(L, C);
        double omega_d = rlc_series_omega_d(R, L, C);

        printf("--- R=%.1f ohm: %s (zeta=%.2f) ---\n", R, regimes[k], zeta);
        printf("omega_n=%.0f rad/s, omega_d=%.0f rad/s\n", omega_n, omega_d);

        double os = overshoot_from_zeta(zeta);
        double tp = peak_time_from_params(omega_n, zeta);
        double ts = settling_time_2pct(zeta, omega_n);

        printf("Overshoot: %.1f%%, Peak time: %.3f ms, Settling: %.3f ms\n",
               os, tp*1000, ts*1000);

        printf("Time(ms)  Vc(V)    IL(mA)\n");
        for (int i = 0; i <= 20; i++) {
            double t = i * ts / 20.0;
            double vc = rlc_series_step_vc(10.0, R, L, C, 0.0, 0.0, t);
            double il = rlc_series_step_il(10.0, R, L, C, 0.0, 0.0, t);
            printf("%8.2f %8.3f %8.3f\n", t*1000, vc, il*1000);
        }
        printf("\n");
    }

    return 0;
}
