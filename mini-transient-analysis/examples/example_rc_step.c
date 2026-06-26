#include <stdio.h>
#include <math.h>
#include "../include/transient_defs.h"
#include "../include/first_order.h"

int main(void) {
    printf("=== RC Step Response Example ===\n");
    double R = 1000.0, C = 1e-6, Vs = 5.0, v0 = 0.0;
    double tau = rc_time_constant(R, C);
    printf("R=%.0f ohm, C=%.1f uF, tau=%.3f ms\n", R, C*1e6, tau*1000);

    printf("\nTime(ms)  Vc(V)    I(mA)    %%Final\n");
    for (double t = 0; t <= 5*tau; t += tau/4) {
        double vc = rc_step_charge_vc(Vs, v0, R, C, t);
        double i = rc_step_charge_i(Vs, v0, R, C, t);
        double pct = vc / Vs * 100;
        printf("%8.2f %8.3f %8.3f %8.1f\n", t*1000, vc, i*1000, pct);
    }

    printf("\nDischarging from Vc=5V:\n");
    printf("Time(ms)  Vc(V)\n");
    for (double t = 0; t <= 5*tau; t += tau/2) {
        double vc = rc_natural_vc(5.0, R, C, t);
        printf("%8.2f %8.3f\n", t*1000, vc);
    }

    double tr = rc_delay_time(R, C, 0.9*Vs, Vs) - rc_delay_time(R, C, 0.1*Vs, Vs);
    printf("\n10-90%% Rise time: %.3f ms (theoretical: %.3f ms)\n", tr*1000, 2.2*tau*1000);

    double v_measure = rc_step_charge_vc(Vs, v0, R, C, tau);
    printf("At t=tau: Vc = %.3f V (%.1f%% of final)\n", v_measure, v_measure/Vs*100);
    return 0;
}
