#include <stdio.h>
#include <math.h>
#include "transient.h"

int main(void) {
  printf("=== RC Transient Analysis ===\n");
  double R = 1000, C = 10e-6, tau = rc_time_constant(R, C);
  printf("RC: R=1kOhm, C=10uF, tau=%.4fs\n", tau);
  printf("Charging 0V->5V:\n");
  for (double t = 0; t <= 5 * tau; t += tau)
    printf("  t=%.4fs vC=%.4fV i=%.4fmA\n", t,
           rc_charging_vc(5, 0, tau, t),
           rc_charging_i(5, 0, R, tau, t) * 1000);

  printf("\nDischarging from 5V:\n");
  for (double t = 0; t <= 3 * tau; t += tau)
    printf("  t=%.4fs vC=%.4fV\n", t, rc_discharging_vc(5, tau, t));

  double f555 = timer_555_astable_freq(1000, 10000, 10e-6);
  double duty = timer_555_astable_duty(1000, 10000);
  printf("\n555 Timer (R1=1k,R2=10k,C=10uF): f=%.2fHz duty=%.1f%%\n",
         f555, duty * 100);
  printf("  t_high=%.3fs t_low=%.3fs\n", duty / f555, (1 - duty) / f555);

  double Trise = rc_rise_time_10_90(R, 10e-6);
  printf("Rise time (10-90%%): %.4fs\n", Trise);

  double Rtemp = temp_dependent_resistance(1000, 25, 85, 0.00393);
  printf("Resistor at 85C: %.2f Ohm\n", Rtemp);
  return 0;
}
