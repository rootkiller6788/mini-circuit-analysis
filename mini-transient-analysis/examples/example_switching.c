#include <stdio.h>
#include <math.h>
#include "../include/transient_defs.h"
#include "../include/switching_transient.h"

int main(void) {
    printf("=== Switching Transient Application Examples ===\n\n");

    /* Buck Converter Example (L7: Toyota Prius DC-DC, 2020) */
    printf("--- Buck Converter: 12V to 5V @ 100kHz ---\n");
    double Vin = 12.0, Vout = 5.0, L = 100e-6, C = 47e-6, R_load = 2.5;
    double f_sw = 100e3;
    double duty = buck_duty_cycle_ideal(Vout, Vin);
    double ripple_I = buck_inductor_current_ripple(Vin, Vout, L, f_sw);
    double ripple_V = buck_output_voltage_ripple(ripple_I, C, 0.05, f_sw);

    printf("Duty cycle: %.1f%%\n", duty*100);
    printf("Inductor ripple: %.0f mA\n", ripple_I*1000);
    printf("Output ripple: %.1f mV\n", ripple_V*1000);
    printf("Critical inductance: %.1f uH\n\n", buck_critical_inductance(Vin, Vout, 0.5, f_sw)*1e6);

    /* MOSFET Gate Driver (L7: Tesla Model 3 inverter, 2017) */
    printf("--- MOSFET Gate Driver ---\n");
    double Qg = 100e-9, Vdrive = 12.0, t_sw = 50e-9;
    double Ig = mosfet_gate_drive_current(Qg, t_sw);
    double Rg = mosfet_gate_resistor(Vdrive, Ig);
    double P_drive = gate_driver_power_dissipation(Qg, Vdrive, f_sw);

    printf("Gate current: %.2f A\n", Ig);
    printf("Gate resistor: %.1f ohm\n", Rg);
    printf("Driver power: %.1f mW\n\n", P_drive*1000);

    /* PLL Example (L8: GPS receiver, 2020) */
    printf("--- PLL Lock Time: GPS L1 (1575.42 MHz) ---\n");
    double Kpd = 0.1, Kvco = 10e6, N = 100;
    double wn = pll_natural_frequency(Kpd, Kvco, N);
    double zeta = pll_damping_factor(Kpd, Kvco, N, 1000, 1e-6);
    double t_lock = pll_lock_time_estimate(wn, zeta, 1e6, 100);

    printf("Natural freq: %.1f kHz\n", wn/1e3);
    printf("Damping: %.2f\n", zeta);
    printf("Lock time: %.1f us\n\n", t_lock*1e6);

    /* Transmission Line (L7: DDR4 memory bus, 2020) */
    printf("--- Transmission Line: DDR4 100mm trace ---\n");
    double length = 0.1, vf = 0.5;
    double Td = tl_propagation_delay(length, vf);
    double tr = 200e-12;
    int is_dist = tl_is_lumped_or_distributed(length, tr, vf);

    printf("Propagation delay: %.1f ps\n", Td*1e12);
    printf("Critical length: %.1f mm\n", tl_critical_length(tr, vf)*1000);
    printf("Is distributed: %s\n", is_dist ? "YES" : "NO");

    return 0;
}
