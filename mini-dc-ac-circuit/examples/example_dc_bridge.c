#include <stdio.h>
#include <math.h>
#include "circuit_elements.h"
#include "dc_analysis.h"

int main(void) {
  printf("=== Wheatstone Bridge DC Analysis ===\n");
  Circuit_t ckt;
  circuit_init(&ckt, "bridge");
  circuit_add_node(&ckt, "A", 0);
  circuit_add_node(&ckt, "B", 0);
  circuit_add_node(&ckt, "C", 0);
  circuit_add_node(&ckt, "D", 0);
  ckt.ground_node = 0;

  Resistor_t R1, R2, R3, Rx;
  resistor_init(&R1, 0, "R1", 0, 1, 1000);
  resistor_init(&R2, 1, "R2", 1, 2, 1000);
  resistor_init(&R3, 2, "R3", 0, 3, 1000);
  resistor_init(&Rx, 3, "Rx", 3, 2, 1100);
  ckt.resistors[0] = R1; ckt.resistors[1] = R2;
  ckt.resistors[2] = R3; ckt.resistors[3] = Rx;
  ckt.n_resistors = 4;

  DCVoltageSource_t Vs;
  dc_voltage_source_init(&Vs, 0, "Vs", 0, 2, 10);
  ckt.dc_vsrc[0] = Vs; ckt.n_dc_vsrc = 1;

  double Vout = wheatstone_bridge_vout(10, 1000, 1000, 1000, 1100);
  double Rx_bal = wheatstone_bridge_rx_balanced(1000, 1000, 1000);
  int is_bal = wheatstone_is_balanced(1000, 1000, 1000, 1000, 1e-6);

  printf("Bridge Vout = %.6f V\n", Vout);
  printf("Rx for balance = %.1f Ohm\n", Rx_bal);
  printf("Is balanced (R1=R2=R3=Rx): %s\n", is_bal ? "yes" : "no");
  printf("Max power transfer efficiency: %.0f%%\n", max_power_efficiency() * 100);

  double R_dv[] = {100.0, 200.0, 300.0};
  double Vout_dv[3];
  voltage_divider_n(10, R_dv, Vout_dv, 3);
  printf("Voltage divider (100+200+300):");
  for (int i = 0; i < 3; i++) printf(" V[%d]=%.2f", i, Vout_dv[i]);
  printf("\n");
  return 0;
}
