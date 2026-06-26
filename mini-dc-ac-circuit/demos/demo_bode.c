#include <stdio.h>
#include <math.h>
#include "ac_analysis.h"

int main(void) {
  printf("=== Bode Plot Demo ===\n");
  double mag[100], ph[100];
  bode_single_pole(1000, 1, 1e6, 100, mag, ph);
  printf("Single pole at 1kHz:\n  Freq(Hz)  Mag(dB)  Phase(deg)\n");
  for (int i = 0; i < 100; i += 10)
    printf("  %8.1f %8.2f %8.1f\n",
           pow(10, log10(1) + log10(1e6) * i / 99), mag[i], ph[i]);

  printf("\nComplex poles (wn=1k, z=0.5):\n");
  bode_complex_poles(2 * M_PI * 1000, 0.5, 1, 1e6, 100, mag, ph);
  for (int i = 0; i < 100; i += 10)
    printf("  %8.1f %8.2f %8.1f\n",
           pow(10, log10(1) + log10(1e6) * i / 99), mag[i], ph[i]);
  return 0;
}
