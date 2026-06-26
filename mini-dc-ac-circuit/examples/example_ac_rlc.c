#include <stdio.h>
#include <math.h>
#include "ac_analysis.h"
#ifndef I
#define I _Complex_I
#endif

int main(void) {
  printf("=== RLC AC Analysis ===\n");
  double R = 50, L = 1e-3, C = 1e-6;
  Resonance_t res = resonance_series_analyze(R, L, C);
  printf("Series RLC (R=50, L=1mH, C=1uF):\n");
  printf("  f0 = %.2f Hz\n", res.f0_hz);
  printf("  Q  = %.2f\n", res.Q);
  printf("  BW = %.2f Hz\n", res.BW_hz);
  printf("  zeta = %.4f\n", res.zeta);
  printf("  Z at res = %.2f Ohm\n", res.Z_at_res);

  FreqPoint_t resp[20];
  freq_response_rlc_series(R, L, C, FR_RLC_ACROSS_R,
                           res.f0_hz / 10, res.f0_hz * 10, 20, resp);
  printf("  Freq(Hz) -> |H|:\n");
  for (int i = 0; i < 20; i += 2)
    printf("  %.1f Hz -> %.4f\n", resp[i].freq_hz, resp[i].magnitude);

  ComplexPower_t cp = complex_power_from_v_z(120, 50 + 0 * I);
  printf("\nPower: P=%.1fW Q=%.1fvar S=%.1fVA pf=%.3f\n",
         cp.P, cp.Q, cp.S, cp.pf);

  double Cpf = pf_correction_unity(1000, 120, 60, 0.65, 1);
  printf("PF correction 1kW 0.65lag->unity: C=%.1fuF\n", Cpf * 1e6);
  return 0;
}
