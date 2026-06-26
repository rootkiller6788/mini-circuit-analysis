#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "circuit_elements.h"
#include "circuit_utils.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %s: ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)

static void test_resistor_init(void) {
  TEST("resistor_init");
  Resistor_t r;
  int ret = resistor_init(&r, 1, "R1", 0, 1, 1000.0);
  assert(ret == 0);
  assert(r.R == 1000.0);
  assert(r.t[0].node_id == 0);
  assert(r.t[1].node_id == 1);
  assert(strcmp(r.name, "R1") == 0);
  PASS();
}

static void test_impedance(void) {
  TEST("impedance");
  Impedance_t Zr = impedance_resistor(50.0);
  assert(fabs(Zr.re - 50.0) < 1e-9);
  assert(fabs(Zr.im) < 1e-9);
  Impedance_t Zc = impedance_capacitor(1e-6, 1000.0);
  assert(Zc.im < 0);
  PASS();
}

static void test_resonance(void) {
  TEST("resonance");
  double f0 = resonance_freq_hz(1e-3, 1e-6);
  double expected = 1.0 / (2.0 * M_PI * sqrt(1e-3 * 1e-6));
  assert(fabs(f0 - expected) < expected * 1e-6);
  PASS();
}

static void test_q_factor(void) {
  TEST("q_factor");
  double Qs = q_factor_series(10.0, 1e-3, 1e-6);
  assert(Qs > 0);
  PASS();
}

static void test_damping(void) {
  TEST("damping");
  double zeta = damping_factor_series(100.0, 1e-3, 1e-6);
  assert(zeta > 0);
  PASS();
}

int main(void) {
  printf("=== test_elements ===\n");
  test_resistor_init();
  test_impedance();
  test_resonance();
  test_q_factor();
  test_damping();
  printf("Results: %d/%d passed\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
