# Course Alignment — mini-power-factor

Mapping to nine global university curricula.

## MIT (Massachusetts Institute of Technology)

| Course | Topic | This Module |
|--------|-------|-------------|
| 6.061 "Electric Power Systems" | Complex power, PF correction, three-phase | L2-L5: complex_power.c, pf_correction.c |
| 6.003 "Signal Processing" | Fourier analysis, Parseval | L3: harmonic_power.c FFT, Goertzel |
| 6.131 "Power Electronics" | Boost PFC, DC-DC converters | L5: pfc_boost_init/step, ACMC |

## Stanford University

| Course | Topic | This Module |
|--------|-------|-------------|
| EE251 "Power Electronics" | PFC control, efficiency | L5: boost PFC with PI cascaded |
| EE292 "Smart Grid" | PQ monitoring, demand response | L7: pq_monitor, demand_interval |

## UC Berkeley

| Course | Topic | This Module |
|--------|-------|-------------|
| EE137A "Power Electronics" | PF correction, harmonic standards | L5-L6: pfc_solve_industrial |
| EE16A/B "Circuits" | AC steady-state, phasors | L1-L3: phasor_operations.c |

## UIUC (Illinois)

| Course | Topic | This Module |
|--------|-------|-------------|
| ECE 431 "Electric Machinery" | Three-phase, dq0 transforms | L3: Clarke/Park transforms |
| ECE 464 "Power Electronics" | Active PFC design | L5: boost PFC |

## Michigan

| Course | Topic | This Module |
|--------|-------|-------------|
| EECS 418 "Power Electronics" | Automotive PFC (EV chargers) | L7: pq_assess_ev_charger |
| EECS 463 "Power Systems" | Grid PQ, harmonics | L6: IEEE 519 compliance |

## Georgia Tech

| Course | Topic | This Module |
|--------|-------|-------------|
| ECE 4330 "Power Electronics" | PFC topologies, efficiency | L5: pfc_boost algorithms |
| ECE 6320 "Power Systems" | PQ, harmonics | L7: data center PQ |

## TU Munich

| Course | Topic | This Module |
|--------|-------|-------------|
| Power Electronics | European PFC standards (IEC) | L7: IEC 61000-3-2/4/12 |
| Electric Power Systems | Three-phase unbalanced analysis | L3: symmetrical components |

## ETH Zurich

| Course | Topic | This Module |
|--------|-------|-------------|
| 227-0517 "Power Electronics" | High-frequency PFC | L5: boost PFC with ACMC |
| 227-0526 "Power Systems" | PQ monitoring, flicker | L4: IEC 61000-4-15 Pst/Plt |

## Tsinghua University (清华)

| Course | Topic | This Module |
|--------|-------|-------------|
| 电力电子技术 (Power Electronics) | PF correction fundamentals | L1-L5: core PF calculations |
| 电力系统分析 (Power System Analysis) | Three-phase, sequence components | L3: Fortescue transform |
| 信号与系统 (Signals and Systems) | Fourier analysis for power | L3: FFT, Goertzel |
