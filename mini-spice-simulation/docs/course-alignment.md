# Course Alignment — mini-spice-simulation

## Nine-School Curriculum Mapping

### MIT
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| 6.003 Signal Processing | Laplace/Fourier in AC analysis | Complex MNA solve |
| 6.450 Digital Comm | Circuit simulation for RF front-end | AC analysis, S-parameter foundation |
| 6.630 EM Waves | Lumped-element circuit models | MNA formulation |

### Stanford
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| EE102A Signal Processing | Time-domain analysis | Transient analysis with companion models |
| EE359 Wireless | RF circuit modeling | AC frequency sweep |
| EE247 Optical | Device modeling | Diode/MOSFET/BJT device models |

### Berkeley
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| EE16A/B Circuits | Nodal analysis, KCL/KVL | MNA with KCL/KVL stamping |
| EE105 Analog IC | Device models, SPICE | Diode (Shockley), MOSFET L1, BJT Ebers-Moll |
| EE117 EM | Numerical methods | LU decomposition, sparse matrices |
| EE123 DSP | Frequency response | AC analysis + Bode export |

### Illinois
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| ECE 310 DSP | Filter analysis | RLC bandpass AC example |
| ECE 459 Comm | System simulation | Full SPICE pipeline |
| ECE 451 EM | Numerical EM | Matrix solver infrastructure |

### Michigan
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| EECS 351 DSP | Transient response | RC/RLC transient analysis |
| EECS 455 Comm | Circuit-system co-design | SPICE netlist interface |
| EECS 411 Microwave | S-parameter analysis | AC complex matrix solve (foundation) |

### Georgia Tech
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| ECE 4270 DSP | Time-frequency analysis | AC + TRAN analysis |
| ECE 6601 Comm | System simulation | Full SPICE engine |
| ECE 6350 EM | Computational EM | LU decomposition, sparse matrices |

### TU Munich
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| Signal Processing | Numerical methods | Newton-Raphson, LU |
| Communications | Circuit modeling | MNA, companion models |
| High-Frequency Eng | Device characterization | Nonlinear device models |

### ETH Zürich
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| 227-0427 Signal Processing | Frequency-domain methods | AC analysis |
| 227-0436 Comm | System simulation | SPICE pipeline |
| 227-0455 EM | Matrix methods | Dense LU, condition number |

### Tsinghua (清华)
| Course | Topic | Module Coverage |
|--------|-------|-----------------|
| 信号与系统 | Time/frequency analysis | TRAN + AC |
| 通信原理 | Circuit simulation | Full SPICE flow |
| 电磁场 | Field-circuit coupling | MNA matrix formulation |
| 数字信号处理 | Filter design | RLC bandpass AC example |

---

## Reference Textbooks

| Textbook | Author(s) | Year | Module Use |
|----------|-----------|------|------------|
| SPICE2: A Computer Program to Simulate Semiconductor Circuits | Nagel | 1975 | Core architecture |
| Inside SPICE | Kielkowski | 1998 | MNA stamping, companion models |
| Computer Methods for Circuit Analysis and Design | Vlach & Singhal | 1994 | AC analysis, complex matrix |
| Computer-Aided Analysis of Electronic Circuits | Chua & Lin | 1975 | Transient analysis, integration |
| Matrix Computations | Golub & Van Loan | 2013 | LU decomposition, condition number |
| Microelectronic Circuits | Sedra & Smith | 2020 | Device equations |
| Semiconductor Device Modeling with SPICE | Massobrio & Antognetti | 1998 | MOSFET, BJT models |
| Iterative Methods for Sparse Linear Systems | Saad | 2003 | CSR sparse matrix |
