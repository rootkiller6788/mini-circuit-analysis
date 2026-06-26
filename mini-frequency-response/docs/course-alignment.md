# Course Alignment — mini-frequency-response

## MIT 6.003 Signal Processing

| Chapter | Topic | Module Coverage |
|---------|-------|----------------|
| Ch. 6 | Frequency Response of LTI Systems | `freq_response_t`, `tf_compute_freq_response()` |
| Ch. 7 | CT Fourier Transform | `bode_compute()`, frequency ↔ time duality |
| Ch. 9 | The Laplace Transform | `tf_polynomial_t`, `tf_evaluate()` |
| Ch. 10 | System Functions | `tf_to_pole_zero()`, `tf_dc_gain()`, `tf_hf_gain()` |
| Ch. 11 | Feedback Systems | `tf_feedback()`, `stability_margins()` |

## Stanford EE102A Signal Processing and Linear Systems

| Topic | Module Coverage |
|-------|----------------|
| Bode plots | `bode_compute()`, `bode_asymptotic()`, `bode_gain_phase_relation()` |
| Stability analysis | `stability_routh_hurwitz()`, `stability_nyquist_check()` |
| Root locus | `stability_root_locus()` |
| Feedback control | `tf_feedback()`, gain/phase margins |

## Berkeley EE16B Designing Information Devices and Systems II

| Topic | Module Coverage |
|-------|----------------|
| RLC circuits | `resonance_series()`, `resonance_parallel()`, step response |
| Frequency response | `freq_response_t`, `bode_plot_t` |
| Transfer functions | `tf_polynomial_t`, pole-zero analysis |
| Bode plots | `bode_compute()`, asymptotic approximations |

## Berkeley EE105 Microelectronic Circuits

| Topic | Module Coverage |
|-------|----------------|
| Amplifier frequency response | Transfer function poles, GBWP via `bode_gain_bandwidth_product()` |
| Feedback and stability | `stability_margins()`, `tf_feedback()`, stability criteria |
| Active filters | `filter_sallen_key_lp()`, `filter_mfb_lp()`, `filter_tow_thomas()` |

## Stanford EE247 Analog-Digital Interface Circuits

| Topic | Module Coverage |
|-------|----------------|
| Continuous-time filters | `filter_butterworth_prototype()`, `filter_chebyshev1_prototype()` |
| Switched-capacitor filters | `stability_sc_filter()` |
| Filter sensitivity | `filter_sensitivity()` |

## ETH 227-0427 Signal Processing

| Topic | Module Coverage |
|-------|----------------|
| Frequency analysis | `freq_response_t`, complete sweep analysis |
| Filter theory | All `filter_*` functions, prototype generation |
| Bode diagrams | `bode_compute()`, `bode_asymptotic()` |

## ETH 227-0455 High-Frequency Engineering

| Topic | Module Coverage |
|-------|----------------|
| S-parameters | `s_params_2port_t`, `network_s_params_from_z()` |
| Impedance matching | `network_vswr()`, `network_input_impedance()` |
| Two-port stability | `network_rollett_k()`, `network_max_gain()` |

## Georgia Tech ECE 6350 Applied Electromagnetics

| Topic | Module Coverage |
|-------|----------------|
| Network parameters | `network_param_convert()`, Z/Y/ABCD/S conversion |
| Microwave network analysis | `s_params_2port_t`, S-parameter computation |

## Illinois ECE 310 Digital Signal Processing

| Topic | Module Coverage |
|-------|----------------|
| Analog filter prototypes | All `filter_*_prototype()` functions |
| Frequency transformations | `filter_lp_to_hp/bp/bs()` |

## Michigan EECS 411 Microwave Circuits

| Topic | Module Coverage |
|-------|----------------|
| S-parameter analysis | `network_s_params_from_z()`, `network_rollett_k()` |
| Filter synthesis | `filter_g_values_butterworth()`, `filter_g_values_chebyshev()` |

## TU Munich Signal Processing

| Topic | Module Coverage |
|-------|----------------|
| System theory in frequency domain | Transfer functions, pole-zero analysis |
| Filter design | All approximation methods (Butterworth, Chebyshev, Bessel, Elliptic) |

## 清华 信号与系统 / 通信原理

| Topic | Module Coverage |
|-------|----------------|
| 频域分析 | `freq_response_t`, frequency sweep and analysis |
| 传递函数 | `tf_polynomial_t`, complete TF manipulation |
| 波特图 | `bode_compute()`, `bode_asymptotic()` |
| 谐振分析 | `resonance_series()`, `resonance_parallel()` |
| 稳定性判据 | `stability_routh_hurwitz()`, `stability_nyquist_check()` |
