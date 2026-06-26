# 1. mini-circuit-analysis

**电路分析 — 从零构建、零依赖的 C 语言实现，覆盖电网理论、直流/交流稳态、
瞬态动力学、频域方法和 SPICE 兼容仿真。**

---

## 子模块

| 子模块 | 主题 | 核心课程 |
|--------|------|----------|
| [mini-circuit-topology](mini-circuit-topology/) | 图论基础、MNA 列式、基尔霍夫定律、特勒根定理、稀疏求解器 | MIT 6.002, Stanford EE101A |
| [mini-dc-ac-circuit](mini-dc-ac-circuit/) | 欧姆定律、节点/网孔分析、复阻抗、相量法、谐振、耦合电感 | MIT 6.002, Berkeley EE16A |
| [mini-frequency-response](mini-frequency-response/) | 波特图、传递函数、模拟滤波器设计（巴特沃斯/切比雪夫）、稳定性判据 | MIT 6.003, Stanford EE102A |
| [mini-network-theorem](mini-network-theorem/) | 叠加定理、戴维南/诺顿等效、米尔曼定理、最大功率传输、互易定理 | MIT 6.002, Berkeley EE16B |
| [mini-power-factor](mini-power-factor/) | 复功率、功率三角形、谐波与 THD、功率因数校正、相量运算、电能质量 | MIT 6.061, Stanford EE253 |
| [mini-spice-simulation](mini-spice-simulation/) | SPICE 网表解析、MNA 矩阵装配、器件模型（二极管/BJT/MOSFET）、DC/AC/TRAN 分析 | MIT 6.002, Berkeley EE105 |
| [mini-transient-analysis](mini-transient-analysis/) | 一阶/二阶瞬态、状态空间方法、数值 ODE 求解器（Euler/RK4/BDF2/Gear）、开关瞬态 | MIT 6.002, Stanford EE101B |
| [mini-two-port-network](mini-two-port-network/) | Z/Y/H/G/ABCD/S 参数、全部 30 种方向转换、互联、网络综合、稳定性 | MIT 6.003, Stanford EE101B |

---

## 设计哲学

1. **零依赖** — 每个子模块仅需 C99/C11 编译器和标准库即可编译。无需外部数学库或线性代数
   包；所有求解器均为手写实现，确保教学透明性。

2. **分层学习** — 每个子模块遵循 L1→L9 递进：定义 → 核心概念 → 数学结构 → 基本定律 →
   算法 → 典型问题 → 应用 → 高级主题 → 研究前沿。这映射了从本科到研究生的完整课程体系。

3. **从方程到代码** — 每个函数直接对应教科书中的方程。源代码中引用了支配性物理定律
   （欧姆、麦克斯韦、法拉第），读者可将每一行 C 代码追溯到第一性原理。

4. **跨模块可组合** — 共享的类型定义（`circuit_elements.h`、`transient_defs.h`、
   `spice_core.h`）使子模块可以组合。拓扑模块的 MNA 矩阵直接输入到 SPICE 仿真器和
   DC/AC 求解器中，形成一个微型但完整的 EDA 流水线。

---

## 构建

一次性构建所有子模块：

```bash
# 在本目录下
for d in mini-*/; do (cd "$d" && make); done
```

或单独构建各模块：

```bash
cd mini-circuit-topology && make && make test
cd mini-dc-ac-circuit && make && make test
cd mini-frequency-response && make && make test
cd mini-network-theorem && make && make test
cd mini-power-factor && make && make test
cd mini-spice-simulation && make && make test
cd mini-transient-analysis && make && make test
cd mini-two-port-network && make && make test
```

---

## 项目结构

```
1. mini-circuit-analysis/
├── mini-circuit-topology/          # 图论、MNA、基尔霍夫定律、特勒根定理、稀疏LU分解
├── mini-dc-ac-circuit/             # 直流/交流稳态、欧姆定律、节点/网孔分析、复功率
├── mini-frequency-response/        # 波特图、传递函数、滤波器设计、稳定性分析
├── mini-network-theorem/           # 叠加定理、戴维南/诺顿等效、最大功率传输
├── mini-power-factor/              # 复功率、谐波分析、THD、功率因数校正、电能质量
├── mini-spice-simulation/          # 网表解析、MNA装配、器件模型、多模式仿真
├── mini-transient-analysis/        # 一阶/二阶系统、状态空间、数值ODE、开关瞬态
├── mini-two-port-network/          # Z/Y/H/G/ABCD/S参数、30种转换、网络综合、稳定性
├── .gitignore                      # 构建产物和IDE排除规则
├── README.md                       # 本文件（英文版）
└── README-CN.md                    # 中文版
```

---

## 许可证

MIT — 详见各子模块源码头部声明。
