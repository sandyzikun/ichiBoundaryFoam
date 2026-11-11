# `ichiBoundaryFoam`：IBM-aware pressure projection 技术预研

日期：2026-11-08

## 本次追加尝试的结论

今晚已经为 IBM-aware pressure projection 建立了一个**默认关闭、双重隔离、
失败时关闭（fail-closed）**的开发框架，但尚未把它写成一个可以求解的数值
后端，也没有编译或运行新代码。

这一区分非常重要：当前代码已经给出了后端的配置入口、PIMPLE 调用位置、
目标块矩阵及未来求解器的边界；它还没有实现 immersed constraint operator、
force-spreading operator 及 pressure--force Schur complement。因此不能把今晚的
改动称作“IBM-aware projection 已经可用”。

现有 binary、正在运行的 V1/R1 及日常 IBM 主线均未受到影响。

## 为什么需要 IBM-aware projection

当前 `ichiBoundaryFoam` 的默认步骤可概括为：

1. 从 momentum equation 得到预测速度 `HbyA`；
2. MDF 把 `HbyA` 拉回 prescribed immersed-body velocity；
3. 常规 PIMPLE pressure projection 再把速度投影到无散空间；
4. pressure correction 没有感知 immersed velocity constraint，因而可能重新
   引入法向或切向 slip。

这正是 V1/R1 中观察到的现象：predictor MDF 后的 cell-centred residual 可以
接近机器零，但 regular pressure projection 后 residual 再次出现。PPMDF 在
每轮 forcing 后重复普通 pressure correction，可以缓解 slip，却没有同时满足
两个约束，而且旧 V1 显示其 wall-clock cost 约为默认路线的 `9.3` 倍。

IBM-aware projection 的目标不是“再多解一次 Poisson equation”，而是把

* 不可压缩条件；
* immersed boundary 的速度条件；

视作需要在同一次 correction 中满足的两个约束。

## 文献中的基本思想

Taira 与 Colonius 把 pressure 和 immersed-boundary force 都解释为约束的
Lagrange multipliers：pressure 消除 velocity divergence，boundary force 消除
immersed surface 上不满足 prescribed velocity 的分量。其方法通过 interpolation
和 regularisation operators 把两种约束加入 projection，并保持与传统
fractional-step method 相近的代数结构。该工作直接说明，真正的 IBM-aware
projection 与“forcing 后再做一次普通压力修正”并不相同。

Kallemov 等人的 rigid-body IBM 进一步采用 unsplit constrained formulation，
联立求解 fluid velocity、pressure 和 Lagrange-multiplier force。他们指出，严格
满足 Lagrangian no-slip constraint 需要一个更大的线性系统；高效求解依赖针对
Schur complement 的近似与预条件，而不是免费获得。

Zhang、Guy 与 Philip 展示了另一种组织方式：先通过 Schur complement 消去
Lagrangian variables，得到 Eulerian saddle-point system，再用 projection
preconditioner 和 Krylov method 求解。这表明“IBM-aware”描述的是约束耦合，
并不强制只能使用某一种完整矩阵存储方式。

Goncharuk、Oshri 与 Feldman 近年把相同思想嵌入 SIMPLE：incompressibility 与
no-slip constraints 被视为 coupled distributed Lagrange multipliers，并组合成
pressure--body-force system，再用 Schur-complement approximation 分解。这项
工作对 OpenFOAM 类 pressure--velocity framework 尤其有启发，但其离散与求解器
不能未经验证直接移植到当前两相 PIMPLE 程序。

## 目标离散系统

今晚在 `src/ibmAwarePEqn.H` 中记录的目标 correction system 为

$$
\begin{bmatrix}
A & G & -S \\
D & 0 & 0 \\
J & 0 & 0
\end{bmatrix}
\begin{bmatrix}
\delta u \\
\delta p \\
\lambda
\end{bmatrix}
=
\begin{bmatrix}
r_m \\
r_c \\
r_b
\end{bmatrix}.
$$

其中：

* $A$ 是 momentum operator，初期也可能采用由 `rAU` 形成的近似；
* $G$ 与 $D$ 分别是 discrete gradient 和 divergence；
* $J$ 把 Eulerian velocity 映射到 immersed constraint space；
* $S$ 把 Lagrange multiplier force 反馈到 Eulerian momentum equation；
* $\lambda$ 是为满足 immersed velocity condition 所需的约束力；
* $r_c$ 对应质量守恒 residual，$r_b$ 对应 immersed velocity residual。

对 no-slip，$J$ 约束三个速度分量；对 free-slip，$J$ 只应约束法向分量。
离散层面还必须定义 $J$ 与 $S$ 在体积加权内积下的关系，否则 coupled system
可能失去对称性、动量一致性或良好的 conditioning。

## 与 PI 提到的连续压力路线是否冲突

不冲突。这里没有把 immersed surface 创建成 pressure equation 的内部
Dirichlet/Neumann patch，也没有要求固体内压力不连续。Pressure 仍然可以作为
整个 Eulerian domain 上的连续未知量；区别在于 pressure correction 和 IBM force
不再彼此无知，而是在同一个约束系统或其 Schur-complement reduction 中协调。

因此，“不把 IBM 当作 pressure internal boundary”与“IBM-aware projection”
可以同时成立。前者描述 pressure field 的空间处理，后者描述 pressure constraint
与 velocity constraint 的代数耦合。

## 为什么当前程序一开始没有采用它

当前实现从 `interFoam`-like PIMPLE 主线逐步加入 SDF、Eulerian cell mask、
direct forcing 和 MDF，优先目标是先得到能够运行的两相 IBM baseline。这个路线
天然复用 OpenFOAM 的 scalar pressure matrix，却没有显式的 Lagrangian force
unknown，也没有独立的 $J$、$S$ operators。

真正的 constrained projection 会引入以下工作：

* 定义 forcing cells/markers 上的 constraint degrees of freedom；
* 构造并验证 $J$、$S$，包括 MPI processor boundaries；
* 把 variable-density VoF pressure operator 纳入 block/Schur system；
* 处理 PIMPLE outer correctors、non-orthogonal corrections、MRF 与 moving body；
* 设计 Krylov solver 和 preconditioner；
* 重新定义 force、slip 与 continuity diagnostics。

这不是在 `pEqn.H` 中增加一个 mask 就能可靠完成的局部修改。此前先采用普通
projection 并保留 PPMDF 作为对照是合理的分阶段开发选择。

## 今晚建立的代码框架

### 编译期隔离

`Makefile` 新增：

```makefile
ENABLE_IBM_AWARE_PROJECTION ?= 0
```

默认值为 `0`。普通 `make` / `make ibm` 不会包含实验调用路径。将来只有显式
使用下列形式重新构建实验 binary，代码才会具备运行期入口：

```bash
make ibm ENABLE_IBM_AWARE_PROJECTION=1
```

今晚没有执行这条命令，也没有执行任何编译命令。

### 运行期隔离

若未来专门编译了实验后端，还需在 `constant/ichiBoundaryGeometry` 中显式设置：

```foam
ibmAwareProjection
{
    enabled                 false;
    formulation             coupledLagrangeMultiplier;
    preconditioner          schurApproximation;
    constraintTolerance     1e-8;
    maxIterations           50;
    writeDiagnostics        false;
}
```

`enabled false` 仍走当前已验证的 regular `pEqn.H`。这形成了“编译期默认关闭 +
运行期默认关闭”的双重保险。

### 新增的源文件

* `src/ibmAwareProjection.H`：读取和校验实验配置，检查 IBM 是否启用；
* `src/ibmAwarePEqn.H`：记录目标 saddle-point system，并占据未来 coupled
  correction 的准确调用位置；
* `src/main.cpp`：只有两个开关均启用时才进入上述实验路径；
* `Makefile`：增加编译开关和 `show-config` 输出。

当前 `ibmAwarePEqn.H` 会明确触发 fatal error，说明 coupled Schur solve 尚未
实现。这样做是有意的：实验路径不得悄悄退化为普通 Poisson correction，也
不得让未完成的代码产生看似合理但物理含义错误的结果。

## 后续开发顺序

### P0：定义与验证 operators

* 明确 current Eulerian forcing cells 是否继续作为 constraint DOFs，或引入
  surface markers；
* 分别定义 no-slip 与 free-slip 的 $J$；
* 定义体积加权 $S$ 并做 discrete-adjoint test；
* 在 serial 与 decomposed mesh 上验证 operator identity。

### P1：只做 diagnostics，不改解

* 在 regular projection 前后计算 $J u-U_b$；
* 输出 divergence residual 与 constraint residual 的空间分布；
* 估计 pressure--force Schur operator 的规模与 conditioning；
* 与现有 detailed SDF-band slip sampling 对照。

### P2：单相静止圆柱 prototype

* 暂时去掉 VoF、wave、turbulence 和 moving geometry；
* 实现 matrix-free $J A^{-1}S$ 与 pressure coupling；
* 比较 regular projection、PPMDF、IBM-aware projection；
* 验收 mass residual、normal/tangential slip、force 与 grid convergence。

### P3：恢复两相与运动

依次恢复 variable density、静止两相 cylinder、Stokes II wave、prescribed body
motion，最后才考虑 `kOmegaSST`。每一步都应保留 body-fitted baseline。

## 主要风险与代价

* coupled system 多出 boundary-force unknowns；
* Schur operator 通常是非局部的，application 可能包含 momentum/Poisson
  solves；
* conditioning 取决于 marker spacing、regularisation 与 operator pairing；
* OpenFOAM collocated finite volume、variable-density VoF 和 non-orthogonal
  correction 与多数 staggered-grid 文献并不完全相同；
* parallel interpolation/spreading 会增加通信；
* 若 preconditioner 不好，成本可能仍显著高于当前 regular projection。

因此不能预先承诺它一定比 PPMDF 更便宜。它的核心优势是数学上同时处理
两个约束，并为系统性降低 slip 提供正确框架；性能收益需要后续实现和测量。

## 当前可以怎样表述

可以说：

> 我们已确认现有 post-projection slip 来自 regular pressure projection 与
> immersed velocity forcing 的非交换性，并依据 constrained/projection IBM
> 文献建立了一个完全隔离的 coupled pressure--force 后端框架。当前主线结果
> 未被改变；下一步先实现和验证 constraint/interpolation/spreading operators，
> 再在单相静止圆柱上开发 Schur-complement prototype。

不能说：

> IBM-aware pressure projection 已经实现或通过验证。

## 参考文献

* Taira, K. & Colonius, T. (2007), [The immersed boundary method: A
  projection approach](https://doi.org/10.1016/j.jcp.2007.03.005),
  *Journal of Computational Physics*, 225, 2118--2137.
* Kallemov, B., Bhalla, A. P. S., Griffith, B. E. & Donev, A. (2016),
  [An Immersed Boundary Method for Rigid
  Bodies](https://doi.org/10.2140/camcos.2016.11.79),
  *Communications in Applied Mathematics and Computational Science*,
  11, 79--141.
* Zhang, Q., Guy, R. D. & Philip, B. (2014),
  [A Projection Preconditioner for Solving the Implicit Immersed Boundary
  Equations](https://doi.org/10.4208/nmtma.2014.1304si),
  *Numerical Mathematics: Theory, Methods and Applications*.
* Goncharuk, K., Oshri, O. & Feldman, Y. (2023),
  [The immersed boundary method: A SIMPLE
  approach](https://doi.org/10.1016/j.jcp.2023.112148),
  *Journal of Computational Physics*, 487, 112148.
