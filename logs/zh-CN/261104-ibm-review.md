# `ichiBoundaryFoam` 当前进展与主流 IBM 路线复审

## 结论先行

目前的 `ichiBoundaryFoam` 已经完成了以下基础工作：

* 在 IBM 关闭时，基本复现 OpenFOAM-com `v2312` 的
  `interFoam` VoF、PIMPLE/PISO 与 Stokes wave 路径；
* 使用解析 Signed Distance Function（SDF）描述静止或预设运动的
  圆柱；
* 使用 SDF 生成平滑的 solid fraction，并在动量预测速度上施加
  direct forcing；
* 在完整背景域上求解连续压力，并增加 multi-direct forcing（MDF）
  与压力增量投影；
* 单相全水域测试中，静止和匀速平移圆柱均达到约
  `slip = O(10^-5) m/s`、`divergence = O(10^-12) 1/s`；
* 移动圆柱测试中，covered/uncovered volume 基本配平。

但当前实现还不能称为经过验证的两相 WSI/FSI solver。最新的
Stokes II 对照实验显示：

| 指标 | 无 IBM 对照 | 静止 IBM 圆柱 |
| --- | ---: | ---: |
| 模拟时间 | 9 s | 9 s |
| 计算期间最大 `alpha.water` | 1.00280 | 1.39537 |
| 末时刻最大 `alpha.water` | 1.00189 | 1.01539 |
| 末时刻 IBM slip | 不适用 | 1.80252e-2 m/s |
| 末时刻最大散度 | 不适用 | 1.77479e-8 1/s |

因此，MDF 已经能把投影后的离散散度压低，但没有同时解决两相界面
附近的相分数、密度梯度和 no-slip 约束。当前版本可以作为研究原型，
尚不能用于可信的 wave load、双向 6DoF FSI 或最终 benchmark。

## 当前算法实际是什么

当前圆柱由解析 SDF 表示。根据到固体表面的有符号距离，在厚度为
`transitionWidth` 的区域内构造平滑 solid fraction：

$$
0\leq\chi_s\leq1.
$$

预测速度通过下式向结构速度靠近：

$$
\mathbf{u}^{*}\leftarrow
\mathbf{u}^{*}+\chi_s
(\mathbf{u}_b-\mathbf{u}^{*}).
$$

随后仍在整个背景域，包括虚构的固体内部流体，求解压力 Poisson
方程。标准压力投影会再次改变已经受约束的速度，因此 MDF 继续重复：

```text
速度约束
→ 从约束速度构造通量
→ 全域压力增量投影
→ 同步修正 U、phi 和 p_rgh
→ 检查 slip
```

这一做法属于 Eulerian、diffuse-interface、discrete/direct forcing
的一种简化实现。它借鉴了 `sdfibm` 的 SDF 分类思想，但并不等同于
完整的 `sdfibm`，也不等同于 Ye 等人的二阶 discrete-forcing IBM：

* 当前没有在真实 IB 表面构造插值点或二阶 Laplacian reconstruction；
* 当前 forcing 由 cell-centred SDF 权重直接混合，不是严格的表面
  Lagrangian forcing；
* 当前没有对 `alpha.water` 向固体内部作 Neumann-compatible extension；
* 当前没有 IB wall function 或 `k`、`omega`、`nut` 的边界处理；
* 当前压力仍跨越固体内部与流体域连续耦合。

这个区别很重要：当前单相测试验证的是“约束—投影数据通路能够稳定
工作”，并没有证明边界位置、边界层、载荷或空间精度已经达到成熟
DF-IBM 的水平。

## 主流 IBM implementations 如何分类

Mittal 与 Iaccarino 的经典综述把 IBM 大体分为 diffuse forcing 与
sharp-interface 两类；后续 FSI 综述仍沿用这一基本谱系。不存在一种
对所有问题都最优的 IBM，所谓 SOTA 取决于更重视几何灵活性、局部
守恒、壁面精度、运动稳定性还是实现成本。

### Regularised-delta / Lagrangian direct forcing

Peskin–Uhlmann 系方法在结构表面布置 Lagrangian markers，通过离散
delta kernel 在 Lagrangian 与 Eulerian 网格之间插值和传播力。Uhlmann
将其发展为适合刚性颗粒的 direct forcing；Breugem 随后加入 MDF、
marker inward retraction 和内部虚构流体惯性修正，使颗粒边界和载荷
更准确。

优点是：

* 固定背景网格，运动几何处理自然；
* 表面 marker 可直接承担力交换；
* 适合大量颗粒和颗粒分辨 DNS。

典型问题是：

* delta kernel 把边界扩散到数个 cells；
* 有效物体尺寸可能偏大；
* 单次 direct forcing 不能精确满足 no-slip；
* 运动时可能出现力振荡；
* 载荷必须修正虚构固体内部流体的动量。

我们的平滑 SDF forcing 与这一类在“diffuse forcing”上相近，但没有
Lagrangian markers、delta kernel、marker retraction 和成熟的内部
流体惯性修正，因此不能直接继承 Breugem 方法的二阶精度结论。

### Eulerian SDF / volume forcing

`sdfibm` 使用 SDF 在 OpenFOAM 中识别流体、边界和固体区域，并通过
discrete forcing 施加结构速度。SDF 对解析几何和运动刚体很方便，
也是当前项目最自然的起点。

近年的 volume-of-solid 或 implicit volume-penalty 路线进一步把固体
体积分数当作守恒变量，并将固体动量隐式嵌入统一动量方程。这类方案
的目标是减轻 moving body 的 fresh/dead cell、质量守恒和显式 penalty
刚性问题。它们并不意味着简单的平滑 `chi_s` 自动具有守恒性；关键
在于 solid fraction 的输运、几何守恒律和动量方程必须一致。

### Ghost-cell / sharp-interface IBM

ghost-cell 方法在固体内建立 ghost values，通过真实边界上的 Dirichlet
或 Neumann 条件重构临近流体 stencil。Tseng 与 Ferziger 的方法是这一
路线的经典实现。

优点是：

* 边界位置较清楚，通常比 diffuse forcing 更锐利；
* 容易得到二阶局部空间精度；
* 压力、速度和湍流变量可分别施加合适边界条件。

困难是：

* moving boundary 会不断改变 ghost stencil；
* fresh/dead cells 会引入质量和压力振荡；
* 基本 ghost-cell 方法并不自动全局守恒；
* 三维复杂几何的鲁棒 reconstruction 成本较高。

### Cut-cell / conservative immersed interface

cut-cell 方法按真实固体表面切割控制体，并在截断面上计算通量。它最
符合有限体积法的局部守恒逻辑。Seo 与 Mittal 指出，moving sharp
interface 的虚假压力振荡主要来自违反 geometric conservation law，
并用局部 conservative cut-cell 处理降低该问题。

优点是局部质量、动量守恒和清晰的壁面应力；缺点是小 cut cells、
时间步限制、拓扑更新和三维实现复杂。对于当前规则圆柱，它可能提供
更好的 benchmark，但若目标是尽快形成可修改的 OpenFOAM 原型，直接
全面转向 cut-cell 会显著扩大开发范围。

## 主流方案如何处理我们正在遇到的问题

### 问题一：压力投影破坏 no-slip

当前连续压力方案允许 pressure Poisson stencil 跨越 fluid–solid 区域。
其好处是简单、稳定，也符合 PI 提到的早期实现经验；但投影会改变
forced velocity，因此出现“散度很小而 slip 仍较大”的现象。

MDF 是传统补救措施。Wang 等提出 combined MDF，通过重复插值与
forcing 改善 no-slip；Breugem 的方法也使用 MDF。不过成熟的 MDF
通常重点迭代 boundary forcing，不一定每轮都完整重复昂贵的全域
Poisson solve。我们的版本每轮增加一次全域压力增量投影，稳定但成本
较高，而且两相工况下仍不能保证 slip 收敛。

2024 年 Yildiran 等专门研究了 IBM 的压力边界。他们指出，仅在动量
方程施加 forcing 会带来两类问题：投影后的 slip，以及 Poisson 方程
使固体和流体域耦合。其方案以 ghost-fluid 思想在 immersed surface
附近施加压力条件并实现 pressure decoupling。

因此建议保留两个可切换后端：

* `continuousPressure`：保留当前全域连续压力，作为 PI 路线和低成本
  baseline；
* `pressureDecoupled`：在 IB 附近修正 pressure stencil 或施加等价
  homogeneous Neumann condition，作为高精度路线。

不建议仅靠无限增加 MDF correctors 掩盖 pressure coupling 问题。

### 问题二：两相 `alpha.water` 超界与密度梯度错误

这是当前最明确的缺项。Ye、Chen 与 Maki 的 OpenFOAM 两相 DF-IBM
指出：固体内部流体是虚构的，若 `alpha` 未合理延拓，空气—水约
`1000:1` 的密度比会在 IB 附近产生错误密度梯度，污染浮力项和动量
forcing。

他们采用逐层 dilation/field extension：从邻近非固体 cells 以逆平方
距离插值 `alpha` 到固体内部，通常延拓约五层。其目的不是把固体内部
体积当成真实水量，而是在 IB 上间接实现

$$
\frac{\partial\alpha}{\partial n}=0,
$$

相当于中性接触角。水量统计必须只计算可接触流体域，不能把 extension
改变的固体内部 `alpha` 计入物理质量。

我们的最新结果——IBM case 的瞬时 `Max(alpha)=1.39537`，而无 IBM
对照峰值仅为 `1.00280`——与这一缺项高度一致。这里是基于代码和文献
作出的推断，仍需在加入 extension 后通过相同对照实验验证。

近期 sharp-interface/VoF 方法通常使用 PLIC-VoF 并在固体边界附近
重构速度、相分数或通量；更保守的 volume-of-solid 方法则让 solid
fraction 与流体输运在统一守恒框架中演化。它们实现更复杂，但特别
适合水线穿过运动物体的 water-entry、floating body 与 breaking-wave
问题。

### 问题三：moving body 的 fresh/dead cells 与质量守恒

当前代码只统计 covered/uncovered volume，并未对新暴露 fluid cells
构造守恒的 `U`、`p_rgh`、`alpha`、`k`、`omega`。单相匀速测试中体积
配平，只能说明 SDF 几何更新平滑，不能证明两相移动边界守恒。

主流处理包括：

* 从流体侧沿法向外推 fresh-cell variables；
* 对相分数作 Neumann-compatible extension；
* 强制满足 geometric conservation law；
* 使用 conservative cut-cell 或 moving control volume；
* 使用守恒输运的 volume-of-solid field。

若后续 monopile 只是小幅、缓慢预设运动，field extension 加局部守恒
修正可能足够；若涉及 water entry/exit、大角度倾倒或双向 FSI，则应
认真考虑 sharp conservative 或 volume-of-solid 路线。

### 问题四：IBM 载荷不可信

当前 `IBM preliminary load` 只是对累计 velocity correction 除以
时间步后的体积分。短时单相测试中，medium IBM 与 body-fitted 的
`Fx` 相差约 4.8 倍，因此还不能用于结构方程。

成熟方法有三类：

* 对 Lagrangian forcing 直接积分，并减去固体内部虚构流体的惯性；
* 在重构的真实表面上积分 pressure 与 viscous stress；
* 使用包围物体的 moving control volume，由动量收支计算总力和力矩。

Ye 等使用包围 IB 的网格面投影 pressure 与 viscous stress；Breugem
路线显式考虑内部虚构流体惯性；moving-control-volume 方法则避免直接
依赖噪声较大的表面梯度。对当前项目，建议首先实现两套相互校验的
载荷：forcing-volume 加 internal-fluid correction，以及 surface-stress
reconstruction。二者与 body-fitted force 达到网格收敛一致后，才可
接入双向 6DoF。

### 问题五：高 Reynolds 数与湍流模型

当前只强制 `U`，未给 `k`、`omega`、`nut` 构造 immersed wall
condition。对于 monopile 的 k-omega SST RANS，这会使近壁湍流生产、
耗散与 wall shear 缺乏物理意义。

Ye 等的 turbulent DF-IBM 使用专门的 IB wall function；其两相工作
采用 Spalart–Allmaras IB wall treatment，并明确指出 DES/LES 在分离流
方面可能更准确。对我们而言，不能直接把 body-fitted SST wall
function 原样套到 cell-centred diffuse layer。需要：

* 由 SDF 提供 wall distance 与法向；
* 重构壁面切向速度和 wall shear；
* 为 `k`、`omega` 或 `nut` 设计 IBM wall condition；
* 用 channel flow、圆柱绕流和 body-fitted monopile 分级验证。

## 当前方案与 SOTA 的相对位置

“SOTA”并不是某个唯一软件，而是以下能力的组合：

| 能力 | 当前实现 | 成熟/SOTA 方向 |
| --- | --- | --- |
| 几何 | 解析圆柱 SDF | triangle surface、level set 或守恒 VoS |
| 速度边界 | 平滑 cell forcing | 二阶 reconstruction、implicit/Lagrange multiplier forcing |
| 压力 | 全域连续 Poisson | pressure decoupling 或 conservative cut-cell |
| 两相相分数 | 无 IB extension | Neumann-compatible extension 或 sharp PLIC flux |
| Moving conservation | 仅诊断 swept volume | GCL、fresh-cell extension、conservative VoS/cut-cell |
| 载荷 | forcing 调试积分 | 内部流体修正、surface stress、moving CV |
| 湍流 | 未处理 | IB wall model、WMLES/DES 或解析壁面 |
| FSI | prescribed one-way motion | added-mass-stable partitioned/monolithic coupling |

当前路线的价值在于：它与 `interFoam`、OpenFOAM wave models 和未来
DPM 的整合成本低，适合快速形成可读、可改的研究框架。它的主要风险
不是 SDF 本身，而是把 SDF 平滑体积分数、连续压力、VoF 相分数和载荷
积分分别实现，却尚未形成一个离散一致且守恒的整体。

## 建议的实施顺序

### 第一优先级：修复 static two-phase IBM

1. 实现 `alpha.water` 向固体域逐层延拓，可配置层数，先复现 Ye 等的
   inverse-distance dilation。
2. 延拓后再更新 `rho` 与两相动量相关量；物理水量只在
   `1-solidFraction` 加权域统计。
3. 增加分区诊断：fluid、forcing layer、solid、free-surface 与 wave
   inlet 分别统计 `alpha`、slip 和 divergence。
4. 重跑完全相同的 `9 s` Stokes II control/IBM 对照。最低通过条件：
   IBM 不显著放大 control 的 `alpha` 超界，accessible phase volume
   无非物理漂移，forcing-layer slip 随迭代单调下降。

### 第二优先级：提高 boundary accuracy

1. 用基于 SDF 法向的 image point / second-order interpolation 代替单纯
   cell blending。
2. 检查 `transitionWidth/delta` 对有效圆柱半径和载荷的影响。
3. 将 MDF 改为基于 slip 与 divergence 的联合停止准则，减少不必要的
   Poisson solves。
4. 增加可切换的 pressure-decoupled 原型，与 continuous-pressure
   baseline 比较。

### 第三优先级：正式载荷与 moving body

1. 实现 forcing load 的内部虚构流体惯性修正；
2. 实现 surface-stress reconstruction 或 moving control volume；
3. 与同分辨率 body-fitted cylinder 比较 `Fx`、`Fy`、相位与涡脱落；
4. 为 fresh cells 实现 `alpha/U/p` extension 与守恒修正；
5. 通过 prescribed translation/rotation 后，再接入双向 6DoF。

### 第四优先级：monopile RANS 与 DPM

1. 增加 SST-compatible IB wall treatment；
2. 比较 body-fitted 与 IBM 的 `yPlus`、wall shear、`k`、`omega`；
3. 待流体—结构载荷验证后，再加入 Lagrangian sediment particles；
4. DPM 的 fluid–particle coupling 应使用可接触流体体积，避免把 IBM
   固体内部虚构流体计入孔隙率或颗粒受力。

## Go / No-Go 判断

当前结论不是停止 IBM，也不是立即改写为 cut-cell。合理判断是：

* **Go**：继续以 SDF + DF-IBM 作为第一条 OpenFOAM 实现路线；
* **Go with conditions**：优先补齐 `alpha` extension、界面 reconstruction
  与正式载荷，不再把单相低 slip 当作两相验证通过；
* **No-Go for now**：不得使用当前 preliminary force 驱动双向 FSI；
* **转向条件**：若完成相分数延拓和二阶界面 reconstruction 后，moving
  two-phase case 仍持续出现质量/压力振荡，或载荷无法随网格收敛到
  body-fitted 结果，则转向 pressure-decoupled ghost-cell、conservative
  cut-cell 或 volume-of-solid backend。

## 主要参考文献

* Mittal, R. & Iaccarino, G. (2005),
  [Immersed Boundary Methods](https://doi.org/10.1146/annurev.fluid.37.061903.175743).
* Uhlmann, M. (2005),
  [An immersed boundary method with direct forcing for the simulation of particulate flows](https://doi.org/10.1016/j.jcp.2005.03.017).
* Tseng, Y.-H. & Ferziger, J. H. (2003),
  [A ghost-cell immersed boundary method for flow in complex geometry](https://doi.org/10.1016/j.jcp.2003.07.024).
* Wang, Z. et al. (2008),
  [Combined multi-direct forcing and immersed boundary method for simulating flows with moving particles](https://doi.org/10.1016/j.ijmultiphaseflow.2007.10.004).
* Seo, J. H. & Mittal, R. (2011),
  [A sharp-interface immersed boundary method with improved mass conservation and reduced spurious pressure oscillations](https://doi.org/10.1016/j.jcp.2011.06.003).
* Kempe, T. & Fröhlich, J. (2012),
  [An improved immersed boundary method with direct forcing for the simulation of particle laden flows](https://doi.org/10.1016/j.jcp.2012.01.021).
* Breugem, W.-P. (2012),
  [A second-order accurate immersed boundary method for fully resolved simulations of particle-laden flows](https://doi.org/10.1016/j.jcp.2012.02.026).
* Ye, H., Chen, Y. & Maki, K. (2020),
  [A Discrete-Forcing Immersed Boundary Method for Moving Bodies in Air–Water Two-Phase Flows](https://doi.org/10.3390/jmse8100809).
* Zhang, C. et al. (2020),
  [sdfibm: a signed distance field based discrete forcing immersed boundary method in OpenFOAM](https://doi.org/10.1016/j.cpc.2020.107370).
* Yildiran, I. N. et al. (2024),
  [Pressure boundary conditions for immersed-boundary methods](https://doi.org/10.1016/j.jcp.2024.113057).
* Tsetoglou, A. et al. (2025),
  [A volume-of-solid implicit volume penalty method for moving-body flows](https://doi.org/10.1002/fld.5334).
* Cao, J. (2025), *CFD-DPM Modelling of Wave-Driven Sediment Transport
  and Scour Around Offshore Structure Foundations*, PhD thesis, University
  of Liverpool（本地项目参考文献）。

