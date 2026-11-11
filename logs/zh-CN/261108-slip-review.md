# Discrete-forcing IBM 中的 slip 控制：文献综述与开发建议

> 日期：2026-11-08  
> 对象：`ichiBoundaryFoam` 的静止／预设运动圆柱 IBM  
> 范围：不可压缩 discrete-forcing IBM，重点讨论 no-slip／no-penetration
> 残差；同时考虑投影法、有限体积、VoF 两相流和运动边界。  
> 资料原则：正文引用均尽量指向原始研究论文，而不是二手综述。

## 结论先行

当前最大的开发难点确实可以概括为 **immersed surface 上 slip 的可控性**，
但必须马上补充一句：这里的 slip 不是一个单一误差，也不能只靠“增加 MDF
次数”统一解决。

文献中至少存在以下几类机制：

* **forcing 算子误差**：插值与铺展并非精确互逆，多个表面 marker 的
  支撑域彼此重叠，一次显式 forcing 后边界速度通常不会恰好等于壁面速度；
* **压力投影造成的 splitting error**：forcing 可令中间速度满足边界条件，
  随后的不可压压力投影却会再次修改速度，重新引入法向或切向残差；
* **diffuse／grid-scale 界面误差**：forcing 作用于有限厚度区域，会产生
  有效直径偏大、近壁梯度被污染和边界位置网格化等问题；
* **局部质量与通量误差**：即便 cell-centred 速度看起来接近目标速度，
  穿过 immersed surface 的面通量仍可能不正确；运动边界还会产生
  fresh/dead cells 和几何守恒误差；
* **时间离散和黏性耦合误差**：显式 forcing、扩散项、时间步长与壁面速度
  约束之间并不总是相容；
* **两相交会误差**：VoF、密度、表面张力和 IBM classification 在
  air–water–solid 交会处不一致时，会把额外压力／通量误差表现成 slip；
* **测量误差**：solid cell 中的最大相对速度、SDF band 中 cell centre 的
  相对速度，以及真实 immersed surface 上插值得到的残差，不是同一个量。

对当前 `ichiBoundaryFoam`，最重要的近期判断是：

1. 默认关闭“每轮 MDF 后重解全域 Poisson 方程”是合理的。经典 MDF 本来
   就是重复 forcing；它不天然要求每次都重做 pressure projection。
2. 预测速度上的 MDF 只能保证 predictor residual 下降，不能保证普通
   PIMPLE projection 后的最终 residual 仍然很小。
3. 如果 predictor residual 随 MDF 次数快速下降、但 projection 后又明显
   回升，主要瓶颈是 **pressure–velocity splitting／flux compatibility**，
   而不是 MDF 次数不足。
4. 如果 predictor residual 本身就形成平台，则应优先修 forcing 算子、
   relaxation、界面位置或局部重构；反复求解压力不是对症手段。
5. `freeSlip` 只要求 no-penetration，因此应以法向残差为主要验收量；总
   slip 和切向 slip 大并不自动表示 free-slip 条件失败。
6. 当前 `maxSlip()` 是固体域 cell-centre 指标，详细 sampling 也是 SDF
   band 中的 cell-centre 指标；在宣布边界条件“通过”以前，还需要真正的
   surface-interpolated residual 和 surface flux 指标。

因此，建议采用的近期 baseline 不是立刻重写整个方法，而是：

> predictor-level MDF + 单次常规 PIMPLE projection + 分阶段 slip 诊断，
> 然后根据误差究竟出现在 forcing 前后还是 projection 前后，选择局部
> forcing correction、pressure-compatible correction 或 sharp-interface
> 重构。

## 1. “Slip”究竟应该怎样定义？

### 1.1 刚体表面相对速度

设真实 immersed surface 为 $\Gamma$，流体速度为 $\boldsymbol{u}$，
刚体在表面位置 $\boldsymbol{X}$ 的速度为 $\boldsymbol{U}_b$，则表面残差为

$$
\boldsymbol{e}_{\Gamma}
=
\mathcal{I}_{\Gamma}\boldsymbol{u}
-
\boldsymbol{U}_b,
$$

其中 $\mathcal{I}_{\Gamma}$ 表示把 Eulerian 速度一致地插值到真实表面，
而不是简单读取最近一个 cell centre。

将其按表面单位法向 $\boldsymbol{n}$ 分解：

$$
e_n
=
\boldsymbol{e}_{\Gamma}\cdot\boldsymbol{n},
$$

$$
\boldsymbol{e}_t
=
\boldsymbol{e}_{\Gamma}
-
e_n\boldsymbol{n}.
$$

于是：

* no-slip/no-penetration 要求 $e_n\rightarrow0$ 且
  $\lVert\boldsymbol{e}_t\rVert\rightarrow0$；
* free-slip/no-penetration 只要求 $e_n\rightarrow0$；切向相对速度可以存在；
* 对不可压流，更直接的渗漏量是表面通量

$$
Q_{\Gamma}
=
\int_{\Gamma}
(\boldsymbol{u}-\boldsymbol{U}_b)
\cdot\boldsymbol{n}\,\mathrm{d}A.
$$

局部 $e_n$ 与整体 $Q_{\Gamma}$ 都应检查：整体积分接近零，仍可能是表面
不同区域的正负渗漏互相抵消。

### 1.2 至少需要四个时间位置的 residual

为了定位误差来源，建议同一时间步记录：

* $e^{*}$：IBM forcing 前的 predictor residual；
* $e^{(k)}$：每一轮 MDF forcing 后的 residual；
* $e^{P}$：普通 pressure projection 后的最终 residual；
* $e^{n+1}$：所有 velocity correction、boundary correction 和 function
  object 执行完毕后，真正写盘前的 residual。

如果只看最后一个 `maxSlip`，无法判断“forcing 没有施加好”还是“forcing
施加好后被投影破坏了”。

### 1.3 最大值不是充分统计量

最大值容易被一个离散位置、相界面 cell 或网格跳变主导。更稳妥的输出应
同时包含：

* $L_{\infty}$：最大值，用于捕捉局部最坏情况；
* $L_2$ 或 area-weighted RMS：用于评估整体边界质量；
* P95/P99：区分“普遍偏差”和“少数离群点”；
* 法向／切向分量；
* liquid／gas／VoF-interface 三个子集；
* 所在表面位置、SDF、$\alpha_\mathrm{water}$ 和局部网格尺度；
* 以特征速度 $U_{\mathrm{ref}}$ 无量纲化后的误差。

## 2. 显式 direct forcing 为什么会留下 slip？

### 2.1 插值—铺展并不互逆

经典 Lagrangian direct-forcing IBM 中，以 $E$ 表示 Eulerian-to-Lagrangian
插值，以 $H$ 表示 Lagrangian-to-Eulerian force spreading。一次显式
forcing 可以抽象写成

$$
\boldsymbol{F}^{(0)}
=
\frac{
\boldsymbol{U}_b-E\boldsymbol{u}^{*}
}{\Delta t},
$$

$$
\widetilde{\boldsymbol{u}}
=
\boldsymbol{u}^{*}
+
\Delta t\,H\boldsymbol{F}^{(0)}.
$$

重新插值到表面后，残差为

$$
\boldsymbol{e}^{(1)}
=
\boldsymbol{U}_b-E\widetilde{\boldsymbol{u}}
=
(I-EH)
(\boldsymbol{U}_b-E\boldsymbol{u}^{*}).
$$

除非 $EH=I$，一次 forcing 不会将误差精确消除。误差来自：

* interpolation 与 spreading 的离散非互易性；
* 相邻 marker kernel 支撑域重叠；
* marker spacing、quadrature weight 和 Eulerian 网格尺度不匹配；
* 曲率、非均匀网格和边界相对网格位置变化。

[Gsell & Favier (2021)](https://doi.org/10.1016/j.jcp.2021.110265)
将这种 boundary-slip error 明确归因于 interpolation/spreading
non-reciprocity，并提出了几乎不增加成本的 force scaling correction。
最新的统一分析则把它命名为 **reciprocity error**，并指出它与
projection-induced splitting error 是两类不同机制
（[2026 年统一分析](https://doi.org/10.1007/s10494-026-00802-7)）。

### 2.2 MDF 在修什么？

MDF 重复计算剩余误差并施加 forcing：

$$
\boldsymbol{e}^{(k+1)}
=
(I-EH)\boldsymbol{e}^{(k)}.
$$

只要相应迭代算子的谱半径小于一，residual 就会下降。Wang、Fan 与 Luo
的 combined MDF 工作证明了它比一次 direct forcing 更好地满足 moving
particle 的 no-slip 条件
（[Wang et al., 2008](https://doi.org/10.1016/j.ijmultiphaseflow.2007.10.004)）。
Breugem 后来把 MDF、marker inward retraction 和有限体积压力修正组合，
并在其粒子测试中发现少量 forcing corrections 就能在精度与成本之间取得
较好平衡
（[Breugem, 2012](https://doi.org/10.1016/j.jcp.2012.02.026)）。

但 MDF 不是万能的：

* 如果 $EH$ 病态或局部 marker 支撑严重重叠，迭代可能很慢或不稳定；
* MDF 通常修的是 forcing/interpolation residual，不会自动修复后续
  pressure projection；
* 单纯增加次数可能在空间离散误差平台以后不再有收益；
* 显式处理黏性扩散时，低 Reynolds number 或较大时间步下仍可能有明显
  residual；
* force 收敛不等于 surface flux、pressure 或 hydrodynamic load 已收敛。

### 2.3 relaxation 与稳定性

[Zhou & Balachandar (2021)](https://doi.org/10.1016/j.jcp.2020.109862)
把 Lagrangian volume weight 解释为控制 velocity error 衰减速率的 relaxation
factor，并指出最有效的选择是稳定条件允许的最大值。这个结果提醒我们：

* “forcing coefficient = 1”不总是最优；
* 更强 forcing 或更多迭代并非越多越好；
* 应画出 residual 对 MDF iteration 的收敛曲线，而不是只比较 1 次和某个
  随意固定次数；
* relaxation／amplification 参数应做稳定性扫描。

## 3. 压力投影为什么会重新制造 slip？

### 3.1 两个约束并不自动相容

不可压求解器通常先得到满足 IBM 条件的中间速度
$\widetilde{\boldsymbol{u}}$，再通过压力投影得到

$$
\boldsymbol{u}^{n+1}
=
P\widetilde{\boldsymbol{u}},
$$

其中 $P$ 把速度投影到离散无散空间。最终边界残差为

$$
\boldsymbol{e}^{P}
=
\boldsymbol{U}_b
-
E P\widetilde{\boldsymbol{u}}.
$$

一般而言，IBM constraint operator 与 pressure projection 不交换。换言之，
“先满足壁面速度，再满足不可压”不等价于“同时满足两者”。

[Domenichini (2008)](https://doi.org/10.1016/j.jcp.2008.03.009)
在尽量排除 interpolation 干扰的简单问题中表明，标准 fractional-step
direct forcing 不必然满足 impenetrability condition。因此，projection 后
法向速度回升是已知的算法问题，并不自动说明圆柱 SDF 或并行计算有 bug。

2026 年的统一分析进一步区分：

* reciprocity error 主要由空间 transfer operators 决定；
* splitting error 来自 pressure 与 IBM multiplier 的顺序解耦；
* 没有 pressure predictor 时，示例中的 splitting residual 呈一阶时间误差；
* 使用上一时刻压力作为 predictor 可把相应阶数提高，但朴素高阶 predictor
  也可能带来稳定性问题；
* 网格加密可以降低空间误差，却不保证消除纯粹的 splitting mismatch。

参见
[On Reciprocity and Splitting Errors in Direct-Forcing Immersed Boundary
Methods (2026)](https://doi.org/10.1007/s10494-026-00802-7)。

### 3.2 为什么“每轮 MDF 后重解 Poisson”有时有效？

交替执行

```text
forcing -> pressure projection -> forcing -> pressure projection -> ...
```

是在迭代追逐两个约束的交集，因此可能同时降低 slip 与 divergence。
[Constant et al. (2017)](https://arxiv.org/abs/1609.04364) 在 OpenFOAM/PISO
中采用 IBM 与 pressure correction sub-iteration，同时施加固体条件和
divergence-free constraint。

但这不等于“经典 MDF 必须每轮重解 Poisson”：

* predictor-level MDF 本来就是一个成立的方法类别；
* 交替 projection 的成本是多次全域 elliptic solve；
* 收敛速度取决于两个 constraint spaces 的夹角与预条件，并不一定高效；
* 它仍不是把 pressure 与 IB force 作为一套 Lagrange multipliers 同时求解
  的 monolithic method；
* 若主要误差是 forcing reciprocity 或有效几何偏差，反复 Poisson 并不
  对症。

因此，`ENABLE_MDF_PRESSURE_PROJECTION=1` 适合保留为诊断／严格对照后端，
不宜在没有证据时成为默认 baseline。

### 3.3 同时施加两个约束的成熟路线

[Taira & Colonius (2007)](https://doi.org/10.1016/j.jcp.2007.03.005)
把 pressure 与 immersed-boundary force 都作为 Lagrange multipliers，构造
immersed boundary projection method。其思想不是反复修补最终速度，而是
在 coupled saddle-point system 中同时满足 incompressibility 与 no-slip。

优点是约束结构清楚，边界 residual 可以做到与线性系统容差一致；缺点是：

* 需要 Schur complement、block factorisation 和专门 preconditioner；
* 与 OpenFOAM collocated finite-volume、变密度 VoF、PIMPLE 的整合远比
  predictor MDF 复杂；
* moving rigid body／two-way FSI 会进一步扩大 multiplier system。

这是长期“严格后端”的合理方向，而不是当前一两次修改就能安全完成的
补丁。

## 4. Diffuse／grid-scale 界面本身造成什么问题？

### 4.1 有效几何偏大

regularized-kernel 或有限厚度 forcing 会把固体作用扩散到真实几何两侧，
从而把物体表现为一个略大的 hydrodynamic body。Breugem 发现 marker
inward retraction 对有效直径与网格收敛阶影响很强，而 MDF 次数主要影响
no-slip/no-penetration residual；二者修的是不同问题
（[Breugem, 2012](https://doi.org/10.1016/j.jcp.2012.02.026)）。

这意味着：即使 slip 很小，drag、added mass 或涡脱落仍可能因为有效直径
偏差而错误。反过来，某个局部 slip 不为零，也不必然导致整体 load 最差。

对当前 cell-centred SDF 方法，`transitionWidth=0` 虽然令 solid fraction
二值化，但边界仍被量化到 cell centre；没有 cut-face area、真实 surface
intersection 或 ghost point。因此它仍具有 grid-scale geometry error，不能
因二值 mask 就称为严格 sharp interface。

### 4.2 流体侧动量和速度梯度被 forcing 污染

传统两侧 diffuse kernel 会把非零 force 铺展到 fluid nodes，使这些位置
不再求解原本的无 forcing Navier–Stokes equations，并污染近壁速度梯度。
[Peng & Wang (2020)](https://doi.org/10.1103/PhysRevE.101.053305) 提出：

* 只向 solid side 分配 forcing，恢复 fluid side 的原方程；
* 使用 force amplification 减少达到 no-slip 所需的迭代次数；
* 在其 laminar／turbulent tests 中，同时改善速度、近壁梯度以及 force／
  torque。

这个方向与当前 `diffuseVelocity` 很相关：如果后续启用平滑
`transitionWidth>0`，应特别检查 correction 是否作用到 fluid side；即使
保持二值 mask，也应检查 interface fluid cells 是否通过重构被过度修改。

### 4.3 Local explicit correction

[Vagnoli et al. (2025)](https://doi.org/10.1016/j.cpc.2025.109741)
从 implicit forcing system 出发，通过近似分解导出每个 marker、每个时间步
动态变化的局部显式 correction。其关键价值是：

* 同时考虑 transfer non-reciprocity 与邻近 marker support overlap；
* 不需要求解完整的昂贵 implicit force linear system；
* 可以嵌入已有 explicit IBM；
* 其 fixed/moving-body tests 显示 no-slip consistency 和 convergence 明显
  改善。

这类 local inverse／diagonal correction 是当前项目很值得优先借鉴的
“中间路线”：比无脑增加 MDF 更有针对性，又比 monolithic projection
容易实现。

### 4.4 Implicit direct forcing

[Park et al. (2016)](https://doi.org/10.1016/j.jcp.2016.03.035) 以
block-LU recurrence 和 preconditioner 构造 implicit direct forcing，指出：

* explicit direct forcing 对 no-slip 的施加不够精确；
* MDF 在某些条件下仍留有不可忽略的 surface velocity error；
* 低 Reynolds number 时，显式处理 diffusion 可能使问题更明显；
* 预条件 implicit recurrence 能在少量迭代内得到更准确 forcing，并允许
  相对更大的时间步。

[Abdol Azis, Evrard & van Wachem (2019)](https://doi.org/10.1016/j.jcp.2018.10.048)
则在 collocated、unstructured finite-volume 背景下发展 regularised discrete
momentum forcing，目标是在每个时间步精确施加 no-slip，并保持局部／全局
质量守恒。这条路线从离散 momentum equation 本身构造 forcing，比简单的
cell velocity overwrite 更接近 OpenFOAM 的自然升级方向。

## 5. Sharp-interface 能否彻底消除问题？

Sharp-interface／ghost-cell／cut-cell 的核心优点是把边界条件施加在真实几何
位置，不把 forcing 平滑到流体侧，因此可显著减少 diffuse geometry error
和近壁梯度污染。Mittal 等人的 ghost-cell 方法是代表性路线
（[Mittal et al., 2008](https://doi.org/10.1016/j.jcp.2008.01.028)）。

但 sharp interface 不是“免费消除所有 slip”：

* 必须稳定地寻找 surface intersection、image point 和 interpolation
  stencil；
* moving body 会产生 fresh/dead cells；
* cut cells 会遇到 small-cell stiffness；
* velocity stencil、pressure Poisson、face flux 和 turbulence quantities
  必须使用相容几何；
* 几何守恒不正确时会产生 spurious pressure／force oscillations。

[Seo & Mittal (2011)](https://doi.org/10.1016/j.jcp.2011.06.003)
把 moving sharp-interface IBM 中的主要压力振荡归因于局部质量／几何守恒
破坏，并只在 pressure Poisson 与 velocity correction 附近引入 cut-cell
treatment，加上 virtual cell merging 处理 small cells；其测试中压力振荡
约降低一个数量级。

因此，当前 `sharpReconstruction` 若要成为正式 sharp backend，至少还需要：

* 真实 surface–grid intersection 和沿法向的 image/ghost construction；
* 对 face flux 与 pressure equation 的一致处理；
* moving geometry 的 GCL／covered–uncovered cell 处理；
* VoF、$k$、$\omega$ 等 scalar fields 的配套 extension／boundary stencil；
* 与 load integration 一致的 surface quadrature。

## 6. Volume conservation 与 slip 的关系

### 6.1 divergence-free Eulerian velocity 不代表表面不渗漏

Eulerian field 在离散意义上满足 $\nabla_h\cdot\boldsymbol{u}=0$，并不保证把
它插值到 Lagrangian／immersed surface 后仍连续无散。反过来，在 surface
marker 上满足速度条件，也不保证 finite-volume face flux 形成正确的局部
质量平衡。

[Bao et al. (2017)](https://doi.org/10.1016/j.jcp.2017.06.041)
构造 continuously divergence-free velocity interpolation，并使用与其伴随的
force spreading；在其问题中，volume conservation 比传统 IB interpolation
改善多个数量级，但代价是非局部算子和额外 Poisson-type operations。

该方法来自 Lagrangian elastic IB，并不能原样移植到当前 SDF mask，但它
给出一个很重要的设计原则：

> velocity interpolation、force application 和 pressure/flux projection
> 应在离散层面彼此相容；只修 cell-centred velocity 不足以证明没有 leakage。

### 6.2 运动边界的额外问题

运动物体扫过固定网格时，cells 会在 fluid 与 solid 身份间切换。若新暴露
fluid cell 的 velocity、pressure 和 VoF 没有物理一致的 extension，或者
solid volume 变化没有进入 pressure/flux balance，会出现：

* 瞬时局部质量源／汇；
* pressure spikes 和 force oscillation；
* slip 指标随 cell classification 跳变；
* 两相界面的非物理撕裂或相变量污染。

这也是为什么 moving IBM 的验证不能只重复 static-cylinder test。

## 7. 两相 VoF + IBM 的特殊困难

### 7.1 `alpha.water` extension 修的不是 slip 本身

[Ye, Chen & Maki (2020)](https://doi.org/10.3390/jmse8100809)
在 OpenFOAM finite-volume、RANS、VoF 和 moving discrete-forcing IBM 中，
对 solid–air–water 交会附近使用 volume-fraction field extension／dilation，
并验证 dam-break obstacle、圆柱出水和船舵等问题。

field extension 的主要目标是让物体内部或刚被几何扫过的 cell 获得与外侧
流体相容的相信息，从而正确计算：

* density 和 viscosity；
* hydrostatic／dynamic pressure coupling；
* surface-tension／VoF advection 附近的物性；
* moving solid 与 free surface 的交会。

它**不是** interpolation/spreading correction，也不直接消除
pressure-projection splitting。因此，打开 alpha extension 后 slip 没有改善、
甚至局部变差，并不能说明 Ye et al. 的思想错误；更可能是当前 forcing、
pressure coupling 与 field extension 采用的离散框架不同，或 extension
范围／时序／stencil 与 MULES 不相容。

### 7.2 两相问题应额外记录什么？

* `min/max(alpha.water)` 与 boundedness；
* 水体总体积及 accessible-fluid volume 中的水体积；
* immersed surface 上按 liquid/gas/interface 分类的 $e_n$、$e_t$；
* free-surface–solid contact 区域的 local mass flux；
* IBM forcing 前后、pressure projection 前后的 $p_{rgh}$ 和 density jump；
* surface tension 与 IBM acceleration 在界面 cell 中的相对量级；
* 开／关 alpha extension 的差值场，而不只看全局最大值。

### 7.3 海工 monopile 中 free-slip 是否足够？

如果 PI 的建模假设是桩面摩擦相对总体波浪载荷较小，可以先将
free-slip/no-penetration 作为低成本 IBM baseline。此时：

* 法向 velocity/flux residual 是硬性指标；
* 切向相对速度不是错误；
* 不能用 free-slip 结果验证 bottom shear、wall shear 或 near-wall turbulence；
* pressure-dominated inline load 可能仍有意义，但必须与 body-fitted free-slip
  baseline 比较；
* 将来若研究 scour、boundary layer、Reynolds stress 或 viscous damping，
  no-slip 与合适 wall treatment 仍然不可替代。

## 8. 主流解决方案对照

| 路线 | 主要修复对象 | 优点 | 主要代价／风险 | 对当前项目的适合度 |
|---|---|---|---|---|
| 一次 explicit direct forcing | 基础速度约束 | 最简单、便宜 | residual 通常不为零；projection 后可回升 | 仅适合作对照 |
| Predictor-level MDF | forcing residual、support overlap | 容易实现；不必重复 Poisson | 不能保证 post-projection slip；次数存在平台 | **当前 baseline** |
| Pressure-projected MDF | slip 与 divergence 的交替约束 | 可降低 projection 后 residual | 多次全域 Poisson；仍非 monolithic | 诊断／严格对照 |
| Global force amplification | reciprocity error | 几乎无额外成本 | 对局部几何和动态变化适应有限 | 值得先试 |
| Local explicit correction | 局部 reciprocity 与 support overlap | 成本低于 implicit；可动态调整 | 需要可靠 local operator/weight | **优先研发候选** |
| Single-sided forcing | fluid-side 方程和梯度污染 | 改善近壁梯度、force/torque | kernel／mask 设计更敏感 | 对 diffuse backend 很相关 |
| Marker retraction | 有效直径偏大、网格收敛 | 方法成熟、成本小 | 最优距离依方法／几何而变 | 若采用 marker backend 很有用 |
| Preconditioned implicit forcing | forcing residual、显式 diffusion error | 少量迭代内更严格 | 修改线性系统、实现复杂 | 中期路线 |
| RDMF／implicit FVM forcing | no-slip 与局部/全局质量守恒 | 适合 collocated unstructured FVM | 必须深入 momentum/flux 离散 | **OpenFOAM 中期重点** |
| Coupled constraint projection | splitting + forcing constraint | 理论最完整、可达线性容差 | saddle-point／Schur solver 成本高 | 长期严格后端 |
| Ghost-cell/sharp reconstruction | diffuse geometry 和梯度污染 | 在真实边界位置施加条件 | stencil、moving GCL、VoF scalar 复杂 | 中长期独立后端 |
| Cut-cell pressure/flux | 局部质量与 moving pressure oscillation | 守恒性强 | small-cell 与几何计算复杂 | moving FSI 后期重要 |

## 9. 当前 `ichiBoundaryFoam` 与文献的对应关系

### 9.1 已实现部分

当前正式源码可概括为：

* OpenFOAM v2312、collocated finite volume；
* `interFoam`-like algebraic VoF／PIMPLE；
* analytic capped-cylinder SDF；
* cell-centred `immersedSolidFraction` 与 `immersedCellType`；
* discrete direct correction of momentum predictor；
* `forcingMethod diffuseVelocity` 或实验性 `sharpReconstruction`；
* runtime `velocityBoundaryCondition noSlip/freeSlip`；
* optional predictor-level MDF；
* ordinary full-domain continuous pressure；
* optional alpha extension；
* optional slip fields／surface-band statistics；
* prescribed rigid-body motion，尚非 load-driven two-way FSI。

### 9.2 默认 MDF 路径

当前默认编译配置为：

```makefile
ENABLE_MDF_PRESSURE_PROJECTION ?= 0
```

运行顺序是：

```text
momentum predictor
    -> initial IBM correction
    -> predictor MDF corrections
    -> one ordinary PIMPLE pressure projection
    -> residual report
```

这与经典 predictor-level MDF 的思想相容，避免了旧版每一轮 forcing 都
重解全域 Poisson 的巨大成本。但它留下一个明确待验证的问题：普通
projection 会把 residual 恢复多少。

### 9.3 当前 `maxSlip()` 的含义与局限

`maxSlip()` 当前计算的是固体域 cell centres 中

$$
\max_c
\chi_s(c)
\left\lVert
\mathcal{C}
[\boldsymbol{u}(c)-\boldsymbol{U}_b(c)]
\right\rVert,
$$

其中 $\mathcal{C}$ 在 no-slip 下保留全部相对速度，在 free-slip 下只保留
法向分量。

它适合快速监控 solver constraint，但不是严格 surface slip，因为：

* sample points 在 cell centres，不在 $\mathrm{SDF}=0$；
* no-slip 时会把深处 solid cells 也纳入最大值；
* 它不显式检查 cut-face／surface flux；
* 它的数值会随 `transitionWidth`、solid mask 和网格相位改变；
* 与 Lagrangian marker IBM 文献报告的 surface residual 不能直接横比。

### 9.4 当前详细 sampling 的含义与局限

启用 `ENABLE_SLIP_SAMPLING=1` 后，程序在
$|\mathrm{SDF}|\leq w_s$ 的 surface band 内输出：

* `immersedSlipVelocity`；
* `immersedSlipMagnitude`；
* `immersedNormalSlip`；
* `immersedTangentialSlip`；
* `immersedSlipSamplingMask`；
* liquid/gas sample counts、max values、位置、alpha 与 signed distance。

这比单个 `maxSlip` 强得多，但仍是 band-cell metric。下一步应在解析圆柱
表面布置独立 quadrature/sample points，沿法向从 fluid side 重构速度，
得到真正 $\mathrm{SDF}=0$ 上的 $e_n$、$e_t$ 和 surface flux。

### 9.5 `diffuseVelocity` 与 Lagrangian diffuse forcing 不完全等价

当前 diffuse correction 是

$$
\Delta\boldsymbol{u}(c)
=
\chi_s(c)
\mathcal{C}
[\boldsymbol{U}_b(c)-\boldsymbol{u}^{*}(c)].
$$

它没有 Lagrangian markers，也没有 $E/H$ regularized-delta operators。
因此，Lagrangian literature 中 $EH$ 的 reciprocity algebra 不能逐项原封不动
套用。对当前实现更准确的说法是：

* `transitionWidth=0` 的 binary mask 属于 Eulerian selection forcing；
* 2026 年统一分析指出理想 Eulerian selection operators 可使传统
  reciprocity term 消失；
* 但真实程序仍有 SDF cell-centre geometry error、neighbor reconstruction、
  face-flux incompatibility 和 pressure splitting；
* `transitionWidth>0` 或 `sharpReconstruction` 引入的 weighted／neighbor
  operations 仍可能产生类似 operator mismatch 的空间残差。

因此，不应把当前所有 slip 都叫作 interpolation/spreading non-reciprocity；
必须用分阶段记录来识别。

## 10. 建议的验证矩阵

### 10.1 第一阶段：在同一个时间步拆开 residual

对静止单相 cylinder，固定网格与 $\Delta t$，输出：

| 时间位置 | no-slip 指标 | free-slip 指标 | 目的 |
|---|---|---|---|
| forcing 前 | $e_n,e_t$ | $e_n$ | 初始误差 |
| MDF-1 后 | $e_n,e_t$ | $e_n$ | 一次 correction 效果 |
| MDF-k 后 | $e_n,e_t$ | $e_n$ | forcing convergence curve |
| PIMPLE projection 后 | $e_n,e_t$ | $e_n$ | splitting error |
| 最终写盘前 | $e_n,e_t,Q_\Gamma$ | $e_n,Q_\Gamma$ | 真正验收量 |

MDF correctors 建议至少扫：

```text
1, 2, 4, 8
```

压力投影后端只做 A/B：

```text
ENABLE_MDF_PRESSURE_PROJECTION=0
ENABLE_MDF_PRESSURE_PROJECTION=1
```

如果 predictor curve 已形成平台，再增加到 16 次通常没有解释价值。

### 10.2 第二阶段：空间与时间收敛

至少进行：

* 三个 cylinder diameter resolution；
* 三个 fixed $\Delta t$；
* 保持相同物理流况、相同 forcing parameters；
* 比较 $L_{\infty}$、RMS、P95、$Q_\Gamma$、drag/lift 和运行时间。

判断规则：

* 误差随 $\Delta t$ 明显变化而随 $\Delta x$ 不敏感：优先怀疑 splitting／
  temporal forcing；
* 误差随 $\Delta x$ 或 cylinder-grid phase 强烈变化：优先怀疑 geometry／
  interface reconstruction；
* MDF-k 期间不下降：优先怀疑 forcing operator／relaxation；
* MDF-k 后很小、projection 后回升：优先修 pressure/flux compatibility；
* static 正常、moving 失效：优先修 GCL、fresh/dead cells 和 field extension；
* monophase 正常、VoF 失效：优先修 alpha/rho/CSF 与 IBM 时序。

### 10.3 第三阶段：物理量敏感性

依次比较：

1. body-fitted free-slip vs IBM free-slip；
2. body-fitted no-slip vs IBM no-slip；
3. laminar monophase static cylinder；
4. laminar VoF + Stokes II + static cylinder；
5. RANS VoF + Stokes II + static cylinder；
6. prescribed translating／rotating cylinder；
7. 最后才是 load-driven motion。

比较量不只包括 slip，还应有：

* $F_x,F_y,F_z$ 及 phase-averaged loads；
* surface elevation 与 probe velocities；
* phase volume、alpha boundedness；
* pressure solver iterations 与 wall time；
* wake/vorticity/$\lambda_2$ structure；
* moving-body pressure/force oscillation。

### 10.4 内部验收阈值如何设？

文献不存在适用于所有 IBM、网格和物理问题的统一“slip 必须小于某个
m/s”标准。建议暂时采用相对指标，并明确标为项目内部工程阈值，而不是
理论定律：

* $e_{n,\mathrm{RMS}}/U_{\mathrm{ref}}$ 目标不高于 $10^{-2}$；
* $e_{n,\mathrm{P95}}/U_{\mathrm{ref}}$ 目标在几个百分点内；
* 最大值必须随网格／时间步收敛，且离群位置可解释；
* body-fitted 与 IBM 的主载荷、相位和自由液面差异应小于项目预定误差；
* free-slip 不以 tangential residual 作为失败判据；
* no-slip 必须另外验证 tangential residual 和 near-wall response。

最终是否“可用”应由目标 observable 的敏感性决定，而不能只由某个单点
最大 slip 决定。

## 11. 面向 `ichiBoundaryFoam` 的实施优先级

### P0：先把诊断做成可信的

* 记录 forcing 前、每轮 MDF 后、projection 后的 residual；
* 输出 max/RMS/P95，而非只输出 max；
* 明确 no-slip 与 free-slip 的不同验收量；
* 增加解析 cylinder surface quadrature，使用 fluid-side reconstruction；
* 增加 $Q_\Gamma$ 和局部 normal-flux 分布；
* 记录最坏 sample 的 alpha、SDF、cell volume、processor 和 surface angle；
* 为 slip sampling 加独立 write interval，避免每个时间步写大场。

### P1：低风险 forcing 改善

* 给 correction 增加 runtime relaxation／amplification coefficient；
* 扫描 MDF convergence curve 并允许 residual-based stopping；
* 实现或试验 local diagonal/inverse correction；
* 检查 correction 是否不必要地污染 fluid-side cells；
* 将 pressure-projected MDF 仅保留为 comparison backend；
* 用 body-fitted baseline 校准 grid-scale effective cylinder radius。

### P2：修正 pressure／flux compatibility

若证据显示 residual 主要在 projection 后回升：

* 检查 `phiHbyA` 与被约束 `HbyA` 的构造时序；
* 在 immersed surface 附近显式检查／修正 normal face flux；
* 研究 pressure predictor 和 IBM-compatible pressure correction；
* 评估只在 interface neighborhood 求解局部 correction，而不是每轮全域
  Poisson；
* 研究 RDMF 或 coupled multiplier Schur complement，避免长期依赖交替
  projection。

### P3：正式 sharp／moving backend

* surface intersection、ghost/image points；
* pressure/velocity/VoF/turbulence 一致 stencil；
* cut-face flux 与 GCL；
* covered/uncovered cell extension；
* small-cell stabilization／virtual merging；
* surface load quadrature；
* prescribed motion 后再进入 two-way 6DoF。

## 12. 对当前项目的推荐技术路线

### 近期 baseline

建议把以下组合定义为开发 baseline：

```text
discrete direct forcing
+ analytic SDF cylinder
+ predictor-level MDF
+ one regular PIMPLE pressure projection
+ continuous pressure
+ freeSlip/no-penetration for the first RANS wave benchmark
+ detailed pre/post-projection slip diagnostics
```

理由是：

* 能维持当前 `interFoam`-like VoF 主体；
* 计算成本比 legacy pressure-projected MDF 小得多；
* 对 pressure-dominated monopile wave load，free-slip 是 PI 已接受的初步物理
  假设；
* 它保留足够清楚的误差分解空间；
* 不会把实验性 sharp reconstruction 误当作已经成熟的 sharp IBM。

### 如果 predictor MDF 仍不收敛

优先顺序建议为：

1. relaxation／force amplification；
2. local explicit correction；
3. single-sided／fluid-side-preserving forcing；
4. implicit momentum forcing/RDMF；
5. 最后才是无上限增加 MDF 次数。

### 如果只有 projection 后 residual 回升

优先顺序建议为：

1. 修正 diagnostic，确认是真实 surface normal slip；
2. 检查 HbyA–phi–pressure correction 的离散一致性；
3. pressure predictor／local flux correction；
4. 以 legacy projected-MDF 做参照，量化“多次 Poisson 换来多少 residual”；
5. 中期转向 coupled constraint/RDMF，而不是永久使用昂贵交替投影。

### 如果只有 VoF case 失效

不要先改 MDF 次数。应检查：

1. alpha extension 的执行时序与范围；
2. solid cells 中的 alpha/rho/mu 是否进入 pressure、MULES 或 CSF；
3. free surface 与 SDF 的交会 cells；
4. liquid/gas/interface 三组 slip；
5. accessible-fluid phase conservation；
6. IBM correction 是否在 density jump 两侧使用一致权重。

## 13. 向 PI 汇报时可以怎样表述？

可以简洁地说：

> 当前的核心问题不是求解器会不会运行，而是 immersed boundary 的速度
> 约束在 pressure projection 后能保持到什么程度。文献表明 direct-forcing
> IBM 的 residual slip 至少由 forcing operator 与 pressure splitting 两类
> 机制产生；MDF 主要修第一类，不一定修第二类。我们已取消默认每轮 MDF
> 都重解 Poisson 的高成本做法，接下来将分别记录 forcing 前后和 projection
> 前后的法向／切向 residual，以决定是做局部 force correction，还是修正
> pressure/flux compatibility。对首个海工 RANS baseline，我们可先采用
> free-slip/no-penetration，只以法向渗漏作为硬性边界指标，并用 body-fitted
> cylinder 对照载荷和波场。

这比笼统说“slip 很大，所以 IBM 失败了”更准确，也比笼统说“多做几次
MDF 就会好”更有文献依据。

## 14. 主要原始文献

* Uhlmann, M. (2005). *An immersed boundary method with direct forcing for
  the simulation of particulate flows*. Journal of Computational Physics,
  209, 448–476.
  [DOI](https://doi.org/10.1016/j.jcp.2005.03.017)
* Taira, K. & Colonius, T. (2007). *The immersed boundary method: A projection
  approach*. Journal of Computational Physics, 225, 2118–2137.
  [DOI](https://doi.org/10.1016/j.jcp.2007.03.005)
* Wang, Z., Fan, J. & Luo, K. (2008). *Combined multi-direct forcing and
  immersed boundary method for simulating flows with moving particles*.
  International Journal of Multiphase Flow, 34, 283–302.
  [DOI](https://doi.org/10.1016/j.ijmultiphaseflow.2007.10.004)
* Domenichini, F. (2008). *On the consistency of the direct forcing method in
  the fractional step solution of the Navier–Stokes equations*. Journal of
  Computational Physics, 227, 6372–6384.
  [DOI](https://doi.org/10.1016/j.jcp.2008.03.009)
* Mittal, R. et al. (2008). *A versatile sharp interface immersed boundary
  method for incompressible flows with complex boundaries*. Journal of
  Computational Physics, 227, 4825–4852.
  [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2834215/)
* Seo, J. H. & Mittal, R. (2011). *A sharp-interface immersed boundary method
  with improved mass conservation and reduced spurious pressure oscillations*.
  Journal of Computational Physics, 230, 7347–7363.
  [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3156558/)
* Breugem, W.-P. (2012). *A second-order accurate immersed boundary method for
  fully resolved simulations of particle-laden flows*. Journal of
  Computational Physics, 231, 4469–4498.
  [DOI](https://doi.org/10.1016/j.jcp.2012.02.026)
* Park, H., Pan, X., Lee, C. & Choi, J.-I. (2016). *A pre-conditioned implicit
  direct forcing based immersed boundary method for incompressible viscous
  flows*. Journal of Computational Physics, 314, 774–799.
  [DOI](https://doi.org/10.1016/j.jcp.2016.03.035)
* Constant, E. et al. (2017). *Implementation of a discrete Immersed Boundary
  Method in OpenFOAM*.
  [Preprint](https://arxiv.org/abs/1609.04364)
* Bao, Y. et al. (2017). *An Immersed Boundary method with divergence-free
  velocity interpolation and force spreading*. Journal of Computational
  Physics.
  [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6781881/)
* Abdol Azis, M. H., Evrard, F. & van Wachem, B. (2019). *An immersed boundary
  method for incompressible flows in complex domains*. Journal of
  Computational Physics, 378, 770–795.
  [DOI](https://doi.org/10.1016/j.jcp.2018.10.048)
* Peng, C. & Wang, L.-P. (2020). *Force-amplified, single-sided
  diffused-interface immersed boundary kernel for correct local velocity
  gradient computation and accurate no-slip boundary enforcement*. Physical
  Review E, 101, 053305.
  [DOI](https://doi.org/10.1103/PhysRevE.101.053305)
* Ye, H., Chen, Y. & Maki, K. (2020). *A Discrete-Forcing Immersed Boundary
  Method for Moving Bodies in Air–Water Two-Phase Flows*. Journal of Marine
  Science and Engineering, 8, 809.
  [Full text](https://doi.org/10.3390/jmse8100809)
* Zhou, K. & Balachandar, S. (2021). *An analysis of the spatio-temporal
  resolution of the immersed boundary method with direct forcing*. Journal of
  Computational Physics, 424, 109862.
  [DOI](https://doi.org/10.1016/j.jcp.2020.109862)
* Gsell, S. & Favier, J. (2021). *Direct-forcing immersed-boundary method: A
  simple correction preventing boundary slip error*. Journal of Computational
  Physics, 435, 110265.
  [DOI](https://doi.org/10.1016/j.jcp.2021.110265)
* Vagnoli, L. et al. (2025). *A local and explicit forcing correction for
  Lagrangian immersed boundary methods*. Computer Physics Communications,
  315, 109741.
  [DOI](https://doi.org/10.1016/j.cpc.2025.109741)
* *On Reciprocity and Splitting Errors in Direct-Forcing Immersed Boundary
  Methods* (2026).
  [DOI](https://doi.org/10.1007/s10494-026-00802-7)

## 15. 最终判断

从主流与近期研究来看，当前项目没有偏离 discrete-forcing IBM 的基本问题
谱系：residual slip、pressure splitting、diffuse geometry、local mass
conservation 和 moving-interface treatment 都是该领域长期存在的核心难点。

但当前方法也还不能因“solver 不崩溃”或“MDF residual 下降”就称为成熟：

* 它尚未在真实 immersed surface 上测量 residual；
* predictor forcing 与 final projection 的误差尚未分离；
* cell velocity 与 face flux 的相容性尚未证明；
* `sharpReconstruction` 尚不是守恒的 sharp IBM；
* VoF、turbulence 和 moving geometry 的边界处理仍需独立验证。

最务实的下一步不是立即放弃现有 backend，而是先用可信诊断确认 dominant
error mechanism。若 dominant error 是 predictor forcing，本项目可沿 local
correction／single-sided／implicit momentum forcing 改进；若 dominant error
是 projection，则应转向 pressure-compatible flux correction 或 coupled
constraint；若只有 moving/two-phase 才出现问题，则应集中处理 GCL 与相场
extension，而不是继续增加 MDF 次数。
