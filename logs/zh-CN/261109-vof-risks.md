# VoF–IBM 在高密度比下的风险与对策

> 日期：2026-11-10  
> 对象：`ichiBoundaryFoam` 的 air–water VoF 与 discrete-forcing IBM  
> 重点：自由液面穿越 immersed boundary 时，密度、黏度、压力、动量与
> IBM forcing 的一致性。

## 结论先行

PI 提出的风险完全成立。VoF–IBM 不是把一个已经成熟的单相 IBM 原样放进
两相 solver 即可。对 air–water 问题，

$$
\rho(\alpha)=\alpha\rho_w+(1-\alpha)\rho_a,
$$

$$
\mu(\alpha)=\alpha\mu_w+(1-\alpha)\mu_a,
$$

且

$$
\frac{\rho_w}{\rho_a}\sim 10^3.
$$

因此，immersed surface 附近很小的 $\alpha$ 误差就可能变成很大的密度误差，
继而污染浮力项、压力修正、动量通量、IBM 载荷和自由液面形状。这里尤其
危险的不是纯水区或纯空气区，而是 **air–water–solid 三者相交的窄带**，
以及运动固体每一步扫过的新旧 cells。

当前 `ichiBoundaryFoam` **已经有一部分正确准备**：

* 继承了 `interFoam` 的 one-fluid VoF 框架，动量方程采用
  `ddt(rho,U)` 与 `div(rhoPhi,U)` 的守恒形式；
* 压力方程使用随 $\rho(\alpha)$ 变化的系数，而非把全域密度视为常量；
* 每次求解 $\alpha$ 后调用 `mixture.correct()`，更新混合物性；
* IBM constraint impulse 的载荷统计使用局部 `rho*correction`；
* 已实现可选的 `alpha.water` 向固体内部逐层延拓，默认五层、反距离平方
  加权，并可把延拓结果限制在 $[0,1]$；
* 已有 `alpha.water` 上下界、phase volume、Courant number、slip 和
  divergence 等诊断。

但目前仍只能把它称为 **具备两相 IBM 开发基础的研究原型**，不能宣称已经
解决高密度比耦合。最重要的缺口是：

* `alpha` 延拓默认关闭，且尚未通过自由液面穿越静止／运动圆柱的系统验证；
* 还没有独立证明 `rhoPhi`、VoF 相通量、IBM velocity correction 与压力
  投影在离散层面保持严格的质量—动量一致性；
* `alpha` 延拓会改变固体内部的“虚构流体体积”，现有总体积诊断若把固体
  内部也计入，就可能误报；必须另行统计 **可达流体域** 内的水体积；
* moving IBM 尚未实现严格的 swept-volume／geometric conservation，
  即固体移动时如何排开水、释放空气并处理 fresh/dead cells；
* 尚无 air–water–solid 接触角／润湿模型；目前的延拓最多对应一种近似中性
  Neumann 条件，不能代表任意物理接触角；
* IBM forcing 虽然按运动学方式修正速度，载荷按局部密度换算，但尚未证明
  forcing、黏性项、压力项和相通量在跨相界面时形成同一个离散平衡；
* immersed turbulence-wall treatment 尚未在三相交会区验证。

所以目前最准确的汇报说法是：

> 我们已经识别并实现了 Ye 等人所强调的相分数延拓机制，也保留了
> variable-density conservative VoF 的基础结构；下一步必须用专门的
> interface-crossing benchmarks 验证离散质量—动量一致性、物性梯度与
> moving-body 排液守恒。短程 V1/R1 没有发散是必要条件，但不是充分验证。

## 1. 为什么 $10^3$ 密度比会放大 IBM 误差？

在 one-fluid VoF 中，空气和水共用一套速度与压力方程，物性由 $\alpha$
插值得到。若某 immersed-boundary stencil 错把本应为水的 cell 赋成空气，
或反过来，则

$$
\delta\rho=(\rho_w-\rho_a)\delta\alpha.
$$

当 $\rho_w-\rho_a$ 很大时，即使 $\delta\alpha$ 很小，$\delta\rho$ 也不可忽略。
随后误差至少通过四条路径传播：

* **重力／浮力路径**：`p_rgh` 形式中包含与 $\nabla\rho$ 有关的项；固体内部
  错误的 $\alpha$ 会在 immersed surface 上制造伪密度梯度；
* **惯性路径**：同一个速度误差在水侧对应的动量误差约为空气侧的
  $10^3$ 倍；
* **压力投影路径**：variable-coefficient Poisson operator 对 $1/\rho$
  敏感，界面两侧的压力修正强度不同；
* **载荷路径**：相同的 velocity correction 换算为 impulse/force 时会乘以
  局部密度；错误的相标签会直接造成数量级错误的 IBM force。

黏度也会跨界面变化，但需要区分动力黏度 $\mu$ 与运动黏度 $\nu=\mu/\rho$。
对于 air–water，密度跳变通常比动力黏度跳变更剧烈；不过近壁剪切、湍流
模型和扩散离散仍会受到 $\mu(\alpha)$ 与 $\nu(\alpha)$ 选择的影响。

## 2. Ye 等人给出的直接答案：延拓 $\alpha$

[Ye、Chen 与 Maki (2020)](https://doi.org/10.3390/jmse8100809)
专门讨论了 moving bodies in air–water VoF。其判断与 PI 的担心几乎完全
一致：固体内部的流体是非物理的，因而其中的 $\alpha$ 和 $\rho$ 没有自然
定义；如果直接让这个未定义／错误的场参与 immersed-boundary stencil，
$\nabla\rho$ 会在水气密度比作用下产生很大误差。

他们的处理是：

1. 找出固体 cells；
2. 从真实流体侧开始，使用共享节点的已知邻居对 $\alpha$ 作反距离平方插值；
3. 一层一层向固体内部延拓；
4. 使用约五层延拓，为 immersed surface 提供近似
   $\partial\alpha/\partial n=0$ 的边界条件；
5. 再由延拓后的 $\alpha$ 计算 $\rho$ 和 $\nu$。

这并不是“只把水向固体内扩散”。若壁外是空气，就延拓空气的 $\alpha=0$；
若壁外是水，就延拓水的 $\alpha=1$；若自由液面正穿过固体，则延拓局部的
混合相分数。因此它的目的，是让物性在 **壁面法向上保持合理延续**，避免
固体内部任意值污染壁外差分。

Ye 等人也明确承认：该操作会改变固体内部计算出来的 fluid volume；但这部分
只是 ghost/fictitious state，不应被解释为真实水量。因为理想 IBM 表面没有
法向通量，真正需要守恒的是 fluid-accessible domain 中的水量。

我们的 `DirectForcingIbm::extendAlpha()` 与这一思路高度接近：

* 以 SDF 区分流体侧与固体侧；
* 使用 point-neighbour stencil；
* 逐层延拓；
* 默认 `nAlphaExtensionLayers = 5`；
* 默认权重幂为 2；
* 可将延拓值限制在 $[0,1]$。

这是目前最直接的准备，但它仍是运行时可选项且默认关闭。此前一次开启后结果
看起来更差，不能据此否定文献方法；更可能的解释是：延拓只修复了物性边界
条件，却暴露了 forcing、flux、pressure projection 或 volume accounting
中的另一处不一致。延拓必须放在完整一致算法中验证，不能作为孤立补丁评判。

## 3. 近期研究强调：相分数正确还不够，质量与动量必须共用通量

[van der Eijk 与 Wellens (2023)](https://doi.org/10.1016/j.jcp.2022.111796)
针对高密度比自由液面与 moving-body FSI 指出，若相质量和动量使用不一致的
离散通量，会出现：

* momentum loss；
* spurious interface velocity；
* 自由液面畸变；
* 甚至数值失稳。

其核心不是简单地“把网格加密”或“减小时间步”，而是让动量对流项使用与
VoF 相输运一致的质量通量。文中两相 moving-body tests 显示，一致方法的动量
守恒误差可比非一致方法改善约一个数量级。

[Ghods 与 Herrmann (2012)](https://doi.org/10.1016/j.compfluid.2012.04.002)
及后续大量高密度比工作也给出同一结论：应使用守恒形式

$$
\frac{\partial(\rho\boldsymbol{u})}{\partial t}
+\nabla\cdot(\rho\boldsymbol{u}\boldsymbol{u}),
$$

并让动量通量中的质量通量来自界面输运，而不是独立重新插值得到一个不一致
的 face density。

当前 `ichiBoundaryFoam` 的有利之处是沿用了

```cpp
fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
```

而不是常密度或非守恒的 `ddt(U)+div(phi,U)`。这说明基础方向正确。
但 IBM correction 是在 predictor 上额外施加的运动学速度修正，之后还要经历
flux reconstruction 和 pressure projection。因此仍需用预算量证明：

$$
\Delta M_w,
\quad
\Delta\int\rho\boldsymbol{u}\,\mathrm{d}V,
\quad
\int\boldsymbol{f}_{IB}\,\mathrm{d}V,
\quad
\int p\boldsymbol{n}\,\mathrm{d}A
$$

在离散层面彼此相容，而不能仅凭源码中出现 `rhoPhi` 就宣布通过。

## 4. 其他成熟路线怎样处理？

### 4.1 守恒 diffuse-interface 路线

[Patel 与 Natarajan (2018)](https://doi.org/10.1016/j.jcp.2018.01.024)
采用 fluid volume fraction 加 solid volume fraction 的全 Eulerian
diffuse-interface 方法：

* 在全域（包括固体内部）求一个混合动量方程；
* 使用 variable-coefficient pressure Poisson equation；
* 令质量与动量输运一致；
* 通过 well-balanced 离散处理高密度比与强体力；
* 将 solid 内虚拟流体物性明确规定，而不是让其任意漂移。

这条路线与我们的 `diffuseVelocity + immersedSolidFraction` 在表面形式上接近，
但他们的重点是 **从方程层面构造离散守恒的混合系统**；我们的当前实现更像
在 `interFoam` predictor 上追加速度约束，尚未完成同等级的守恒证明。

### 4.2 Sharp embedded/cut-cell 路线

[Tavares 等 (2024)](https://doi.org/10.1016/j.compfluid.2024.106317)
把 geometric VoF 与 sharp embedded boundary/cut-cell 结合，并专门处理
复杂固体表面的接触角与 moving contact line。其优点是固体边界位置、流体
体积和接触几何更明确；代价是 cut-cell 几何、small-cell 稳定性、运动几何
更新和三维实现复杂得多。

这类工作说明：若最终问题非常依赖 wetted surface、瞬态水线、接触角与冲击
载荷，单靠 diffuse forcing 和简单 $\alpha$ 延拓可能不够；但对当前“先建立
可运行的 prescribed-motion monopile IBM baseline”，直接切换 cut-cell
并非最低风险路线。

### 4.3 Consistent level-set／DLM 路线

[Nangia 等人的高密度比 WSI 实现](https://arxiv.org/abs/1901.07892)
把 conservative variable-coefficient Navier–Stokes、consistent mass–momentum
transport 与 distributed Lagrange multiplier IBM 结合。这再次说明，主流
高密度比方法通常把以下三件事视为一个整体：

* fluid-interface transport；
* rigid-body constraint；
* variable-density pressure projection。

它们不是三个可以完全独立验证的插件。

## 5. 对当前实现的逐项审计

| 问题 | 当前准备 | 判断 |
|---|---|---|
| $\rho,\mu$ 随 $\alpha$ 更新 | `mixture.correct()` | 已具备 |
| 守恒变量密度动量式 | `ddt(rho,U)+div(rhoPhi,U)` | 方向正确，需预算验证 |
| variable-density pressure equation | 沿用 `interFoam` 式 `rAUf/rho` 路径 | 已具备基础 |
| 固体内 $\alpha$ 的定义 | 可选五层 SDF 延拓 | 已实现，未验证，默认关闭 |
| 延拓后的 boundedness | `boundExtendedAlpha` | 已具备 |
| IBM load 的局部密度 | `rho*correction` | 已具备基础，尚未与 BF load 对照 |
| 水／气／界面分区 slip | 有 liquid threshold 等诊断框架 | 部分具备 |
| 流体可达域水量守恒 | 暂无专门积分 | 缺失 |
| IBM force 前后动量预算 | 暂无闭合预算 | 缺失 |
| moving-solid swept volume | 暂无严格 GCL/排液方案 | 关键缺失 |
| 接触角／润湿 | 无 immersed contact-angle model | 缺失 |
| 表面张力与 IBM 联合平衡 | 沿用 VoF CSF，未做交会区验证 | 未验证 |
| immersed RANS 三相交会 | turbwall 初稿存在 | 未验证 |

## 6. 现在最值得做的验证，而不是继续堆功能

### 6.1 静止 cylinder 穿过静止自由液面

设置一个部分浸没、静止的 IBM cylinder，不造波，分别运行：

* body-fitted reference；
* IBM，`alpha extension off`；
* IBM，`alpha extension on`；
* extension layer 数 $N=1,3,5,7$。

验收量：

* hydrostatic pressure；
* 静水自由液面是否产生 parasitic motion；
* 可达流体域的水体积；
* cylinder 的浮力与理论排水体积；
* liquid/gas/interface 三类 IBM slip；
* $\max|\nabla\rho|$ 的位置与幅值；
* 网格与时间步收敛。

这是最直接回答 PI 问题的 baseline。

### 6.2 cylinder 匀速穿越自由液面

在无波、低速条件下令 cylinder 做 prescribed vertical translation，检查：

* 固体扫入水中时是否把水“凭空删除”或“凭空制造”；
* 出水时是否留下不合理水膜、空穴或 $\alpha$ 尖峰；
* 总水量、总动量和 IBM impulse 的预算闭合；
* force history 是否在 cell crossing 时出现强烈网格频率振荡；
* extension on/off 对结果的影响。

这一步通过后，再上简谐倾倒与 Stokes II waves。

### 6.3 高密度比敏感性

不要只跑真实 $1000:1$。建议依次测试

$$
\rho_w/\rho_a=1,10,100,1000.
$$

若误差随密度比近似线性或突然放大，就能把问题定位到物性／通量耦合，而
不是把所有异常都归因于 IBM slip 或波浪非线性。

## 7. 推荐的近期开发顺序

1. **保留当前 conservative VoF baseline**，不要先重写成另一套界面方法。
2. **把 alpha extension 纳入正式 A/B 验证**，但在通过前仍默认关闭。
3. **新增 fluid-accessible volume 统计**，积分时排除 solid fraction；同时
   分别报告 raw domain volume 与 physical fluid volume。
4. **新增质量—动量—IBM impulse 预算**，至少做到每个 write interval 可核对。
5. **先完成静水部分浸没和匀速穿越 benchmark**，再研究 M1/V1/R1 长程波浪。
6. 若 moving-body 排液误差明显，再实现 swept-volume/GCL correction；这比
   单纯增加 MDF 或 pressure projections 更对症。
7. 只有当研究问题确实需要动态接触线与湿润角时，再决定实现 immersed
   contact-angle model，或转向 sharp embedded/cut-cell 路线。

## 8. 对目前 V1/R1 结果的谨慎解释

V1 和 R1 在真实 air–water 物性、静止 IBM cylinder 下完成前 $20\,\mathrm{s}$，
自由液面与 phase volume 没有发散，这是积极结果，说明：

* variable-density VoF 与当前 static IBM 至少能够共同运行；
* 取消 PPMDF 没有立即造成两相 PDE 失稳；
* $\alpha$ 小幅超界随后衰减，而非持续增长。

但它们尚不能证明：

* 自由液面穿越 immersed surface 时的 $\nabla\rho$ 正确；
* moving cylinder 能守恒排液；
* IBM load 在水线附近可靠；
* extension on/off 已经收敛到同一物理解；
* 1000:1 密度比下的离散质量—动量预算闭合。

因此，当前不是“VoF–IBM 已经失败”，也不是“已经完成”；更准确地说，我们
已经进入 **可以做针对性 verification 的阶段**，而 PI 提出的高密度比问题，
正应成为下一轮 baseline tests 的核心。

## 参考文献

* H. Ye, Y. Chen, K. Maki, “A Discrete-Forcing Immersed Boundary Method
  for Moving Bodies in Air–Water Two-Phase Flows,” *Journal of Marine
  Science and Engineering*, 2020.
  <https://doi.org/10.3390/jmse8100809>
* M. van der Eijk, P. Wellens, “Two-phase free-surface flow interaction
  with moving bodies using a consistent, momentum preserving method,”
  *Journal of Computational Physics*, 2023.
  <https://doi.org/10.1016/j.jcp.2022.111796>
* M. Ghods, M. Herrmann, “Consistent mass and momentum transport for
  simulating incompressible interfacial flows with large density ratios
  using the level set method,” *Computers & Fluids*, 2012.
  <https://doi.org/10.1016/j.compfluid.2012.04.002>
* J. K. Patel, G. Natarajan, “Diffuse interface immersed boundary method
  for multi-fluid flows with arbitrarily moving rigid bodies,” *Journal
  of Computational Physics*, 2018.
  <https://doi.org/10.1016/j.jcp.2018.01.024>
* M. Tavares et al., “A coupled VOF/embedded boundary method to model
  two-phase flows on arbitrary solid surfaces,” *Computers & Fluids*,
  2024.
  <https://doi.org/10.1016/j.compfluid.2024.106317>
* N. Nangia et al., “A DLM immersed boundary method based wave-structure
  interaction solver for high density ratio multiphase flows,” 2019.
  <https://arxiv.org/abs/1901.07892>
