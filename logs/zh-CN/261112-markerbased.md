# Marker-based IBM 顶层配置与第一阶段路线

## 本次改动

IBM 配置由原先单层的 `forcingMethod`，调整为两层选择：

```text
ibm
{
    enabled          true;
    couplingMethod   cellBased;
    forcingMethod    diffuseVelocity;
}
```

* `couplingMethod` 先决定固体几何与 Eulerian 流场如何耦合；
* `forcingMethod` 再决定该耦合路线内部采用哪一种强迫实现；
* 若旧算例没有显式给出 `couplingMethod`，程序默认使用
  `cellBased`，因而既有的 `diffuseVelocity` 和
  `sharpReconstruction` 算例不会改变行为；
* 不兼容的组合会在初始化阶段直接报错，而不会静默退回另一种
  算法。

当前允许的组合为：

| `couplingMethod` | `forcingMethod` | 当前状态 |
|---|---|---|
| `cellBased` | `diffuseVelocity` | 既有基线 |
| `cellBased` | `sharpReconstruction` | 实验性 sharp 路线 |
| `markerBased` | `markerDirectForcing` | 第一阶段原型 |

## 旧的 cell-based 路线

旧路线以 SDF 对 Eulerian cells 进行区域划分，并直接在 cells 上
构造目标速度与修正量。

* `diffuseVelocity` 使用平滑的固相体积分数，在有限厚度的过渡带
  内混合流体预测速度与结构速度；
* `sharpReconstruction` 使用较清晰的 inside / forcing layer /
  fluid 分类，在靠近零等值面的 Eulerian cells 上重构目标速度；
* 速度约束、通量和压力修正仍嵌入现有的 `interFoam`-like
  PIMPLE 流程；
* 网格 cell 本身既承担几何分类，也承担强迫自由度，因此实现
  简单，但边界位置、作用宽度和网格方向之间存在较强联系。

## 新的 marker-based 路线

marker-based 路线将 immersed surface 上的 Lagrangian markers 与
Eulerian 流场自由度分开处理。第一阶段原型的主要过程为：

1. 根据 SDF 零等值面建立或更新 surface markers；
2. 从周围 Eulerian cells 向 marker 插值当前流体速度；
3. 在 marker 上计算与结构目标速度之间的 residual；
4. 将 marker correction 通过紧支撑核函数 spread 回 Eulerian
   cells；
5. 继续执行现有的速度、通量与压力修正流程。

这使结构边界的表示不再等同于某一层 Eulerian cells，并为后续
实现 Constant--Pinelli 风格的插值—强迫—扩散算子、伴随算子及
速度—通量—压力迭代耦合留出位置。

## 两条路线的核心区别

| 方面 | `cellBased` | `markerBased` |
|---|---|---|
| 几何自由度 | Eulerian cells | surface markers |
| 速度取样 | cell/邻域重构 | Eulerian-to-marker 插值 |
| 强迫回写 | 直接修改 cells | marker-to-Eulerian spreading |
| 界面厚度 | 由 cell 分类或 SDF 带宽控制 | 由离散 delta 核的支撑宽度控制 |
| 与网格的绑定 | 较强 | 较弱，但仍受插值核和分辨率影响 |
| 当前成熟度 | 已有长程算例 | 第一阶段 operator prototype |

`markerBased` 并不意味着程序要另开一个 solver 项目。两条路线
可以共享 `ichiBoundaryFoam` 的 VoF、几何运动、SDF、载荷诊断、
算例字典以及绝大部分时间推进框架；真正需要逐步分离的是 IBM
约束算子的构造与速度—通量—压力耦合方式。

## 当前边界与下一步测试

目前的 `markerDirectForcing` 仍是第一阶段原型：

* marker stencil 暂时是 processor-local，跨 processor 的插值与
  spreading 尚未完善；
* 插值和 spreading 尚未证明为严格互为伴随，也没有完成离散功
  与载荷守恒验证；
* 压力与通量仍主要沿用旧的 PIMPLE 路线，尚不能等同于完整的
  Constant--Pinelli 迭代耦合；
* SDF 仍用于生成 marker、法向量和邻域尺度，但不再作为直接
  选择全部强迫 cells 的唯一机制。

## 验证方案：两层测试

marker-based 路线不应直接从完整的 Stokes II + VoF +
`kOmegaSST` 算例开始验证。否则即使结果出现问题，也很难判断错误
来自 marker 几何、插值、spreading、压力耦合、VoF，还是湍流模型。

现有 480 s 的 cell-based R1、BF-1、载荷和 slip 结果可冻结为旧路线
baseline；新的 BF-1/R1 长程运行有助于最终比较，但不必等待它们完成
后才继续开发 marker-based 路线。

### 第一层：不求解 N--S 方程的离散算子测试

这里的“不求解 N--S 方程”并不表示不能使用 OpenFOAM。更合适的
实现是编写一个小型 OpenFOAM test utility：

1. 使用 `Time` 和 `fvMesh` 读取普通 OpenFOAM case 与网格；
2. 在内存中构造具有解析表达式的 `volVectorField`；
3. 调用与正式 solver 完全相同的 marker generation、interpolation
   和 spreading 代码；
4. 计算并输出误差、合力和力矩；
5. 随后退出，不建立时间推进循环，也不求解动量方程或压力方程。

因此它仍然可以利用 OpenFOAM 的网格、场、并行通信和维度系统，
只是把 IBM 的离散算子当成可单独验证的软件组件。`potentialFoam`
会额外求解势流方程，既不是这些测试所必需的，也会把 PDE 离散误差
混入 IBM 算子误差；故不建议用它替代这一层测试。

第一层至少应包括：

#### Marker 几何

* marker 是否位于目标圆柱表面；
* 法向量是否为单位外法向量；
* marker 面积权重之和是否等于相应的圆柱表面积；
* 几何平移或旋转后，marker 位置、法向量和刚体速度
  是否正确；
* marker spacing 与局部 Eulerian 网格尺度之比是否处于预期范围。

#### 插值一致性

先在 Eulerian cells 上人为设置常量场

$$
\boldsymbol{u}(\boldsymbol{x})=\boldsymbol{u}_0,
$$

再执行 Eulerian-to-marker 插值。所有 marker 都应恢复同一速度；
该测试直接检查离散核是否满足 partition of unity。若后续采用 MLS
或其他高阶重构，还应测试线性场的再现能力。

#### Spreading 的合力与力矩守恒

给 marker 人为指定已知的力，检查 spread 后是否满足

$$
\sum_i \boldsymbol{f}_i V_i
\simeq
\sum_l \boldsymbol{F}_l\Delta A_l,
$$

以及

$$
\sum_i
(\boldsymbol{x}_i-\boldsymbol{x}_c)
\times\boldsymbol{f}_iV_i
\simeq
\sum_l
(\boldsymbol{X}_l-\boldsymbol{x}_c)
\times\boldsymbol{F}_l\Delta A_l.
$$

前者检查总载荷，后者检查一阶矩和总力矩。

#### 插值与展布的伴随性

对任意测试速度和 marker force，检查

$$
\langle\boldsymbol{u},S\boldsymbol{F}\rangle_E
\simeq
\langle J\boldsymbol{u},\boldsymbol{F}\rangle_L.
$$

这里 $J$ 是 Eulerian-to-Lagrangian interpolation，$S$ 是
Lagrangian-to-Eulerian spreading。该测试用于判断两个算子是否具有
一致的离散功率交换，而不是仅凭最终流场看起来是否合理。

#### 静止零场与并行一致性

* 当 $\boldsymbol{U}=\boldsymbol{U}_b=0$ 时，应得到零 correction、
  零载荷和零 slip；
* 同一算子测试分别串行运行和并行运行，结果应在舍入误差范围内
  一致；
* 必须特别测试 marker stencil 跨越 processor boundary 的情形。

第一层的目的不是验证真实流动物理，而是在进入 CFD 前确认 marker
算子本身没有明显 coding 或离散错误。

### 第二层：逐级增加物理复杂度的 CFD 测试

第一层通过后，再使用 `ichiBoundaryFoam` 求解控制方程。建议按以下
顺序进行：

1. **单相静水、静止圆柱**：初始流体和结构速度均为零，检查是否
   产生寄生速度、虚假载荷、slip 或 divergence；
2. **单相层流绕静止圆柱**：先采用低 Reynolds 数稳态算例，再用
   非稳态算例检查平均阻力、升力和 Strouhal 数，并与 body-fitted
   结果比较；
3. **单相预设运动圆柱**：采用小振幅简谐平移或旋转，检查 marker
   跨越 cells 时载荷是否连续、结构速度是否正确以及 slip 是否
   有界；
4. **两相静水 + 静止圆柱**：让圆柱穿越自由液面但不造波，检查
   密度比约 $10^3$ 时的水量守恒、alpha 有界性、寄生流和界面附近
   forcing；
5. **Stokes II + laminar + 静止圆柱**：与相同网格、固定时间步和
   波浪配置下的 cell-based V1/BF baseline 比较；
6. **Stokes II + `kOmegaSST`**：最后进入 R1；第一轮关闭实验性的
   turbwall treatment，以免同时引入两个未稳定因素。

每一级至少记录：

* marker slip 的 maximum、RMS、normal 和 tangential 分量；
* constraint residual 与 maximum divergence；
* `min/max(alpha.water)` 和总体水量；
* IBM 合力、力矩及其时间连续性；
* 压力迭代次数与计算成本；
* 至少两档 Eulerian 分辨率和两档 marker spacing 下的敏感性。

### 近期执行顺序

1. 验证旧字典省略 `couplingMethod` 时仍使用 `cellBased`；
2. 验证非法 coupling/forcing 组合在启动时明确失败；
3. 将 marker/interpolation/spreading 从 solver 流程中整理为可复用
   组件；
4. 建立上述 OpenFOAM operator-test utility，并先通过常量插值、
   force/torque 守恒和伴随性测试；
5. 再运行最小的单相静止零场 CFD 算例；
6. 依次进入单相绕流、运动结构、两相静水和波浪算例。
