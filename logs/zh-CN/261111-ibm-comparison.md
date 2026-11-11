# IBM 实现路径比较与 `ichiBoundaryFoam` 的位置

## 核查范围

本次对以下三类路线进行了区分和比较：

* 当前 `ichiBoundaryFoam` 采用的 SDF-based discrete forcing IBM；
* foam-extend 4.1 官方源码中的 Immersed Boundary Surface（IBS）；
* 文献和研究项目中的 `interIbFoam` / `interDyMIbFoam` 两相实现。

需要特别说明：在当前 foam-extend 4.1 官方源码树中，正式提供的
IBM solver 是单相不可压缩的 `pimpleDyMIbFoam`。本次没有在官方源码
树中找到 `interIbFoam` 或 `interDyMIbFoam`；后二者存在于相关论文、
博士项目和研究代码中，不能直接当作 foam-extend 4.1 官方维护的
开箱即用 solver。

## 三类实现路径

### 当前 `ichiBoundaryFoam`

目前的主要路线为：

* 以解析 Signed Distance Function（SDF）表示 capped cylinder；
* 根据 SDF 区分固体内部、界面带与外部流体；
* 在 SDF forcing band 内施加 discrete direct forcing；
* 通过 Multi-Direct Forcing（MDF）反复降低结构表面的速度误差；
* 在标准两相 VoF/PIMPLE 求解流程中进行压力修正；
* 支持 `noSlip` / `freeSlip`、静止结构及预设主动运动；
* 已有 constraint-load、slip 和其他诊断；
* 已实验性加入 PPMDF、IBM-aware projection、sharp reconstruction、
  alpha extension 和 immersed turbulence-wall treatment 等可选功能。

其基础思想仍然是：

```text
标准 VoF 动量预测
        ↓
SDF forcing band 内施加 MDF
        ↓
压力/通量修正
        ↓
得到满足连续性但可能重新产生 slip 的速度场
```

它的优势是容易与 OpenFOAM-com v2312、Stokes waves、VoF 和现有水槽
算例结合，也便于快速测试主动运动结构。当前最明显的问题是速度约束
与压力投影之间并非同一个耦合离散系统，因此压力修正可能重新产生
immersed-boundary slip。

### foam-extend 4.1 官方 IBS

foam-extend 4.1 官方实现更接近 modified cut-cell / Immersed Boundary
Surface（IBS），不是我们的 MDF 路线。其主要过程为：

* 读取任意三角面 `triSurface`；
* 分类 fluid cells、solid/dead cells 和 cut cells；
* 计算 cut cell 中的湿润体积及其中心；
* 计算被结构截断后的流体 face 面积；
* 在结构表面构造具有面心、面积矢量和法向的虚拟 immersed patch；
* 通过专用 patch fields 修改有限体积矩阵系数；
* 让速度、压力、通量和湍流壁面条件共同识别该内部边界。

背景网格的拓扑连接不需要改变，但求解器使用的 cell volume、face
area、几何中心、边界贡献和矩阵系数会识别 immersed surface。因此，
结构不只是一个“事后修改速度的区域”，而是离散方程中的实际边界。

官方源码同时包含：

* `movingImmersedBoundaryVelocity`；
* moving immersed surface 和动态局部加密；
* `immersedBoundaryOmegaWallFunction`；
* `immersedBoundaryEpsilonWallFunction`；
* `immersedBoundaryNutWallFunction`；
* `immersedBoundaryKqRWallFunction`。

### 研究项目中的 `interIbFoam` / `interDyMIbFoam`

相关研究把 IBM 与自由液面 VoF 结合，并进一步加入动态 immersed body：

* `interIbFoam`：两相不可压缩 VoF + IBM；
* `interDyMIbFoam`：在前者基础上加入动态结构、运动求解和 IB mask
  更新，并用于 iceberg-calving wave 等问题。

这些工作证明了 foam-extend 系 IBM 可以扩展到两相运动结构，但它们
不是当前官方源码树中的标准 solver，也不能由其存在直接推断它们已经
解决了水气密度比、alpha extension、结构穿越自由液面及全部载荷振荡
问题。

## 核心差异

| 层面 | `ichiBoundaryFoam` | foam-extend 4.1 IBS |
|---|---|---|
| 几何 | 解析 capped-cylinder SDF | 任意三角面 `triSurface` |
| 固液分类 | SDF 符号与 forcing band | fluid/dead/cut cells |
| 边界表示 | 若干层 cell 内的强迫区 | 显式虚拟 immersed patch |
| 速度条件 | direct forcing + MDF | patch BC + 矩阵系数修改 |
| 压力与通量 | 标准投影；实验性 aware projection | cut faces 和 patch 已进入离散算子 |
| 湍流壁面 | 实验性 SDF-band treatment | 专用 immersed wall functions |
| 载荷 | forcing/constraint impulse 诊断 | 可在虚拟表面积分压力和剪切应力 |
| 运动结构 | 解析 SDF prescribed motion | 移动三角面、mask 更新和动态加密 |
| 双向 6DoF | 尚未完成 | 基础设施更成熟，具体能力依 solver |
| 两相 VoF | 已集成，是本项目核心 | 官方 IBM solver 目前为单相 |
| Stokes waves | 已集成 OpenFOAM-com 波浪功能 | 官方 IBM solver 未覆盖本项目配置 |

## 当前实现尚未覆盖的 foam-extend 能力

### 矩阵级内部边界

这是最根本的差距。当前 MDF 在预测速度上施加结构约束，之后压力
投影为了恢复不可压缩性又会改变速度。IBS 则把结构边界写入有限体积
矩阵和边界条件，使压力、速度和通量从离散层面共同知道此处存在固壁。

这可以解释当前长期存在的现象：

* MDF 能明显降低约束前的速度误差；
* 标准压力修正后 slip 再次增大；
* PPMDF 和当前 IBM-aware projection 只能在现有 forcing 架构内折衷，
  尚未等价于完整的虚拟 patch 离散。

### 成熟 immersed turbulence-wall conditions

foam-extend 的 `k`、`omega`、`epsilon` 和涡粘性条件以 immersed patch
边界条件参与矩阵。当前试验性 turbwall treatment 主要在 SDF band
邻近 cells 上近似恢复壁面变量，缺少对应的显式壁面几何和 patch
离散，这也是该功能暂未达到预期效果的重要原因。

### 表面载荷积分

当前 constraint load 主要由实现速度约束所需的动量改变量得到。它是
有用的 IBM 载荷诊断，但不完全等价于严格的表面应力积分：

$$
\mathbf F
=
\int_{\Gamma_{\mathrm{IB}}}
\left(-p\mathbf I+\boldsymbol\tau\right)
\mathbf n\,\mathrm dS.
$$

IBS 具有显式虚拟 surface faces，因此更容易分别计算压力载荷与粘性
剪切载荷。这对以后进行双向 WSI/FSI 尤其重要。

### 任意几何

解析 SDF 对规则圆柱及其主动运动非常高效，但 foam-extend 的三角面
方案更适合不规则结构、倾倒后的复杂姿态或破坏几何。相应代价是 cut
geometry、并行搜索、退化交切和移动表面更新都更加复杂。

## 两相 VoF 方面的判断

foam-extend 在内部边界离散方面明显更加成熟，但不能因此认为它已经
全面解决本项目的两相问题。官方 `pimpleDyMIbFoam` 是单相 solver；研究
项目中的两相实现也仍需面对：

* 水气密度比约为 1000 时的压力与动量敏感性；
* immersed surface 附近的相分数定义与延拓；
* VoF face flux 与固体边界通量的一致性；
* 结构运动穿越自由液面时的相体积守恒；
* 自由液面、immersed interface 和湍流壁面的三重耦合。

因此，`ichiBoundaryFoam` 已完成的 OpenFOAM-com v2312 两相 VoF、Stokes
wave 和现有水槽兼容性仍是实际优势，不能简单用官方单相 IBS 替代。

## 对后续技术路线的建议

目前不宜整体移植 foam-extend IBS。完整移植需要深入修改 OpenFOAM-com
的 mesh geometry、patch field 和有限体积离散层，工作量和回归风险都
很大。更现实的渐进路线为：

1. 由 SDF 与背景 cell/face 的交点构造显式 immersed surface faces；
2. 先利用这些 faces 分别积分压力和剪切载荷；
3. 让速度边界条件逐步通过矩阵系数施加，而非完全依靠事后 forcing；
4. 让压力和通量修正识别 solid-cut faces；
5. 在显式 immersed surface 稳定后，再实现 `k`、`omega` 和 `nut`
   immersed-wall conditions；
6. 最后用可信的表面载荷驱动结构 6DoF，实现双向 WSI/FSI。

这一比较得到的核心认识是：

> foam-extend 最值得借鉴的并不是增加 MDF 次数，而是把 immersed
> boundary 从“速度修正区域”提升为有限体积离散中的真实内部边界。

同时也不应把 IBS 当作完美答案。移动 immersed surface 仍可能产生
spurious force oscillations；对于较轻或较小的运动物体，这些载荷振荡
仍可能影响运动轨迹。因此后续仍需要 body-fitted baseline、网格/时间步
收敛及载荷频谱验证。

## 参考资料

* [foam-extend 4.1 官方源码树](https://sourceforge.net/p/foam-extend/foam-extend-4.1/ci/master/tree/)
* [Verification and validation of a finite volume immersed boundary method for static and moving geometries](https://doi.org/10.1016/j.jnnfm.2021.104510)
* [A critical assessment of the immersed boundary method for modeling flow around fixed and moving bodies](https://www.sciencedirect.com/science/article/pii/S004579302300066X)
* [Benchmark simulations of flow past rigid bodies using an open-source, sharp interface immersed boundary method](https://collaborate.princeton.edu/en/publications/benchmark-simulations-of-flow-past-rigid-bodies-using-an-open-sou/)
* [A two-phase immersed-boundary application to iceberg-calving waves](https://doi.org/10.1016/j.advwatres.2020.103647)
