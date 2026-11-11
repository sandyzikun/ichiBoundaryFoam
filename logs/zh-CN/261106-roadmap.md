# `ichiBoundaryFoam` 的验证与对照路线

## 基本原则

后续验证应逐次只增加一类新的物理或数值过程。这样一旦结果发生变化，
才能判断它来自 VoF、immersed-boundary forcing、湍流近壁处理，还是结构
运动本身。

body-fitted 与 IBM 的桩周网格不可能逐 cell 相同。公平对照应建立在：

* 相同的计算域、远场传播网格和垂向分层；
* 相同的波浪、物性、时间步和离散格式；
* 可比的桩周有效分辨率；
* 相同的测点、统计量和收敛标准；
* 两类网格各自的分辨率敏感性分析；

之上，而不是要求两套网格的节点一一对应。

## 阶段 BF-1: body-fitted VoF control

配置：

* 使用 `gen-piled-blockmesh-benchmark.py` 生成的 body-fitted monopile；
* 使用包含 IBM 代码的常规 `ichiBoundaryFoam` binary；
* 运行时设置 `ibm.enabled false`；
* Stokes II wave；
* `kOmegaSST`；
* 固定时间步与既定的 `limitedLinearV 1` 动量对流格式。

理由：

这一阶段完全不施加 IBM forcing，用于确认开发版 solver 在真实水槽中能否
复现稳定的波浪传播、VoF 界面、贴体单桩受力和湍流近壁行为。它是后续
所有 IBM 结果的 control case。

进入下一阶段前至少应确认：

* Courant number 与 `alpha.water` 界限没有持续恶化；
* 水量、测波点水位和压力没有非物理长期漂移；
* 单桩载荷和主要流场统计量可以重复分析。

## 阶段 IBM-L1: laminar static IBM cylinder

配置：

* 使用无圆柱孔洞、无 `monopile` patch 的背景网格；
* 通过 SDF 定义与 BF-1 相同的静止圆柱；
* 在桩周采用与几何无关的局部 Cartesian refinement；
* 单相或 laminar 流动；
* 启用 discrete forcing 与 multi-direct forcing。

理由：

先移除自由液面和湍流模型，单独验证 IBM 最核心的几何识别、速度 forcing、
不可穿透条件与载荷计算。如果这一阶段的 slip、质量守恒或受力仍不正确，
加入 VoF 与 k-omega SST 只会掩盖根因。

重点检查：

* immersed surface 上的法向穿透和切向 slip；
* continuity error；
* forcing 的收敛及其对 MDF corrector 数量的敏感性；
* forcing 体积分载荷与重建表面载荷的一致性；
* coarse/medium/fine 背景网格下的收敛趋势。

## 阶段 IBM-V1: 2-phase VoF with static IBM

配置：

* 保持 IBM-L1 的静止圆柱和背景网格；
* 加入 air-water VoF；
* 先采用已通过 BF-1 验证的波浪和离散设置；
* 暂不加入 immersed turbulence wall treatment。

理由：

这一阶段专门检验 VoF 通量与 IBM forcing 的耦合，包括固体内部的
`alpha.water`、界面附近的通量修正以及 alpha extension 是否必要。它不能
同时承担湍流近壁模型验证，否则无法判断 alpha 异常究竟来自哪一层。

重点检查：

* `min/max(alpha.water)` 及异常位置；
* Phase-1 总体积及其长期漂移；
* 圆柱附近的界面形态与非物理通量；
* 连续压力场在 immersed interface 两侧的行为；
* diffuse/sharp reconstruction 和 alpha extension 的独立影响。

## 阶段 IBM-R1: RANS IBM monopile

配置：

* 在 IBM-V1 通过后加入 `kOmegaSST`；
* 增加 SDF-aware wall distance 或适当的 immersed wall model；
* 为 immersed surface 定义可验证的 `k`、`omega` 和 `nut` 近壁处理；
* 与 BF-1 的压力、速度、尾流和载荷进行比较。

理由：

原生 OpenFOAM wall distance 只识别真实 mesh patches，并不会自动把 SDF
零等值面当作壁面。若直接在 IBM 背景网格上启用原生 `kOmegaSST`，所得
近壁湍流状态与 body-fitted case 并不等价，因此不能据此评价 IBM 精度。

这一阶段通过后，才可把 Stokes II + k-omega SST + static IBM 称为正式的
monopile WSI baseline。

## 阶段 IBM-M1: prescribed moving IBM

配置：

* 从已通过验证的 IBM-R1 或 IBM-V1 出发；
* 先加入保守的匀速平移或单轴旋转；
* refinement zone 覆盖结构的完整运动包络；
* 仍采用单向耦合，即结构运动预先给定。

理由：

运动结构会额外引入时变 SDF、moving-interface forcing、几何守恒与载荷功率
平衡问题。只有静止 IBM 已被验证后，才能把新的误差明确归因于运动算法。

重点检查：

* SDF 与结构轨迹是否严格一致；
* 新进入或离开固体区域的 cells 是否产生 alpha/pressure 瞬变；
* 表面 slip 与载荷是否随时间步和结构速度收敛；
* 流体对结构做功与 IBM forcing 功率是否一致。

## 暂不提前进行的工作

在前一阶段尚未通过时，不应：

* 用 RANS 结果掩盖基础 forcing 问题；
* 用运动圆柱诊断静止 IBM 的 alpha 上溢；
* 仅凭一次计算选择 diffuse 或 sharp interface；
* 将 body-fitted patch force 与未经验证的 IBM volume force 直接等同；
* 把数值稳定等同于物理与载荷精度已经通过验证。

这条路线的目标不是把所有功能一次加入，而是建立一条可以定位误差、能够
回归测试，并最终支持 moving monopile WSI 的开发链。
