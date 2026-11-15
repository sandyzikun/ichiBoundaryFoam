# `ichiBoundaryFoam` marker-based IBM 路线阶段性封存记录

日期：2026-11-15

## 封存结论

本轮开发尝试沿 Constant--Pinelli 思路，将旧 cell-based IBM 后端替换为
marker-based interpolation/spreading，并进一步实现 forcing--pressure
sub-iteration，以同时控制 immersed no-slip residual 与离散连续性误差。

这条路线完成了从 marker geometry、离散算子、独立 tester、marker-space
constraint solve，到小水槽 CFD、sub-iteration 和 120-rank L1-MB preflight
的完整探索。串行和小规模测试证明若干基础构件可以工作，但最终的大规模测试
连续暴露 processor-coupled compact stencil、并行 spreading/interpolation、
stencil conditioning 和 reproducibility 问题。当前实现尚不能作为可靠的生产
后端，继续局部修补的成本已经高于现阶段收益，因此决定阶段性封存。

封存不是否定 marker IBM 或 Constant--Pinelli 方法本身，而是承认当前实现从
cell-based/SDF 原型逐步演化而来的并行数据结构不够稳健。后续应先回到文献，
重新确认适合 OpenFOAM distributed finite-volume mesh 的 marker ownership、
kernel support、pressure coupling 和验证路线，再决定是否重启。

## 开发起点与动机

旧 cell-based 路线使用 SDF、forcing cells、diffuse/sharp velocity forcing 与
MDF。它已能稳定进行 VoF、Stokes waves、静止及 prescribed-motion cylinder
计算，但 immersed surface 附近仍存在明显 slip，尤其是 tangential slip。

问题的核心被归纳为：IBM velocity forcing 与 PIMPLE pressure projection 并非
同一离散约束系统。forcing 可以暂时把速度推向结构速度，随后 pressure
projection 为恢复 divergence-free condition 又可能重新引入 slip。

Constant--Pinelli 路线因此被选作候选方案：

* 使用 Lagrangian surface markers 表示 immersed boundary；
* 以插值算子 $J$ 获得 marker velocity；
* 以展布算子 $S$ 将 marker correction/force 返回 Eulerian mesh；
* 求解 marker constraint；
* 在 forcing 与 pressure correction 之间进行 sub-iteration，同时满足 no-slip
  与 continuity。

SDF 在该路线中没有被完全废弃，仍用于结构几何、signed distance、法向、
表面投影和初始 marker 生成。

## 第一阶段：marker geometry 与离散算子

第一版从 SDF forcing cells 向解析 cylinder surface 投影生成 markers，并采用
compact cosine kernel 构造 Eulerian-to-marker interpolation。初始测试发现，
简单 partition-of-unity normalization 只能精确重构常量场，无法保证线性场。

随后加入 affine first-moment correction，使理想测试网格上的权重满足：

$$
\sum_i w_{li}=1,
\qquad
\sum_iw_{li}(\boldsymbol x_i-\boldsymbol X_l)=\boldsymbol0.
$$

由此，串行规则测试网格上的常量和 affine-linear vector field 可以在浮点误差
范围内重构。二次场被保留为高阶误差表征，不要求精确。

为了避免 processor boundary 截断 local stencil，算子改用 OpenFOAM v2312
的 `centredCPCCellToCellStencilObject` 和 `mapDistribute`。同时实现：

* processor-coupled Eulerian-to-marker interpolation；
* raw algebraic transpose spreading；
* normalized overlap spreading；
* marker area weights；
* cell-volume/marker-area weighted physical force-density spreading；
* total force、torque 与 weighted-adjoint diagnostics。

## 独立 marker operator tester

顶层 `Makefile` 增加第三种构建形态：

```bash
make marker-operator-test
```

对应 executable 为：

```text
ichiMarkerOperatorTest
```

tester 复用生产 `ImmersedGeometry` 与 `DirectForcingIbm`，但不求解 VoF、动量
或压力方程，用于隔离验证：

* marker surface projection 与 normal；
* constant、linear、quadratic interpolation；
* zero interpolation/spreading；
* algebraic transpose；
* surface area；
* force/torque conservation；
* area/volume weighted adjointness；
* matrix-free marker constraint solve。

该 tester 的价值应保留：它确实多次在进入 CFD 前暴露离散算子缺陷，也说明
未来若重启 marker 路线，必须先通过 operator-level verification，不能直接以
完整 CFD 是否运行作为唯一判断。

## Marker-space constraint solve

生产 correction 最初只是把 marker residual 做一次 normalized spreading；这不
等价于真正满足离散 no-slip constraint。随后实现 matrix-free system：

$$
JS\boldsymbol\lambda
=\boldsymbol U_b-J\boldsymbol U,
$$

并以

$$
\delta\boldsymbol U=S\boldsymbol\lambda
$$

形成 Eulerian velocity correction。

求解器采用预条件 CG，初版 Jacobi preconditioner 为：

$$
M_{ll}^{-1}
=\left(\sum_iw_{li}^2\right)^{-1}.
$$

在规则小网格的独立测试中，constraint residual 可以降到严格 tolerance，但常需
约 $O(10^2)$ iterations。该结果说明算子可以形成可求解系统，同时也提前暴露
了 $JS$ conditioning 与预条件性能问题。

## Forcing--pressure sub-iteration

在 marker constraint solve 之后，新增隔离 binary：

```bash
make marker-subiteration
```

输出：

```text
ichiBoundaryFoam-markerSubIter
```

该版本复用此前的 IBM-aware continuous-pressure projection，并实现外层循环：

1. marker constraint solve；
2. Eulerian velocity correction；
3. continuous-pressure correction；
4. velocity/flux update；
5. 同时检查 constraint residual 与 continuity residual；
6. 未同时收敛时继续 sub-iteration。

小水槽测试 `dampless-marker-cfd1/` 先验证零场保持，随后加入速度扰动。串行
短 CFD 能运行至 `End`，marker constraint 与 continuity residual 均可达到当前
配置阈值。这证明 sub-iteration 初稿能够执行，但不能证明其大规模并行实现、
长期稳定性或 decomposition independence。

## Reproducibility 异常

相同配置的 marker-based 小规模并行测试很快出现 run-to-run 差异。初期曾将
其解释为 decomposition、MPI collective reduction 和浮点加法非结合性。这些
因素确实会造成末位差异，例如：

* marker CG 的 global dot products；
* 多 ranks 对同一 owner cell 的贡献求和；
* MPI reduction tree；
* residual 接近 tolerance 时相差一次迭代。

但这种一般解释不足以说明：长期 `interFoam` 与 cell-based IBM 计算没有表现出
类似的快速分叉，而 marker-based 在两三步内就出现可观察差异。

后续 120-rank 失败表明，早期差异更可能受到未初始化 compact slots 污染：

* 未启用 `FOAM_SETNAN` 时，空槽可能读取堆内存残值；
* 启用 `FOAM_SETNAN=true` 后，同一缺陷转化为 FPE；
* 加入 availability 检查后，它稳定报告相同 ranks、markers 与 compact labels。

因此，早期 marker reproducibility 结果不应被视为该方法固有的随机性，也不应
使用 random seed 解释。当前算法没有 stochastic sampling；真正的问题是并行
数据语义和浮点 reduction/停止分支。旧的 reproducibility 数据应作废。

## L1-MB 与 compact-halo 问题

`dampless-ibm-l1-mb/` 继承单相 laminar L1 的大水槽、cylinder、物理模型和
120-rank hierarchical decomposition：

```foam
numberOfSubdomains 120;
method hierarchical;
n (40 3 1);
```

最初两次运行在时间循环前确定性 FPE，调用栈落在 marker stencil rebuild。
第一轮修补扩展 compact arrays 至 `map.constructSize()` 并加入 availability
检查。再次运行后，两份作业在 ranks 57--62 上报告相同 unavailable compact
labels，证明问题是确定性 map/layout 错误，而非 stochastic failure。

根因是把 CPC extended-stencil compact addressing 误当成纯 internal-cell map。
实际 compact field 还包含 physical-boundary、processor-patch 和 remote-halo
slots。只扩大数组而不按 patch compact offset 写入 boundary field，仍然错误。

随后将 gather 重构为 volField-aware 路径：

1. 填充 internal values；
2. 按 patch offset 填充 boundary values；
3. processor patch 交换 neighbour values；
4. `mapDistribute` 填充 remote halo；
5. centres、velocity interpolation、cell volumes 与 force spreading 使用统一语义。

重构后，串行完整 operator suite 和小水槽 CFD 再次通过。为了避免继续用正式
CFD 充当 map tester，L1-MB submission 又增加了同一 120-rank decomposition
上的 `ichiMarkerOperatorTest` preflight。

## 最终 120-rank preflight 结果

2026-11-15 的两份 preflight 日志结果一致：

```text
markers                     = 2880
fallback stencils           = 0
affine fallbacks            = 2784
negative-weight markers     = 96
max first-moment error      = 0.002669010895
```

两份运行均得到：

* geometry surface distance：PASS；
* unit/stored normal：PASS；
* constant interpolation：PASS，误差约 $2.0\times10^{-15}$；
* linear interpolation：在当前宽松阈值下 PASS；
* linear RMS error：约 $3.19\times10^{-3}$；
* linear maximum error：约 $6.42\times10^{-3}$；
* 随后在下一段并行 operator 测试中触发 FPE，未完成 suite。

两次结果数值一致，说明该失败仍是确定性的。更重要的是，2880 个 markers 中
有 2784 个发生 affine correction fallback；这意味着即便继续修复当前 FPE，
L1 网格上的 stencil geometry 也没有稳健支持设计中的 3D affine correction。
规则小网格测试所得的近机器精度线性重构，不能外推到实际 L1 decomposition。

因此当前障碍不再只是单个越界或漏填槽位，而包括：

* marker ownership 与跨 processor 支持域；
* compact gather/reverse accumulation 的完整语义；
* 大量 rank 上没有 local markers 时的 collective/operator 行为；
* affine correction 的秩与 conditioning；
* negative weights 与 $JS$ conditioning；
* transpose/adjoint 与 conservation 在实际 decomposition 下的稳健性；
* reproducible reductions 与 tolerance branch amplification。

## 本轮可保留成果

尽管生产路线失败，本轮工作仍留下可复用成果：

* marker operator tester 及其测试分类；
* affine moment consistency 的需求与失败案例；
* interpolation/spreading、algebraic transpose 和 physical adjoint 的区分；
* marker-area/cell-volume weighted force spreading；
* matrix-free marker constraint formulation；
* forcing--pressure sub-iteration 的控制框架；
* compact-halo diagnostics 与 availability checks；
* 小水槽与 L1-MB 分层验证算例；
* 对 reproducibility、MPI reduction 和未初始化数据污染的经验；
* 一套明确证明“小规模 PASS 不等于大规模并行正确”的失败记录。

这些代码和日志应放入独立 archive branch，而不是删除。未来若在文献支持下采用
新的 marker mesh、kernel、distributed ownership 或 monolithic/Schur coupling，
仍可从 tester、diagnostics 和失败样例中复用部分基础设施。

## 封存后的开发方向

主线应回到 marker 路线开始前的稳定 VoF--cell-based IBM baseline。近期不再把
“彻底消灭 slip”作为阻塞其他工作的前置条件，而是：

1. 保留旧 cell-based IBM 作为可运行研究基线；
2. 明确报告其 diffuse constraint layer 和 residual/slip 局限；
3. review immersed-boundary forcing/pressure coupling 文献；
4. 对比 direct forcing、continuous forcing、ghost-cell、cut-cell、BDIM、IBPM、
   distributed Lagrange multiplier 和 fictitious-domain 等路线；
5. 特别关注适用于 collocated finite-volume、parallel unstructured mesh、VoF
   free surface 与 moving rigid body 的方法；
6. 在确定新的离散结构和验证标准后，再决定是否重启 marker backend。

本轮最终状态应概括为：

> Constant--Pinelli marker-based IBM 的离散算子、constraint solve 和
> forcing--pressure sub-iteration 初稿已经完成；规则小网格与串行短 CFD
> 可以工作，但实际 120-rank L1 preflight 暴露了尚未解决的并行 operator
> 与 stencil-conditioning 问题。当前实现不进入生产主线，完整代码、tester、
> 算例设计与失败日志转入 archive branch，主线回到稳定 baseline，并先以
> 文献 review 寻找下一条更可靠的技术路线。
