# `ichiBoundaryFoam` 开发后记：marker sub-iteration、并行 compact halo 与可重复性

日期：2026-11-14（晚间）

## 后记结论

在 `261114.md` 完成后，Constant--Pinelli marker-based 路线继续推进到了
forcing--pressure sub-iteration、小水槽 CFD 和 L1-MB 大水槽并行测试。晚间
最重要的结论不是“L1-MB 已通过”，而是：**串行及小规模 operator/CFD 测试
尚不足以证明 marker 并行后端正确；120-rank L1-MB 暴露了 CPC compact
addressing 中 boundary/processor slots 没有被完整装配的问题。**

该问题同时为此前 marker-based 短测试中异常快速出现的 run-to-run 差异提供了
更可信的解释。它们不能再简单归因于“一般 MPI 浮点误差”：在未启用
`FOAM_SETNAN` 时，未初始化 compact slots 可能带入不同的堆内存残值；启用
`FOAM_SETNAN=true` 后，同一问题转化为确定性的 FPE；加入显式 availability
检查后，则稳定报告相同 ranks 和 compact labels。

截至本文结束时，compact-field 通路已按 OpenFOAM extended-stencil 的
internal、boundary、processor-patch 和 remote-halo 语义重构，独立 binary
已重新编译，串行完整 operator suite 与小水槽 CFD 回归均通过。下一关不是
直接宣布 L1-MB 成功，而是在与 L1-MB 完全相同的 120-rank `(40 3 1)`
decomposition 上先通过新增的 marker-operator preflight，再进入正式 CFD。

## Forcing--pressure sub-iteration 的实现

在此前独立 marker constraint solve 的基础上，新增了
Constant--Pinelli 风格的交替约束--压力修正循环。每次外层 iteration 执行：

1. 对当前速度计算 marker residual；
2. matrix-free 求解

   $$
   JS\boldsymbol\lambda
   =\boldsymbol U_b-J\boldsymbol U;
   $$

3. 以

   $$
   \delta\boldsymbol U=S\boldsymbol\lambda
   $$

   修正 Eulerian velocity；
4. 求解连续压力修正并更新 velocity/flux；
5. 同时检查 marker constraint residual 与 continuity residual；
6. 两者均达到配置阈值后退出，否则继续交替迭代。

当前控制项包括：

```foam
algorithm             constantPinelliSubIteration;
constraintTolerance   1e-5;
continuityTolerance   1e-9;
minimumIterations     1;
maximumIterations     40;
requireConvergence    true;
```

该实现复用了此前的 IBM-aware continuous-pressure projection，而没有重新复制
一套压力方程。marker forcing 和 pressure correction 通过外层 sub-iteration
交替调用，从而避免旧 predictor-only 路线中“forcing 满足 no-slip、随后压力
投影又重新引入 slip”的单向覆盖关系。

为了隔离该开发形态，顶层 `Makefile` 提供：

```bash
make marker-subiteration
```

其 binary 为：

```text
platforms/$WM_OPTIONS/bin/ichiBoundaryFoam-markerSubIter
```

它不覆盖常规 `ichiBoundaryFoam`。随后建立的 `dampless-marker-cfd1/` 用于
小水槽、单相、静止 cylinder 的短 CFD 测试；`dampless-ibm-l1-mb/` 则继承
L1 的大水槽、物理模型和 120-rank hierarchical decomposition，用于正式的
marker-based L1 层验证。

## 小水槽 CFD 阶段

`dampless-marker-cfd1/` 先验证零场保持，后改为带初始扰动的测试。该阶段确认：

* marker constraint CG 可以收敛；
* forcing--pressure sub-iteration 能同时降低 constraint 与 continuity residual；
* 串行短 CFD 能够推进到 `End`；
* 当前 Jacobi 预条件器可以工作，但每次 constraint solve 仍需约
  $O(10^2)$ iterations，性能与 conditioning 仍待改进。

两次相同配置的小规模并行运行产生了细小但可见的 run-to-run 差异。当时曾将
其主要解释为 decomposition、MPI collective reduction 和浮点非结合性。
这一解释只能说明末位差异的来源，却不能充分说明为什么 marker-based 路线
在两三步内就出现可观察分叉，而长期运行的 `interFoam` 与 cell-based IBM
没有表现出类似现象。因此该判断现已修正：普通并行舍入可能是次要来源，但
不能作为此前快速分叉的首要解释。

## L1-MB 的两轮确定性失败

`dampless-ibm-l1-mb/` 保留了 L1 的单相 laminar cylinder 配置，并使用：

```foam
numberOfSubdomains 120;
method hierarchical;
n (40 3 1);
```

最初两次运行在进入时间循环前确定性 FPE。调用栈均落在
`DirectForcingIbm::rebuildMarkerStencils()`，主要涉及 ranks 57--62。第一轮
修补将 compact arrays 从 local `nCells` 扩展到 `map.constructSize()`，并加入
availability 检查。重新运行后不再得到含糊的 FPE，而是两次都报告相同错误，
例如：

```text
Marker stencil references unavailable compact cell data:
marker = 0, seed cell = ..., compact label = ...,
construct size = 16336, availability = 0
```

两次运行中的 ranks、marker 和 compact labels 一致。这证明失败不是随机采样、
CFD 发散或 stochastic solver 行为，而是确定性的并行数据布局错误。

第一轮修补仍然不完整，因为它错误假定 compact map 的前 `nCells` 个元素就是
全部有效输入。OpenFOAM 的 `centredCPCCellToCellStencilObject` 使用
extended cell stencil；其 compact field 包含：

* internal-cell values；
* physical-boundary values；
* processor-patch values；
* 经 `mapDistribute` 收集的 remote values。

只扩大数组而不按 patch compact offset 写入 boundary values，仍会留下空槽。
因此不能通过关闭 availability 检查、把空槽设零或跳过该 stencil 来“让它跑
起来”；那样即使不崩溃，插值、展布、Schur operator 和 force conservation
仍可能是错误的。

## 异常“随机性”的重新归因

当前没有任何基于 probability distribution 的 marker sampling，也不存在可由
`random-state` 或 seed 控制的随机数。此前现象应分为两类。

### 未初始化 compact slots

这是目前最能解释“marker-based 两三步即分叉”的主因：

* 未开启 `FOAM_SETNAN` 时，未写入的 slots 可能保留堆内存残值；
* 不同 process layout、allocation history 或运行时序可使残值不同；
* marker interpolation、CG 和 sub-iteration 会快速放大该污染；
* 开启 `FOAM_SETNAN=true` 后，污染变成可捕获的 NaN/FPE；
* availability 检查进一步把它转成带 rank/marker/label 的确定性错误。

cell-based IBM 不使用该 marker compact interpolation/spreading backend，因此
长期 `interFoam` 与 cell-based 经验不能用来排除这个 marker 专属缺陷。反过来，
也不应再用“旧路线其实也有，只是不明显”来解释此前结果。

### 正常的并行浮点非 bitwise reproducibility

即使 compact halo 完全正确，以下操作仍可能造成末位差异：

* 多个 ranks 对同一 owner cell 的浮点贡献求和；
* marker CG 中的 global dot products；
* MPI reduction tree 的加法顺序；
* CG 与外层 sub-iteration 在 tolerance 附近相差一次迭代。

这是浮点加法不满足结合律导致的，不是随机采样，固定 seed 对它无效。稳定的
CFD 系统也可能耗散而非单调累积舍入扰动，所以总 steps 很多不意味着必然出现
明显偏差。但是，正常舍入应从机器精度附近开始；若修复后的 marker 路线仍在
两三步内形成明显物理解差异，就必须继续按 ill-conditioning、marker ownership、
非伴随 spreading 或 tolerance branch amplification 排查，不能再笼统归因于
“MPI 本来就不确定”。

## Compact-halo 重构

晚间已将原来的 internal-only helper 替换为 volField-aware compact gather。
新通路遵循 OpenFOAM extended-stencil 的装配方式：

1. 写入 internal field；
2. 按

   ```text
   patch.start - nInternalFaces + nCells
   ```

   计算每个 boundary patch 的 compact offset；
3. 写入 physical-boundary 与 processor-patch field values；
4. processor patches 先完成相邻 rank 数据交换；
5. 再由 `mapDistribute` 填充 remote compact slots。

该规则已统一用于：

* cell-centre gather 与 marker stencil construction；
* Eulerian-to-marker velocity interpolation；
* matrix-free $JS$ operator；
* marker-to-Eulerian reverse accumulation；
* remote cell-volume gather 与 physical force-density spreading。

availability 检查被保留为防御性诊断；它不再负责“修复”空槽，而是在 map
语义再次被破坏时立即报告具体 marker、seed、compact label 和 construct size。

重构后重新编译了：

```text
ichiBoundaryFoam-markerSubIter
ichiMarkerOperatorTest
```

没有覆盖常规 `ichiBoundaryFoam`。

## 重构后的本地验证

重构完成后，串行 `dampless-marker-cfd1` 扰动测试重新运行至 `End`。最终一组
diagnostics 中：

* marker constraint solve 达到配置 tolerance；
* maximum no-slip constraint residual 约为 $9.8\times10^{-6}$；
* maximum divergence 约为 $1.7\times10^{-11}$。

完整 `ichiMarkerOperatorTest -test all -strict` 也全部通过：

| 测试 | 结果 |
| --- | ---: |
| constant interpolation | PASS，误差约 $2.5\times10^{-15}$ |
| linear interpolation | PASS，maximum error 约 $5.6\times10^{-16}$ |
| zero interpolation/spreading | PASS |
| marker constraint solve | PASS，maximum residual 约 $1.5\times10^{-11}$ |
| algebraic transpose | PASS |
| surface area | PASS |
| force conservation | PASS，误差约 $1.6\times10^{-15}$ |
| torque conservation | PASS，误差约 $2.6\times10^{-16}$ |
| weighted adjoint | PASS，误差约 $2.2\times10^{-16}$ |

这些结果只证明重构没有破坏串行算子与短 CFD。由于登录节点的 PMIx listener
不能启动，本地环境无法执行真正的多 rank 回归；因此本文没有把 120-rank
后端标记为“已通过”。

## L1-MB 新增 120-rank preflight

`dampless-ibm-l1-mb/sbatch-main.sh` 已改为在 `decomposePar` 与 `checkMesh`
之后、`setFields` 与正式 solver 之前执行：

```bash
mpirunApp ichiMarkerOperatorTest \
    -fileHandler collated \
    -test all \
    -strict
```

它使用与正式 L1-MB 完全相同的 mesh、120 ranks 和 `(40 3 1)` decomposition。
输出为：

```text
logmpi-ichiMarkerOperatorTest.txt
```

任何 compact-halo availability、interpolation、spreading、transpose/adjoint、
conservation 或 constraint-solve 失败都会在进入 CFD 前终止作业。这样不再用
昂贵的 L1 时间推进充当最基本的并行 map tester。

## Reproducibility 的目标与测试方法

当前应明确区分：

1. **数值可重复**：物理解与 diagnostics 在规定 tolerance 内一致；
2. **同 decomposition bitwise reproducible**：复用同一 processor mesh、相同
   ranks、MPI 与 rank mapping 时逐位一致；
3. **跨 decomposition/硬件 bitwise reproducible**：不同 ranks 或设备仍逐位
   一致，该目标需要 reproducible reductions 等额外代价。

当前开发阶段至少应要求第 1 项，并争取第 2 项。未来可加入类似 DL 框架
deterministic algorithms 的调试模式，但它不是 random seed：

```foam
deterministicMarkerOperations true;
```

候选实现包括：

* 为 markers/stencils 建立稳定 global IDs 并排序；
* 跨 rank contributions 按 source rank/global cell ID 的固定顺序累加；
* 为 CG dot products 使用固定 rank 顺序或 compensated/reproducible reduction；
* 调试模式下固定 CG 与 projection iteration 数，隔离 tolerance branch；
* 每个 operator stage 输出 field/checksum diagnostics。

该模式会增加通信、排序与计算开销，主要用于验证和回归，不必默认用于所有
正式长算例。

## 下一步

下一阶段应严格按以下顺序进行。

1. **提交更新后的 L1-MB，但先只看 preflight。** 只有
   `logmpi-ichiMarkerOperatorTest.txt` 全部 PASS，才解释后续 CFD 日志。
2. **检查 120-rank operator invariants。** 重点检查 marker count、fallbacks、
   negative weights、first-moment error、transpose/adjoint、force/torque conservation
   与 constraint residual。
3. **执行短 reproducibility pair。** 复用同一 decomposed mesh，固定 nodes、
   ranks 与 mapping，各运行 3--5 steps，不重新 `decomposePar`。
4. **定位首个差异阶段。** 比较每步的 marker interpolation、CG initial/final
   residual、CG iterations、outer sub-iterations、pressure residual、field checksum。
5. **隔离 stopping-branch amplification。** 若仍分叉，临时固定 CG 与外层
   sub-iteration 次数；若差异消失，则主因是 tolerance 分支，若仍存在则继续查
   scatter/reduction/ownership。
6. **评估 $JS$ conditioning。** 记录 iteration count 和 wall time，再决定是否
   优化 Jacobi 预条件器；不能在并行算子正确性尚未确认前先做性能优化。
7. **最后才延长 L1-MB。** 短测试可重复、约束与连续性同时收敛后，再进行
   L1 长时间积分，并与 cell-based L1 比较 slip、loads 和成本。

在这些关口通过前，最准确的项目状态是：

> Constant--Pinelli forcing--pressure sub-iteration 初稿和串行 marker operator
> 已实现；此前 120-rank 失败揭示的 compact-halo 数据语义已重构；新的
> 120-rank operator preflight 与 reproducibility pair 仍是进入正式 L1-MB
> 之前必须通过的验证节点。
