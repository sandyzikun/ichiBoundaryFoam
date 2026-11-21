# Roulund Test 1 body-fitted 基线：结果与推进判断

日期：2026-11-18

## 1. 本次工作的定位

本次完成了 `roulund-t1-bf/` 的首轮计算与 `ichiroulundplot t1` 后处理，目的
是先建立 Roulund, Sumer et al. (2005) Test 1 的 body-fitted 参照，再判断
是否可以进入 `roulund-t1-ib/` 的 cell-based immersed-boundary 对照。

这里必须区分两种结论：

* **数值运行是否正常、主要流动结构是否形成**；
* **是否已经逐项、定量复现了实验 benchmark**。

当前 T1-BF 对前者的回答是肯定的，但尚不能把后者写成已经完成。

## 2. 论文中下游环流的原始事实

论文 §5.1 对 Fig. 10 的说明非常关键：

* Fig. 10(a) 是作者在真实线下水槽中对 Test 1（smooth rigid bed）得到的
  **实验测量结果**，其中桩后平均环流为 anticlockwise；
* Fig. 10(b) 是论文自身数值模型的 steady-state 结果；该模型**未能捕捉**
  Test 1 实验中的 anticlockwise circulation；
* 在 Test 2 的 rough-bed 实验中，桩后平均环流改为 clockwise，而论文模型
  能够捕捉该方向；
* 作者对 Test 1 的这一差异明确写道：`No clear explanation has been found
  for this behaviour.`

因此，Test 1 的 anticlockwise lee-side circulation 不是概念图或作者主观推测，
而是实体水槽实验观测；但连论文用于提出该 benchmark 的 RANS 数值模型也没有
将它复现出来。当前 T1-BF 的桩后流态更接近论文自己的数值结果，并不能据此
判定我们的 solver 或算例实现错误。

## 3. 当前 T1-BF 得到的结果

计算顺利推进至约 `119.996 s`，后段平均场变化已经较小。当前图像显示：

### 3.1 Undisturbed inflow profile

入口剖面已经形成合理的近床边界层，上部速度约为 `U/U0 = 1`。由现有剖面
估计的 friction velocity 约为：

```text
u_tau ~= 1.12 cm/s
```

论文 Test 1 的 log-law 结果为 `1.30 cm/s`，当前值约低 `14%`。相应远场床面
剪切也偏低约 `20%`。这说明入口基线尚未严格校准到实验，但不是发散或错误
流场。

### 3.2 Cylinder 附近速度场

当前结果已经产生 benchmark 所关心的主要结构：

* 桩前近床减速、分离和逆流；
* 桩前 horseshoe-vortex 对应的 downwash/upwash；
* 桩后明显的 recirculation；
* 随高度变化的水平速度亏损与垂向速度结构。

Fig. 11/12 的剖面连续且物理趋势合理；其中上游近床逆流、桩后速度亏损和
垂向输运都清晰可辨。它们还没有叠加论文实验数据，因此目前只支持定性和
内部一致性判断，不能据此宣称定量 validation 已完成。

### 3.3 Bed shear

`ichiroulundplot` 已修正 OpenFOAM `wallShearStress` 的作用方向：OpenFOAM 输出
wall-on-fluid shear，而论文使用 fluid-on-bed shear，故绘图时需要反号。修正后：

* 远上游床面剪切为正；
* 在桩前分离点附近过零；
* horseshoe-vortex 区域出现负剪切；
* 整体符号与 Fig. 16 的物理定义相符。

### 3.4 Fig. 10 的独立采样问题

现有 `roulundSymmetryPlane` 取在精确的 `y=0`，与 body-fitted 网格分块面重合，
导致已有 VTP 仅提取到圆柱两侧的狭窄区域，未覆盖目标
`-2 <= (x-xc)/D <= 2`。因此 Fig. 10 中间的空白是 cutting-plane 采样缺失，
不是 CFD 流场缺失。

模板已把平面偏移到 `y=1e-4 m`，避免后续计算重复这一问题。旧计算若要补全
Fig. 10，只需从最终场离线重新采样，不需要重跑 120 s CFD；本阶段暂不执行。

## 4. 当前结论：T1-BF 是否“没问题”

可以作如下分层判断：

1. **运行与物理结构层面：通过。** 算例稳定完成，入口边界层、horseshoe
   vortex、桩前逆流、桩后回流和床面剪切换向均已出现，没有证据表明 BF
   求解流程存在致命错误。
2. **作为 IBM 对照基线：可以使用。** 当前 BF 结果足以进入同一几何、网格
   尺度、湍流模型、入口条件与采样定义下的 T1-IB 计算，以比较 body-fitted
   与 cell-based IBM 的差异。
3. **作为严格实验 validation：尚未完成。** 入口 `u_tau` 仍偏低，Fig. 10
   采样需要补全，且 Fig. 11/12/16 尚未数字化叠加实验曲线与计算误差指标。
4. **不能把 Test 1 桩后环流方向单独作为否决条件。** 论文自身模型也未能
   复现实验中的 anticlockwise circulation，作者亦未给出明确原因。要求我们
   的 RANS BF 基线必须先解决原论文未解决的这一点，才允许测试 IBM，并不
   合理。

一句话概括：

> T1-BF 是一个运行正常、主要物理结构可信、足以充当 IBM 对照的首轮数值
> 基线；它不是已经完成全部实验定量复现的 benchmark validation。

## 5. 下一步

下一步可以推进 `roulund-t1-ib/`，但比较时应保持控制变量，并优先回答：

* IBM 是否得到与 BF 相同方向和相近尺度的 horseshoe/lee-side circulation；
* Fig. 11/12 的速度剖面差异集中在何处；
* Fig. 16 的分离位置、负剪切区和远场基线差异多大；
* IBM diffuse forcing、近壁 slip 与 turbulence treatment 对这些差异分别有
  多大贡献。

后续若要形成正式 benchmark validation，再补做：

1. 用偏移平面离线补采 Fig. 10；
2. 数字化并叠加论文实验/模型曲线；
3. 校准入口边界层，使 `u_tau` 更接近 `1.30 cm/s`；
4. 统一统计平均窗口并报告误差指标。
