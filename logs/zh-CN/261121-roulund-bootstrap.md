# Roulund Tests 1 & 2：IBM 基线验证状态与后续动态试验准入

日期：2026-11-21

## 目的

使用 Roulund et al. (2005) 的刚性床 Tests 1 & 2，检查 `ichiBoundaryFoam` 中 VoF–IBM 耦合在固定圆柱、稳态来流条件下能否复现圆柱附近的三维近床流动，并判断该方案是否可以进入 prescribed force oscillation / vibration 阶段。

参考文献：`bibs/Roulund et al. - 2005 - Numerical and experimental investigation of flow a.pdf`

比较对象：

- `T1-BF`：光滑床、boundary-fitted 圆柱；
- `T1-IB`：光滑床、immersed-boundary 圆柱；
- `T2-BF`：粗糙床、boundary-fitted 圆柱；
- `T2-IB`：粗糙床、immersed-boundary 圆柱。

结果图位于两个 `~/fastscratch/waves2-3d-ibm-*` workspace 下的：

- `figs-roulund-t1-bf/roulund-t1/`；
- `figs-roulund-t1-ib/roulund-t1/`；
- `figs-roulund-t2-bf/roulund-t2/`；
- `figs-roulund-t2-ib/roulund-t2/`。

## 论文规定的试验条件

Tests 1 & 2 都是 steady-current rigid-bed experiments，而不是振荡流试验。因此论文没有为它们规定 KC 数；在不存在振荡周期 `T` 时，`KC = U_m T / D` 不是这两个试验的适用控制参数。

论文明确给出的主要参数为：

| 参数 | Tests 1 & 2 |
| --- | ---: |
| 水深 `h` | 0.54 m |
| 边界层厚度 `delta` | 0.54 m |
| 圆柱直径 `D` | 0.536 m |
| 积分平均速度 `V` | 0.326 m/s |
| `Re_D = VD/nu` | `1.7e5` |
| `Re_delta = V delta/nu` | `1.8e5` |
| `Fr = V/sqrt(gh)` | 0.14 |

Test 1 为光滑床。对数律拟合得到 `u_tau = 0.013 m/s`，hot-film 测量约为 `0.015 m/s`。论文给出的来流拟合为

```text
u/u_tau = 2.5 ln(z+) + 4.9
```

Test 2 为粗糙床。物理粗糙高度约为 0.007 m，等效 Nikuradse roughness 为 `k_s = 0.010 m`，`u_tau = 0.023 m/s`，`k_s+` 约为 230。论文给出的来流拟合为

```text
u/u_tau = 2.5 ln(30 z/k_s)
```

圆柱表面在两个试验中均为 hydraulically smooth。实验圆柱位于水槽入口下游约 20 m，即约 `37.3D`。

## 当前算例的参数符合程度

当前四个算例均使用 `h = 0.54 m`、`D = 0.536 m`、名义入口速度 `U = 0.326 m/s` 和水的 `nu = 1e-6 m2/s`。由此得到：

```text
Re_D = 0.326 * 0.536 / 1e-6 = 1.747e5
Re_delta = 0.326 * 0.54 / 1e-6 = 1.760e5
Fr = 0.326 / sqrt(9.81 * 0.54) = 0.142
```

这些数值与论文的舍入值一致。T2-BF 和 T2-IB 的床面使用 `k_s = 0.010 m` 的 rough-wall treatment，圆柱保持光滑，亦符合论文设定。

但应区分 nominal parameter matching 与 actual inflow matching：

- 入口 `U = 0.326 m/s` 并不自动证明测量断面的积分平均速度仍为 0.326 m/s；
- 水深为 0.54 m 并不自动证明实际边界层厚度 `delta` 已达到 0.54 m；
- T2 设置了正确的 `k_s`，但实际 `k_s+` 仍取决于计算得到的 `u_tau`；
- 当前入口 `U`、`k`、`omega` 仍属于 provisional uniform placeholders；
- 当前入口至圆柱中心约为 `10D`，而 `x/D = -6` 的来流采样面距入口仅约 `4D`，明显短于实验约 `37D` 的发展距离。

所以当前算例已经满足名义 Reynolds/Froude similarity，但尚未严格证明来流边界层、`u_tau`、`delta` 和 `k_s+` 与实验一致。

## 图像比较结论

### Fig. 6：无扰动来流剖面

T1/T2 的数值剖面均已形成合理的边界层形状，T2 也表现出比 T1 更强的近床速度亏损，趋势正确。不过当前图只给出 `Ux/U0` 数值曲线，尚未叠加论文实验点及对数律，也没有报告拟合得到的 `u_tau`、`delta` 和积分平均速度。因此 Fig. 6 尚未完成定量闭合。

### Figs. 10 与 13：中心平面流动

IBM 已复现圆柱前方流线下弯、迎流减速、近床分离、马蹄涡以及圆柱后方尾流等主要流动拓扑。T2 中更强的粗糙床近床旋转也能辨认。因此，固定圆柱 IBM 的主要流动物理在定性意义上得到支持。

需要注意：当前 BF 的 Fig. 10/13 中心平面图存在大面积采样缺失，只剩圆柱前后窄条，较可能属于并行 surface sampling、merge 或 mask 后处理问题。这两张 BF 图在修复前不能用于判断 BF 与 IBM 的物理优劣。

论文指出，Test 1 实验尾流存在 mean anticlockwise circulation，而论文自身的数值模型也未能正确捕获；Test 2 的 downstream clockwise circulation 则被捕获。当前 T1-IB 的尾流循环同样偏弱或不够明确，应在叠加实验矢量场后专门检查。

### Figs. 11、12、14 与 15：`Ux`、`Uz` 剖面

四个算例均捕捉到：

- 圆柱前方随接近圆柱而产生的减速；
- 近床局部反向流；
- 圆柱后方速度亏损和回流；
- 随高度增加而减弱的圆柱扰动。

IBM 与 BF 的外流场总体一致，但 IBM 的近壁反向速度、竖向速度和局部极值通常更平滑、更弱；差异主要集中在 `z = 0.005--0.02 m`、迎流分离区和尾流重新附着区。这个现象可能反映 IBM forcing、插值/extension 或 pressure correction 带来的局部耗散，后续动态载荷试验中应重点监测。

论文自身也报告了低高度、近圆柱区域的速度低估，因此这种偏弱趋势不是完全反常；但尚需通过实验数据叠加判断 IBM 是否比论文误差范围更差。

### Fig. 16：床面剪切

论文 Fig. 16 只给出了 smooth-bed Test 1 的 hot-film bed-shear validation。T1-BF/T1-IB 均复现了迎流侧剪切变化、马蹄涡附近符号变化及尾流负剪切区，但当前图使用 kinematic shear，尚未按论文的 bed-shear amplification 定义叠加实验数据。

论文原模型在 `-0.75 < x/D < -0.5` 区域低估床面剪切超过约 30%，该区间应作为后续判断 IBM 是否改善近壁载荷的重点区域。

T2 没有对应的 rough-bed hot-film shear validation；当前 T2 bed-shear 图只能作为诊断结果，不能称为对论文 Fig. 16 的复现。

## 当前验收结论

当前 IBM 方案可以认为已经**基本通过 Roulund Tests 1 & 2 的定性/工程级固定圆柱基线验证**，具体包括：

- 几何尺度、名义 `Re_D`、`Re_delta` 和 `Fr` 匹配；
- Test 2 的等效粗糙度参数匹配；
- 主要近床流动拓扑得到复现；
- IBM 与 BF 的外流场总体一致；
- MPI 修复后，两套 workspace 中已检查的 T2-BF 图逐字节一致，支持 plot-level reproducibility。

当前结果尚不足以宣称完成严格的 quantitative validation，原因包括：

- 来流剖面、`u_tau`、`delta`、实际 Reynolds 数和 T2 `k_s+` 尚未闭合；
- 论文实验曲线尚未直接叠加，尚无 RMSE、峰值误差或位置误差；
- IBM 近壁 `Uz`、分离区和局部极值存在偏弱迹象；
- Test 1 bed shear 尚未按论文定义比较；
- BF 中心平面图仍有后处理缺失。

## 对后续 force oscillation / vibration 试验的准入判断

当前方案可以进入 prescribed force oscillation / vibration 的方法开发和探索性试验，包括：

- moving/oscillating cylinder 的实现与 MPI 调试；
- IBM forcing、alpha extension 和 pressure correction 顺序的比较；
- 振幅、频率、KC 和 reduced velocity 的趋势扫描；
- 不同 IBM/BF 或算法组件之间的相对比较。

但 Roulund Tests 1 & 2 只验证了 fixed cylinder + steady current，并没有验证：

- moving-boundary kinematics；
- 圆柱速度和加速度项；
- 非定常水动力的幅值与相位；
- added mass、damping 和 vortex-induced force；
- VIV lock-in 区间；
- 长时间 fluid–structure coupling 稳定性。

因此，在完成以下第二阶段动态验证之前，forced-oscillation 结果应标为 exploratory，不应用于最终工程结论或正式的 VIV 定量预测。

## 第二阶段动态验证门槛

1. **补齐静止圆柱基线**
   - 在 `x/D = -6` 拟合 T1/T2 对数律；
   - 输出实际 `V`、`u_tau`、`delta`、`Re_D`、`Re_delta`、`Fr` 和 T2 `k_s+`；
   - 将论文 Figs. 6、11、12、14、15、16 的实验数据数字化并叠加；
   - 分别报告 upstream、horseshoe-vortex region 和 downstream wake 的误差。

2. **静水 prescribed oscillation**
   - 规定正弦位移；
   - 核对圆柱位置、速度、加速度和 IBM forcing 的幅值与相位；
   - 检查完整周期后的周期重复性和无物理漂移。

3. **均匀流 forced oscillation benchmark**
   - 选择具有公开实验或高可信数值数据的振荡圆柱 benchmark；
   - 至少测试三组 frequency/KC 或 reduced velocity；
   - 比较 `CD(t)`、`CL(t)` 的均值、谐波幅值及相位。

4. **时间步与空间收敛**
   - 至少两级 `Delta t`；
   - 至少两级 IBM 附近空间分辨率；
   - 对位移、载荷幅值和载荷相位分别报告变化。

5. **力计算交叉检查**
   - IBM forcing 积分力；
   - 控制体动量平衡力；
   - 条件允许时，与压力/剪切表面积分进行交叉验证。

## 推荐的对外表述

在动态验证完成前，建议采用如下措辞：

> The IBM simulations reproduce the principal flow topology of Roulund Tests 1 and 2 under nominally matched Reynolds and Froude numbers. The results provide a fixed-cylinder baseline for subsequent prescribed-oscillation studies. Quantitative validation remains conditional on matching the measured approach-flow boundary layers and validating the unsteady force amplitude and phase for a moving cylinder.

不宜写成“Roulund Tests 1 & 2 已经完全验证了 moving-body IBM 或 VIV prediction”。
