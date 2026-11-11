# 圆柱 IBM 初步验证记录

记录日期：2026-11-04。

## 验证目的与边界

本轮测试首先排查低级实现问题，而不是建立网格无关的正式
benchmark。当前程序只支持预设刚体运动对流体的单向作用；它尚未由
流体载荷反求结构运动，因此这里准确的称呼是 prescribed one-way WSI，
而不是双向 FSI。

基础工况为全水域、零重力、均匀来流圆柱绕流：

* 圆柱半径为 0.06 m；
* 来流速度为 0.012984 m/s；
* 水的运动黏度为 3.8952e-6 m2/s，对应 Re 约为 200；
* 仿真时长为 0.5 s；
* 贴体网格包含 34,764 个 cells；
* IBM 背景网格包含 20,760 个 cells。

因为两种离散使用了不同网格，本轮只能严格比较原生 `interFoam` 与
关闭 IBM 后的 `ichiBoundaryFoam`；贴体网格与 IBM 的差异只能作为早期
诊断，不能解释为 IBM 的最终精度。

## 贴体网格 solver 等价性

原生 `interFoam` 与关闭 IBM 的 `ichiBoundaryFoam` 使用完全相同的
贴体网格、初边值条件和数值设置。二者均正常运行至 0.5 s。末时刻全场
差异为：

| 场 | 最大绝对差 | RMS 差 |
|---|---:|---:|
| U | 3.53e-9 m/s | 1.40e-10 m/s |
| alpha.water | 1.00e-11 | 2.36e-13 |
| p_rgh | 5.09e-9 m2/s2 | 4.03e-10 m2/s2 |

圆柱末时刻总阻力分别为 3.28920598815e-4 N 和
3.28920597149e-4 N，相对差约 5.1e-9。由此可认为：在本工况下，关闭
IBM 后 `ichiBoundaryFoam` 已基本复现原生 `interFoam` 的流体路径。

## 静止 IBM 圆柱

静止 IBM 圆柱正常运行至 0.5 s：

* alpha.water 始终位于约 [0.999999999994, 1.00000000001]；
* 最大 Courant No. 为 0.0143；
* 启动阶段最大 slip 为 7.08e-3 m/s，末时刻降至 1.30e-4 m/s；
* 最大散度为 2.72e-9 1/s，末时刻为 2.83e-11 1/s；
* 圆柱上下游出现了预期的阻塞、侧向加速和尾流速度亏损。

末时刻 IBM forcing 的体积分给出约 5.65e-4 N 的反作用力估计，比贴体
网格结果高约 72%。该数值暂时不能作为正式载荷结果：当前保存的
`immersedAcceleration` 是末次 predictor correction，并不是经过严格
时间与压力/黏性分解的结构载荷；同时 IBM 与贴体网格的分辨率也不同。

## 预设刚体运动

为便于检查 SDF 与壁面速度是否随时间更新，新增了
`constant6DoF` 运动类型。

匀角速度算例令圆柱绕自身 z 轴以 0.20 rad/s 转动。0.5 s 时记录的
转角为 0.10 rad。它实际表示二维圆柱自转，不应称为倾倒。末时刻 slip
为 1.66e-4 m/s，散度为 9.88e-12 1/s，计算稳定完成。

匀速平移算例令圆柱沿 -y 方向以 0.002 m/s 移动。0.5 s 时圆心由
(0, 0, 0) 移至 (0, -0.001, 0)，SDF 与 forcing-cell 分类同步更新。
末时刻 slip 为 1.35e-4 m/s，散度为 1.57e-11 1/s，alpha.water 未出现
有意义的越界。该结果说明预设平移的基本数据通路可以工作。

## Stokes II + IBM 烟雾测试

在已有 Stokes II 回归算例中加入静止 IBM 圆柱，并运行半个波周期
（1.5 s）。计算没有崩溃，最大 Courant No. 为 0.234；水相下界约为
-3.53e-19，上界达到 1.00065。

启动阶段记录的最大 slip 为 0.0365 m/s，末时刻仍为 0.0177 m/s，说明
波浪工况下的速度约束仍明显弱于温和单相算例。全域最大散度在启动时
达到 16.3 1/s，末时刻为 1.30 1/s；但随后进行的空间定位表明，末时刻
峰值位于入口附近的空气区 `(x,z)=(0.06,0.425) m`，并不位于 IBM。
同一时刻固体内部、IBM transition layer 和自由液面 cells 中的最大
绝对散度分别约为 2.89e-2、2.88e-3 和 5.77e-2 1/s。此前无 IBM 的
Stokes II 回归算例也在入口空气区出现约 0.825 1/s 的全域峰值。因此
`IBM projection residuals` 当前打印的全域最大散度容易被 wave inlet
主导，不能单独据此判定 IBM 压力耦合失效。

这一定位缓和了“IBM 产生 1.30 1/s 局部散度”的判断，但尚不能视为
自由液面 WSI 已通过：波浪算例的 slip 仍偏大，且不同区域的散度应与
完全相同时间和网格下的无 IBM 基线逐项比较。

## 本轮发现的低级问题

* 自编译 executable 位于项目的 `platforms/.../bin/`，并未自动加入
  Slurm OpenFOAM module 的 PATH。验证脚本必须显式调用该路径，或在
  安装阶段设置 FOAM_USER_APPBIN；否则可能只生成日志而没有真正运行。
* 原先的旋转工况是绕二维圆柱自身轴线自转，不是三维 monopile 的
  倾倒。真正的倾倒必须使用三维有限长度结构与绕水平轴的转动。
* `immersedAcceleration` 的简单体积分只适合调试 forcing 的量级，
  当前还不能替代正式的 IBM hydrodynamic force reconstruction。
* 单相温和绕流通过不代表两相 WSI 已通过；当前全域散度日志会被入口
  主导，后续诊断必须按 fluid/IB/free-surface 区域分别统计。

## 下一步建议

1. 将分区散度诊断正式加入 solver 日志，避免把入口峰值误判为 IBM
   transition-layer 峰值。
2. 对 Stokes II 工况分别运行：无 IBM、远离自由液面的 IBM、减半网格
   尺度和减半时间步，以分离波浪基线误差、IBM 误差与离散误差。
3. 实现守恒且时间一致的 IBM 载荷积分，再与相同局部分辨率的贴体网格
   比较阻力、升力及相位。
4. 在此之前，不把当前 IBM 结果称为正式 benchmark，也不开展双向
   6DoF FSI。

## Multi-direct forcing 追加实验

2026-11-04 增加了可选的 multi-direct-forcing pressure projection。
它不改变 SDF 几何表示，也不在 immersed surface 上设置压力边界。
每轮 correction 执行：

```text
SDF-weighted velocity constraint
→ global continuous pressure-increment solve
→ consistent U/phi/p_rgh update
```

与原单次 forcing 路线相比，末时刻结果如下：

| 工况 | 原 slip | MDF slip | 原散度 | MDF 散度 |
|---|---:|---:|---:|---:|
| 静止、全水域 | 1.30e-4 | 3.68e-5 | 2.83e-11 | 5.27e-12 |
| 匀速平移、全水域 | 1.35e-4 | 8.48e-5 | 1.57e-11 | 1.16e-11 |
| 静止、Stokes II | 1.77e-2 | 2.91e-3 | 1.30 | 1.42e-9 |

Stokes II 工况的最大 alpha.water 同时由 1.00065 降至 1.0。静止
全水域工况使用 4 轮额外 correction；平移和 Stokes II 使用 1 轮。
静止工况中第一轮已经取得绝大部分 slip 改善，后续三轮的边际收益很
小，因此当前建议从 1 轮开始，而不是默认使用 4 轮。

严格的压力增量收敛是必要的。开发中的第一版误用了普通 p_rgh 的相对
容差，虽然降低 slip，却把散度放大到 3.18e-4 1/s；改用
`p_rghFinal` 的零相对容差后，slip 和散度才同时改善。额外投影会增加
计算成本，本轮记录的 CPU 时间约增至原来的 2.3–6.7 倍，
后续仍需优化线性求解器和停止条件。
