# `ichiBoundaryFoam` 的 IBM 开发 baseline

## 目标

先把两类误差源分开：

* VoF 与波浪传播模块本身；
* 启用 immersed boundary 后新增的几何、forcing 与通量耦合。

因此，首个正式 water-tank baseline 不启用 IBM forcing，但仍使用
包含 IBM 代码的常规 `ichiBoundaryFoam` binary。

## Baseline BF-0：贴体网格对照

* Solver：常规、支持 IBM 的 `ichiBoundaryFoam`；
* Runtime：`ibm.enabled false`；
* Geometry：原有 body-fitted `blockMesh` monopile wall；
* Wave：Stokes II，周期 `T = 2 s`；
* Turbulence：`kOmegaSST`；
* Time step：固定 `deltaT = 0.001 s`；
* Momentum convection：

```text
div(rhoPhi,U) Gauss limitedLinearV 1;
```

这个算例首先回答：在完全没有 IBM forcing 的情况下，开发版 solver
能否复现原有水槽的波浪传播、桩体受力和长期稳定性。

## Baseline IBM-0：静止 immersed cylinder

只有 BF-0 通过后，才在相同波浪与湍流配置下建立背景网格版本：

* static SDF cylinder：中心 `(10 0 0.6) m`；
* radius：`0.045 m`；
* length：`1.2 m`；
* discrete forcing；
* `forcingMethod diffuseVelocity`；
* `multiDirectForcing true`；
* 压力仍采用连续压力场，不把 immersed boundary 当作压力方程中的
  内部边界；
* `extendAlphaIntoSolid false`；
* `sharpReconstruction` 与 alpha extension 暂列为实验选项，不作为
  第一个 baseline。

该组合不是最终方法的声明，而是当前代码中测试历史最完整、最适合
继续定位问题的起点。

## 验收顺序

1. 网格、初始场和 solver 启动无误；
2. Courant number 与 `alpha.water` 上下界不随时间恶化；
3. 总水量、测波点水位和压力不存在长期非物理漂移；
4. `F_x`、`F_y`、`F_z` 及桩周流场可重复分析；
5. BF-0 通过后，再比较 IBM-0 与 BF-0 的相位、载荷和局部流场；
6. 最后才比较 diffuse/sharp、MDF 次数和 alpha extension。

## 关于“修缮 VoF”

开发版可以用于审查和改进 VoF，但不能先假定官方 `interFoam` 存在
缺陷。任何改动都应在相同网格、固定时间步、相同离散格式和初始场下，
用守恒性、界面界限、相位误差与长期漂移给出可复现证据。纯
`interFoam` 路径已经完成固定时间步的短程回归；BF-0 则负责验证真实
Stokes II monopile water tank 中的表现。
