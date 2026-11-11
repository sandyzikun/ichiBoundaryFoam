# MDF 迭代优化路线初步评估

## 核心判断

Multi-direct forcing（MDF）在形式上确实类似于迭代优化，但更准确地说，它是
对 immersed-boundary 离散约束方程进行的 stationary/fixed-point iteration，
而不是含有 mini-batch noise 的机器学习训练。

因此，可以借鉴 ML/DL optimizer 的“历史信息、动量、自适应步长”等思想，
但落实时应优先采用适合 immersed-boundary operator 的数值线性代数方法：

* relaxation parameter；
* Aitken dynamic relaxation；
* Anderson acceleration；
* Krylov method；
* preconditioned implicit/direct forcing；
* residual-based stopping criterion。

Adam、AdaDelta 等可以保留为实验性对照，但不宜直接成为默认 MDF backend。

## MDF 对应的离散约束问题

记：

* `J` 为 Eulerian velocity 到 immersed markers 的 interpolation operator；
* `S` 为 marker force 到 Eulerian cells 的 spreading operator；
* `M^-1` 为离散动量算子的近似逆；
* `f` 为待求的 immersed-boundary force；
* `Ub` 为结构速度。

经过 forcing 后，immersed boundary 上的速度约束可以抽象为：

$$
A\boldsymbol f=\boldsymbol b,
\qquad
A=JM^{-1}S.
$$

boundary slip residual 可记为：

$$
\boldsymbol r_k
=
\boldsymbol U_b-J\boldsymbol U_k
=
\boldsymbol b-A\boldsymbol f_k.
$$

朴素 MDF 的主要形式相当于 Richardson iteration：

$$
\boldsymbol f_{k+1}
=
\boldsymbol f_k
+
\omega_k\boldsymbol r_k.
$$

当 `A` 具有适当的对称正定性质时，该过程也可以解释为对下列二次目标函数作
梯度下降：

$$
\mathcal L(\boldsymbol f)
=
\frac12\boldsymbol f^T A\boldsymbol f
-
\boldsymbol b^T\boldsymbol f.
$$

不过，这一解释依赖 interpolation/spreading、质量矩阵和离散内积之间的相容性。
若 `S` 并不是 `J` 在适当加权内积下的伴随，MDF 就不能简单当作普通 Euclidean
gradient descent，盲目套用 ML optimizer 可能破坏守恒性和稳定性。

## 文献中已有的 MDF 加速思路

### Relaxed MDF

已有研究通过估计 relaxation parameter 加快 MDF 收敛。在某些固定圆柱算例
中，合适的 relaxation parameter 可显著减少满足 no-slip 所需的迭代次数。

这与 optimizer 中的 learning rate 最相近，但该参数应由 IBM operator 的谱性质
或残差历史确定，而不是照搬机器学习中的经验学习率。

参考：[Jiang and Liu, relaxed MDF IB-CLBM (2020)](https://pureportal.strath.ac.uk/en/publications/a-relaxed-multi-direct-forcing-immersed-boundary-cascaded-lattice/)

### Preconditioned implicit direct forcing

Park 等人在 implicit direct forcing recurrence 中引入 preconditioner matrix，
以改善 forcing iteration 的收敛。这比 Adam 式逐 marker 缩放更符合当前问题的
线性代数结构：

$$
\boldsymbol f_{k+1}
=
\boldsymbol f_k
+
P^{-1}
\left(
\boldsymbol b-A\boldsymbol f_k
\right).
$$

理想的 `P^-1` 应廉价近似 `A^-1`，同时保留 marker/cell coupling 和离散守恒
性质。

参考：[Park et al. (2016)](https://doi.org/10.1016/j.jcp.2016.03.035)

### Residual-based stopping

固定 MDF 次数容易出现两种问题：

* 简单时间步中做了不必要的 forcing corrections；
* 困难时间步中固定次数仍不足以满足 boundary condition。

已有工作以相邻 forcing cycles 的 force difference 作为停止条件，例如：

$$
\left\|
\boldsymbol F_{k+1}-\boldsymbol F_k
\right\|
\leq
\varepsilon_F.
$$

Fang 与 Tan 报告 residual/stopping criterion 能在保持精度的同时减少固定次数
MDF 的冗余计算。

参考：

* [Fang and Tan (2023)](https://doi.org/10.2139/ssrn.4538797)
* [Improved MDF and Fourier pseudo-spectral method (2022)](https://doi.org/10.1007/s40430-022-03679-5)

### Operator-level improvement

有些误差并不是因为 MDF 迭代次数不足，而是 interpolation 与 spreading operator
本身不相容。Bao 等人构造了 divergence-free interpolation 和与之伴随的 force
spreading，从 operator 层面改善体积守恒和 Lagrangian force 光滑性。

这说明“增加更聪明的 optimizer”不能替代对 `J`、`S` 和 pressure projection
离散相容性的检查。

参考：[Bao et al. (2017)](https://doi.org/10.1016/j.jcp.2017.06.041)

## 对 Adam/AdaDelta 的评估

Adam 式 MDF 可以形式化为：

$$
\boldsymbol m_k
=
\beta_1\boldsymbol m_{k-1}
+
(1-\beta_1)\boldsymbol r_k,
$$

$$
\boldsymbol v_k
=
\beta_2\boldsymbol v_{k-1}
+
(1-\beta_2)\boldsymbol r_k^2,
$$

$$
\boldsymbol f_{k+1}
=
\boldsymbol f_k
+
\eta
\frac{\widehat{\boldsymbol m}_k}
{\sqrt{\widehat{\boldsymbol v}_k}+\epsilon}.
$$

它并非绝对不可用，但存在明显风险：

* MDF residual 基本是确定性的，Adam 对 stochastic gradient 的优势并不明显；
* marker-wise adaptive scaling 可能破坏旋转与镜像对称性；
* 不同 marker 使用不同有效步长，可能破坏 interpolation/spreading 的伴随关系；
* moving IBM 中，上一时刻或上一位置的历史量可能污染当前几何状态；
* `eta`、`epsilon`、`beta1`、`beta2` 引入额外尺度和调参问题；
* 即使 predictor-level slip 下降，也不能保证 pressure projection 后的 slip、
  divergence 和载荷同时改善。

因此，如果以后实现 `adamExperimental`，应默认关闭，并且：

* 每个新时间步默认清空历史；
* moving geometry 的 marker topology/position 显著变化时强制清空历史；
* 优先使用全局或分组尺度，而不是任意逐 marker 缩放；
* 同时检查 boundary work、总 constraint impulse 和镜像对称性；
* 仅把它作为与 Aitken/Anderson 的研究性对照。

AdaDelta、RMSProp 等也面临类似问题。它们的 adaptive scaling 思想可用于设计
preconditioner，但不宜不经修改地直接替换 MDF recurrence。

## 更适合当前 IBM 的 acceleration backends

### Fixed relaxation

$$
\boldsymbol f_{k+1}
=
\boldsymbol f_k
+
\omega\boldsymbol r_k.
$$

优点是实现简单、内存开销小，适合先建立基准。需要扫描或估计稳定的 `omega`
范围。

### Aitken dynamic relaxation

Aitken 根据连续两次 residual/correction 自动更新一个全局 relaxation parameter。
它保留 MDF 的基本结构，并避免 Adam 式逐 marker 自适应缩放，因此建议作为第一
个正式 acceleration backend。

实施时必须：

* 限制 `omega` 的上下界；
* 当 residual 方向突变或分母过小时回退到安全值；
* moving geometry 后重置历史；
* 分别记录 predictor-level 与 post-projection slip。

### Anderson acceleration

Anderson acceleration 保存最近若干次 force/residual 历史，并通过小型 least-
squares problem 组合 correction directions。它比 Aitken 使用更多历史信息，
适合 fixed-point form 的 MDF：

$$
\boldsymbol f_{k+1}=G(\boldsymbol f_k).
$$

建议采用较小 history depth，例如 `m=2--5`，并提供：

* regularization；
* ill-conditioned least-squares fallback；
* residual 增大时自动 restart；
* moving geometry 时 history reset。

### Krylov/GMRES

若可以显式或 matrix-free 地计算：

$$
A\boldsymbol p=JM^{-1}S\boldsymbol p,
$$

则可直接以 GMRES 等 Krylov method 求解 IBM force constraint。它在数学上通常
比 Adam 更适合非对称线性系统，但每个 matrix-vector product 都可能需要一次
spreading、momentum response 和 interpolation，因此必须评估实际成本。

### Preconditioned MDF

长期来看，这是最值得和 Krylov/Anderson 结合的路线。可从廉价的局部近似开始：

* diagonal marker self-coupling；
* block diagonal coupling；
* 基于 marker spacing、kernel overlap 和 local fluid mass 的近似；
* 对平移/旋转圆柱可复用的局部 operator information。

必须确保 preconditioner 不因 marker ordering、镜像网格或圆柱转动而产生非物理
偏置。

## 与压力投影优化的关系

MDF acceleration 主要解决：

$$
J\boldsymbol U
\rightarrow
\boldsymbol U_b.
$$

普通或 IBM-aware pressure projection 主要解决：

$$
\nabla\cdot\boldsymbol U
\rightarrow
0,
$$

并在 coupled formulation 中同时维持 boundary constraint。两者属于不同维度，
可以分别配置：

```text
forcing accelerator:
    plain
    fixedRelaxation
    aitken
    anderson
    krylov
    preconditioned
    adamExperimental

projection scheduling/backend:
    finalOnly
    everyIteration
    everyM
    residualTriggered
    ibmAware
```

更快的 MDF 不能从根本上防止普通 pressure projection 重新引入 slip；更多压力
投影也不能修复 interpolation/spreading operator 本身的 reciprocity error。

## 建议的软件接口

后续可在 `constant/ichiBoundaryGeometry` 中增加类似配置：

```text
forcingSolver
{
    accelerator         plain;
    maxIterations       4;
    slipTolerance       1e-5;
    forceTolerance      1e-6;
    relaxation          1.0;

    aitken
    {
        omegaMin        0.05;
        omegaMax        1.5;
        resetEachStep   true;
    }

    anderson
    {
        depth           3;
        regularisation  1e-10;
        resetEachStep   true;
    }
}
```

具体字段和默认值需要通过验证决定；上述内容仅用于规划接口，尚未写入 solver。

建议使用运行时选项选择 acceleration backend，而不是为每种方法生成不同 binary。
只有涉及额外依赖或显著内存结构的 backend，才考虑使用预编译开关控制是否编译。

## 推荐开发顺序

1. 把当前 MDF residual、force change 和 iteration count 记录完整；
2. 把固定 MDF 次数改成 tolerance 加 safety cap；
3. 加入 bounded fixed relaxation；
4. 实现 Aitken，并在静止单相圆柱上验证；
5. 实现 Anderson acceleration，和 Aitken 作等成本比较；
6. 研究 `JM^-1S` 的 matrix-free application 和局部 preconditioner；
7. 最后才加入默认关闭的 `adamExperimental`；
8. 在单相静止、单相运动、V1、R1 的顺序上逐级验证；
9. 再与 batched/adaptive pressure projection 和 IBM-aware projection 组合。

## 验证指标

不能只看 MDF 迭代末尾的 slip。每种 accelerator 都应记录：

* 每轮 predictor-level slip；
* final post-projection slip；
* immersed band 和全域 divergence；
* force residual 与 constraint impulse；
* IBM force/moment 的平滑性；
* MDF iteration count；
* pressure projection 次数与 GAMG iterations；
* wall-clock time；
* `alpha.water` 有界性与 phase-volume conservation；
* 对镜像网格、marker ordering 和 cylinder rotation 的敏感性。

## 当前结论

可以借鉴 ML optimizer 的思想，但实现应以 IBM constraint operator 为中心。近期
最合理的路线不是直接加入 Adam，而是：

```text
residual-based stopping
    -> bounded relaxation
    -> Aitken
    -> Anderson
    -> Krylov/preconditioning
```

这条路线既保留“利用历史信息加快迭代”的核心思想，也更容易维持 IBM 的守恒、
对称性、量纲一致性和可解释性。Adam/AdaDelta 可在上述基准建立后作为有趣的
研究性分支，而不应抢占当前 slip、projection compatibility 和 moving IBM 的
主线开发。
