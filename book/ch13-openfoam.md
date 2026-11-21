# 第 13 章 用 OpenFOAM 再实现一次 DFN

> ★★ 选学 · 前置：第 3、4、7、9 章

> 第 9 章我们从零搭出了 `dfn-sim`：自己写有限体积、自己组装 Jacobian、把时间推进交给 IDA。这一章换一条路——不自己写矩阵，而是把同一组方程翻译成 OpenFOAM 的“场”和 `fvScalarMatrix`，让一个成熟的有限体积框架替我们做离散和线性求解。代码会短很多，但那些必须自己想清楚的东西——非线性怎么线性化、边界条件怎么进矩阵、Picard 迭代为什么会不收敛——一个都少不了。

## 本章你会得到

- 一张“第 3 章方程 → OpenFOAM 代码”的逐行映射表，知道每个 `fvm::` 项对应哪一项物理。
- 一个可运行的求解器 `dfnFoam`：4 个 OpenFOAM 场，加上每个 x 单元一列的颗粒径向浓度。
- 一份与 `dfn-sim` 的 1C 对照：全程电压差小于 2 mV，3.5 V 截止时间相差 0.2%。
- 三个最容易踩的坑：旧时间层被重复推进、矩阵欠松弛用错符号、液相电位整体平移不定。

## 13.1 为什么可以把 DFN “搬”进 OpenFOAM

OpenFOAM 处理的是三维非结构网格上的守恒方程。DFN 是一维的，但一维只是三维的一个特例：把 y、z 方向各划一格，x 方向按负极/隔膜/正极分成三段，就得到一个 20 个单元的小网格。面积因子由网格自己算，方程写法和三维完全一样。

OpenFOAM 把“离散 + 解线性方程”封装成三个对象：

- `volScalarField`：网格中心上的未知量场，例如电解液浓度 `c_e`；
- `surfaceScalarField`：面上的系数，例如有效扩散系数 `dEff`；
- `fvScalarMatrix`：一张已经离散好的线性方程，用 `fvm::laplacian`、`fvm::ddt`、`fvm::Sp` 这些算子拼出来，最后 `.solve()`。

第 3 章的方程可以几乎逐字翻译：

| 连续方程 | `dfnFoam` 中的写法 |
| --- | --- |
| $\varepsilon\partial_t c_e=\nabla\cdot(D_{\rm eff}\nabla c_e)+(1-t_+)a_s j$ | `fvm::ddt(eps, cE) - fvm::laplacian(dEff, cE) == srcCE` |
| $\nabla\cdot(\kappa_{\rm eff}\nabla\phi_e)+a_sFj=0$ | `fvm::laplacian(kappaEff, phiE) - fvm::Sp(coefE, phiE) == srcPhiE` |
| $\nabla\cdot(\sigma_{\rm eff}\nabla\phi_s)=a_sFj$ | `fvm::laplacian(sigmaPosf, phiSPos) - fvm::Sp(coefPos, phiSPos) == srcPhiPos` |
| 球形颗粒扩散 | 自定义径向有限体积 + 三对角求解 |

固相电位有两个场：`phi_s_neg` 和 `phi_s_pos`。它们各自只在自己的电极里有物理意义，在隔膜和另一个电极里只是占位未知量。`dfnFoam` 的做法是给每个场加一个掩码 `mask`：电极内为 0，电极外为 1，方程里多一项 `- fvm::Sp(mask, phiS)`，把无关区域的电位钉在 0；同时把电极/隔膜界面上的面电导率置零，防止电子通路穿过隔膜。

颗粒内部的锂浓度不是 OpenFOAM 场，而是每个 x 单元一列的数组 `cNeg[ix*nR + ir]`。每个单元独立解一个球对称扩散方程，用有限体积离散成三对角矩阵，再用 Thomas 算法回代。它的边界条件是 $r=R$ 处 $-D_s\partial c_s/\partial r=j$，离散后在 rhs 的最后一行减去面通量。

## 13.2 非线性从哪里来：把 $j$ 线性化

四个方程里，固相和液相的电荷守恒只对电位是线性的；非线性全部藏在 Butler–Volmer 反应

$$
j=j_0(c_e,c_s,\eta)\left[\exp\left(\frac{\alpha_aF\eta}{RT}\right)
-\exp\left(-\frac{\alpha_cF\eta}{RT}\right)\right],
\qquad \eta=\phi_s-\phi_e-U(c_s)
\tag{13.1}
$$

里。第 6 章讲过：Newton 法需要 $\partial j/\partial\eta$。这里不组装全局 Newton，而是在 Picard 迭代中把 $j$ 对**本方程的主角**做隐式线性化。

对固相电位，把 $j$ 在上一轮的 $\phi_s^{*}$ 附近展开：

$$
j\approx j^{*}+\frac{\partial j}{\partial\eta}(\phi_s-\phi_s^{*}),
\qquad \phi_s^{*}=\phi_e+U+\eta^{*}
\tag{13.2}
$$

代入 $\nabla\cdot(\sigma\nabla\phi_s)=a_sFj$，得到

$$
\nabla\cdot(\sigma\nabla\phi_s)-a_sF\frac{\partial j}{\partial\eta}\phi_s
=a_sF\left(j^{*}-\frac{\partial j}{\partial\eta}\phi_s^{*}\right)
\tag{13.3}
$$

这正是代码里的 `- fvm::Sp(coefPos, phiSPos)` 和 `srcPhiPos`。对液相电位则要注意符号：$\eta$ 里 $\phi_e$ 前面是负号，所以

$$
\nabla\cdot(\kappa\nabla\phi_e)-a_sF\frac{\partial j}{\partial\eta}\phi_e
=-a_sF\left(j^{*}+\frac{\partial j}{\partial\eta}\phi_e^{*}\right),
\qquad \phi_e^{*}=\phi_s-U-\eta^{*}
\tag{13.4}
$$

式 (13.3) 和 (13.4) 的右端括号里一个是减号、一个是加号。第一次实现时把 (13.4) 也写成减号，平衡态会被驱动偏离，几个 Picard 迭代后数值爆炸——这是本章最常见的符号错误。

还有一处细节：表面浓度不是独立未知量，而是由 $c_{\rm surf}=c_{\rm last}-\alpha j$ 外推出来的。所以 $j$ 对 $\phi_s$ 的“有效”导数要再把这条关系代进去：

$$
\left.\frac{\partial j}{\partial\eta}\right|_{c_{\rm last}}
=\frac{\partial j/\partial\eta}{1+\alpha\,\partial j/\partial c_s}
\tag{13.5}
$$

代码里的 `dj_dEtaEff` 就是它。少了这个约化，1C 下反应项会刚硬大约四个数量级，线性求解器解出来的电位会剧烈振荡。

## 13.3 恒流边界：电压是一个要解出来的标量

恒流放电时，正极集流体上的总电流是已知的，端电压是未知的。`dfnFoam` 用 Dirichlet 边界 `fixedValue` 施加正极集流体电位，然后反解这个电位。

第一步是“均匀电位假设”：假设正极内 $\phi_s$ 处处等于集流体电位 $V$，对每个单元解局部界面方程，累加得到总电流 $I(V)$。$I(V)$ 随 $V$ 单调变化，用二分法就能找到使 $I(V)=I_{\rm target}$ 的初值。

但均匀假设只在薄电极、高电导时成立。真正解出 $\phi_s$ 剖面后，实际电流会与目标差几个百分点。于是每个 Picard 迭代末尾再做一次标量 Newton：

$$
\delta V=\frac{I_{\rm target}-I_{\rm actual}}
{\sum_{\rm pos} a_sF\,\dfrac{\partial j}{\partial\eta}\Big|_{\rm eff}\Delta V}
\tag{13.6}
$$

其中 $I_{\rm actual}$ 是用当前解出的 $\phi_s$、$\phi_e$、$c_{\rm last}$ 重新算出的正极反应电流。加上这一步后，1C 下正极积分电流与 36 A/m² 的偏差从约 5% 降到舍入误差量级。

符号约定要一次说清：本章沿用第 3 章，$j>0$ 表示锂**离开**固相。放电时负极脱锂（$j>0$）、正极嵌锂（$j<0$），所以正极的积分电流 $\sum a_sFj$ 是负的。代码里用 `current -= ...` 把它翻成正的外加电流；如果忘了这个负号，二分法会收敛到开路电位以下的“充电”根，电压会在放电过程中不降反升。

## 13.4 液相电位的整体平移

$\phi_e$ 的方程

$$
\nabla\cdot(\kappa_{\rm eff}\nabla\phi_e)+a_sFj=0
\tag{13.7}
$$

在两端都是零通量边界，只有梯度被确定，整体加一个常数仍是解。这个常数不是任意的：它要正好让正负极的总反应电流互相抵消，

$$
I_{\rm neg}+I_{\rm pos}=0
\tag{13.8}
$$

因为电解液里不能凭空多出或减少锂。数值上，这个条件对应一个标量方程。`dfnFoam` 在每次 Picard 迭代里对 $\phi_e$ 的整体平移做内层 Newton：

$$
\phi_e\leftarrow\phi_e-\frac{I_{\rm neg}+I_{\rm pos}}
{\partial(I_{\rm neg}+I_{\rm pos})/\partial\phi_e},
\qquad
\frac{\partial(I_{\rm neg}+I_{\rm pos})}{\partial\phi_e}
=-\sum a_sF\frac{\partial j}{\partial\eta}\Big|_{\rm eff}\Delta V
\tag{13.9}
$$

（求和覆盖两个电极，正极的导数带负号，因为 $\eta$ 里 $\phi_e$ 前是负号。）迭代到相对残差 $10^{-10}$ 以下。没有这一步，Picard 会以每步约 0.85 的因子慢慢爬向平衡，一个时间步要几十次迭代；加上之后通常 2～6 次就收敛。

## 13.5 时间推进的三个坑

### 坑一：旧时间层被重复推进

时间项用后向欧拉。直觉写法是每个 Picard 迭代都解一次

```
fvm::ddt(eps, cE) - fvm::laplacian(dEff, cE) == src
```

但 OpenFOAM 的 `fvMatrix::solve()` 在末尾会调用 `correctBoundaryConditions()`，而后者会把当前场存成“旧时间层”。第二次迭代再装配 `fvm::ddt` 时，旧时间层已经变成了第一轮的解，于是同一个 `dt` 又被推进了一次。Picard 迭代 20 次，等于把时间推进了 20 个 `dt`，正极表面浓度会被迅速抽干。

解决办法是在每次装配前显式恢复步初场：

```cpp
cE.oldTime().primitiveFieldRef() = cEOld;
syncBoundary(cE.oldTime());
```

颗粒径向方程更容易犯这个错，因为它没有 `oldTime` 机制。`solveRadial` 必须同时接收旧场和新场，每次都从旧场重建：

```cpp
solveRadial(negGrid, cNegOld, cNeg, jNeg, dt);
```

### 坑二：矩阵欠松弛假定对角为正

三个电位方程写成 `+ fvm::laplacian(...)`，OpenFOAM 的矩阵对角是**负**的。`fvMatrix::relax(alpha)` 的实现会先取对角绝对值、强行保证对角占优，再加松弛源项；对负对角方程，这套操作会把源项放大成非物理量。第一版用 `phiSPosEqn.relax(0.5)`，结果正极电位每步减半、最后变成 0。

正确做法是把欠松弛作用在场上，而不是矩阵上：

```cpp
phiSPos.storePrevIter();
phiSPosEqn.solve();
phiSPos.correctBoundaryConditions();
phiSPos.relax(params.relaxPhi);   // phi = prev + alpha*(phi - prev)
```

`c_e` 的方程是标准的 `ddt - laplacian`，对角为正，矩阵松弛本身没有问题；但见坑三。

### 坑三：欠松弛会破坏守恒

`c_e` 方程对给定的 $j$ 是线性的，从旧时间层解一次就得到精确的后向欧拉解。如果用 `relax_ce=0.5` 做矩阵欠松弛，而 Picard 又只迭代有限次，最终场并不是精确的后向欧拉解，固相和液相各自差一点，全锂量就会缓慢漂移。实测 600 s 后相对漂移约 $2.4\times10^{-4}$；把 `relax_ce` 设为 1 后立刻降到 $10^{-13}$。

| `relax_ce` | 600 s 后全锂相对漂移 |
| ---: | ---: |
| 0.5 | $2.4\times10^{-4}$ |
| 1.0 | $8.1\times10^{-14}$ |

## 13.6 跑起来，和 `dfn-sim` 对照

```bash
cd openfoam/cases/demo-1c
./run.sh
```

`run.sh` 会自动构建求解器（如果还没构建）、清理旧时间目录、运行 `blockMesh`，再运行 `dfnFoam`。1C 放电到 3.5 V 大约十几秒，结果在

```
postProcessing/dfnFoam/terminalVoltage.csv
```

参考解由第 9 章的 `dfn-sim` 在相同网格（8/4/8，径向 6/6）和相同参数下给出：

| t (s) | dfnFoam V (V) | dfn-sim V (V) | 电压差 (mV) |
| ---: | ---: | ---: | ---: |
| 600 | 3.93953 | 3.93885 | +0.67 |
| 1200 | 3.91254 | 3.91186 | +0.68 |
| 1800 | 3.87903 | 3.87834 | +0.69 |
| 2400 | 3.82749 | 3.82675 | +0.73 |
| 3000 | 3.72348 | 3.72261 | +0.87 |
| 3.5 V 截止 | 3536 s | 3533.7 s | — |

化学计量数在 3000 s 时为 $0.2149/0.6712$（参考 $0.2147/0.6713$）。电压差随放电加深缓慢增大，主要来自一阶后向欧拉与 IDA 可变阶格式的差别；截止时间差 0.07%。全锂相对漂移约 $4\times10^{-13}$，与参考解同一量级。

## 13.7 调试清单

| 现象 | 最可能的原因 | 处理 |
| --- | --- | --- |
| 表面浓度几步内掉到 0 | 径向方程每个 Picard 都用更新后的场再推进 `dt` | 从步初场重建 |
| 电位出现棋盘振荡并指数增长 | 电位方程用了矩阵 `relax()`，或 $\phi_e$ 线性化符号写反 | 改场松弛；核对式 (13.4) 的 `+ dEta` |
| 放电时电压反而上升 | 恒流反解的电流符号反了 | 放电时正极 $j<0$，`current -= ...` |
| 端电压卡在 6 V 或 0 | 二分法括号/比较方向反了 | 先打印 $I(V)$ 的几个点核对单调性 |
| 电压、浓度都“收敛”但全锂漂移 | `relax_ce < 1` 或正负极电流未平衡 | 设 `relax_ce=1`；内层 Newton 平衡电流 |
| Picard 每步都要几十次 | $\phi_e$ 的整体平移模没有单独求解 | 用式 (13.9) 做标量 Newton |
| `phi_s` 在隔膜里出现大数 | 掩码符号或面电导率没有在界面置零 | 检查 `mask` 与 `sigmaPosf` |

## 13.8 练习

1. **换成电流边界。** 把正极集流体的 `fixedValue` 改为 `fixedGradient`，梯度取 $-J/\sigma_{\rm eff}$，直接用 Neumann 条件施加上电流，端电压由 $\phi_s$ 外推半步得到。与现在的电压反解比较收敛速度。
2. **自适应时间步。** Picard 达到上限仍未收敛时，把 `dt` 减半并从步初场重试；收敛快时再放大。注意保存所有场和径向浓度的步初副本。
3. **网格收敛。** 把 `blockMeshDict` 的三段单元数加密一倍，观察 1C 电压曲线与截止时间的变化，估计空间离散误差。
4. **加 CV 阶段。** 电压到 3.5 V 后把边界切换成 `fixedValue` 3.5 V，电流自由，观察电流衰减。
5. **换参数。** 修改 `constant/dfnProperties` 里的颗粒半径或固相扩散系数，比较电压平台与截止时间，和 `configs/demo.yaml` + `dfn-sim` 对照。

## 本章文件地图

| 文件 | 作用 |
| --- | --- |
| `openfoam/dfnFoam/dfnFoam.C` | 求解器：场、方程、径向求解、恒流控制、输出 |
| `openfoam/dfnFoam/build.sh` | Debian 运行时环境下的直接编译脚本 |
| `openfoam/dfnFoam/Make/` | 标准 OpenFOAM 环境下的 `wmake` 文件 |
| `openfoam/cases/demo-1c/` | 1C 算例：网格、物性、初始场、控制参数 |
| `openfoam/cases/demo-1c/run.sh` | 一键构建 + `blockMesh` + 求解 |
| `openfoam/README.md` | 构建、运行、验证与限制的完整说明 |

到这里，同一组 DFN 方程已经有了两种完全不同的实现：一种把每个矩阵元素都握在自己手里，另一种把离散交给框架、自己只负责物理和耦合。两条路得到的电压曲线在 2 mV 以内一致——这比任何单条曲线都更能说明：方程和离散是对的，差别只在实现细节。
