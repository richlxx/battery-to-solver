# OpenFOAM 版 DFN 求解器

`openfoam/` 用 OpenFOAM 的有限体积框架重新实现本书第 3 章的 DFN（P2D）模型。
它与仓库里的 `dfn_core`（作者自研的 IDA/KLU 求解器）共享同一套参数与方程，
因此可以直接对照验证。

- `dfnFoam/dfnFoam.C`：求解器源码。四个 OpenFOAM 场（`c_e`、`phi_e`、
  `phi_s_neg`、`phi_s_pos`）加每个 x 单元一列的颗粒径向浓度数组。
- `dfnFoam/build.sh`：在 Debian 的运行时安装（只有 OpenFOAM 库和头文件源树、
  没有完整 `wmake` 环境）下直接编译。
- `dfnFoam/Make/`：标准 OpenFOAM 安装下用 `wmake` 构建所需文件。
- `cases/demo-1c/`：与 `configs/demo.yaml` 一致的 1C 恒流放电算例。

## 快速开始

```bash
cd openfoam/cases/demo-1c
./run.sh
```

`run.sh` 会在需要时调用 `build.sh`，删除旧的时间目录，运行 `blockMesh`，
再运行 `dfnFoam`。结果写入
`postProcessing/dfnFoam/terminalVoltage.csv`，列为

```
t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,li_drift,picard_iter
```

1C 放电到 3.5 V 截止约需 15 秒墙钟时间（3600 个 1 s 步，每步 5～45 次 Picard）。

与仓库根目录的 `dfn-sim`（需要先构建出 `build-ref/dfn-sim`）自动对照：

```bash
openfoam/validation/compare_1c.sh
```

脚本会重跑 1C 算例、运行参考求解器，并打印每 600 s 的电压差。也可以用
`DFN_SIM=/path/to/dfn-sim` 指定参考解。

## 方程与场的对应

| 连续方程 | OpenFOAM 实现 | 说明 |
| --- | --- | --- |
| $\varepsilon\partial_t c_e=\nabla\cdot(D_{\rm eff}\nabla c_e)+(1-t_+)a_s j$ | `fvm::ddt(eps,cE) - fvm::laplacian(dEff,cE)` | 后向欧拉；源项按电极分区写入 |
| $\nabla\cdot(\kappa_{\rm eff}\nabla\phi_e)+a_sFj=0$ | `fvm::laplacian(kappaEff,phiE) - fvm::Sp(coefE,phiE)` | $j$ 对 $\phi_e$ 隐式线性化 |
| $\nabla\cdot(\sigma_{\rm eff}\nabla\phi_s)=a_sFj$ | `fvm::laplacian(sigmaNegf/sigmaPosf,phiS)` | 隔膜与另一电极内用单位阵钉在 0 |
| $\partial_t c_s=\nabla\cdot(D_s\nabla c_s)$（球） | 自定义径向有限体积 + Thomas 三对角求解 | 每个 x 单元独立求解 |

反应采用与 `src/models/reaction.cpp` 相同的 Butler–Volmer 形式，并保留
$c_{\rm surf}=c_{\rm last}-\alpha j$ 的表面浓度修正。电位方程中的 $j$ 用
$\partial j/\partial\eta$ 隐式线性化；考虑表面浓度反馈后使用约化导数
$\partial j/\partial\eta\big|_{c_{\rm last}}
=(\partial j/\partial\eta)/(1+\alpha\,\partial j/\partial c_s)$。

## 几个实现要点

1. **一维网格。** `blockMeshDict` 在 x 方向划分负极/隔膜/正极，y、z 各一格，
   面积的单位厚度由网格自动给出，方程与三维有限体积完全一致。
2. **恒流控制。** 正极集流体用 `fixedValue` 电位。先由均匀电位假设反解一个初值，
   之后每个 Picard 迭代用“解出的实际正极反应电流”对集流体电压做标量 Newton
   修正，使总电流严格等于 `current_density`。
3. **液相电位参考点。** $\phi_e$ 的方程只有 Neumann 边界，整体平移量不唯一。
   它的值由正负极总反应电流为零这一条件确定；求解器对平移量做内层 Newton
   迭代，残差收敛到约 $10^{-10}$ 的相对量级。
4. **时间推进。** 后向欧拉 + Picard。颗粒径向方程和 $c_e$ 方程必须从**时间步
   开始时的旧场**重新求解；如果直接拿上一轮 Picard 更新过的场再推进一个
   `dt`，一个时间步内会重复推进多次，表面浓度会被抽干。$c_e$ 的 `oldTime`
   因此在每次装配前恢复；`relax_ce` 应为 1，欠松弛会让后向欧拉解不满足
   守恒。
5. **场的欠松弛。** 三个电位方程写成 `+fvm::laplacian` 形式，矩阵对角为负；
   OpenFOAM 的 `fvMatrix::relax()` 假定对角为正。因此这里不用矩阵级松弛，
   而是先 `storePrevIter()`、求解后对场调用 `relax(alpha)`。

## 验证

参考解由仓库根目录的 `dfn-sim`（IDA + KLU，与本书正文同一实现）在相同网格
（8/4/8，径向 6/6）与相同参数下给出。

### 0.1C（3.6 A/m²，100 s）

| t (s) | dfnFoam V (V) | dfn-sim V (V) | ΔV (mV) |
| ---: | ---: | ---: | ---: |
| 20 | 3.984342 | 3.984271 | +0.07 |
| 60 | 3.983992 | 3.983920 | +0.07 |
| 100 | 3.983712 | 3.983640 | +0.07 |

### 1C（36 A/m²，放电到 3.5 V）

| t (s) | dfnFoam V (V) | dfn-sim V (V) | ΔV (mV) | dfnFoam $\theta_-/\theta_+$ | dfn-sim $\theta_-/\theta_+$ |
| ---: | ---: | ---: | ---: | --- | --- |
| 600 | 3.93953 | 3.93885 | +0.67 | 0.68298/0.37423 | 0.68294/0.37425 |
| 1200 | 3.91254 | 3.91186 | +0.68 | 0.56595/0.44846 | 0.56589/0.44850 |
| 1800 | 3.87903 | 3.87834 | +0.69 | 0.44893/0.52269 | 0.44883/0.52275 |
| 2400 | 3.82749 | 3.82675 | +0.73 | 0.33191/0.59692 | 0.33178/0.59701 |
| 3000 | 3.72348 | 3.72261 | +0.87 | 0.21488/0.67116 | 0.21472/0.67126 |
| 截止 | 3536 s | 3533.7 s | — | 0.11034/0.73747 | 0.11060/0.73731 |

全锂量相对漂移约 $4\times10^{-13}$（参考解约 $4\times10^{-15}$）。电压偏差随放电
加深缓慢增大，来自一阶后向欧拉与参考解可变阶 IDA 的差异；截止时间相差 0.07%。

## 构建

### 标准 OpenFOAM 环境

已经 `source` 过 OpenFOAM 环境（`wmake` 可用）时：

```bash
cd openfoam/dfnFoam
wmake
```

### 本仓库的 Debian 运行时环境

Debian 的 `openfoam` 包只提供运行库，没有完整的开发环境。`build.sh` 会：

- 从 `$FOAM_SRC`、`$WM_PROJECT_DIR` 或 `~/openfoam-src/openfoam-OpenFOAM-*`
  找到 OpenFOAM 头文件源树；
- 用源树里的 `wmakeLnInclude` 生成 `lnInclude`；
- 用 `CXX`（默认 `g++`，本机为 conda 的 g++ 12）把求解器与
  `/usr/lib` 下的 `libOpenFOAM`、`libfiniteVolume`、`libmeshTools` 链接。

```bash
cd openfoam/dfnFoam
./build.sh
```

OpenFOAM 运行时需要 `WM_PROJECT_DIR=/usr/share/openfoam`（`run.sh` 已设置）。

## 输出与诊断

除了 CSV，算例目录下的时间目录（600、1200、……）保存 `c_e`、`phi_e`、
`phi_s_neg`、`phi_s_pos`、`c_surf_neg`、`c_surf_pos` 场，可直接用 ParaView
打开。把 `constant/dfnProperties` 里的 `debug` 改成 `true`，求解器会打印
每个 Picard 迭代的电位范围、正极电位剖面和电流平衡量，便于教学时观察。

## 当前限制

- 一维等温模型，不含热耦合、SEI、析锂和老化。
- 时间推进是一阶后向欧拉，时间步固定；参考求解器使用可变阶 IDA。
- 恒流通过集流体电位反解实现，适合本算例的薄电极；更厚的电极可以改用
  集流体电流（Neumann）边界。
- 目前只演示 CC 放电。CV 阶段、充电工况与从中间时间目录重启都还没有实现；
  每次运行都从 `0/` 的一致初值开始。
- 参数与 `configs/demo.yaml` 一致，是教学用的简化拟合，不是真实电芯标定值。
