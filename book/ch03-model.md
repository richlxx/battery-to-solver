# 第 3 章 案例模型：P2D/DFN

> ★ 核心 · 前置：第 0–2 章

> 把电池拆成三样东西：两颗颗粒、一条电解液通道、一个界面反应。本章不求解，只把问题说清楚。

## 本章你会得到

- 一张 DFN 方程总表：11 个编号公式，后续章节直接引用。
- 知道每个未知量由哪个方程决定、边界在哪里。
- 一套可加载的电池参数与状态布局。

## 3.1 电池是怎么工作的（三分钟版）

充电时锂从正极脱出、穿过电解液、嵌入负极；放电反过来。电子走外电路，离子走电解液。

## 3.2 五个方程与编号总表

DFN 的未知量是四个场：固相浓度 $c_s(x,r,t)$、液相浓度 $c_e(x,t)$、固相电位 $\phi_s(x,t)$、液相电位 $\phi_e(x,t)$，外加界面通量 $j(x,t)$。它们由五个物理过程决定。

本章把所有独立公式统一编号为 (3.1)–(3.11)；后续章节只引用编号，不重复推导。第一遍只需看懂下表，公式可以后面再回查。

| 物理过程 | 内容 | 编号 | 未知量 | 类型 |
| --- | --- | --- | --- | --- |
| 固相扩散 | 球坐标 Fick + 边界通量 | (3.1)(3.2) | $c_s$ | 抛物型 PDE |
| 液相输运 | 扩散-迁移 + 界面源项 | (3.3)(3.4) | $c_e$ | 抛物型 PDE |
| 界面反应 | Butler–Volmer + 交换电流 | (3.5)(3.6) | $j$ | 非线性代数 |
| 固相电荷 | 欧姆定律 + 电荷守恒 | (3.7) | $\phi_s$ | 椭圆型 PDE |
| 液相电荷 | 欧姆 + 浓溶液修正 | (3.8) | $\phi_e$ | 椭圆型 PDE |
| 约束 | 总电流、参考电位 | (3.9)(3.10) | — | 代数 |
| 初值 | 初始浓度 | (3.11) | — | 初值 |

**符号与约定**

- 下标 $s$/$e$ 表示固相/液相；上标 $\mathrm{eff}$ 表示多孔电极中的有效值。
- $j>0$ 表示锂离开固相进入电解液。放电时负极 $j>0$、正极 $j<0$，全书统一此约定。
- $c_{s,\mathrm{surf}}=c_s(R_p,t)$ 为颗粒表面浓度；$U(\cdot)$ 为开路电位函数。
- $F$、$R$、$T$ 为法拉第常数、气体常数、温度；$\alpha_a$、$\alpha_c$ 为传递系数；$k$ 为反应速率常数。
- $t_+^0$ 为迁移数；$\varepsilon_e$ 为孔隙率；$a_s$ 为比表面积；$\sigma^{\mathrm{eff}}$、$\kappa^{\mathrm{eff}}$ 为有效电导率。

## 3.3 固相：锂在颗粒里怎么走

颗粒内锂的扩散由球坐标 Fick 方程描述：

$$
\frac{\partial c_s}{\partial t}
= \frac{1}{r^2}\frac{\partial}{\partial r}
\left(r^2 D_s \frac{\partial c_s}{\partial r}\right),
\qquad 0<r<R_p
\tag{3.1}
$$

中心对称，表面通量由界面反应决定：

$$
\left.\frac{\partial c_s}{\partial r}\right|_{r=0}=0,
\qquad
-D_s\left.\frac{\partial c_s}{\partial r}\right|_{r=R_p}=j
\tag{3.2}
$$

其中 $j$ 由 3.5 节给出。扩散时间 $\tau_s=R_p^2/D_s$：半径减半，时间变为四分之一。这是电压滞后的主要来源。

## 3.4 液相：离子在电解液里怎么走

电解液中的锂离子同时扩散和迁移。对多孔电极做体积平均后：

$$
\varepsilon_e \frac{\partial c_e}{\partial t}
= \frac{\partial}{\partial x}
\left(D_e^{\mathrm{eff}}\frac{\partial c_e}{\partial x}\right)
+ (1-t_+^0)a_s j
\tag{3.3}
$$

右端第二项是界面反应产生的锂离子源；隔膜中 $a_s=0$、$\varepsilon_e=1$。集流体处没有离子进出：

$$
\left.\frac{\partial c_e}{\partial x}\right|_{x=0}=0,
\qquad
\left.\frac{\partial c_e}{\partial x}\right|_{x=L}=0
\tag{3.4}
$$

有效扩散系数用 Bruggeman 关系 $D_e^{\mathrm{eff}}=D_e\varepsilon_e^{b}$，$b$ 通常取 1.5。

## 3.5 界面：反应有多快

界面摩尔通量由 Butler–Volmer 方程给出：

$$
\begin{aligned}
j &= \frac{i_0}{F}
\left[
\exp\left(\frac{\alpha_a F\eta}{RT}\right)
-\exp\left(-\frac{\alpha_c F\eta}{RT}\right)
\right],\\
\eta &= \phi_s-\phi_e-U(c_{s,\mathrm{surf}}).
\end{aligned}
\tag{3.5}
$$

交换电流密度：

$$
i_0=Fk\,c_e^{\alpha_a}
\left(c_{s,\max}-c_{s,\mathrm{surf}}\right)^{\alpha_a}
c_{s,\mathrm{surf}}^{\alpha_c}
\tag{3.6}
$$

$\eta>0$ 时 $j>0$，锂离开固相；指数关系意味着强非线性——这是第 6 章必须用 Newton 法的原因。本书默认 $\alpha_a=\alpha_c=0.5$。

## 3.6 电荷守恒与"海平面"

固相中电子电流满足欧姆定律，电荷守恒给出：

$$
\frac{\partial}{\partial x}
\left(\sigma^{\mathrm{eff}}\frac{\partial \phi_s}{\partial x}\right)
= a_s F j
\tag{3.7}
$$

电解液中离子电流同时含扩散与迁移，浓溶液理论给出：

$$
\frac{\partial}{\partial x}
\left(\kappa^{\mathrm{eff}}\frac{\partial \phi_e}{\partial x}\right)
= \frac{\partial}{\partial x}
\left(\kappa_D^{\mathrm{eff}}\frac{\partial \ln c_e}{\partial x}\right)
- a_s F j
\tag{3.8}
$$

其中 $\kappa_D^{\mathrm{eff}}=\dfrac{2RT\kappa^{\mathrm{eff}}}{F}(1-t_+^0)\left(1+\dfrac{\partial \ln f_\pm}{\partial \ln c_e}\right)$。第一遍可令 $\kappa_D^{\mathrm{eff}}=0$，先理解扩散-迁移的主体。

两相电流之和等于外电流：

$$
\begin{aligned}
i_s+i_e &= i_{\mathrm{app}}(t),\\
i_s &= -\sigma^{\mathrm{eff}}\frac{\partial \phi_s}{\partial x},\\
i_e &= -\kappa^{\mathrm{eff}}\frac{\partial \phi_e}{\partial x}
+\kappa_D^{\mathrm{eff}}\frac{\partial \ln c_e}{\partial x}.
\end{aligned}
\tag{3.9}
$$

电位只有差值有意义，固定负集流体为参考：

$$
\phi_s(0,t)=0,
\qquad
V(t)=\phi_s(L,t)-\phi_s(0,t)
\tag{3.10}
$$

式 (3.7)–(3.10) 都不含时间导数，离散后是代数约束——这就是第 7 章 DAE 的来源。

## 3.7 边界、初值与"账要平"

初始浓度：

$$
c_s(r,0)=c_{s,0},
\qquad
c_e(x,0)=c_{e,0}
\tag{3.11}
$$

负极与正极的 $c_{s,0}$ 不同。电位初值不能随便给：必须满足 (3.9) 与 (3.10)，否则第一步残差就会爆表（第 7.6 节）。

边界条件汇总：固相中心对称与表面通量见 (3.2)；电解液在集流体零通量见 (3.4)；固相电流在隔膜界面为零，在集流体等于外电流；电解液电位在集流体满足 $i_e=0$。

## 3.8 时间尺度：为什么难算

用参数算出 $\tau_s$、$\tau_e$、双电层时间的数量级，直观看到"快慢相差几个数量级"。这是第 7 章"刚性"的定量依据。

## 3.9 模型阶梯：SPM / SPMe / DFN

- **SPM**：保留 (3.1)(3.2)(3.5)(3.6)(3.7)(3.9)(3.10)，令 $\phi_e=0$，去掉 (3.3)(3.4)(3.8)。
- **SPMe**：保留全部方程，但对 (3.3) 做渐近简化。
- **DFN**：全部保留。

自由度由网格设置决定，第 4 章会用同样的状态布局统计未知量。

## 3.10 代码：参数与状态布局

本书先用 `configs/demo.yaml` 打通结构；其中的数值是示意值，换真实材料（Chen2020）的路径见第 9 章与第 10 章。

`include/dfn/params.hpp` 把参数分成四组：电极（几何、输运、动力学、电导）、隔膜、电解液、电芯。`validate()` 在加载时检查量纲和范围：厚度、半径、扩散系数必须为正，孔隙率和迁移数必须在 $(0,1)$，初始浓度不能超过最大浓度。参数不合法时直接抛异常，避免"参数错了但程序照常跑完"。

状态向量按固定顺序打包，`include/dfn/state.hpp` 给出索引：

| 块 | 未知量 | 长度 |
| --- | --- | --- |
| $c_s$ 负极 | 每个 x 节点 × 每个 r 节点 | $N_{x,-}N_{r,-}$ |
| $c_s$ 正极 | 每个 x 节点 × 每个 r 节点 | $N_{x,+}N_{r,+}$ |
| $c_e$ | 全电池 x 节点 | $N_x$ |
| $\phi_s$ 负极 | 负极 x 节点 | $N_{x,-}$ |
| $\phi_s$ 正极 | 正极 x 节点 | $N_{x,+}$ |
| $\phi_e$ | 全电池 x 节点 | $N_x$ |

顺序不是随便定的：把同一物理量的未知量放在一起，第 4 章装配出的矩阵会呈现清晰的块结构，第 5 章可以据此选择 KLU 的排序。

```cpp
auto layout = dfn::make_layout(3, 2, 4, 5, 6);  // 64 unknowns
```

`dfn_state_tests` 会检查索引不重不漏、初值正确，以及非法参数会被拒绝：

```bash
ctest --test-dir build -C Release --output-on-failure
```

## 陷阱

- 符号：$j$ 的正方向定义为"锂离开固相"，换一篇文献可能相反。混用会导致充放电曲线互换。
- 把 (3.7)(3.8) 当成时间演化方程：它们没有时间导数，是代数约束。
- 忘记 (3.10) 的参考电位：Jacobian 会出现零空间，线性求解器报奇异。

## 检查点

1. (3.2) 中 $j$ 的符号如何影响颗粒表面浓度？
2. 为什么 (3.7)(3.8) 离散后是代数约束而不是微分方程？
3. SPM 去掉了哪些编号方程？为什么可以去掉？

## 代码

`configs/chen2020.yaml`、`include/dfn/params.hpp`、`include/dfn/state.hpp`

## 练习

1. 用 Chen2020 参数计算正负极的 $\tau_s$，比较差异。
2. 令 (3.3) 中 $j=0$，解释此时电解液浓度会怎样演化。
3. 写出 (3.9) 在 3 个 x 节点上的离散形式，并说明未知量个数。
