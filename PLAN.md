# 《从一块电池到求解器》立项计划

## 1. 定位

面向**完全没有数值计算经验**的读者，用锂离子电池 P2D/DFN 模型作为唯一案例，从"什么是数值解"一路讲到能独立运行、修改并验证一个 C++ 求解器（手写 BDF + SUNDIALS IDA/KLU）。

不是把专业书写薄，而是把入口放低：**先看到现象，再理解原理，最后写出代码。**

## 2. 读者假设

| 读者具备 | 读者不需要具备 |
| --- | --- |
| 高中数学（函数、斜率、面积） | 数值分析、微积分熟练、线性代数熟练 |
| 能读懂简单代码（变量、循环、函数） | 写过数值程序、用过求解器 |
| 对电池有常识性认识 | 电化学背景、PDE 知识 |
| 一台能装 CMake 的电脑 | C++ 项目经验 |

补充路径：

- 数学不熟 → 第 1 章按需查阅，附录 A 速查。
- 编程不熟 → 附录 B 的 C++ 最小入门；每章代码先运行再阅读。
- 电化学零基础 → 第 3 章从"锂离子怎么进出颗粒"讲起。

## 3. 教学原则

1. **现象先行**：每章用一个可观察的现象开场（电压掉得快、曲线爆炸、迭代不动），再引入方法。
2. **最小例子**：新方法先在 1 个未知量、3 个格点上讲透，再放到 DFN 上。
3. **一个模型贯穿**：DFN 只在第 3 章定义一次，后续章节反复使用同一份残差。
4. **代码四步法**：先运行 → 改参数 → 读关键函数 → 自己写（填空练习）。
5. **数学分级**：正文只保留必要推导；长推导放进"公式"提示框或附录，第一遍可跳过。
6. **错误驱动**：每章至少一个失败案例（显式法爆炸、Newton 发散、KLU 奇异），先看失败再学解法。
7. **可验证**：每个定量结论都有测试或图，读者能亲手复现。
8. **难度标记**：★ 核心必读、★★ 需要耐心、★★★ 选学，第一遍可跳过。

## 4. 章节地图

| 部分 | 章 | 标题 | 难度 | 读者产出 |
| --- | --- | --- | --- | --- |
| 0 起步 | 0 | 从一块电池说起 | ★ | 跑通第一个程序，画出 OCV 曲线 |
| 0 起步 | 1 | 数学预备：只讲用得上的 | ★ | 看懂后续公式，能手算欧拉一步 |
| 0 起步 | 2 | 数值计算第一课：近似、误差与收敛 | ★ | 画出误差-步长图，理解"阶" |
| 1 模型 | 3 | 案例模型：P2D/DFN | ★ | 参数集与状态布局 |
| 2 离散 | 4 | 空间离散：把连续世界切成小块 | ★ | 装配器 + 守恒测试 |
| 2 线性 | 5 | 解线性方程组：从消元到稀疏 | ★ | Eigen/KLU 双后端 |
| 3 非线性 | 6 | 解非线性方程：Newton 法 | ★★ | 手写 Newton + 收敛历史 |
| 3 时间 | 7 | 随时间前进：显式、隐式与刚性 | ★★ | 最小自适应 BDF2 |
| 3 工业级 | 8 | 站在巨人肩上：IDA + KLU | ★★ | IDAKLU 后端 |
| 4 集成 | 9 | 组装完整求解器：跑真实工况 | ★ | CLI + 三层验证报告 |
| 5 延伸 | 10 | 走得更远（选学） | ★★★ | SPM/SPMe 对照与扩展地图 |

附录：A 数学与参数速查、B 环境构建与 C++ 最小入门、C 调试清单、D 术语表、E API 对照、F 习题提示与答案、G 参考文献。

## 5. 方法 → 模型 → 代码 映射

| 方法 | 在 DFN 中的落点 | 章 | 代码 |
| --- | --- | --- | --- |
| 泰勒展开与差分 | 误差分析、通量近似 | 2, 4 | `labs/diff_error.cpp` |
| 有限体积法 | 球壳固相扩散、电解液 1D 输运 | 4 | `include/dfn/grid.hpp`, `operators.hpp` |
| 稀疏直接法 | 块三对角 + 边界箭头 Jacobian | 5 | `include/dfn/linalg/{eigen,klu}.hpp` |
| Newton 法 | 稳态与每步内的非线性方程 | 6 | `include/dfn/solvers/newton.hpp` |
| BDF/自适应步长 | 固相扩散与双电层造成的刚性 | 7 | `include/dfn/solvers/bdf2.hpp` |
| IDA + KLU | 生产级 DAE 求解 | 8 | `include/dfn/solvers/idaklu.hpp` |
| 事件检测 | 电压截止、CC→CV | 7, 9 | `include/dfn/sim/events.hpp` |
| 降阶模型 | SPM/SPMe | 10 | `include/dfn/reduction/` |

## 6. 代码架构

```text
numerical-p2d/
├─ labs/                    # 单文件小实验，零依赖，先于库代码出现
│  ├─ ocv_curve.cpp         # 第 0 章：读参数、打印曲线
│  ├─ diff_error.cpp        # 第 2 章：差分误差实验
│  ├─ euler_stiff.cpp       # 第 7 章：显式法爆炸与隐式法对比
│  └─ newton_scalar.cpp     # 第 6 章：单变量 Newton
├─ include/dfn/
│  ├─ params.hpp            # 参数集、单位、缩放因子
│  ├─ grid.hpp              # x 与 r 两套网格
│  ├─ state.hpp             # 未知量打包/解包、索引映射
│  ├─ operators.hpp         # FV 算子与 CSR 装配
│  ├─ model.hpp             # DFN 残差/Jacobian 装配
│  ├─ linalg/               # eigen.hpp, klu.hpp
│  ├─ solvers/              # newton, bdf2, idaklu, dae_solver
│  ├─ sim/                  # 工况、事件、输出、CLI
│  └─ reduction/            # SPM/SPMe
├─ src/                     # 实现与 main
├─ bindings/python.cpp      # pybind11
├─ tests/                   # 守恒、收敛阶、回归、一致性
├─ bench/                   # N 扫描、后端对比
├─ configs/                 # chen2020.yaml
├─ scripts/                 # 绘图、PyBaMM 对比
└─ CMakeLists.txt
```

统一接口（第 8 章对照实验的基础）：

```cpp
class DaeModel {                       // 残差形式 F(t, y, y') = 0
public:
  virtual int size() const = 0;
  virtual void residual(double t, const double* y, const double* yp,
                        double* out) const = 0;
  virtual CscPattern jacobian_pattern() const = 0;
  virtual void jacobian(double t, const double* y, const double* yp, double cj,
                        CscBuilder& out) const = 0;
  virtual std::vector<int> variable_types() const;   // 1 = 微分, 0 = 代数
};

struct IdaOptions {
  double rtol = 1e-6;
  std::vector<double> atol;
  double max_step = 0;
  bool supply_jacobian = true;
  bool calculate_consistent_initial = true;
  LinearSolverKind linear_solver = LinearSolverKind::klu;
};

SolveReport ida_integrate(const DaeModel& model, const std::vector<double>& y0,
                          const std::vector<double>& yp0, double t0,
                          double t_end, const IdaOptions& options = {});
```

手写 BDF2-DAE 使用同一个 `DaeModel`，只是时间推进换成固定步长，因此两套后端可以逐点对照。

## 7. 构建与依赖

- C++20、CMake ≥ 3.24、vcpkg manifest（`vcpkg.json`）。
- 依赖：`eigen3`、`sundials[klu]`（使用 `cmake/vcpkg-overlay` 的 overlay 端口）、`pybind11`、`yaml-cpp`、`catch2`。
- `labs/` 下的单文件程序只依赖标准库，可独立编译，保证第 0–2 章不卡在环境上。
- Python 侧用 `uv` 建独立 venv（numpy/matplotlib；pybamm 可选）。
- 平台：Windows 主开发 + Linux 用于基准；无 SUNDIALS 时不构建 IDA 后端，IDA 内部可在 KLU 与稠密求解器之间切换。

## 8. 质量门禁

### 8.1 教学可读性（每章自检）

- [ ] 开头有一个**不需要前置知识**的现象或问题。
- [ ] 至少一个可运行的 `labs/` 程序或命令。
- [ ] 至少一个失败案例，并解释失败原因。
- [ ] 首次出现的术语有一句"人话解释"和中英对照。
- [ ] 每个公式前有文字版说明；长推导可跳过且不影响主线。
- [ ] 公式统一用 `$...$` / `$$...$$`，块级公式前后空行，编号用 `\tag{章.序号}`。
- [ ] 章末 3 个检查点，答案在附录 F。
- [ ] 全章不依赖"读者应该知道"的隐含知识。

### 8.2 数值正确性

| 验证项 | 方法 | 通过标准 |
| --- | --- | --- |
| 锂守恒 | 总量随时间漂移 | 相对误差 < 1e-10/步 |
| 电荷守恒 | 离散电流恒等式 | 残差 < 1e-12 |
| 空间收敛阶 | N 加倍 | 观测阶 ≥ 1.9（当前表面单侧差分约 1.4，二阶重构后达标） |
| 时间收敛阶 | BDF2 vs 参考解 | 观测阶 ≥ 1.9 |
| SPM 精度 | 与球扩散级数解对比 | 电压误差 < 1 mV |
| DFN 回归 | 与 PyBaMM 同参数对比 | 1C 电压 RMSE < 10 mV；≤C/4 时 < 3 mV（省略 $\kappa_D$） |
| 后端一致性 | BDF2 vs IDAKLU | 状态差 < 10×tol |

### 8.3 试读门槛

每完成一部分，请一位没有数值计算经验的读者按章操作；记录卡住的位置，作为修订清单。没有试读反馈的章节不标记"完成"。

## 9. 里程碑

| 批次 | 内容 | 交付物 | 验收 | 预估 |
| --- | --- | --- | --- | --- |
| M0 | 仓库骨架、构建、`labs/` 框架 | 可编译的 `dfn-sim --help` 与 3 个 lab | Windows/Linux 构建通过 | 1–2 天 |
| M1 | 第 0–2 章 | 起步三章 + 2 个 lab + 误差图 | 新手能独立跑通 | 5–6 天 |
| M2 | 第 3–4 章 | 模型定义 + FV 装配 | 守恒与二阶收敛 | 5–6 天 |
| M3 | 第 5–6 章 | 稀疏后端 + Newton | 稳态 SPM/DFN 收敛 | 5–6 天 |
| M4 | 第 7–8 章 | BDF2 + IDAKLU | 两后端逐点一致 | 5–6 天 |
| M5 | 第 9 章 | CLI、验证矩阵、性能 | 1C 电压 RMSE < 10 mV（PyBaMM，省略 $\kappa_D$） | 3–4 天 |
| M6 | 第 10 章 + 附录 + 试读修订 | 全书与速查 | 命令可复现、试读通过 | 4–5 天 |

## 10. 风险与对策

| 风险 | 影响 | 对策 |
| --- | --- | --- |
| 内容过载，新手中途放弃 | 完读率低 | 难度标记、选学章节、每章一个可见成果 |
| C++ 门槛 | 第 3 章后掉队 | labs 单文件零依赖；完整参考代码；填空式练习 |
| 数学焦虑 | 跳过推导也读不懂 | 每章"直觉"框；公式可跳过；第 1 章按需查阅 |
| 篇幅膨胀 | 写作周期失控 | 每章 3500–6000 字、代码 200–500 行；超出的入附录 |
| 只跑代码不理解 | 无法迁移 | 每章检查点 + 失败案例 + 手算练习 |
| KLU/SuiteSparse 构建失败 | 第 5/8 章受阻 | vcpkg 优先；banded 降级并在附录 B 说明 |
| PyBaMM 版本漂移 | 对比失效 | 锁定版本并记录参数集来源 |

## 11. 约定

- `IDAKLU` 指 SUNDIALS **IDA** + **KLU**（PyBaMM 命名）。
- 代码标识符英文，讲解中文；代码注释只解释非显然的原因。
- 缩放变量后缀 `_s`，物理量保留原单位，转换集中在 `params.hpp`。
- 公式统一 SI 单位；正文首次出现处给出符号含义。
- 每章源码可独立构建，对应 `tests/` 同名测试。

## 12. 下一步

1. 确认本计划与 [book/SUMMARY.md](book/SUMMARY.md) 大纲。
2. M0：搭 CMake/vcpkg 骨架与 `labs/`，验证 Windows 上 `sundials[klu]` 可编译。
3. M1：写第 0 章与 `labs/ocv_curve.cpp`，让读者第一小时就有可见成果。
