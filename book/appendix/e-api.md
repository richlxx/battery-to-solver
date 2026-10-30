# 附录 E API 对照表

> 第 8 章的接口已经接到 SUNDIALS；查表时不必记住 SUNDIALS 的全部名字。

## E.1 求解器接口

| 本书接口 | SUNDIALS | PyBaMM | 说明 |
| --- | --- | --- | --- |
| `DaeModel::residual(t, y, yp, r)` | `IDAResFn` | `BaseModel.rhs/algebraic` | 残差形式 |
| `DaeModel::jacobian(t, y, yp, cj, J)` | `IDALsJacFn` | 由 CasADi 生成 | 稀疏 Jacobian，`cj` 由 IDA 传入 |
| `DaeModel::jacobian_pattern()` | `SUNSparseMatrix` 的结构 | 稀疏模式 | 只需构造一次 |
| `DaeModel::variable_types()` | `IDASetId` | 变量分类 | 1 = 微分，0 = 代数 |
| `ida_integrate(model, y0, yp0, t0, t_end, opt)` | `IDASolve` | `Solver.solve` | 时间推进入口 |
| `IdaOptions::rtol/atol` | `IDASStolerances` / `IDASVtolerances` | `Solver(atol, rtol)` | 容差 |
| `LinearSolverKind::klu` | `SUNLinSol_KLU` | 默认后端 | 稀疏直接求解 |
| `LinearSolverKind::dense` | `SUNLinSol_Dense` | — | 小规模对照/降级 |
| `ida_error_message(flag)` | `IDAProcessError` 的返回值 | — | 错误码翻译 |

## E.2 状态与向量

| 本书 | SUNDIALS | 说明 |
| --- | --- | --- |
| `std::vector<double>` | `N_Vector` | 状态向量（`N_VNew_Serial`） |
| `CscPattern` | `SUNSparseMatrix`（CSC） | 稀疏结构 |
| `CscBuilder` | `SM_DATA_S` 写入 | 按固定模式装配 |
| `state.layout()` | 无对应 | 未知量索引映射 |

## E.3 常用配置项

| 本书 | SUNDIALS | 建议 |
| --- | --- | --- |
| `max_step` | `IDASetMaxStep` | 事件前限制步长 |
| `calculate_consistent_initial` | `IDACalcIC` | 代数约束不满足时打开 |
| `supply_jacobian=false` | 不设 `IDASetJacFn` | 仅稠密矩阵可用（内部差商） |
| 终点精确停止 | `IDASetStopTime` | 返回 `IDA_TSTOP_RETURN` |
| 规模较大 | `SUNLinSol_KLU` | 稀疏装配 + 符号复用 |

## E.4 第 9 章的仿真接口

| 本书接口 | 作用 | 对应命令行 |
| --- | --- | --- |
| `SimulationConfig` | 参数、网格、工况、容差 | `--config`、`--nx`、`--nr`、`--rtol` |
| `SimulationMode::cc/cv/cccv/pulse` | 工况类型 | `--mode` |
| `SimulationReport::samples` | 采样表（时间、电压、电流、化学计量数、分段） | CSV 主体 |
| `SimulationReport::capacity` | 放电容量，由负极锂存量差算出（与 (9.7) 等价，且不依赖采样间隔） | 摘要 `cap=` |
| `SimulationReport::lithium_drift` | 总锂量相对漂移 | 摘要 `Li_drift=` |
| `run_simulation(config)` | 分段、事件、重启的驱动 | 整个程序 |
| `load_params_yaml(path)` | YAML 子集加载 | `--config` |
| `write_csv(report, path)` | CSV 输出 | `--out` |

与 SUNDIALS 的对应关系：每个工况分段是一次 `ida_integrate`；`SimulationReport` 里的 `steps/jacobian_evals` 是各分段 `IDAGetNum*` 的累加；`event_index` 只用于分段内部判断是否切换，不进入 CSV。

## E.5 降阶模型（第 10 章）

| 本书接口 | 作用 | 对应方程 |
| --- | --- | --- |
| `SpmModel(params, n_r_neg, n_r_pos)` | 单颗粒模型，状态只有两个颗粒 + 两个表面浓度 | (10.1)–(10.3) |
| `SpmModel::set_current_density(J)` | 反应电流由外加电流直接决定 | (10.1) |
| `SpmModel::reaction_flux(positive)` | 均匀反应电流 $j$ | (10.1) |
| `SpmModel::terminal_voltage(y)` | $U_{\mathrm{pos}}-U_{\mathrm{neg}}+\eta_{\mathrm{pos}}-\eta_{\mathrm{neg}}$，过电位用标量 Newton 反解 | (10.2)(10.3) |
| `particle_lithium / average_stoichiometry` | 颗粒锂存量与平均化学计量数 | — |

SPM 与 DFN 共用 `DaeModel`、`ida_integrate` 和 `CscBuilder`，所以对比实验只需要换模型对象；`labs/spm_compare.cpp` 就是同一个积分循环跑两遍。
