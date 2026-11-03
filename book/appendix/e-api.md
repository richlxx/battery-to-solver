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
| `SimulationConfig` | 参数、模型、网格、工况、容差 | `--config`、`--model`、`--nx`、`--nr`、`--rtol` |
| `ModelKind::spm/spme/dfn` | 模型阶梯 | `--model` |
| `SimulationMode::cc/cv/cccv/pulse` | 工况类型 | `--mode` |
| `SimulationReport::samples` | 采样表（时间、电压、电流、化学计量数、分段） | CSV 主体 |
| `SimulationReport::capacity` | 放电容量，由负极锂存量差算出（与 (9.7) 等价，且不依赖采样间隔） | 摘要 `cap=` |
| `SimulationReport::lithium_drift` | 总锂量相对漂移 | 摘要 `Li_drift=` |
| `run_simulation(config)` | 分段、事件、重启的驱动 | 整个程序 |
| `load_params_yaml(path)` | YAML 子集加载 | `--config` |
| `write_csv(report, path)` | CSV 输出 | `--out` |

与 SUNDIALS 的对应关系：每个工况分段是一次 `ida_integrate`；`SimulationReport` 里的 `steps/jacobian_evals` 是各分段 `IDAGetNum*` 的累加；`event_index` 只用于分段内部判断是否切换，不进入 CSV。

## E.5 模型阶梯（第 3/10 章）

| 本书接口 | 作用 | 对应方程 |
| --- | --- | --- |
| `CellModel` | SPM/SPMe/DFN 公共接口：在 `DaeModel` 上补初值、端电压、化学计量数等查询 | — |
| `SpmModel(params, n_r_neg, n_r_pos)` | 单颗粒模型，状态只有两个颗粒 + 两个表面浓度 | (3.1)–(3.6) |
| `SpmeModel(params, nx_neg, nx_sep, nx_pos, nr_neg, nr_pos)` | 加回电解液浓度与电位，反应仍均匀 | (3.5)、(3.7)–(3.10) |
| `DfnModel(params, layout, options)` | 全阶模型，局部反应与固相电位 | (3.1)–(3.15) |
| `set_current_density(J)` | 恒流控制；SPM/SPMe 的反应通量由此直接决定 | (3.5) |
| `set_applied_voltage(V)` | 恒压控制，只有 DFN 支持 | (9.5) |
| `terminal_voltage(y)` | 端电压；SPM/SPMe 反解 Butler–Volmer，DFN 读固相电位 | (3.6)、(3.10)、(3.14) |
| `particle_lithium / average_stoichiometry / total_lithium` | 颗粒锂存量、平均化学计量数、守恒检查 | — |

`--model spm|spme|dfn` 选择具体实现，默认 DFN；`cv`/`cccv` 只有 DFN 支持，SPM/SPMe 会在入口处报错。三个模型共用 `DaeModel`、`ida_integrate` 和 `CscBuilder`，所以对比实验只需要换模型对象；`labs/model_compare.cpp` 是同一个积分循环跑三遍。

## E.6 参数拟合与 CUDA（第 11 章）

| 本书接口 | 作用 | 对应公式 |
| --- | --- | --- |
| `dfn::fit::SpmBatch(params, nr, dt)` | 固定步长后向欧拉 SPM，CPU 与 GPU 的共同参考 | (11.2) 的正演 |
| `SpmBatch::step(current)` | 前进一步，返回端电压 | (3.6) |
| `rows_from_protocol(params, segments, nr, dt)` | 按分段电流生成时间/电压表 | 协议 C |
| `simulate_rows(params, rows, nr, dt)` | 按 CSV 中的电流序列重放 | 拟合残差 |
| `dfn::levenberg_marquardt(theta, residual, options)` | 通用 LM，阻尼正规方程 | (11.5) |
| `param_fit --generate/--sensitivity/--fit/--evaluate` | 数据生成、敏感度、拟合、交叉评估 | (11.2)–(11.6) |
| `spm_batch --samples N` | NVRTC CUDA 批量正演与计时 | (11.7) |

`SpmBatch` 的状态只有两个颗粒的径向浓度；表面浓度由 $c_{\mathrm{surf}}=c_{\mathrm{last}}-\alpha j$ 代数消去，所以每步只需要一次三对角求解，没有 Newton 迭代。`LmOptions` 的 `step` 是前向差分步长，`initial_lambda` 是初始阻尼；`LmResult` 返回最优参数、残差与迭代次数。
