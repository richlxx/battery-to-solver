# 附录 E API 对照表

> 第 8 章会把本项目的接口接到 SUNDIALS；查表时不必记住 SUNDIALS 的全部名字。

## E.1 求解器接口

| 本书接口 | SUNDIALS | PyBaMM | 说明 |
| --- | --- | --- | --- |
| `IModel::residual(t, y, ydot, r)` | `IDAResFn` | `BaseModel.rhs/algebraic` | 残差形式 |
| `IModel::jacobian(...)` | `IDASetJacFn` | 由 CasADi 生成 | 稀疏 Jacobian |
| `IModel::root(...)` | `IDARootFn` | `events` | 事件函数 |
| `DaeSolver::integrate(...)` | `IDASolve` | `Solver.solve` | 时间推进入口 |
| `SolverOptions::rtol/atol` | `IDASetTolerances` | `Solver(atol, rtol)` | 容差 |
| `LinearSolver` | `SUNLinearSolver` | `Solver` 内部 | 稀疏线性求解 |

## E.2 状态与向量

| 本书 | SUNDIALS | 说明 |
| --- | --- | --- |
| `Vec`（Eigen） | `N_Vector` | 状态向量 |
| `SparseMatrix`（Eigen/CSR） | `SUNMatrix`（`SUNSparseMatrix`） | 稀疏矩阵 |
| `state.layout()` | 无对应 | 未知量索引映射 |

## E.3 常用配置项

| 本书 | SUNDIALS | 建议 |
| --- | --- | --- |
| `max_step` | `IDASetMaxStep` | 事件前限制步长 |
| `jac_reuse` | `IDASetJacUpdates` | 减少装配次数 |
| `linear_solver=KLU` | `SUNLinSol_KLU` | DAE 默认 |
| `linear_solver=Banded` | `SUNLinSol_Band` | 降级路径 |
