# 附录 D 术语表（中英对照）

| 中文 | 英文 | 一句话解释 |
| --- | --- | --- |
| 数值解 | numerical solution | 用有限步近似算出的解，不保证精确 |
| 解析解 | analytical solution | 能写成公式直接计算的解 |
| 离散化 | discretization | 把连续问题切成有限个未知数 |
| 网格 | mesh / grid | 离散后的小块划分 |
| 控制体 | control volume | 有限体积法里记账的一个小格子 |
| 通量 | flux | 单位时间穿过单位面积的量 |
| 边界条件 | boundary condition | 区域边界上的约束 |
| 初始条件 | initial condition | 时间起点的状态 |
| 截断误差 | truncation error | 近似公式丢掉的高阶项 |
| 舍入误差 | round-off error | 浮点数精度造成的误差 |
| 收敛阶 | order of convergence | 步长减半时误差下降的速度 |
| 观测阶 | observed order | 从数值实验测出的收敛阶 |
| 刚性 | stiffness | 快慢过程混杂，显式方法被迫用极小步长 |
| 显式方法 | explicit method | 用当前信息直接算下一步 |
| 隐式方法 | implicit method | 下一步出现在方程两边，需要求解 |
| 稳定性 | stability | 误差不会被算法放大 |
| 稀疏矩阵 | sparse matrix | 绝大多数元素为零的矩阵 |
| 填充 | fill-in | 消元过程中新产生的非零元 |
| Jacobian | Jacobian matrix | 多维函数的偏导数矩阵 |
| 残差 | residual | 把当前解代回方程后的偏差 |
| 阻尼 | damping | 缩小 Newton 步长以避免发散 |
| 缩放 | scaling | 把不同量级的变量化成相近数量级 |
| 质量矩阵 | mass matrix | DAE 中乘以时间导数的那部分 |
| 微分代数方程 | DAE | 同时含微分方程和代数约束的系统 |
| 事件 | event | 满足某个条件时暂停或切换工况 |
| 降阶模型 | reduced-order model | 用更少未知量近似原模型 |
| 开路电压 | OCV | 无电流时的端电压 |
| 过电位 | overpotential | 实际电位与平衡电位之差 |
| 化学计量比 | stoichiometry | 电极中锂占可容纳位的比例 |
| 比表面积 | specific surface area | 单位体积电极的活性面积 |
| 正问题 | forward problem | 已知参数算电压/浓度 |
| 反问题 | inverse problem | 已知测量反推参数 |
| 可辨识性 | identifiability | 数据能否唯一确定参数 |
| 敏感度 | sensitivity | 参数小幅变化引起的输出变化 |
| 条件数 | condition number | 矩阵对扰动的放大倍数；越大越难反演 |
| 置信区间 | confidence interval | 参数估计的不确定范围 |
| Levenberg–Marquardt | LM | 带阻尼的 Gauss–Newton，用于非线性最小二乘 |
| GITT | galvanostatic intermittent titration | 恒流脉冲 + 静置，分离动力学与扩散 |
| NVRTC | NVIDIA Runtime Compilation | 运行时把 CUDA C++ 编译成 GPU 代码 |
