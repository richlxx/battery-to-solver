#pragma once

#include <functional>
#include <vector>

namespace dfn {

// 常微分方程组 y' = f(t, y) 的定步长 BDF2 时间推进，第一步用 BDF1 启动。
// 面向教学：固定步长便于观察二阶精度与 L-稳定性。

using RhsFunction =
    std::function<std::vector<double>(double, const std::vector<double>&)>;
using RhsJacobian = std::function<std::vector<std::vector<double>>(
    double, const std::vector<double>&)>;

struct Bdf2Options {
  double newton_tolerance = 1e-10;  // 每步牛顿迭代的收敛阈值
  int newton_max_iterations = 20;   // 每步牛顿迭代次数上限
};

struct Bdf2Result {
  std::vector<double> t;               // 采样时间点
  std::vector<std::vector<double>> y;  // 各时间点对应的状态
  bool ok = true;                      // 是否所有步都成功
  int steps = 0;                       // 实际执行的时间步数
};

// 从 (t0, y0) 以步长 h 积分到 t_end，返回完整时间序列。
Bdf2Result bdf2_integrate(const RhsFunction& rhs,
                          const RhsJacobian& jacobian,
                          std::vector<double> y0,
                          double t0,
                          double t_end,
                          double h,
                          const Bdf2Options& options = {});

}  // namespace dfn
