#pragma once

#include <functional>
#include <vector>

namespace dfn {

// 非线性方程组的阻尼牛顿法（可选回溯线搜索与变量缩放），
// 用于教学演示以及不依赖外部库时的小规模求解。

struct NewtonOptions {
  int max_iterations = 50;        // 最大迭代次数
  double tolerance = 1e-10;       // 残差范数收敛阈值
  double step_tolerance = 1e-14;  // 步长阈值，低于它认为已停滞
  bool use_line_search = true;    // 是否启用回溯线搜索
  std::vector<double> variable_scale;  // 逐分量变量尺度（典型量级）
};

struct NewtonResult {
  std::vector<double> x;                 // 最终解
  int iterations = 0;                    // 实际迭代次数
  bool converged = false;                // 是否满足收敛判据
  double residual_norm = 0.0;            // 最终残差范数
  std::vector<double> residual_history;  // 每次迭代的残差范数
  std::vector<double> step_history;      // 每次迭代的更新步长
};

// 回调类型：输入当前 x，返回 f(x)；Jacobian 回调返回逐行的 df/dx。
using ResidualFunction =
    std::function<std::vector<double>(const std::vector<double>&)>;
using JacobianFunction =
    std::function<std::vector<std::vector<double>>(const std::vector<double>&)>;

// 以 x0 为初值求解 f(x) = 0，返回最终解与迭代历史。
NewtonResult newton_solve(const ResidualFunction& residual,
                          const JacobianFunction& jacobian,
                          std::vector<double> x0,
                          const NewtonOptions& options = {});

}  // namespace dfn
