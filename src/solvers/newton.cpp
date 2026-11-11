// 阻尼 Newton 求解器：求解非线性方程组 F(x) = 0。
// 提供变量缩放（右预条件）与回溯线搜索，供 BDF2 时间步的内层迭代调用。

#include "dfn/solvers/newton.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "dfn/linalg/dense.hpp"

namespace dfn {

namespace {

// 无穷范数 ‖v‖_∞：Newton 的残差范数与步长范数都用它衡量。
double max_abs(const std::vector<double>& values) {
  double norm = 0.0;
  for (double value : values) {
    norm = std::max(norm, std::abs(value));
  }
  return norm;
}

}  // namespace

// 用阻尼 Newton 法求解 F(x) = 0。
// 参数：
//   residual —— 残差函数 F(x)，返回值长度必须为 n；
//   jacobian —— Jacobian 函数，返回 n×n 嵌套 vector，J[i][j] = dF_i/dx_j；
//   x0       —— 初始猜测（按值传入，函数内部直接把它当解向量 x 迭代）；
//   options  —— 最大迭代数、残差容差、步长容差、线搜索开关与变量缩放。
// 返回：NewtonResult，含末次迭代点 x、迭代次数、收敛标志与残差/步长历史。
NewtonResult newton_solve(const ResidualFunction& residual,
                          const JacobianFunction& jacobian,
                          std::vector<double> x0,
                          const NewtonOptions& options) {
  const int n = static_cast<int>(x0.size());
  if (n == 0) {
    throw std::invalid_argument("newton_solve: empty initial guess");
  }

  // 变量缩放：默认全 1；非空时长度必须为 n 且每个分量严格为正。
  std::vector<double> scale = options.variable_scale;
  if (scale.empty()) {
    scale.assign(static_cast<std::size_t>(n), 1.0);
  }
  if (static_cast<int>(scale.size()) != n) {
    throw std::invalid_argument("newton_solve: scale size mismatch");
  }
  for (double value : scale) {
    if (value <= 0.0) {
      throw std::invalid_argument("newton_solve: scale must be positive");
    }
  }

  NewtonResult result;
  std::vector<double> x = std::move(x0);

  // 主迭代：每轮先算残差，再解线性化方程 J·δ = -F，最后线搜索更新 x。
  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    const std::vector<double> residual_values = residual(x);
    if (static_cast<int>(residual_values.size()) != n) {
      throw std::invalid_argument("newton_solve: residual size mismatch");
    }
    const double residual_norm = max_abs(residual_values);
    result.residual_history.push_back(residual_norm);

    // 收敛判据：残差无穷范数低于 tolerance 即接受当前点。
    if (residual_norm < options.tolerance) {
      result.x = x;
      result.iterations = iteration;
      result.converged = true;
      result.residual_norm = residual_norm;
      return result;
    }

    const std::vector<std::vector<double>> J = jacobian(x);
    if (static_cast<int>(J.size()) != n) {
      throw std::invalid_argument("newton_solve: Jacobian size mismatch");
    }
    for (const auto& row : J) {
      if (static_cast<int>(row.size()) != n) {
        throw std::invalid_argument("newton_solve: Jacobian must be square");
      }
    }

    // 右预条件：令 x = S·z，则缩放后的 Jacobian 是 J·S（每列乘 scale[j]）；
    // 先解出 z 的增量，再换算回真实步长 step_j = scale_j * scaled_step_j。
    std::vector<std::vector<double>> scaled_J = J;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        scaled_J[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *=
            scale[static_cast<std::size_t>(j)];
      }
    }

    // 线性求解 J·δ = -F（稠密 LU 分解，内部含部分选主元），得到 Newton 增量。
    std::vector<double> rhs(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
      rhs[static_cast<std::size_t>(i)] = -residual_values[static_cast<std::size_t>(i)];
    }
    const std::vector<double> scaled_step = solve_dense(scaled_J, rhs);

    std::vector<double> step(static_cast<std::size_t>(n), 0.0);
    for (int j = 0; j < n; ++j) {
      step[static_cast<std::size_t>(j)] =
          scale[static_cast<std::size_t>(j)] * scaled_step[static_cast<std::size_t>(j)];
    }
    const double step_norm = max_abs(scaled_step);
    result.step_history.push_back(step_norm);

    // 回溯线搜索（Armijo 型充分下降条件）：
    // 从 α = 1 开始折半，直到 ‖F(x + α·δ)‖_∞ <= (1 - 1e-4·α)·‖F(x)‖_∞，
    // 下限 α > 1e-8；条件始终不满足时就用最后一个很小的步长继续迭代。
    double alpha = 1.0;
    if (options.use_line_search) {
      const auto trial_norm = [&](double factor) {
        std::vector<double> trial(static_cast<std::size_t>(n), 0.0);
        for (int j = 0; j < n; ++j) {
          trial[static_cast<std::size_t>(j)] =
              x[static_cast<std::size_t>(j)] +
              factor * step[static_cast<std::size_t>(j)];
        }
        return max_abs(residual(trial));
      };

      double norm = trial_norm(alpha);
      while (norm > (1.0 - 1e-4 * alpha) * residual_norm && alpha > 1e-8) {
        alpha *= 0.5;
        norm = trial_norm(alpha);
      }
    }

    for (int j = 0; j < n; ++j) {
      x[static_cast<std::size_t>(j)] +=
          alpha * step[static_cast<std::size_t>(j)];
    }

    // 步长已足够小：迭代停滞，以当前点结束，是否收敛由最终残差决定。
    if (step_norm < options.step_tolerance) {
      const double final_norm = max_abs(residual(x));
      result.x = x;
      result.iterations = iteration + 1;
      result.converged = final_norm < options.tolerance;
      result.residual_norm = final_norm;
      return result;
    }
  }

  // 迭代次数用尽仍未收敛：返回最后迭代点并标记 converged = false。
  result.x = x;
  result.iterations = options.max_iterations;
  result.converged = false;
  result.residual_norm = max_abs(residual(x));
  return result;
}

}  // namespace dfn
