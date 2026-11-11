// 定步长 BDF2 积分器，用于一般 DAE 残差形式 F(t, y, y') = 0。
// 质量矩阵信息隐含在 model.residual 里（代数行不依赖 y'）；
// Jacobian 由 dense_jacobian 组装 dF/dy + a0·dF/dy'，a0 为 BDF 的 y' 系数。

#include "dfn/solvers/bdf2_dae.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "dfn/solvers/newton.hpp"

namespace dfn {

// 从 t0 积分到 t_end：区间均分成 n_steps 步，实际步长 dt = (t_end-t0)/n_steps。
// 参数：model 为 DAE 模型（residual / jacobian / size）；y0 为初始状态（长度 n）；
//       h 为期望步长；options 为内层 Newton 的容差与最大迭代数。
// 返回：Bdf2DaeResult，t/y 为输出序列，另含总 Newton 迭代数、ok 与失败信息。
Bdf2DaeResult bdf2_dae_integrate(const DaeModel& model, std::vector<double> y0,
                                 double t0, double t_end, double h,
                                 const Bdf2DaeOptions& options) {
  if (h <= 0.0 || t_end <= t0) {
    throw std::invalid_argument("bdf2_dae_integrate: invalid time interval or step");
  }
  const int n = model.size();
  if (static_cast<int>(y0.size()) != n) {
    throw std::invalid_argument("bdf2_dae_integrate: initial state size mismatch");
  }

  // 步数用 ceil 并留 1e-12 余量，避免浮点误差多算一步；再把区间均分成等步长。
  const int n_steps = std::max(
      1, static_cast<int>(std::ceil((t_end - t0) / h - 1e-12)));
  const double dt = (t_end - t0) / static_cast<double>(n_steps);

  Bdf2DaeResult result;
  result.t.push_back(t0);
  result.y.push_back(y0);

  std::vector<double> y_prev = y0;
  std::vector<double> y_curr = y0;
  std::vector<double> yp(static_cast<std::size_t>(n), 0.0);
  std::vector<double> residual_buffer(static_cast<std::size_t>(n), 0.0);
  double t = t0;
  bool have_previous_step = false;

  for (int step = 0; step < n_steps; ++step) {
    const double t_new = t + dt;
    // 启动步没有 y_{n-1}，只能用 BDF1；第二步起用 BDF2。
    const bool use_bdf2 = have_previous_step;
    // a0 是 y' 对 y 的系数：BDF2 为 3/(2dt)，BDF1 为 1/dt；
    // Jacobian 的 dF/dy + a0·dF/dy' 必须与残差使用同一个 a0。
    const double a0 = use_bdf2 ? 3.0 / (2.0 * dt) : 1.0 / dt;

    // Newton 初值：BDF2 用线性外推 y_{n+1} ≈ 2·y_n - y_{n-1}，BDF1 用 y_n。
    std::vector<double> guess = y_curr;
    if (use_bdf2) {
      for (int j = 0; j < n; ++j) {
        guess[static_cast<std::size_t>(j)] =
            2.0 * y_curr[static_cast<std::size_t>(j)] -
            y_prev[static_cast<std::size_t>(j)];
      }
    }

    // 先用 BDF 公式把 y' 表示成已知历史值与待求 y 的组合，再代入残差 F(t, y, y')：
    //   BDF2: y' = (3y - 4y_n + y_{n-1}) / (2dt)
    //   BDF1: y' = (y - y_n) / dt
    // y' 写到共享缓冲区 yp 中，避免每次求值都分配内存。
    const auto residual = [&](const std::vector<double>& y) {
      for (int j = 0; j < n; ++j) {
        yp[static_cast<std::size_t>(j)] =
            use_bdf2
                ? (3.0 * y[static_cast<std::size_t>(j)] -
                   4.0 * y_curr[static_cast<std::size_t>(j)] +
                   y_prev[static_cast<std::size_t>(j)]) /
                      (2.0 * dt)
                : (y[static_cast<std::size_t>(j)] -
                   y_curr[static_cast<std::size_t>(j)]) /
                      dt;
      }
      model.residual(t_new, y.data(), yp.data(), residual_buffer.data());
      return residual_buffer;
    };

    // 用同一套 y' 计算组装 Jacobian dF/dy + a0·dF/dy'，
    // 保证线性化与残差完全一致（否则 Newton 会失去二阶收敛性）。
    const auto jacobian = [&](const std::vector<double>& y) {
      for (int j = 0; j < n; ++j) {
        yp[static_cast<std::size_t>(j)] =
            use_bdf2
                ? (3.0 * y[static_cast<std::size_t>(j)] -
                   4.0 * y_curr[static_cast<std::size_t>(j)] +
                   y_prev[static_cast<std::size_t>(j)]) /
                      (2.0 * dt)
                : (y[static_cast<std::size_t>(j)] -
                   y_curr[static_cast<std::size_t>(j)]) /
                      dt;
      }
      return dense_jacobian(model, t_new, y.data(), yp.data(), a0);
    };

    // 每个时间步解一次非线性方程组，并把迭代数累加用于诊断非线性难度；
    // 不收敛则停止积分，保留已完成部分的解，由调用方决定后续处理。
    NewtonOptions newton_options;
    newton_options.tolerance = options.newton_tolerance;
    newton_options.max_iterations = options.newton_max_iterations;
    const NewtonResult newton =
        newton_solve(residual, jacobian, guess, newton_options);
    result.newton_iterations += newton.iterations;

    if (!newton.converged) {
      result.ok = false;
      result.message = "Newton did not converge at t = " + std::to_string(t_new);
      break;
    }

    // 接受该步：滚动历史（y_prev/y_curr），推进时间并记录输出。
    y_prev = y_curr;
    y_curr = newton.x;
    t = t_new;
    result.t.push_back(t);
    result.y.push_back(y_curr);
    result.steps = step + 1;
    have_previous_step = true;
  }

  return result;
}

}  // namespace dfn
