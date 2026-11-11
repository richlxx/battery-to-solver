// 定步长 BDF2 积分器，用于显式给出右端项的常微分方程 y' = f(t, y)。
// 第一个时间步用 BDF1（隐式 Euler）启动，之后切换到二阶 BDF2；
// 每一步的隐式方程都用 Newton 求解，主要供教学与对照使用。

#include "dfn/solvers/bdf2.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "dfn/solvers/newton.hpp"

namespace dfn {

// 从 t0 积分到 t_end：把区间均分成 n_steps 步，实际步长 dt = (t_end-t0)/n_steps。
// 参数：rhs/jacobian 为右端项 f(t, y) 与 ∂f/∂y；y0 为初始状态（长度 n）；
//       h 为期望步长（向上取整成整数步）；options 为 Newton 容差与迭代上限。
// 返回：Bdf2Result，t[k] 与 y[k] 分别是第 k 个输出时刻及其状态向量。
// 索引约定：Jacobian 采用 J[i][j] = ∂f_i/∂y_j，与 newton_solve 保持一致。
Bdf2Result bdf2_integrate(const RhsFunction& rhs,
                          const RhsJacobian& jacobian,
                          std::vector<double> y0,
                          double t0,
                          double t_end,
                          double h,
                          const Bdf2Options& options) {
  if (h <= 0.0 || t_end <= t0) {
    throw std::invalid_argument("bdf2_integrate: invalid time interval or step");
  }
  if (y0.empty()) {
    throw std::invalid_argument("bdf2_integrate: empty initial state");
  }

  const int n = static_cast<int>(y0.size());
  // 步数用 ceil 并留 1e-12 余量，避免浮点误差多算一步；再把区间均分成等步长。
  const int n_steps = std::max(
      1, static_cast<int>(std::ceil((t_end - t0) / h - 1e-12)));
  const double dt = (t_end - t0) / static_cast<double>(n_steps);

  Bdf2Result result;
  result.t.push_back(t0);
  result.y.push_back(y0);

  std::vector<double> y_prev = y0;
  std::vector<double> y_curr = y0;
  double t = t0;
  bool have_previous_step = false;

  for (int step = 0; step < n_steps; ++step) {
    const double t_new = t + dt;
    // 启动步没有 y_{n-1} 可用，只能用一阶 BDF1；第二步起才切换到 BDF2。
    const bool use_bdf2 = have_previous_step;

    // Newton 初值：BDF2 用线性外推 y_{n+1} ≈ 2·y_n - y_{n-1}，
    // BDF1 直接用 y_n；好的初值能明显减少内层迭代次数。
    std::vector<double> guess = y_curr;
    if (use_bdf2) {
      for (int j = 0; j < n; ++j) {
        guess[static_cast<std::size_t>(j)] =
            2.0 * y_curr[static_cast<std::size_t>(j)] -
            y_prev[static_cast<std::size_t>(j)];
      }
    }

    // 残差就是 BDF 公式移项的结果：
    //   BDF2: y - (4/3)·y_n + (1/3)·y_{n-1} - (2/3)·dt·f(t_{n+1}, y) = 0
    //   BDF1: y - y_n - dt·f(t_{n+1}, y) = 0
    // f 在 t_{n+1} 与未知的 y 处取值，所以是隐式方程，需要 Newton 迭代求解。
    const auto residual = [&](const std::vector<double>& y) {
      const std::vector<double> f = rhs(t_new, y);
      if (static_cast<int>(f.size()) != n) {
        throw std::invalid_argument("bdf2_integrate: RHS size mismatch");
      }
      std::vector<double> r(static_cast<std::size_t>(n), 0.0);
      for (int j = 0; j < n; ++j) {
        if (use_bdf2) {
          r[static_cast<std::size_t>(j)] =
              y[static_cast<std::size_t>(j)] -
              (4.0 / 3.0) * y_curr[static_cast<std::size_t>(j)] +
              (1.0 / 3.0) * y_prev[static_cast<std::size_t>(j)] -
              (2.0 / 3.0) * dt * f[static_cast<std::size_t>(j)];
        } else {
          r[static_cast<std::size_t>(j)] =
              y[static_cast<std::size_t>(j)] -
              y_curr[static_cast<std::size_t>(j)] -
              dt * f[static_cast<std::size_t>(j)];
        }
      }
      return r;
    };

    // 隐式方程对 y 的 Jacobian：I - factor·∂f/∂y，
    // 其中 factor = (2/3)·dt（BDF2）或 dt（BDF1），与残差中的系数一一对应。
    const auto jac = [&](const std::vector<double>& y) {
      const std::vector<std::vector<double>> dfdy = jacobian(t_new, y);
      if (static_cast<int>(dfdy.size()) != n) {
        throw std::invalid_argument("bdf2_integrate: Jacobian size mismatch");
      }
      std::vector<std::vector<double>> J(
          static_cast<std::size_t>(n),
          std::vector<double>(static_cast<std::size_t>(n), 0.0));
      const double factor = use_bdf2 ? (2.0 / 3.0) * dt : dt;
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          J[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
              (i == j ? 1.0 : 0.0) -
              factor * dfdy[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
      }
      return J;
    };

    // 每个时间步解一次非线性方程组；Newton 不收敛就放弃后续步并标记失败。
    NewtonOptions newton_options;
    newton_options.tolerance = options.newton_tolerance;
    newton_options.max_iterations = options.newton_max_iterations;
    const NewtonResult newton =
        newton_solve(residual, jac, guess, newton_options);

    if (!newton.converged) {
      result.ok = false;
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
