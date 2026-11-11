#pragma once

#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <vector>

#include "dfn/linalg/dense.hpp"

namespace dfn {

// Levenberg-Marquardt 最小二乘：用正演模型残差对参数做拟合。
// 头文件内为模板实现，残差回调由调用方提供（例如 SPM 批量模型）。

struct LmOptions {
  int max_iterations = 50;       // 外层迭代次数上限
  int max_backtracks = 12;       // 每次迭代内放大阻尼重试的次数上限
  double initial_lambda = 1e-3;  // 初始阻尼系数
  double step = 1e-6;            // 参数空间前向差分步长
  double xtol = 1e-9;            // 参数步长收敛阈值（相对）
  double ftol = 1e-12;           // 代价下降量收敛阈值（相对）
};

struct LmResult {
  bool ok = false;  // 是否正常结束（收敛或到达迭代上限）
  int iterations = 0;    // 外层迭代次数
  int evaluations = 0;   // 残差求值次数（含差分扰动）
  double cost = std::numeric_limits<double>::quiet_NaN();  // 0.5 * ||r||^2
  std::vector<double> theta;     // 最优参数
  std::vector<double> residual;  // 最优参数处的残差向量
};

// 最小化 0.5 * ||r(theta)||^2。回调返回 false 表示正演模型失败，
// 该试探点按被拒绝的步处理。
template <typename Residual>
LmResult levenberg_marquardt(std::vector<double> theta, Residual&& residual,
                             const LmOptions& options = {}) {
  LmResult result;
  const int n = static_cast<int>(theta.size());

  std::vector<double> r;
  auto evaluate = [&](const std::vector<double>& x, std::vector<double>& out) {
    if (!residual(x, out)) {
      return false;
    }
    ++result.evaluations;
    return true;
  };

  if (!evaluate(theta, r)) {
    return result;
  }
  const std::size_t m = r.size();
  auto cost_of = [](const std::vector<double>& values) {
    double sum = 0.0;
    for (const double value : values) {
      sum += value * value;
    }
    return 0.5 * sum;
  };

  // cost 为当前代价，lambda 为阻尼；result 始终保存当前最优解。
  double cost = cost_of(r);
  double lambda = options.initial_lambda;
  result.theta = theta;
  result.residual = r;
  result.cost = cost;

  std::vector<double> r_trial(m, 0.0);
  std::vector<double> r_perturbed(m, 0.0);
  std::vector<std::vector<double>> jacobian(
      m, std::vector<double>(static_cast<std::size_t>(n), 0.0));
  std::vector<std::vector<double>> normal(
      static_cast<std::size_t>(n),
      std::vector<double>(static_cast<std::size_t>(n), 0.0));
  std::vector<double> gradient(static_cast<std::size_t>(n), 0.0);
  std::vector<double> delta(static_cast<std::size_t>(n), 0.0);
  std::vector<double> trial = theta;

  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    // 前向差分 Jacobian：每个参数一次正演。
    for (int j = 0; j < n; ++j) {
      trial = theta;
      trial[static_cast<std::size_t>(j)] += options.step;
      if (!evaluate(trial, r_perturbed)) {
        return result;
      }
      const double scale = 1.0 / options.step;
      for (std::size_t i = 0; i < m; ++i) {
        jacobian[i][static_cast<std::size_t>(j)] =
            (r_perturbed[i] - r[i]) * scale;
      }
    }

    // 组装正规方程：normal = J^T J，gradient = J^T r。
    for (int i = 0; i < n; ++i) {
      gradient[static_cast<std::size_t>(i)] = 0.0;
      for (int j = 0; j < n; ++j) {
        double sum = 0.0;
        for (std::size_t k = 0; k < m; ++k) {
          sum += jacobian[k][static_cast<std::size_t>(i)] *
                 jacobian[k][static_cast<std::size_t>(j)];
        }
        normal[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = sum;
      }
      for (std::size_t k = 0; k < m; ++k) {
        gradient[static_cast<std::size_t>(i)] +=
            jacobian[k][static_cast<std::size_t>(i)] * r[k];
      }
    }

    // 解 (J^T J + lambda * diag(J^T J)) delta = -J^T r；求解失败或代价不降
    // 就放大 lambda（更接近梯度下降）后重试。
    bool accepted = false;
    for (int backtrack = 0; backtrack < options.max_backtracks; ++backtrack) {
      std::vector<std::vector<double>> system = normal;
      for (int i = 0; i < n; ++i) {
        system[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] +=
            lambda * (normal[static_cast<std::size_t>(i)]
                            [static_cast<std::size_t>(i)] +
                      1e-12);
      }
      std::vector<double> rhs(static_cast<std::size_t>(n), 0.0);
      for (int i = 0; i < n; ++i) {
        rhs[static_cast<std::size_t>(i)] = -gradient[static_cast<std::size_t>(i)];
      }
      bool solved = true;
      try {
        delta = solve_dense(system, rhs);
      } catch (const std::exception&) {
        solved = false;
      }
      if (!solved) {
        lambda *= 5.0;
        continue;
      }
      bool finite = true;
      for (const double value : delta) {
        finite = finite && std::isfinite(value);
      }
      if (!finite) {
        lambda *= 5.0;
        continue;
      }

      trial = theta;
      for (int i = 0; i < n; ++i) {
        trial[static_cast<std::size_t>(i)] += delta[static_cast<std::size_t>(i)];
      }
      // 代价下降则接受新点并缩小阻尼；下降量小于 ftol 视为收敛。
      if (evaluate(trial, r_trial)) {
        const double trial_cost = cost_of(r_trial);
        if (trial_cost < cost) {
          const double change = cost - trial_cost;
          theta = trial;
          r = r_trial;
          cost = trial_cost;
          lambda = std::max(lambda / 3.0, 1e-12);
          result.theta = theta;
          result.residual = r;
          result.cost = cost;
          accepted = true;
          if (change < options.ftol * (1.0 + cost)) {
            result.iterations = iteration + 1;
            result.ok = true;
            return result;
          }
          break;
        }
      }
      lambda *= 5.0;
    }

    // 未接受任何步，或参数步长已足够小，则终止迭代。
    result.iterations = iteration + 1;
    double step_norm = 0.0;
    double theta_norm = 0.0;
    for (int i = 0; i < n; ++i) {
      step_norm = std::max(step_norm, std::abs(delta[static_cast<std::size_t>(i)]));
      theta_norm = std::max(theta_norm, std::abs(theta[static_cast<std::size_t>(i)]));
    }
    if (!accepted || step_norm < options.xtol * (theta_norm + options.xtol)) {
      break;
    }
  }

  result.ok = true;
  result.theta = theta;
  result.residual = r;
  result.cost = cost;
  return result;
}

}  // namespace dfn
