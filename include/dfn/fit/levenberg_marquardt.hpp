#pragma once

#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <vector>

#include "dfn/linalg/dense.hpp"

namespace dfn {

struct LmOptions {
  int max_iterations = 50;
  int max_backtracks = 12;
  double initial_lambda = 1e-3;
  double step = 1e-6;       // forward-difference step in parameter space
  double xtol = 1e-9;
  double ftol = 1e-12;
};

struct LmResult {
  bool ok = false;
  int iterations = 0;
  int evaluations = 0;
  double cost = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> theta;
  std::vector<double> residual;
};

// Minimizes 0.5 * ||r(theta)||^2. The callback returns false when the
// forward model fails; a failed trial is treated like a rejected step.
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
    // Forward-difference Jacobian: one forward solve per parameter.
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
