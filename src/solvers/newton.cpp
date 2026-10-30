#include "dfn/solvers/newton.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "dfn/linalg/dense.hpp"

namespace dfn {

namespace {

double max_abs(const std::vector<double>& values) {
  double norm = 0.0;
  for (double value : values) {
    norm = std::max(norm, std::abs(value));
  }
  return norm;
}

}  // namespace

NewtonResult newton_solve(const ResidualFunction& residual,
                          const JacobianFunction& jacobian,
                          std::vector<double> x0,
                          const NewtonOptions& options) {
  const int n = static_cast<int>(x0.size());
  if (n == 0) {
    throw std::invalid_argument("newton_solve: empty initial guess");
  }

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

  for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
    const std::vector<double> residual_values = residual(x);
    if (static_cast<int>(residual_values.size()) != n) {
      throw std::invalid_argument("newton_solve: residual size mismatch");
    }
    const double residual_norm = max_abs(residual_values);
    result.residual_history.push_back(residual_norm);

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

    std::vector<std::vector<double>> scaled_J = J;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        scaled_J[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *=
            scale[static_cast<std::size_t>(j)];
      }
    }

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

    if (step_norm < options.step_tolerance) {
      const double final_norm = max_abs(residual(x));
      result.x = x;
      result.iterations = iteration + 1;
      result.converged = final_norm < options.tolerance;
      result.residual_norm = final_norm;
      return result;
    }
  }

  result.x = x;
  result.iterations = options.max_iterations;
  result.converged = false;
  result.residual_norm = max_abs(residual(x));
  return result;
}

}  // namespace dfn
