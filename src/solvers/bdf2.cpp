#include "dfn/solvers/bdf2.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "dfn/solvers/newton.hpp"

namespace dfn {

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
    const bool use_bdf2 = have_previous_step;

    std::vector<double> guess = y_curr;
    if (use_bdf2) {
      for (int j = 0; j < n; ++j) {
        guess[static_cast<std::size_t>(j)] =
            2.0 * y_curr[static_cast<std::size_t>(j)] -
            y_prev[static_cast<std::size_t>(j)];
      }
    }

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

    NewtonOptions newton_options;
    newton_options.tolerance = options.newton_tolerance;
    newton_options.max_iterations = options.newton_max_iterations;
    const NewtonResult newton =
        newton_solve(residual, jac, guess, newton_options);

    if (!newton.converged) {
      result.ok = false;
      break;
    }

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
