#include "dfn/solvers/bdf2_dae.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "dfn/solvers/newton.hpp"

namespace dfn {

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
    const bool use_bdf2 = have_previous_step;
    const double a0 = use_bdf2 ? 3.0 / (2.0 * dt) : 1.0 / dt;

    std::vector<double> guess = y_curr;
    if (use_bdf2) {
      for (int j = 0; j < n; ++j) {
        guess[static_cast<std::size_t>(j)] =
            2.0 * y_curr[static_cast<std::size_t>(j)] -
            y_prev[static_cast<std::size_t>(j)];
      }
    }

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
