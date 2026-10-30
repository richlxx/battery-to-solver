#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dfn/solvers/newton.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

bool close(double a, double b, double tolerance = 1e-8) {
  return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

}  // namespace

int main() {
  const dfn::ResidualFunction residual =
      [](const std::vector<double>& x) {
        return std::vector<double>{x[0] * x[0] + x[1] - 3.0,
                                   x[0] + x[1] * x[1] - 5.0};
      };
  const dfn::JacobianFunction jacobian =
      [](const std::vector<double>& x) {
        return std::vector<std::vector<double>>{{2.0 * x[0], 1.0},
                                                {1.0, 2.0 * x[1]}};
      };

  dfn::NewtonOptions options;
  options.tolerance = 1e-12;
  options.max_iterations = 30;
  const dfn::NewtonResult result =
      dfn::newton_solve(residual, jacobian, {0.5, 1.5}, options);

  if (!result.converged) {
    return fail("coupled system did not converge");
  }
  if (!close(result.x[0], 1.0) || !close(result.x[1], 2.0)) {
    return fail("coupled system solution");
  }
  if (result.residual_history.empty() ||
      result.residual_history.back() > 1e-12) {
    return fail("residual history");
  }
  if (result.iterations > 15) {
    return fail("Newton took too many iterations");
  }

  const dfn::ResidualFunction scaled_residual =
      [](const std::vector<double>& x) {
        return std::vector<double>{x[0] / 1.0e4 + x[1] - 2.0,
                                   x[0] / 1.0e4 - x[1]};
      };
  const dfn::JacobianFunction scaled_jacobian =
      [](const std::vector<double>&) {
        return std::vector<std::vector<double>>{{1.0e-4, 1.0},
                                                {1.0e-4, -1.0}};
      };
  dfn::NewtonOptions scaled_options;
  scaled_options.tolerance = 1e-12;
  scaled_options.variable_scale = {1.0e4, 1.0};
  const dfn::NewtonResult scaled_result = dfn::newton_solve(
      scaled_residual, scaled_jacobian, {0.5e4, 0.5}, scaled_options);
  if (!scaled_result.converged ||
      !close(scaled_result.x[0], 1.0e4, 1e-8) ||
      !close(scaled_result.x[1], 1.0, 1e-8)) {
    return fail("scaled system");
  }

  bool singular_threw = false;
  try {
    (void)dfn::newton_solve(
        [](const std::vector<double>& x) {
          return std::vector<double>{x[0] - 1.0};
        },
        [](const std::vector<double>&) {
          return std::vector<std::vector<double>>{{0.0}};
        },
        {0.0});
  } catch (const std::exception&) {
    singular_threw = true;
  }
  if (!singular_threw) {
    return fail("singular Jacobian must throw");
  }

  std::cout << "newton ok (" << result.iterations << " iterations)\n";
  return 0;
}
