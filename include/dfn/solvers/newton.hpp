#pragma once

#include <functional>
#include <vector>

namespace dfn {

struct NewtonOptions {
  int max_iterations = 50;
  double tolerance = 1e-10;
  double step_tolerance = 1e-14;
  bool use_line_search = true;
  std::vector<double> variable_scale;
};

struct NewtonResult {
  std::vector<double> x;
  int iterations = 0;
  bool converged = false;
  double residual_norm = 0.0;
  std::vector<double> residual_history;
  std::vector<double> step_history;
};

using ResidualFunction =
    std::function<std::vector<double>(const std::vector<double>&)>;
using JacobianFunction =
    std::function<std::vector<std::vector<double>>(const std::vector<double>&)>;

NewtonResult newton_solve(const ResidualFunction& residual,
                          const JacobianFunction& jacobian,
                          std::vector<double> x0,
                          const NewtonOptions& options = {});

}  // namespace dfn
