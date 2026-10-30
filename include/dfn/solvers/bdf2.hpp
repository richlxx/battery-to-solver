#pragma once

#include <functional>
#include <vector>

namespace dfn {

using RhsFunction =
    std::function<std::vector<double>(double, const std::vector<double>&)>;
using RhsJacobian = std::function<std::vector<std::vector<double>>(
    double, const std::vector<double>&)>;

struct Bdf2Options {
  double newton_tolerance = 1e-10;
  int newton_max_iterations = 20;
};

struct Bdf2Result {
  std::vector<double> t;
  std::vector<std::vector<double>> y;
  bool ok = true;
  int steps = 0;
};

Bdf2Result bdf2_integrate(const RhsFunction& rhs,
                          const RhsJacobian& jacobian,
                          std::vector<double> y0,
                          double t0,
                          double t_end,
                          double h,
                          const Bdf2Options& options = {});

}  // namespace dfn
