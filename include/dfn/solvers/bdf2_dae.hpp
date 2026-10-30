#pragma once

#include <string>
#include <vector>

#include "dfn/solvers/dae.hpp"

namespace dfn {

struct Bdf2DaeOptions {
  double newton_tolerance = 1e-10;
  int newton_max_iterations = 20;
};

struct Bdf2DaeResult {
  std::vector<double> t;
  std::vector<std::vector<double>> y;
  bool ok = true;
  int steps = 0;
  int newton_iterations = 0;
  std::string message;
};

// Fixed-step BDF2 for F(t, y, y') = 0. The first step uses BDF1.
Bdf2DaeResult bdf2_dae_integrate(const DaeModel& model, std::vector<double> y0,
                                 double t0, double t_end, double h,
                                 const Bdf2DaeOptions& options = {});

}  // namespace dfn
