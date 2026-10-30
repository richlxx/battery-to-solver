#pragma once

#include <functional>
#include <string>
#include <vector>

#include "dfn/solvers/dae.hpp"

namespace dfn {

enum class LinearSolverKind { klu, dense };

// Root function g(t, y, y'). IDA reports an event when g changes sign.
using EventFunction =
    std::function<double(double, const double*, const double*)>;

struct IdaOptions {
  double rtol = 1e-6;
  // Empty: 1e-8 for every component. Size 1: same absolute tolerance for all.
  // Size n: one tolerance per component.
  std::vector<double> atol;
  double max_step = 0.0;  // <= 0 means "let IDA decide"
  bool supply_jacobian = true;
  bool calculate_consistent_initial = true;
  LinearSolverKind linear_solver = LinearSolverKind::klu;
  std::vector<EventFunction> events;
  std::vector<double> output_times;  // sampled in addition to t0 and the end
};

struct SolveReport {
  bool ok = false;
  std::string message;
  int error_code = 0;
  int steps = 0;
  int residual_evals = 0;
  int jacobian_evals = 0;
  int linear_setups = 0;
  int nonlinear_iterations = 0;
  double t_final = 0.0;
  std::vector<double> y;
  std::vector<double> yp;
  int event_index = -1;   // index into IdaOptions::events, -1 if none
  double t_event = 0.0;
  std::vector<double> output_t;
  std::vector<std::vector<double>> output_y;
};

// Integrates F(t, y, y') = 0 with SUNDIALS IDA. Requires DFN_USE_KLU.
SolveReport ida_integrate(const DaeModel& model, const std::vector<double>& y0,
                          const std::vector<double>& yp0, double t0,
                          double t_end, const IdaOptions& options = {});

// Human-readable description of an IDA return flag.
std::string ida_error_message(int flag);

}  // namespace dfn
