#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/models/chain_dae.hpp"
#include "dfn/solvers/bdf2_dae.hpp"
#include "dfn/solvers/idaklu.hpp"

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

void tolerance_sweep() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);
  const std::vector<double> yp0 = dfn::chain_dae_consistent_derivative(model, 0.0);
  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);

  std::cout << "== single cell DAE: tolerance sweep (IDA + KLU) ==\n";
  std::cout << "rtol\tatol\tsteps\tres_evals\tjac_evals\terror\n";
  for (const double rtol : {1e-3, 1e-5, 1e-7, 1e-9}) {
    dfn::IdaOptions ida;
    ida.rtol = rtol;
    ida.atol = {1e-3 * rtol};
    ida.calculate_consistent_initial = false;
    const dfn::SolveReport report =
        dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
    if (!report.ok) {
      std::cout << rtol << "\t" << ida.atol[0] << "\tfailed: " << report.message
                << "\n";
      continue;
    }
    const double error = std::abs(report.y[0] - exact);
    std::cout << rtol << "\t" << ida.atol[0] << "\t" << report.steps << "\t"
              << report.residual_evals << "\t\t" << report.jacobian_evals << "\t"
              << error << "\n";
  }
  std::cout << "\n";
}

void bdf2_versus_ida() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);
  const std::vector<double> yp0 = dfn::chain_dae_consistent_derivative(model, 0.0);
  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);

  std::cout << "== hand-written BDF2 versus IDA (single cell) ==\n";
  std::cout << "solver\tparams\tsteps\tnewton\tjac\terror\tms\n";

  for (const double h : {0.01, 0.005, 0.0025}) {
    dfn::Bdf2DaeOptions bdf;
    bdf.newton_tolerance = 1e-12;
    bdf.newton_max_iterations = 30;
    const auto start = Clock::now();
    const dfn::Bdf2DaeResult result =
        dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, h, bdf);
    const double ms = elapsed_ms(start);
    if (!result.ok) {
      std::cout << "BDF2\th=" << h << "\tfailed: " << result.message << "\n";
      continue;
    }
    std::cout << "BDF2\th=" << h << "\t" << result.steps << "\t"
              << result.newton_iterations << "\t-\t"
              << std::abs(result.y.back()[0] - exact) << "\t" << ms << "\n";
  }

  for (const double rtol : {1e-5, 1e-7, 1e-9}) {
    dfn::IdaOptions ida;
    ida.rtol = rtol;
    ida.atol = {1e-3 * rtol};
    ida.calculate_consistent_initial = false;
    const auto start = Clock::now();
    const dfn::SolveReport report =
        dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
    const double ms = elapsed_ms(start);
    if (!report.ok) {
      std::cout << "IDA\trtol=" << rtol << "\tfailed: " << report.message << "\n";
      continue;
    }
    std::cout << "IDA+KLU\trtol=" << rtol << "\t" << report.steps << "\t"
              << report.nonlinear_iterations << "\t" << report.jacobian_evals
              << "\t" << std::abs(report.y[0] - exact) << "\t" << ms << "\n";
  }
  std::cout << "\n";
}

void chain_scale_up() {
  std::cout << "== chain DAE: KLU versus dense linear solver ==\n";
  std::cout << "cells\tn\tlin_solver\tsteps\tjac\tsetups\tms\n";

  for (const int cells : {80, 320, 1280}) {
    dfn::ChainDaeOptions options;
    options.cells = cells;
    options.lambda = 1e4;
    options.diffusion = 100.0;
    const dfn::ChainDae model(options);

    const double t_end = 0.02;
    const std::vector<double> y0 = dfn::chain_dae_mode_state(model, 0.0, 0.1);
    const std::vector<double> yp0 =
        dfn::chain_dae_derivative_from_state(model, 0.0, y0);

    for (const auto kind :
         {dfn::LinearSolverKind::klu, dfn::LinearSolverKind::dense}) {
      dfn::IdaOptions ida;
      ida.rtol = 1e-6;
      ida.atol = {1e-8};
      ida.calculate_consistent_initial = false;
      ida.linear_solver = kind;
      const auto start = Clock::now();
      const dfn::SolveReport report =
          dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
      const double ms = elapsed_ms(start);
      const char* name =
          kind == dfn::LinearSolverKind::klu ? "KLU" : "dense";
      if (!report.ok) {
        std::cout << cells << "\t" << model.size() << "\t" << name
                  << "\tfailed: " << report.message << "\n";
        continue;
      }
      std::cout << cells << "\t" << model.size() << "\t" << name << "\t\t"
                << report.steps << "\t" << report.jacobian_evals << "\t"
                << report.linear_setups << "\t" << ms << "\n";
    }
  }
  std::cout << "\n";
}

}  // namespace

int main() {
  std::cout << std::setprecision(4);
  tolerance_sweep();
  bdf2_versus_ida();
  chain_scale_up();
  return 0;
}
