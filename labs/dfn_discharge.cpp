#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/solvers/idaklu.hpp"

namespace {

using Clock = std::chrono::steady_clock;

dfn::SolveReport run_cc(dfn::DfnModel& model, const std::vector<double>& y0,
                        double t_end, double rtol, double* ms,
                        const std::vector<double>& output_times = {}) {
  const std::vector<double> yp0(static_cast<std::size_t>(model.size()), 0.0);
  dfn::IdaOptions options;
  options.rtol = rtol;
  options.atol = model.suggested_atol();
  options.calculate_consistent_initial = true;
  options.output_times = output_times;
  const auto start = Clock::now();
  dfn::SolveReport report =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, options);
  *ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  return report;
}

void grid_convergence(const dfn::Params& params) {
  struct LayoutSpec {
    int nx_neg;
    int nx_sep;
    int nx_pos;
    int nr_neg;
    int nr_pos;
  };

  const auto study = [&params](const char* name,
                               const std::vector<LayoutSpec>& specs) {
    std::cout << "== " << name << ": V at t = 60 s, 8 A/m^2 ==\n";
    std::cout << "n_x(n/s/p)\tn_r(n/p)\tn\tnnz\tsteps\tjac\tV(60)\tms\n";
    std::vector<double> voltages;
    for (const auto& spec : specs) {
      dfn::DfnOptions control;
      control.applied_current_density = 8.0;
      dfn::DfnModel model(
          params,
          dfn::make_layout(spec.nx_neg, spec.nx_sep, spec.nx_pos, spec.nr_neg,
                           spec.nr_pos),
          control);
      double ms = 0.0;
      const dfn::SolveReport report =
          run_cc(model, model.initial_state(), 60.0, 1e-9, &ms);
      if (!report.ok) {
        std::cout << "failed: " << report.message << "\n";
        return;
      }
      voltages.push_back(model.terminal_voltage(report.y.data()));
      std::cout << spec.nx_neg << "/" << spec.nx_sep << "/" << spec.nx_pos
                << "\t" << spec.nr_neg << "/" << spec.nr_pos << "\t"
                << model.size() << "\t" << model.jacobian_pattern().nnz() << "\t"
                << report.steps << "\t" << report.jacobian_evals << "\t"
                << model.terminal_voltage(report.y.data()) << "\t" << ms << "\n";
    }
    const double reference = voltages.back();
    std::cout << "differences vs finest:\n";
    for (std::size_t k = 0; k + 1 < voltages.size(); ++k) {
      std::cout << "  level " << k << ": " << std::abs(voltages[k] - reference)
                << "\n";
    }
    for (std::size_t k = 0; k + 2 < voltages.size(); ++k) {
      const double e1 = std::abs(voltages[k] - reference);
      const double e2 = std::abs(voltages[k + 1] - reference);
      if (e2 > 0.0) {
        std::cout << "  observed order " << k << "->" << k + 1 << ": "
                  << std::log2(e1 / e2) << "\n";
      }
    }
    std::cout << "\n";
  };

  study("x refinement (n_r fixed at 8)",
        {{4, 2, 4, 8, 8}, {8, 4, 8, 8, 8}, {16, 8, 16, 8, 8}, {32, 16, 32, 8, 8}});
  study("radial refinement (n_x fixed at 8/4/8)",
        {{8, 4, 8, 4, 4}, {8, 4, 8, 6, 6}, {8, 4, 8, 8, 8}, {8, 4, 8, 12, 12}});
  study("both refined",
        {{4, 2, 4, 4, 4}, {8, 4, 8, 6, 6}, {16, 8, 16, 8, 8}, {32, 16, 32, 12, 12}});
}

void cc_cv_curve(const dfn::Params& params) {
  dfn::DfnOptions control;
  control.applied_current_density = 8.0;
  dfn::DfnModel model(params, dfn::make_layout(8, 4, 8, 6, 6), control);
  const std::vector<double> y0 = model.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(model.size()), 0.0);
  const double cutoff = 3.5;

  dfn::IdaOptions cc;
  cc.rtol = 1e-6;
  cc.atol = model.suggested_atol();
  cc.calculate_consistent_initial = true;
  cc.events.push_back([&model, cutoff](double, const double* y, const double*) {
    return model.terminal_voltage(y) - cutoff;
  });
  for (int k = 0; k <= 8; ++k) {
    cc.output_times.push_back(2000.0 * static_cast<double>(k));
  }
  const dfn::SolveReport cc_report =
      dfn::ida_integrate(model, y0, yp0, 0.0, 20000.0, cc);
  if (!cc_report.ok) {
    std::cout << "CC failed: " << cc_report.message << "\n";
    return;
  }

  std::cout << "== CC discharge, cutoff " << cutoff << " V ==\n";
  std::cout << "t\tV\tsteps(cumulative)\n";
  for (std::size_t k = 0; k < cc_report.output_t.size(); ++k) {
    std::cout << cc_report.output_t[k] << "\t"
              << model.terminal_voltage(cc_report.output_y[k].data()) << "\n";
  }
  std::cout << "event at " << cc_report.t_event << " s, steps "
            << cc_report.steps << ", jac " << cc_report.jacobian_evals << "\n\n";

  model.set_applied_voltage(cutoff);
  dfn::IdaOptions cv;
  cv.rtol = 1e-6;
  cv.atol = model.suggested_atol();
  cv.calculate_consistent_initial = true;
  for (int k = 1; k <= 10; ++k) {
    cv.output_times.push_back(cc_report.t_event + 60.0 * static_cast<double>(k));
  }
  const dfn::SolveReport cv_report = dfn::ida_integrate(
      model, cc_report.y, cc_report.yp, cc_report.t_event,
      cc_report.t_event + 600.0, cv);
  if (!cv_report.ok) {
    std::cout << "CV failed: " << cv_report.message << "\n";
    return;
  }
  std::cout << "== CV hold at " << cutoff << " V ==\n";
  std::cout << "t since switch\tI (A/m^2)\n";
  for (std::size_t k = 0; k < cv_report.output_t.size(); ++k) {
    std::cout << cv_report.output_t[k] - cc_report.t_event << "\t"
              << model.terminal_current(cv_report.output_y[k].data()) << "\n";
  }
  std::cout << "CV steps " << cv_report.steps << "\n\n";
}

}  // namespace

int main() {
  std::cout << std::setprecision(6);
  const dfn::Params params = dfn::make_demo_params();
  grid_convergence(params);
  cc_cv_curve(params);
  return 0;
}
