#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "dfn/models/cell_model.hpp"
#include "dfn/models/dfn.hpp"
#include "dfn/models/spm.hpp"
#include "dfn/models/spme.hpp"
#include "dfn/solvers/idaklu.hpp"

namespace {

struct Curve {
  std::vector<double> t;
  std::vector<double> v;
  double t_end = 0.0;
  double ms = 0.0;
  double rmse = 0.0;
  double worst = 0.0;
  int steps = 0;
  int n = 0;
  int nnz = 0;
};

double interpolate(const Curve& curve, double t) {
  if (t <= curve.t.front()) {
    return curve.v.front();
  }
  if (t >= curve.t.back()) {
    return curve.v.back();
  }
  std::size_t k = 1;
  while (k < curve.t.size() && curve.t[k] < t) {
    ++k;
  }
  const double t0 = curve.t[k - 1];
  const double t1 = curve.t[k];
  const double w = (t - t0) / (t1 - t0);
  return curve.v[k - 1] * (1.0 - w) + curve.v[k] * w;
}

Curve run(dfn::CellModel& model, double current, double t_max) {
  model.set_current_density(current);
  const std::vector<double> y0 = model.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(model.size()), 0.0);

  dfn::IdaOptions options;
  options.rtol = 1e-6;
  options.atol = model.suggested_atol();
  options.calculate_consistent_initial = true;
  options.events.push_back(
      [&model](double, const double* y, const double*) {
        return model.terminal_voltage(y) - 3.5;
      });
  for (double t = 60.0; t < t_max; t += 60.0) {
    options.output_times.push_back(t);
  }

  const auto start = std::chrono::steady_clock::now();
  const dfn::SolveReport report =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_max, options);
  const double ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - start)
                        .count();
  Curve curve;
  if (!report.ok) {
    std::cerr << "integration failed: " << report.message << "\n";
    return curve;
  }
  for (std::size_t k = 0; k < report.output_t.size(); ++k) {
    curve.t.push_back(report.output_t[k]);
    curve.v.push_back(model.terminal_voltage(report.output_y[k].data()));
  }
  curve.t_end = report.output_t.back();
  curve.ms = ms;
  curve.steps = report.steps;
  curve.n = model.size();
  curve.nnz = model.jacobian_pattern().nnz();
  return curve;
}

void compare(Curve& curve, const Curve& reference) {
  if (curve.t.empty() || reference.t.empty()) {
    return;
  }
  const double limit = std::min(curve.t_end, reference.t_end);
  double sum = 0.0;
  double worst = 0.0;
  int count = 0;
  for (std::size_t k = 0; k < reference.t.size(); ++k) {
    const double t = reference.t[k];
    if (t > limit) {
      break;
    }
    const double diff = reference.v[k] - interpolate(curve, t);
    sum += diff * diff;
    worst = std::max(worst, std::abs(diff));
    ++count;
  }
  curve.rmse = std::sqrt(sum / static_cast<double>(count));
  curve.worst = worst;
}

struct Row {
  double current = 0.0;
  const char* name = "";
  Curve curve;
};

}  // namespace

int main() {
  const dfn::Params params = dfn::make_demo_params();
  const double i_1c = 36.0;
  const double t_max = 20000.0;

  std::cout << std::setprecision(6);
  std::cout << "SPM / SPMe / DFN, discharge to 3.5 V (DFN 8/4/8, nr 6/6; "
               "SPMe 8/4/8, nr 6/6; SPM nr 20/20)\n";
  std::cout << "I(A/m^2)\tC-rate\tmodel\tn\tnnz\tt_cut(s)\tRMSE(mV)"
               "\tmax(mV)\twall(ms)\tsteps\n";

  for (double current : {2.0, 8.0, 16.0, 36.0, 72.0}) {
    dfn::DfnModel dfn(
        params, dfn::make_layout(8, 4, 8, 6, 6),
        dfn::DfnOptions{dfn::ControlMode::current, current, 0.0});
    dfn::SpmeModel spme(params, 8, 4, 8, 6, 6);
    dfn::SpmModel spm(params, 20, 20);

    Row rows[] = {{current, "dfn", run(dfn, current, t_max)},
                  {current, "spme", run(spme, current, t_max)},
                  {current, "spm", run(spm, current, t_max)}};
    if (rows[0].curve.t.empty() || rows[1].curve.t.empty() ||
        rows[2].curve.t.empty()) {
      return 1;
    }

    for (Row& row : rows) {
      compare(row.curve, rows[0].curve);
      std::cout << row.current << '\t' << row.current / i_1c << '\t'
                << row.name << '\t' << row.curve.n << '\t' << row.curve.nnz
                << '\t' << row.curve.t_end << '\t' << 1000.0 * row.curve.rmse
                << '\t' << 1000.0 * row.curve.worst << '\t' << row.curve.ms
                << '\t' << row.curve.steps << '\n';
    }
  }
  return 0;
}
