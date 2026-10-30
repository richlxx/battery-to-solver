#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/fit/levenberg_marquardt.hpp"
#include "dfn/fit/spm_batch.hpp"
#include "dfn/reduction/spm.hpp"

#ifdef DFN_USE_KLU
#include "dfn/solvers/idaklu.hpp"
#endif

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

dfn::Params fit_params() {
  dfn::Params params = dfn::make_demo_params();
  params.neg.k = 3.0e-11;
  params.pos.k = 5.0e-11;
  return params;
}

std::vector<dfn::fit::BatchSegment> protocol_c() {
  return {{600.0, 0.0},   {1200.0, 18.0}, {1800.0, 0.0},
          {1200.0, -18.0}, {1800.0, 0.0},  {600.0, 36.0},
          {1800.0, 0.0},  {300.0, 72.0},   {1800.0, 0.0}};
}

}  // namespace

int main() {
  const dfn::Params params = fit_params();
  const int nr = 20;

  // The fixed-step SPM must approach the IDA solution as dt is reduced.
  const std::vector<dfn::fit::BatchRow> discharge =
      dfn::fit::rows_from_protocol(params, {{600.0, 18.0}}, nr, 5.0);
  const std::vector<dfn::fit::BatchRow> coarse =
      dfn::fit::rows_from_protocol(params, {{600.0, 18.0}}, nr, 10.0);

#ifdef DFN_USE_KLU
  dfn::SpmModel spm(params, nr, nr);
  spm.set_current_density(18.0);
  const std::vector<double> y0 = spm.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(spm.size()), 0.0);
  dfn::IdaOptions options;
  options.rtol = 1e-8;
  options.atol = spm.suggested_atol();
  options.calculate_consistent_initial = true;
  const dfn::SolveReport reference =
      dfn::ida_integrate(spm, y0, yp0, 0.0, 600.0, options);
  if (!reference.ok) {
    std::cerr << reference.message << "\n";
    return fail("IDA reference run failed");
  }
  const double reference_voltage = spm.terminal_voltage(reference.y.data());
  const double fine_error =
      std::abs(discharge.back().voltage - reference_voltage);
  const double coarse_error =
      std::abs(coarse.back().voltage - reference_voltage);
  if (fine_error > 2e-4 || coarse_error > 2e-3) {
    std::cerr << "fixed-step error " << fine_error << " / " << coarse_error
              << " V\n";
    return fail("fixed-step SPM disagrees with IDA");
  }

  // Charge must be the mirror image in the sign of the current.
  dfn::SpmModel charge_model(params, nr, nr);
  charge_model.set_current_density(-18.0);
  const dfn::SolveReport charge_ref = dfn::ida_integrate(
      charge_model, y0, yp0, 0.0, 600.0, options);
  if (!charge_ref.ok) {
    return fail("IDA charge reference failed");
  }
  const std::vector<dfn::fit::BatchRow> charge_rows =
      dfn::fit::rows_from_protocol(params, {{600.0, -18.0}}, nr, 5.0);
  const double charge_error =
      std::abs(charge_rows.back().voltage -
               charge_model.terminal_voltage(charge_ref.y.data()));
  if (charge_error > 2e-4) {
    std::cerr << "charge error " << charge_error << " V\n";
    return fail("fixed-step charge disagrees with IDA");
  }
#endif

  // Noise-free parameter recovery: the LM must return the truth.
  const std::vector<dfn::fit::BatchRow> rows =
      dfn::fit::rows_from_protocol(params, protocol_c(), nr, 5.0);
  const std::vector<double> truth = {0.0, std::log(params.neg.d_s),
                                     std::log(params.pos.d_s)};
  const std::vector<double> reproduced =
      dfn::fit::simulate_rows(params, rows, nr, 5.0);
  double reproduce_error = 0.0;
  for (std::size_t i = 0; i < rows.size(); ++i) {
    reproduce_error =
        std::max(reproduce_error, std::abs(reproduced[i] - rows[i].voltage));
  }
  if (reproduce_error > 1e-10) {
    std::cerr << "reproduction error " << reproduce_error << " V\n";
    return fail("simulate_rows does not reproduce the generated rows");
  }
  auto apply = [&](const std::vector<double>& theta) {
    dfn::Params p = params;
    const double scale = std::exp(theta[0]);
    p.neg.k *= scale;
    p.pos.k *= scale;
    p.neg.d_s = std::exp(theta[1]);
    p.pos.d_s = std::exp(theta[2]);
    return p;
  };
  const auto residual = [&](const std::vector<double>& theta,
                            std::vector<double>& values) {
    const std::vector<double> predicted =
        dfn::fit::simulate_rows(apply(theta), rows, nr, 5.0);
    values.resize(predicted.size());
    for (std::size_t i = 0; i < predicted.size(); ++i) {
      values[i] = predicted[i] - rows[i].voltage;
    }
    return true;
  };

  std::vector<double> start = truth;
  start[0] += 0.3;
  start[1] += 0.25;
  start[2] -= 0.25;
  dfn::LmOptions lm;
  lm.max_iterations = 40;
  const dfn::LmResult fit = dfn::levenberg_marquardt(start, residual, lm);
  if (!fit.ok) {
    return fail("LM did not converge");
  }
  double worst = 0.0;
  for (std::size_t j = 0; j < truth.size(); ++j) {
    worst = std::max(worst, std::abs(fit.theta[j] - truth[j]));
  }
  if (worst > 0.02 || fit.cost > 1e-16) {
    std::cerr << "worst parameter error " << worst << ", cost " << fit.cost
              << "\n";
    return fail("LM did not recover the noise-free parameters");
  }

  std::cout << "fit ok (samples=" << rows.size() << ", parameter error "
            << worst << ", cost " << fit.cost << ")\n";
  return 0;
}
