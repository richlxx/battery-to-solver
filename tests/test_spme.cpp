#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/models/spme.hpp"
#include "dfn/solvers/dae.hpp"

#ifdef DFN_USE_KLU
#include "dfn/solvers/idaklu.hpp"
#endif

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

}  // namespace

int main() {
  dfn::Params params = dfn::make_demo_params();
  dfn::SpmeModel spme(params, 8, 4, 8, 6, 6);

  if (spme.size() != 54) {
    std::cerr << "size " << spme.size() << "\n";
    return fail("SPMe state size does not match 2N_r + 2 + 2N_x");
  }

  const std::vector<double> y0 = spme.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(spme.size()), 0.0);
  const double voltage0 = spme.terminal_voltage(y0.data());
  if (std::abs(voltage0 - 3.98572) > 1e-3) {
    std::cerr << "V0 " << voltage0 << "\n";
    return fail("SPMe initial voltage does not match the OCV");
  }

  // j is fixed by the applied current, so the residual is linear in the
  // state; the analytic Jacobian must agree with finite differences.
  spme.set_current_density(8.0);
  std::vector<double> y = y0;
  y[static_cast<std::size_t>(spme.cs_neg(0))] *= 0.97;
  y[static_cast<std::size_t>(spme.cs_pos(2))] *= 1.03;
  y[static_cast<std::size_t>(spme.cs_surf_neg())] =
      y[static_cast<std::size_t>(spme.cs_neg(5))];
  y[static_cast<std::size_t>(spme.cs_surf_pos())] =
      y[static_cast<std::size_t>(spme.cs_pos(5))];
  for (int i = 0; i < 20; ++i) {
    y[static_cast<std::size_t>(spme.ce(i))] *= 1.01;
  }

  const int n = spme.size();
  const double cj = 1.0e3;
  dfn::CscBuilder analytic(dfn::dense_pattern(n));
  spme.jacobian(0.0, y.data(), yp0.data(), cj, analytic);
  dfn::CscBuilder numeric(dfn::dense_pattern(n));
  dfn::finite_difference_jacobian(spme, 0.0, y.data(), yp0.data(), cj, numeric);

  double worst = 0.0;
  for (std::size_t k = 0; k < analytic.values().size(); ++k) {
    const double a = analytic.values()[k];
    const double b = numeric.values()[k];
    worst = std::max(worst, std::abs(a - b) / (1.0 + std::abs(a) + std::abs(b)));
  }
  if (worst > 1e-6) {
    std::cerr << "worst Jacobian diff " << worst << "\n";
    return fail("SPMe analytic Jacobian disagrees with finite differences");
  }

#ifdef DFN_USE_KLU
  dfn::SpmeModel resting(params, 8, 4, 8, 6, 6);
  dfn::IdaOptions rest_options;
  rest_options.rtol = 1e-8;
  rest_options.atol = resting.suggested_atol();
  rest_options.calculate_consistent_initial = true;
  const dfn::SolveReport rest =
      dfn::ida_integrate(resting, resting.initial_state(), yp0, 0.0, 100.0,
                         rest_options);
  if (!rest.ok) {
    std::cerr << rest.message << "\n";
    return fail("SPMe rest integration failed");
  }
  if (std::abs(resting.terminal_voltage(rest.y.data()) - voltage0) > 1e-9) {
    return fail("SPMe voltage drifted at zero current");
  }

  dfn::SpmeModel driven(params, 8, 4, 8, 6, 6);
  driven.set_current_density(8.0);
  const std::vector<double> yd = driven.initial_state();
  const double lithium0 = driven.particle_lithium(yd.data(), false);
  const double pos0 = driven.particle_lithium(yd.data(), true);
  const double total0 = driven.total_lithium(yd.data());
  dfn::IdaOptions options;
  options.rtol = 1e-7;
  options.atol = driven.suggested_atol();
  options.calculate_consistent_initial = true;
  const dfn::SolveReport run =
      dfn::ida_integrate(driven, yd, yp0, 0.0, 600.0, options);
  if (!run.ok) {
    std::cerr << run.message << "\n";
    return fail("SPMe discharge failed");
  }
  const double lithium1 = driven.particle_lithium(run.y.data(), false);
  const double pos1 = driven.particle_lithium(run.y.data(), true);
  const double total1 = driven.total_lithium(run.y.data());
  const double charge = dfn::kFaraday * (lithium0 - lithium1);
  if (std::abs(charge - 8.0 * 600.0) > 1e-4 * 8.0 * 600.0) {
    std::cerr << "charge " << charge << " C/m^2\n";
    return fail("SPMe charge passed does not match the applied current");
  }
  if (std::abs((pos1 - pos0) + (lithium1 - lithium0)) >
      1e-8 * std::abs(lithium0)) {
    return fail("SPMe lithium is not transferred between electrodes");
  }
  if (std::abs(total1 - total0) > 1e-8 * std::abs(total0)) {
    std::cerr << "Li drift " << (total1 - total0) / total0 << "\n";
    return fail("SPMe total lithium is not conserved");
  }

  // Low rate: SPMe and DFN must agree closely.
  dfn::DfnModel dfn(
      params, dfn::make_layout(8, 4, 8, 6, 6),
      dfn::DfnOptions{dfn::ControlMode::current, 0.5, 0.0});
  const std::vector<double> ydfn0 = dfn.initial_state();
  const std::vector<double> ypdfn(static_cast<std::size_t>(dfn.size()), 0.0);
  dfn::IdaOptions low;
  low.rtol = 1e-8;
  low.atol = dfn.suggested_atol();
  low.calculate_consistent_initial = true;
  const dfn::SolveReport dfn_run =
      dfn::ida_integrate(dfn, ydfn0, ypdfn, 0.0, 60.0, low);

  dfn::SpmeModel slow(params, 8, 4, 8, 6, 6);
  slow.set_current_density(0.5);
  const std::vector<double> yslow0 = slow.initial_state();
  const std::vector<double> ypslow(static_cast<std::size_t>(slow.size()), 0.0);
  dfn::IdaOptions low_spme = low;
  low_spme.atol = slow.suggested_atol();
  const dfn::SolveReport spme_run =
      dfn::ida_integrate(slow, yslow0, ypslow, 0.0, 60.0, low_spme);
  if (!dfn_run.ok || !spme_run.ok) {
    return fail("low-rate comparison failed to integrate");
  }
  const double difference =
      std::abs(slow.terminal_voltage(spme_run.y.data()) -
               dfn.terminal_voltage(dfn_run.y.data()));
  if (difference > 2e-3) {
    std::cerr << "SPMe-DFN difference " << difference << " V\n";
    return fail("SPMe deviates from DFN at low rate");
  }

  std::cout << "spme ok (n=" << n << ", V0 " << voltage0 << ", Jacobian diff "
            << worst << ", 8 A/m^2 charge " << charge
            << " C/m^2, low-rate |dV| " << 1000.0 * difference << " mV)\n";
#else
  std::cout << "spme ok (Jacobian only, no SUNDIALS)\n";
#endif
  return 0;
}
