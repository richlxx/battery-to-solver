#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/reduction/spm.hpp"
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
  dfn::SpmModel spm(params, 8, 8);

  const std::vector<double> y0 = spm.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(spm.size()), 0.0);
  const double voltage0 = spm.terminal_voltage(y0.data());
  if (std::abs(voltage0 - 3.98572) > 1e-3) {
    std::cerr << "V0 " << voltage0 << "\n";
    return fail("SPM initial voltage does not match the OCV");
  }

  // The equations are linear in the state, so the analytic Jacobian must
  // agree with finite differences to round-off.
  spm.set_current_density(8.0);
  std::vector<double> y = y0;
  y[static_cast<std::size_t>(spm.cs_neg(0))] *= 0.97;
  y[static_cast<std::size_t>(spm.cs_pos(2))] *= 1.03;
  y[static_cast<std::size_t>(spm.cs_surf_neg())] =
      y[static_cast<std::size_t>(spm.cs_neg(7))];
  y[static_cast<std::size_t>(spm.cs_surf_pos())] =
      y[static_cast<std::size_t>(spm.cs_pos(7))];

  const int n = spm.size();
  const double cj = 1.0e3;
  dfn::CscBuilder analytic(dfn::dense_pattern(n));
  spm.jacobian(0.0, y.data(), yp0.data(), cj, analytic);
  dfn::CscBuilder numeric(dfn::dense_pattern(n));
  dfn::finite_difference_jacobian(spm, 0.0, y.data(), yp0.data(), cj, numeric);

  double worst = 0.0;
  for (std::size_t k = 0; k < analytic.values().size(); ++k) {
    const double a = analytic.values()[k];
    const double b = numeric.values()[k];
    worst = std::max(worst, std::abs(a - b) / (1.0 + std::abs(a) + std::abs(b)));
  }
  if (worst > 1e-6) {
    std::cerr << "worst Jacobian diff " << worst << "\n";
    return fail("SPM analytic Jacobian disagrees with finite differences");
  }

#ifdef DFN_USE_KLU
  // Zero current: the state must stay put and the voltage must stay at OCV.
  dfn::SpmModel resting(params, 8, 8);
  dfn::IdaOptions rest_options;
  rest_options.rtol = 1e-8;
  rest_options.atol = resting.suggested_atol();
  rest_options.calculate_consistent_initial = true;
  const dfn::SolveReport rest =
      dfn::ida_integrate(resting, resting.initial_state(), yp0, 0.0, 100.0,
                         rest_options);
  if (!rest.ok) {
    std::cerr << rest.message << "\n";
    return fail("SPM rest integration failed");
  }
  if (std::abs(resting.terminal_voltage(rest.y.data()) - voltage0) > 1e-9) {
    return fail("SPM voltage drifted at zero current");
  }

  // 8 A/m^2 for 600 s: the charge passed must equal the lithium removed
  // from the negative particle.
  dfn::SpmModel driven(params, 8, 8);
  driven.set_current_density(8.0);
  const std::vector<double> yd = driven.initial_state();
  const double lithium0 = driven.particle_lithium(yd.data(), false);
  const double pos0 = driven.particle_lithium(yd.data(), true);
  dfn::IdaOptions options;
  options.rtol = 1e-7;
  options.atol = driven.suggested_atol();
  options.calculate_consistent_initial = true;
  const dfn::SolveReport run =
      dfn::ida_integrate(driven, yd, yp0, 0.0, 600.0, options);
  if (!run.ok) {
    std::cerr << run.message << "\n";
    return fail("SPM discharge failed");
  }
  const double lithium1 = driven.particle_lithium(run.y.data(), false);
  const double pos1 = driven.particle_lithium(run.y.data(), true);
  const double charge = dfn::kFaraday * (lithium0 - lithium1);
  if (std::abs(charge - 8.0 * 600.0) > 1e-4 * 8.0 * 600.0) {
    std::cerr << "charge " << charge << " C/m^2\n";
    return fail("SPM charge passed does not match the applied current");
  }
  if (std::abs((pos1 - pos0) + (lithium1 - lithium0)) >
      1e-8 * std::abs(lithium0)) {
    return fail("SPM lithium is not transferred between electrodes");
  }

  // Low rate: SPM and DFN must agree closely.
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

  dfn::SpmModel slow(params, 20, 20);
  slow.set_current_density(0.5);
  const std::vector<double> yslow0 = slow.initial_state();
  const std::vector<double> ypslow(static_cast<std::size_t>(slow.size()), 0.0);
  dfn::IdaOptions low_spm = low;
  low_spm.atol = slow.suggested_atol();
  const dfn::SolveReport spm_run =
      dfn::ida_integrate(slow, yslow0, ypslow, 0.0, 60.0, low_spm);
  if (!dfn_run.ok || !spm_run.ok) {
    return fail("low-rate comparison failed to integrate");
  }
  const double difference =
      std::abs(slow.terminal_voltage(spm_run.y.data()) -
               dfn.terminal_voltage(dfn_run.y.data()));
  if (difference > 2e-3) {
    std::cerr << "SPM-DFN difference " << difference << " V\n";
    return fail("SPM deviates from DFN at low rate");
  }

  std::cout << "spm ok (n=" << n << ", V0 " << voltage0 << ", Jacobian diff "
            << worst << ", 8 A/m^2 charge " << charge
            << " C/m^2, low-rate |dV| " << 1000.0 * difference << " mV)\n";
#else
  std::cout << "spm ok (Jacobian only, no SUNDIALS)\n";
#endif
  return 0;
}
