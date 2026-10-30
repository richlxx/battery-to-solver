#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/solvers/dae.hpp"
#include "dfn/thermo.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

}  // namespace

int main() {
  dfn::Params params = dfn::make_demo_params();
  // A slower reaction keeps the finite-difference reference well conditioned;
  // the analytic derivatives are linear in the rate constant.
  params.neg.k = 1.0e-9;
  params.pos.k = 1.0e-9;
  const dfn::StateLayout layout = dfn::make_layout(2, 1, 2, 3, 3);
  dfn::DfnOptions options;
  options.applied_current_density = 8.0;
  dfn::DfnModel model(params, layout, options);

  const int n = model.size();
  std::vector<double> y = model.initial_state();
  std::vector<double> yp(static_cast<std::size_t>(n), 0.0);

  const double i_neg0 = model.electrode_current(y.data(), false);
  const double i_pos0 = model.electrode_current(y.data(), true);
  if (std::abs(i_neg0) > 1e-6 || std::abs(i_pos0) > 1e-6) {
    std::cerr << "I_neg=" << i_neg0 << " I_pos=" << i_pos0 << "\n";
    return fail("initial state is not at equilibrium");
  }

  // Move the state away from the initial guess so that every block is active.
  y[layout.cs_neg(0, 0)] *= 0.9;
  y[layout.cs_neg(1, 2)] *= 1.05;
  y[layout.cs_pos(0, 1)] *= 1.1;
  y[layout.ce(0)] *= 1.2;
  y[layout.ce(2)] *= 0.85;
  y[layout.ce(4)] *= 1.05;
  for (int i = 0; i < layout.n_x_total(); ++i) {
    y[layout.phie(i)] = -0.01 * static_cast<double>(i + 1);
  }
  for (int ix = 0; ix < layout.n_x_neg; ++ix) {
    y[model.cs_surf_neg(ix)] = y[layout.cs_neg(ix, layout.n_r_neg - 1)];
  }
  for (int ix = 0; ix < layout.n_x_pos; ++ix) {
    y[model.cs_surf_pos(ix)] = y[layout.cs_pos(ix, layout.n_r_pos - 1)];
  }
  // Put the potentials at local equilibrium so that the surface equations are
  // well conditioned; this keeps the finite-difference reference clean.
  for (int ix = 0; ix < layout.n_x_neg; ++ix) {
    const double stoich = y[model.cs_surf_neg(ix)] / params.neg.c_max;
    y[layout.phis_neg(ix)] =
        y[layout.phie(ix)] + dfn::negative_ocp(stoich).value + 0.005;
  }
  for (int ix = 0; ix < layout.n_x_pos; ++ix) {
    const int gx = layout.n_x_neg + layout.n_x_sep + ix;
    const double stoich = y[model.cs_surf_pos(ix)] / params.pos.c_max;
    y[layout.phis_pos(ix)] =
        y[layout.phie(gx)] + dfn::positive_ocp(stoich).value - 0.005;
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[static_cast<std::size_t>(i)])) {
      return fail("initial state is not finite");
    }
  }

  const double cj = 1.0e3;
  dfn::CscBuilder analytic(dfn::dense_pattern(n));
  model.jacobian(0.0, y.data(), yp.data(), cj, analytic);

  dfn::CscBuilder numeric(dfn::dense_pattern(n));
  dfn::finite_difference_jacobian(model, 0.0, y.data(), yp.data(), cj, numeric);

  double worst = 0.0;
  std::size_t worst_index = 0;
  for (std::size_t k = 0; k < analytic.values().size(); ++k) {
    const double a = analytic.values()[k];
    const double b = numeric.values()[k];
    const double diff = std::abs(a - b) / (1.0 + std::abs(a) + std::abs(b));
    if (diff > worst) {
      worst = diff;
      worst_index = k;
    }
  }
  if (worst > 1e-4) {
    const int col = static_cast<int>(worst_index) / n;
    const int row = static_cast<int>(worst_index) % n;
    std::cerr << "worst mismatch " << worst << " at (" << row << "," << col
              << ") analytic " << analytic.values()[worst_index] << " numeric "
              << numeric.values()[worst_index] << "\n";
    return fail("analytic Jacobian disagrees with finite differences");
  }

  const double total = model.total_lithium(y.data());
  const double solid = model.solid_lithium(y.data());
  const double electrolyte = model.electrolyte_lithium(y.data());
  if (std::abs(total - solid - electrolyte) > 1e-9 * std::abs(total)) {
    return fail("total lithium decomposition");
  }

  const double ineg = model.electrode_current(y.data(), false);
  const double ipos = model.electrode_current(y.data(), true);
  if (!std::isfinite(ineg) || !std::isfinite(ipos)) {
    return fail("electrode current is not finite");
  }

  std::cout << "dfn ok (n=" << n << ", worst Jacobian diff " << worst
            << ", Li=" << total << ", I_neg=" << ineg << ", I_pos=" << ipos
            << ")\n";
  return 0;
}
