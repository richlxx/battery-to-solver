#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dfn/models/chain_dae.hpp"
#include "dfn/solvers/bdf2_dae.hpp"
#include "dfn/solvers/dae.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

bool close(double a, double b, double tolerance) {
  return std::abs(a - b) <= tolerance * (1.0 + std::abs(a) + std::abs(b));
}

int test_pattern() {
  const dfn::CscPattern tridiagonal = dfn::build_pattern(3, {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 1}, {1, 2}, {2, 2}});
  if (tridiagonal.nnz() != 7) {
    return fail("tridiagonal nnz");
  }
  if (tridiagonal.index_of(1, 0) < 0 || tridiagonal.index_of(0, 2) >= 0) {
    return fail("tridiagonal lookup");
  }

  dfn::CscBuilder builder(tridiagonal);
  bool threw = false;
  try {
    builder.set(0, 2, 1.0);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  if (!threw) {
    return fail("missing pattern entry must throw");
  }

  const dfn::CscPattern dense = dfn::dense_pattern(4);
  if (dense.nnz() != 16) {
    return fail("dense nnz");
  }
  return 0;
}

int test_chain_model() {
  dfn::ChainDaeOptions options;
  options.cells = 3;
  options.lambda = 50.0;
  options.diffusion = 4.0;
  const dfn::ChainDae model(options);

  if (model.size() != 6) {
    return fail("chain size");
  }
  if (model.jacobian_pattern().nnz() != 13) {
    return fail("chain pattern nnz");
  }

  double t = 0.3;
  std::vector<double> y = dfn::chain_dae_consistent_state(model, t);
  std::vector<double> yp = dfn::chain_dae_consistent_derivative(model, t);
  std::vector<double> residual(static_cast<std::size_t>(model.size()), 0.0);
  model.residual(t, y.data(), yp.data(), residual.data());
  for (double value : residual) {
    if (std::abs(value) > 1e-12) {
      return fail("consistent state does not zero the residual");
    }
  }

  // The same must hold for a non-uniform profile, where the diffusion term is
  // actually active.
  y = dfn::chain_dae_mode_state(model, t, 0.1);
  yp = dfn::chain_dae_derivative_from_state(model, t, y);
  model.residual(t, y.data(), yp.data(), residual.data());
  for (double value : residual) {
    if (std::abs(value) > 1e-12) {
      return fail("non-uniform consistent state does not zero the residual");
    }
  }

  // Analytic Jacobian must match the finite-difference one, including cj.
  y[0] = 0.4;
  y[2] = -0.2;
  yp[0] = 1.5;
  yp[4] = -0.7;
  const double cj = 120.0;

  dfn::CscBuilder analytic(dfn::dense_pattern(model.size()));
  model.jacobian(t, y.data(), yp.data(), cj, analytic);

  dfn::CscBuilder numeric(dfn::dense_pattern(model.size()));
  dfn::finite_difference_jacobian(model, t, y.data(), yp.data(), cj, numeric);

  for (std::size_t k = 0; k < analytic.values().size(); ++k) {
    if (!close(analytic.values()[k], numeric.values()[k], 1e-6)) {
      std::cerr << "analytic " << analytic.values()[k] << " vs numeric "
                << numeric.values()[k] << " at " << k << "\n";
      return fail("Jacobian mismatch");
    }
  }
  return 0;
}

int test_single_cell_solution() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);

  const std::vector<double> yp_at_start =
      dfn::chain_dae_consistent_derivative(model, 0.4);
  if (!close(yp_at_start[0],
             dfn::chain_dae_exact_concentration_rate(100.0, 0.4), 1e-12)) {
    return fail("consistent derivative disagrees with the closed form");
  }

  dfn::Bdf2DaeOptions bdf;
  bdf.newton_tolerance = 1e-12;
  bdf.newton_max_iterations = 30;
  const dfn::Bdf2DaeResult result =
      dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, 0.005, bdf);
  if (!result.ok) {
    return fail(result.message.c_str());
  }

  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);
  const double computed = result.y.back()[0];
  if (!close(computed, exact, 2e-4)) {
    std::cerr << "computed " << computed << " exact " << exact << "\n";
    return fail("single cell solution");
  }
  const double w = result.y.back()[1];
  if (std::abs(w - std::sin(t_end)) > 1e-10) {
    return fail("algebraic constraint at t_end");
  }

  // Second-order convergence in h.
  const auto error = [&](double h) {
    const dfn::Bdf2DaeResult run =
        dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, h, bdf);
    if (!run.ok) {
      throw std::runtime_error(run.message);
    }
    return std::abs(run.y.back()[0] - exact);
  };
  const double e1 = error(0.04);
  const double e2 = error(0.02);
  const double e3 = error(0.01);
  const double ratio1 = e1 / e2;
  const double ratio2 = e2 / e3;
  if (ratio1 < 3.0 || ratio1 > 5.5 || ratio2 < 3.0 || ratio2 > 5.5) {
    std::cerr << "ratios " << ratio1 << ", " << ratio2 << "\n";
    return fail("BDF2-DAE observed order");
  }

  // Without diffusion the cells are independent: every cell must match the
  // single-cell closed form.
  dfn::ChainDaeOptions chain_options;
  chain_options.cells = 4;
  chain_options.lambda = 100.0;
  chain_options.diffusion = 0.0;
  const dfn::ChainDae chain(chain_options);
  const dfn::Bdf2DaeResult chain_result = dfn::bdf2_dae_integrate(
      chain, dfn::chain_dae_consistent_state(chain, 0.0), 0.0, t_end, 0.005, bdf);
  if (!chain_result.ok) {
    return fail(chain_result.message.c_str());
  }
  for (int i = 0; i < chain.cells(); ++i) {
    if (!close(chain_result.y.back()[static_cast<std::size_t>(2 * i)], exact, 2e-4)) {
      return fail("decoupled chain cells disagree with closed form");
    }
  }

  std::cout << "dae ok (BDF2-DAE order ratios " << ratio1 << ", " << ratio2
            << ")\n";
  return 0;
}

}  // namespace

int main() {
  if (const int rc = test_pattern(); rc != 0) {
    return rc;
  }
  if (const int rc = test_chain_model(); rc != 0) {
    return rc;
  }
  return test_single_cell_solution();
}
