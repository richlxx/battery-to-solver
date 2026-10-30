#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "dfn/models/chain_dae.hpp"
#include "dfn/solvers/bdf2_dae.hpp"
#include "dfn/solvers/idaklu.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

bool close(double a, double b, double tolerance) {
  return std::abs(a - b) <= tolerance * (1.0 + std::abs(a) + std::abs(b));
}

// Pure ODE: y' + y = cos(t), so every variable is differential.
class ScalarOdeModel final : public dfn::DaeModel {
 public:
  int size() const override { return 1; }
  void residual(double t, const double* y, const double* yp,
                double* out) const override {
    out[0] = yp[0] + y[0] - std::cos(t);
  }
  dfn::CscPattern jacobian_pattern() const override {
    return dfn::dense_pattern(1);
  }
  void jacobian(double /*t*/, const double* /*y*/, const double* /*yp*/,
                double cj, dfn::CscBuilder& out) const override {
    out.add(0, 0, cj + 1.0);
  }
};

// Pattern deliberately missing the (0,1) entry, as if assembly forgot it.
class BrokenPatternModel final : public dfn::DaeModel {
 public:
  explicit BrokenPatternModel(dfn::ChainDaeOptions options) : inner_(options) {}

  int size() const override { return inner_.size(); }
  void residual(double t, const double* y, const double* yp,
                double* out) const override {
    inner_.residual(t, y, yp, out);
  }
  dfn::CscPattern jacobian_pattern() const override {
    return dfn::build_pattern(2, {{0, 0}, {1, 1}});
  }
  void jacobian(double t, const double* y, const double* yp, double cj,
                dfn::CscBuilder& out) const override {
    inner_.jacobian(t, y, yp, cj, out);
  }

 private:
  dfn::ChainDae inner_;
};

int test_single_cell() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);
  const std::vector<double> yp0 = dfn::chain_dae_consistent_derivative(model, 0.0);

  dfn::IdaOptions ida;
  ida.rtol = 1e-9;
  ida.atol = {1e-11};
  ida.calculate_consistent_initial = false;

  const dfn::SolveReport klu =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
  if (!klu.ok) {
    std::cerr << klu.message << "\n";
    return fail("IDA+KLU integration failed");
  }

  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);
  if (!close(klu.y[0], exact, 1e-6)) {
    std::cerr << "KLU " << klu.y[0] << " exact " << exact << "\n";
    return fail("IDA+KLU single-cell accuracy");
  }
  if (std::abs(klu.y[1] - std::sin(t_end)) > 1e-8) {
    return fail("IDA+KLU algebraic constraint");
  }
  if (klu.steps <= 0 || klu.steps > 500) {
    std::cerr << "steps " << klu.steps << "\n";
    return fail("IDA+KLU step count");
  }

  ida.linear_solver = dfn::LinearSolverKind::dense;
  const dfn::SolveReport dense =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
  if (!dense.ok) {
    std::cerr << dense.message << "\n";
    return fail("IDA+dense integration failed");
  }
  if (!close(dense.y[0], klu.y[0], 1e-6)) {
    return fail("KLU and dense linear solvers disagree");
  }

  // Dense matrices can use IDA's internal difference-quotient Jacobian.
  ida.supply_jacobian = false;
  const dfn::SolveReport internal_jac =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
  if (!internal_jac.ok || !close(internal_jac.y[0], klu.y[0], 1e-6)) {
    std::cerr << internal_jac.message << "\n";
    return fail("IDA internal Jacobian on dense matrix");
  }

  // KLU cannot use the internal Jacobian because new nonzeros are unknown.
  ida.linear_solver = dfn::LinearSolverKind::klu;
  const dfn::SolveReport klu_no_jac =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
  if (klu_no_jac.ok || klu_no_jac.message.find("explicit Jacobian") == std::string::npos) {
    return fail("KLU without a Jacobian must be rejected");
  }

  std::cout << "ida single cell: klu steps " << klu.steps << ", jac "
            << klu.jacobian_evals << ", dense steps " << dense.steps << "\n";
  return 0;
}

int test_chain_against_bdf2() {
  dfn::ChainDaeOptions options;
  options.cells = 8;
  options.lambda = 1e4;
  options.diffusion = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 0.2;
  const std::vector<double> y0 = dfn::chain_dae_mode_state(model, 0.0, 0.1);
  const std::vector<double> yp0 =
      dfn::chain_dae_derivative_from_state(model, 0.0, y0);

  dfn::IdaOptions ida;
  ida.rtol = 1e-8;
  ida.atol = {1e-10};
  ida.calculate_consistent_initial = false;
  const dfn::SolveReport report =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
  if (!report.ok) {
    std::cerr << report.message << "\n";
    return fail("IDA+KLU chain integration failed");
  }

  dfn::Bdf2DaeOptions bdf;
  bdf.newton_tolerance = 1e-11;
  bdf.newton_max_iterations = 30;
  const dfn::Bdf2DaeResult bdf_result =
      dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, 5e-4, bdf);
  if (!bdf_result.ok) {
    return fail(bdf_result.message.c_str());
  }

  double worst = 0.0;
  for (int i = 0; i < model.cells(); ++i) {
    const double a = report.y[static_cast<std::size_t>(2 * i)];
    const double b = bdf_result.y.back()[static_cast<std::size_t>(2 * i)];
    worst = std::max(worst, std::abs(a - b) / (1.0 + std::abs(b)));
  }
  if (worst > 1e-3) {
    std::cerr << "worst relative difference " << worst << "\n";
    return fail("IDA and BDF2-DAE disagree on the chain");
  }

  std::cout << "ida chain: steps " << report.steps << ", res " << report.residual_evals
            << ", jac " << report.jacobian_evals << ", setups "
            << report.linear_setups << ", bdf2 steps " << bdf_result.steps
            << ", worst diff " << worst << "\n";
  return 0;
}

int test_missing_jacobian_entry() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const BrokenPatternModel model(options);

  const std::vector<double> y0(2, 0.0);
  const std::vector<double> yp0 = {0.0, 1.0};
  dfn::IdaOptions ida;
  ida.calculate_consistent_initial = false;
  const dfn::SolveReport report = dfn::ida_integrate(model, y0, yp0, 0.0, 1.0, ida);
  if (report.ok) {
    return fail("broken Jacobian pattern must fail");
  }
  if (report.message.find("pattern") == std::string::npos) {
    std::cerr << report.message << "\n";
    return fail("failure message must mention the pattern");
  }
  return 0;
}

int test_pure_ode_default_options() {
  const ScalarOdeModel model;
  dfn::IdaOptions ida;
  ida.rtol = 1e-9;
  ida.atol = {1e-11};
  // yp0 is deliberately wrong; IDACalcIC must fix it.
  const dfn::SolveReport report =
      dfn::ida_integrate(model, {1.0}, {0.0}, 0.0, 1.0, ida);
  if (!report.ok) {
    std::cerr << report.message << "\n";
    return fail("pure ODE with IDACalcIC failed");
  }
  const double exact = (std::sin(1.0) + std::cos(1.0) + std::exp(-1.0)) / 2.0;
  if (!close(report.y[0], exact, 1e-6)) {
    std::cerr << "computed " << report.y[0] << " exact " << exact << "\n";
    return fail("pure ODE accuracy");
  }
  return 0;
}

int test_consistent_initial_conditions() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t0 = 0.3;
  const double t_end = 0.6;
  std::vector<double> y0 = dfn::chain_dae_consistent_state(model, t0);
  y0[1] = 0.0;  // wrong algebraic value
  const std::vector<double> yp0(2, 0.0);  // wrong derivative as well

  dfn::IdaOptions ida;
  ida.rtol = 1e-9;
  ida.atol = {1e-11};
  ida.calculate_consistent_initial = true;
  const dfn::SolveReport report =
      dfn::ida_integrate(model, y0, yp0, t0, t_end, ida);
  if (!report.ok) {
    std::cerr << report.message << "\n";
    return fail("IDACalcIC integration failed");
  }

  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);
  if (!close(report.y[0], exact, 1e-6)) {
    std::cerr << "computed " << report.y[0] << " exact " << exact << "\n";
    return fail("IDACalcIC result accuracy");
  }

  std::cout << "idacalcic recovered from inconsistent y0/yp0\n";
  return 0;
}

}  // namespace

int main() {
  if (const int rc = test_single_cell(); rc != 0) {
    return rc;
  }
  if (const int rc = test_chain_against_bdf2(); rc != 0) {
    return rc;
  }
  if (const int rc = test_missing_jacobian_entry(); rc != 0) {
    return rc;
  }
  if (const int rc = test_pure_ode_default_options(); rc != 0) {
    return rc;
  }
  return test_consistent_initial_conditions();
}
