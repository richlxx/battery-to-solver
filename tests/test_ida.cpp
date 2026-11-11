// 测试文件：SUNDIALS IDA 封装与 DAE 初值/错误处理（对应 solvers/idaklu、
// models/chain_dae 模块）。覆盖 KLU/dense 后端、一致初值、缺失 Jacobian 项等。

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

// 混合绝对/相对容差，适配 O(1) 状态量与不同求解器的精度差异。
bool close(double a, double b, double tolerance) {
  return std::abs(a - b) <= tolerance * (1.0 + std::abs(a) + std::abs(b));
}

// 纯 ODE 模型：y' + y = cos(t)，所有变量都是微分变量（variable_types 默认 1）。
// 用于验证默认 IDACalcIC 路径以及纯 ODE 的积分精度。
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
  // dF/dy + cj*dF/dy' = 1 + cj，与 cj 成线性关系。
  void jacobian(double /*t*/, const double* /*y*/, const double* /*yp*/,
                double cj, dfn::CscBuilder& out) const override {
    out.add(0, 0, cj + 1.0);
  }
};

// 故意漏掉 (0,1) 非零元的模式，模拟装配阶段漏登记耦合项；
// 用 ChainDae 的残差/Jacobian 驱动，期望 IDA 封装给出明确错误。
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

// 单单元 DAE：比较 KLU 与 dense 后端、内部差分 Jacobian 与显式 Jacobian。
int test_single_cell() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  // 传入严格一致的 y0/yp0，并关闭 IDACalcIC，使后续差异只来自线性求解器。
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);
  const std::vector<double> yp0 = dfn::chain_dae_consistent_derivative(model, 0.0);

  dfn::IdaOptions ida;
  ida.rtol = 1e-9;    // 收紧容差，暴露稀疏与稠密求解的舍入差异
  ida.atol = {1e-11};
  ida.calculate_consistent_initial = false;

  const dfn::SolveReport klu =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
  if (!klu.ok) {
    std::cerr << klu.message << "\n";
    return fail("IDA+KLU integration failed");
  }

  // 与闭式解比较：rtol = 1e-9 下 1e-6 的容差足以覆盖 1 s 积分的累积误差。
  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);
  if (!close(klu.y[0], exact, 1e-6)) {
    std::cerr << "KLU " << klu.y[0] << " exact " << exact << "\n";
    return fail("IDA+KLU single-cell accuracy");
  }
  // 代数约束 w = sin(t) 由求解器强制满足，1e-8 对应非线性收敛精度。
  if (std::abs(klu.y[1] - std::sin(t_end)) > 1e-8) {
    return fail("IDA+KLU algebraic constraint");
  }
  // 1 s 内步数应在合理范围：0 步说明未积分，异常偏多则提示步长控制问题。
  if (klu.steps <= 0 || klu.steps > 500) {
    std::cerr << "steps " << klu.steps << "\n";
    return fail("IDA+KLU step count");
  }

  // 同一解析 Jacobian 换用稠密后端，解应在 1e-6 内一致。
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

  // 稠密矩阵允许 IDA 用内部差商近似 Jacobian；结果仍须与解析版本一致。
  ida.supply_jacobian = false;
  const dfn::SolveReport internal_jac =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
  if (!internal_jac.ok || !close(internal_jac.y[0], klu.y[0], 1e-6)) {
    std::cerr << internal_jac.message << "\n";
    return fail("IDA internal Jacobian on dense matrix");
  }

  // KLU 依赖固定的 CSC 结构，无法预知差分产生的新非零元，
  // 因此必须拒绝“无显式 Jacobian”的组合并返回明确提示。
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

// 交叉验证：同一 8 单元刚性链分别用 IDA（自适应步长）和
// 自研 BDF2-DAE（定步长 5e-4）求解，浓度轨迹的相对差应很小。
int test_chain_against_bdf2() {
  dfn::ChainDaeOptions options;
  options.cells = 8;
  options.lambda = 1e4;     // 强刚性，显式方法需要极小步长
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

  // BDF2 步长 5e-4、Newton 容差 1e-11，期望时间离散误差明显低于 1e-3 阈值。
  dfn::Bdf2DaeOptions bdf;
  bdf.newton_tolerance = 1e-11;
  bdf.newton_max_iterations = 30;
  const dfn::Bdf2DaeResult bdf_result =
      dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, 5e-4, bdf);
  if (!bdf_result.ok) {
    return fail(bdf_result.message.c_str());
  }

  // 只比较微分变量 c_i（偶数下标）；1e-3 的相对容差容纳两种方法的
  // 离散差异，同时能发现量级或符号错误。
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

// 错误路径：Jacobian 模式缺少真实耦合项时，封装层必须拒绝运行，
// 且错误信息要指出 pattern，便于用户定位装配问题。
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

// 纯 ODE 默认选项：验证默认 IDACalcIC 路径可运行，且终点精度达标。
int test_pure_ode_default_options() {
  const ScalarOdeModel model;
  dfn::IdaOptions ida;
  ida.rtol = 1e-9;
  ida.atol = {1e-11};
  // 初值 y0 = 1 时真实导数为 y'(0) = cos(0) - 1 = 0，与给定的 yp0 一致；
  // 本用例验证默认 IDACalcIC 路径可正常运行且不破坏精度。
  const dfn::SolveReport report =
      dfn::ida_integrate(model, {1.0}, {0.0}, 0.0, 1.0, ida);
  if (!report.ok) {
    std::cerr << report.message << "\n";
    return fail("pure ODE with IDACalcIC failed");
  }
  // 解析解：y(t) = (sin t + cos t)/2 + C e^{-t}，由 y(0)=1 得 C=1/2。
  const double exact = (std::sin(1.0) + std::cos(1.0) + std::exp(-1.0)) / 2.0;
  if (!close(report.y[0], exact, 1e-6)) {
    std::cerr << "computed " << report.y[0] << " exact " << exact << "\n";
    return fail("pure ODE accuracy");
  }
  return 0;
}

// IDACalcIC 的 DAE 版：初值同时包含错误的代数值 (w=0) 与零导数，
// 求解器应重建一致初值，并保持后续积分精度。
int test_consistent_initial_conditions() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t0 = 0.3;
  const double t_end = 0.6;
  std::vector<double> y0 = dfn::chain_dae_consistent_state(model, t0);
  y0[1] = 0.0;  // 故意给错的代数值，正确值应为 sin(t0)
  const std::vector<double> yp0(2, 0.0);  // 导数也全部置零，迫使 IDACalcIC 修正

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

  // 修正初值后，0.3 s 积分仍应与闭式解在 1e-6 内一致。
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
