// 测试文件：带电解液的单颗粒模型 SPMe（对应 models/spme、solvers/dae 模块）。
// 验证状态规模、初始 OCV、解析 Jacobian、电荷/锂守恒与低倍率下与 DFN 一致。

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
  // x 网格 8/4/8（共 20 节点），两颗粒径向各 6 层。
  dfn::SpmeModel spme(params, 8, 4, 8, 6, 6);

  // 未知量 = 两颗粒浓度 6+6 + 两个表面浓度 2 + c_e 20 + phi_e 20 = 54。
  if (spme.size() != 54) {
    std::cerr << "size " << spme.size() << "\n";
    return fail("SPMe state size does not match 2N_r + 2 + 2N_x");
  }

  const std::vector<double> y0 = spme.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(spme.size()), 0.0);
  const double voltage0 = spme.terminal_voltage(y0.data());
  // 平衡初态端电压应等于 OCV；1 mV 容差覆盖 OCP 与初始一致化误差。
  if (std::abs(voltage0 - 3.98572) > 1e-3) {
    std::cerr << "V0 " << voltage0 << "\n";
    return fail("SPMe initial voltage does not match the OCV");
  }

  // j 由外加电流直接决定，残差对状态线性；解析 Jacobian 应与差分一致。
  // 扰动颗粒内部、表面与电解液浓度，确保各耦合列非零。
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
  const double cj = 1.0e3;  // 对应 h ~ 1 ms 的隐式时间步
  dfn::CscBuilder analytic(dfn::dense_pattern(n));
  spme.jacobian(0.0, y.data(), yp0.data(), cj, analytic);
  dfn::CscBuilder numeric(dfn::dense_pattern(n));
  dfn::finite_difference_jacobian(spme, 0.0, y.data(), yp0.data(), cj, numeric);

  // 线性问题的差分误差仅为舍入；1e-6 超差说明解析导数漏项或缩放不一致。
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
  // 零电流静止 100 s：无反应通量，端电压漂移应小于 1e-9 V。
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
  // 8 A/m^2 放电 600 s 应通过 4800 C/m^2 电荷，等于负极损失锂量 * F；
  // 1e-4 相对容差覆盖积分误差与浓度求和舍入。
  const double charge = dfn::kFaraday * (lithium0 - lithium1);
  if (std::abs(charge - 8.0 * 600.0) > 1e-4 * 8.0 * 600.0) {
    std::cerr << "charge " << charge << " C/m^2\n";
    return fail("SPMe charge passed does not match the applied current");
  }
  // 负极失去的锂全部转移到正极，1e-8 只允许求和舍入。
  if (std::abs((pos1 - pos0) + (lithium1 - lithium0)) >
      1e-8 * std::abs(lithium0)) {
    return fail("SPMe lithium is not transferred between electrodes");
  }
  // 含电解液的全电池锂总量还必须守恒，验证 c_e 方程的源项与化学计量一致。
  if (std::abs(total1 - total0) > 1e-8 * std::abs(total0)) {
    std::cerr << "Li drift " << (total1 - total0) / total0 << "\n";
    return fail("SPMe total lithium is not conserved");
  }

  // 低倍率 0.5 A/m^2：x 方向浓度梯度很小，SPMe 应逼近 DFN；
  // 60 s 后两者端电压差应小于 2 mV。
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
  // 2 mV 是低倍率下 SPMe 忽略反应沿 x 分布差异的合理误差上限；
  // 超出说明两模型参数或方程不一致，而非模型层级差异。
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
