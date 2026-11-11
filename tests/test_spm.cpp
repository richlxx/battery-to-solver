// 测试文件：单颗粒模型 SPM（对应 models/spm、solvers/dae、solvers/idaklu 模块）。
// 验证初始 OCV、解析 Jacobian、零电流静止、电荷守恒与低倍率下与 DFN 的一致性。

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/models/spm.hpp"
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
  dfn::SpmModel spm(params, 8, 8);  // 两电极各 8 个径向单元

  const std::vector<double> y0 = spm.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(spm.size()), 0.0);
  const double voltage0 = spm.terminal_voltage(y0.data());
  // 平衡初态端电压应等于 OCV（约 3.98572 V）；1 mV 容差覆盖 OCP 拟合误差。
  if (std::abs(voltage0 - 3.98572) > 1e-3) {
    std::cerr << "V0 " << voltage0 << "\n";
    return fail("SPM initial voltage does not match the OCV");
  }

  // 方程对状态是线性的，解析 Jacobian 应与前向差分吻合到舍入量级；
  // 扰动两个颗粒内部节点与表面值，确保所有耦合列都非零。
  spm.set_current_density(8.0);
  std::vector<double> y = y0;
  y[static_cast<std::size_t>(spm.cs_neg(0))] *= 0.97;
  y[static_cast<std::size_t>(spm.cs_pos(2))] *= 1.03;
  y[static_cast<std::size_t>(spm.cs_surf_neg())] =
      y[static_cast<std::size_t>(spm.cs_neg(7))];
  y[static_cast<std::size_t>(spm.cs_surf_pos())] =
      y[static_cast<std::size_t>(spm.cs_pos(7))];

  const int n = spm.size();
  const double cj = 1.0e3;  // 对应 h ~ 1 ms 的隐式时间步
  dfn::CscBuilder analytic(dfn::dense_pattern(n));
  spm.jacobian(0.0, y.data(), yp0.data(), cj, analytic);
  dfn::CscBuilder numeric(dfn::dense_pattern(n));
  dfn::finite_difference_jacobian(spm, 0.0, y.data(), yp0.data(), cj, numeric);

  // 线性方程组的差分误差应仅为舍入；1e-6 已很宽松，超差说明解析导数错误。
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
  // 零电流静止 100 s：没有锂通量，状态应保持不动，端电压漂移 < 1e-9 V。
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

  // 8 A/m^2 放电 600 s：通过的总电荷量（C/m^2）必须等于负极颗粒损失的
  // 锂量乘以 Faraday 常数，这是电荷守恒的直接检验。
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
  // 期望电荷 = 8 A/m^2 * 600 s = 4800 C/m^2；1e-4 的相对容差覆盖
  // IDA 积分误差与浓度求和的舍入。
  const double charge = dfn::kFaraday * (lithium0 - lithium1);
  if (std::abs(charge - 8.0 * 600.0) > 1e-4 * 8.0 * 600.0) {
    std::cerr << "charge " << charge << " C/m^2\n";
    return fail("SPM charge passed does not match the applied current");
  }
  // 负极失去的锂应全部转移到正极；1e-8 的相对容差只允许求和误差。
  if (std::abs((pos1 - pos0) + (lithium1 - lithium0)) >
      1e-8 * std::abs(lithium0)) {
    return fail("SPM lithium is not transferred between electrodes");
  }

  // 低倍率 0.5 A/m^2：电解液与固相浓度梯度很小，SPM 应逼近 DFN；
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
  // 2 mV 是低倍率下忽略电解液/固相 x 方向梯度的合理误差上限；
  // 超出说明 SPM 与 DFN 的参数或方程不一致，而非模型层级差异。
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
