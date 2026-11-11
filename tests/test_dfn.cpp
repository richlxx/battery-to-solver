// 测试文件：完整 DFN/P2D 模型的 Jacobian 与锂守恒记账（对应 models/dfn 模块）。
// 验证平衡态电流为零、解析 Jacobian 与有限差分一致、锂总量可分解且有限。

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
  // 把反应速率常数调慢到 1e-9，使有限差分参考值良态：
  // 反应项对 k 线性，导数尺度可控，同时避免强非线性放大差分误差。
  params.neg.k = 1.0e-9;
  params.pos.k = 1.0e-9;
  // 小网格：负极 2 个 x、1 个隔膜、正极 2 个 x，径向各 3 层，
  // 足以激活所有方程块，又让稠密有限差分 Jacobian 保持廉价。
  const dfn::StateLayout layout = dfn::make_layout(2, 1, 2, 3, 3);
  dfn::DfnOptions options;
  options.applied_current_density = 8.0;  // A/m^2，典型 1C 量级
  dfn::DfnModel model(params, layout, options);

  const int n = model.size();
  std::vector<double> y = model.initial_state();
  std::vector<double> yp(static_cast<std::size_t>(n), 0.0);

  // 平衡态下两侧反应电流都应为 0；1e-6 A 的绝对阈值用于排除舍入噪声。
  const double i_neg0 = model.electrode_current(y.data(), false);
  const double i_pos0 = model.electrode_current(y.data(), true);
  if (std::abs(i_neg0) > 1e-6 || std::abs(i_pos0) > 1e-6) {
    std::cerr << "I_neg=" << i_neg0 << " I_pos=" << i_pos0 << "\n";
    return fail("initial state is not at equilibrium");
  }

  // 人为扰动各状态块，使固相、电解液、电势与反应耦合全部进入非平凡状态；
  // 否则大量 Jacobian 元素为 0，无法发现漏项。
  y[layout.cs_neg(0, 0)] *= 0.9;
  y[layout.cs_neg(1, 2)] *= 1.05;
  y[layout.cs_pos(0, 1)] *= 1.1;
  y[layout.ce(0)] *= 1.2;
  y[layout.ce(2)] *= 0.85;
  y[layout.ce(4)] *= 1.05;
  // 沿 x 方向施加线性电解液电位，使电势方程与欧姆项都有非零贡献。
  for (int i = 0; i < layout.n_x_total(); ++i) {
    y[layout.phie(i)] = -0.01 * static_cast<double>(i + 1);
  }
  // 表面浓度辅助未知量先与最外层径向节点对齐，作为代数变量的初始猜测。
  for (int ix = 0; ix < layout.n_x_neg; ++ix) {
    y[model.cs_surf_neg(ix)] = y[layout.cs_neg(ix, layout.n_r_neg - 1)];
  }
  for (int ix = 0; ix < layout.n_x_pos; ++ix) {
    y[model.cs_surf_pos(ix)] = y[layout.cs_pos(ix, layout.n_r_pos - 1)];
  }
  // 把固相电位放在局部平衡附近（OCV 加减 5 mV 过电位），
  // 使 Butler-Volmer 表面方程良态，有限差分参考不被指数项放大。
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
  // 扰动后仍须全部为有限值，否则后续 Jacobian 比较没有意义。
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[static_cast<std::size_t>(i)])) {
      return fail("initial state is not finite");
    }
  }

  // cj = 1e3 对应 h ~ 1 ms，覆盖隐式时间步对 dF/dy' 的放大；
  // 解析 Jacobian 与稠密前向差分逐元素比较。
  const double cj = 1.0e3;
  dfn::CscBuilder analytic(dfn::dense_pattern(n));
  model.jacobian(0.0, y.data(), yp.data(), cj, analytic);

  dfn::CscBuilder numeric(dfn::dense_pattern(n));
  dfn::finite_difference_jacobian(model, 0.0, y.data(), yp.data(), cj, numeric);

  // 归一化误差 |a-b|/(1+|a|+|b|)：避免大元素掩盖小元素错误，
  // 也避免接近 0 的元素因绝对误差被误判。
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
  // 1e-4 容纳差分截断与 Butler-Volmer 指数项的抵消误差；
  // 超差通常意味着解析导数漏项、缩放因子不一致或符号错误。
  if (worst > 1e-4) {
    const int col = static_cast<int>(worst_index) / n;
    const int row = static_cast<int>(worst_index) % n;
    std::cerr << "worst mismatch " << worst << " at (" << row << "," << col
              << ") analytic " << analytic.values()[worst_index] << " numeric "
              << numeric.values()[worst_index] << "\n";
    return fail("analytic Jacobian disagrees with finite differences");
  }

  // 锂记账恒等式：总量 = 固相 + 电解液；1e-9 的相对容差仅允许求和舍入。
  const double total = model.total_lithium(y.data());
  const double solid = model.solid_lithium(y.data());
  const double electrolyte = model.electrolyte_lithium(y.data());
  if (std::abs(total - solid - electrolyte) > 1e-9 * std::abs(total)) {
    return fail("total lithium decomposition");
  }

  // 扰动状态下反应电流应为有限值；出现 NaN/Inf 说明动力学或化学计量数越界。
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
