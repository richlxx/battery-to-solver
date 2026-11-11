// Butler-Volmer 反应动力学实现：由局部过电位求界面摩尔通量及其偏导，
// 并提供“给定通量反解过电位”的牛顿迭代，供 SPM/SPMe 等降阶模型复用。
// 符号约定：j > 0 表示锂离开固相颗粒（脱锂，阳极方向）。
#include "dfn/models/reaction.hpp"

#include <algorithm>
#include <cmath>

#include "dfn/thermo.hpp"

namespace dfn {

namespace {

// 交换电流密度 i0 = F k (c_e)^alpha_a (c_max - c_s)^alpha_a (c_s)^alpha_c，
// 单位 A/m^2。固相表面浓度先夹紧到物理区间，避免幂运算出现零或负底数。
double exchange_current(const Params& params, const ElectrodeParams& electrode,
                        double ce, double cs_surf) {
  const double cmax = electrode.c_max;
  // 上、下界各留 1e-6 c_max 的余量，防止 OCV 与幂次在端点发散。
  const double cs = std::clamp(cs_surf, 1e-6 * cmax, cmax * (1.0 - 1e-6));
  // 电解液浓度下限取参考浓度的相对值，保证正底数。
  const double ce_clamped = std::max(ce, 1e-6 * params.elyte.c0);
  return kFaraday * electrode.k * std::pow(ce_clamped, electrode.alpha_a) *
         std::pow(cmax - cs, electrode.alpha_a) *
         std::pow(cs, electrode.alpha_c);
}

}  // namespace

// 计算单个界面的 Butler-Volmer 反应量。
// 输入：ce 电解液浓度 [mol/m^3]、cs_surf 固相表面浓度 [mol/m^3]、
//       phis/phie 固相与电解液电位 [V]。
// 返回 j [mol/(m^2 s)] 以及对 eta、ce、cs 的偏导，供雅可比装配使用。
ReactionValue reaction_at(const Params& params, const ElectrodeParams& electrode,
                          bool positive, double ce, double cs_surf, double phis,
                          double phie) {
  const double cmax = electrode.c_max;
  const double cs = std::clamp(cs_surf, 1e-6 * cmax, cmax * (1.0 - 1e-6));
  const double ce_clamped = std::max(ce, 1e-6 * params.elyte.c0);
  // 化学计量比用 c_s / c_max 表示，作为开路电位 OCV 的自变量。
  const double stoich = cs / cmax;
  const OcpValue ocp =
      positive ? positive_ocp(stoich) : negative_ocp(stoich);

  const double i0 = exchange_current(params, electrode, ce_clamped, cs);
  // 过电位 eta = phi_s - phi_e - U，为正表示偏阳极（有利于脱锂）。
  const double eta = phis - phie - ocp.value;
  // a、c 分别是 alpha_a F/(RT) 与 alpha_c F/(RT)，即两支 Tafel 斜率的倒数。
  const double a =
      electrode.alpha_a * kFaraday / (kGasConstant * params.temperature);
  const double c =
      electrode.alpha_c * kFaraday / (kGasConstant * params.temperature);
  // 指数参数夹紧到 ±60，防止 exp 上溢；夹紧处函数已饱和，不影响根。
  const double ea = std::exp(std::clamp(a * eta, -60.0, 60.0));
  const double ec = std::exp(std::clamp(-c * eta, -60.0, 60.0));

  ReactionValue result;
  // Butler-Volmer：j = i0/F [exp(a eta) - exp(-c eta)]，除以 F 得摩尔通量。
  result.j = i0 / kFaraday * (ea - ec);
  // dj/deta = i0/F (a exp(a eta) + c exp(-c eta))，恒为正（单调反应）。
  result.dj_deta = i0 / kFaraday * (a * ea + c * ec);
  // 由 i0 正比于 ce^alpha_a 得 dj/dce = j * alpha_a / ce。
  result.dj_dce = result.j * electrode.alpha_a / ce_clamped;
  // 两项贡献：i0 对 cs 的幂次求导，以及 U(cs) 经 eta 的间接影响。
  // 后者系数为 -dj_deta * dU/dcs，其中 dU/dcs = slope / c_max。
  result.dj_dcs =
      result.j * (-electrode.alpha_a / (cmax - cs) +
                  electrode.alpha_c / cs) -
      result.dj_deta * ocp.slope / cmax;
  return result;
}

// Butler-Volmer 反问题：已知通量 j，用牛顿迭代求解过电位 eta。
// 降阶模型中外加电流固定了 j，因此需要反过来求 eta 再算端电压。
// ce、cs_surf 为局部浓度（决定 i0）；j = 0 或 i0 <= 0 时直接返回 0。
double overpotential_from_flux(const Params& params,
                               const ElectrodeParams& electrode,
                               bool /*positive*/, double ce, double cs_surf,
                               double j) {
  if (j == 0.0) {
    return 0.0;
  }
  const double i0 = exchange_current(params, electrode, ce, cs_surf);
  if (i0 <= 0.0) {
    return 0.0;
  }
  const double a =
      electrode.alpha_a * kFaraday / (kGasConstant * params.temperature);
  const double c =
      electrode.alpha_c * kFaraday / (kGasConstant * params.temperature);

  double eta = 0.0;
  for (int iteration = 0; iteration < 60; ++iteration) {
    const double ea = std::exp(std::clamp(a * eta, -60.0, 60.0));
    const double ec = std::exp(std::clamp(-c * eta, -60.0, 60.0));
    // 残差为 Butler-Volmer 通量与目标通量之差，导数用解析式。
    const double residual = i0 / kFaraday * (ea - ec) - j;
    const double derivative = i0 / kFaraday * (a * ea + c * ec);
    const double step = residual / derivative;
    eta -= step;
    // 牛顿步长足够小即认为收敛；否则继续迭代直到上限。
    if (std::abs(step) < 1e-12) {
      break;
    }
  }
  return eta;
}

}  // namespace dfn
