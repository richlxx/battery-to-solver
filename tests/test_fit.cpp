// 测试文件：SPM 批量正演与 Levenberg-Marquardt 参数拟合
// （对应 fit/spm_batch、fit/levenberg_marquardt、models/spm 模块）。
// 验证定步长近似、充放电对称性以及无噪声参数恢复。

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/fit/levenberg_marquardt.hpp"
#include "dfn/fit/spm_batch.hpp"
#include "dfn/models/spm.hpp"

#ifdef DFN_USE_KLU
#include "dfn/solvers/idaklu.hpp"
#endif

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

// 拟合实验用参数：降低两侧反应速率常数，使电压对 k 和 D_s 的敏感度适中，
// LM 的数值 Jacobian 不易病态。
dfn::Params fit_params() {
  dfn::Params params = dfn::make_demo_params();
  params.neg.k = 3.0e-11;
  params.pos.k = 5.0e-11;
  return params;
}

// 工况 C：包含静置、放电、充电与高倍率脉冲的 9 段协议（duration s，A/m^2）。
// 多段激励让 k 与 D_s 的不同时间尺度特征都能被拟合识别。
std::vector<dfn::fit::BatchSegment> protocol_c() {
  return {{600.0, 0.0},   {1200.0, 18.0}, {1800.0, 0.0},
          {1200.0, -18.0}, {1800.0, 0.0},  {600.0, 36.0},
          {1800.0, 0.0},  {300.0, 72.0},   {1800.0, 0.0}};
}

}  // namespace

int main() {
  const dfn::Params params = fit_params();
  const int nr = 20;  // 每个颗粒 20 个径向单元，兼顾精度与批量速度

  // 定步长后向欧拉 SPM 应随 dt 减小收敛到 IDA 参考解；
  // 先分别用 dt = 5 s 与 10 s 生成 600 s 恒流放电的末态电压。
  const std::vector<dfn::fit::BatchRow> discharge =
      dfn::fit::rows_from_protocol(params, {{600.0, 18.0}}, nr, 5.0);
  const std::vector<dfn::fit::BatchRow> coarse =
      dfn::fit::rows_from_protocol(params, {{600.0, 18.0}}, nr, 10.0);

#ifdef DFN_USE_KLU
  // IDA 参考解：rtol = 1e-8，且由 IDA 计算一致初值，误差远小于时间离散。
  dfn::SpmModel spm(params, nr, nr);
  spm.set_current_density(18.0);
  const std::vector<double> y0 = spm.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(spm.size()), 0.0);
  dfn::IdaOptions options;
  options.rtol = 1e-8;
  options.atol = spm.suggested_atol();
  options.calculate_consistent_initial = true;
  const dfn::SolveReport reference =
      dfn::ida_integrate(spm, y0, yp0, 0.0, 600.0, options);
  if (!reference.ok) {
    std::cerr << reference.message << "\n";
    return fail("IDA reference run failed");
  }
  const double reference_voltage = spm.terminal_voltage(reference.y.data());
  const double fine_error =
      std::abs(discharge.back().voltage - reference_voltage);
  const double coarse_error =
      std::abs(coarse.back().voltage - reference_voltage);
  // 后向欧拉是一阶方法，误差约与 dt 成正比：dt = 5 s 允许 0.2 mV，
  // dt = 10 s 允许 2 mV；超差说明批量模型与 SPM 方程不一致。
  if (fine_error > 2e-4 || coarse_error > 2e-3) {
    std::cerr << "fixed-step error " << fine_error << " / " << coarse_error
              << " V\n";
    return fail("fixed-step SPM disagrees with IDA");
  }

  // 充电是放电的镜像：只改变电流符号，电压偏差应与放电同量级。
  dfn::SpmModel charge_model(params, nr, nr);
  charge_model.set_current_density(-18.0);
  const dfn::SolveReport charge_ref = dfn::ida_integrate(
      charge_model, y0, yp0, 0.0, 600.0, options);
  if (!charge_ref.ok) {
    return fail("IDA charge reference failed");
  }
  const std::vector<dfn::fit::BatchRow> charge_rows =
      dfn::fit::rows_from_protocol(params, {{600.0, -18.0}}, nr, 5.0);
  const double charge_error =
      std::abs(charge_rows.back().voltage -
               charge_model.terminal_voltage(charge_ref.y.data()));
  // 充电方向采用相同的 0.2 mV 阈值，验证模型关于电流符号的对称性。
  if (charge_error > 2e-4) {
    std::cerr << "charge error " << charge_error << " V\n";
    return fail("fixed-step charge disagrees with IDA");
  }
#endif

  // 无噪声参数恢复：用真值生成工况 C 的目标电压，再从扰动初值拟合。
  const std::vector<dfn::fit::BatchRow> rows =
      dfn::fit::rows_from_protocol(params, protocol_c(), nr, 5.0);
  // 待拟合参数 theta = [log k 缩放, log D_s,neg, log D_s,pos]；
  // 用对数参数化保证 k、D_s 始终为正，并改善量级差异。
  const std::vector<double> truth = {0.0, std::log(params.neg.d_s),
                                     std::log(params.pos.d_s)};
  // 正演必须可复现：同一参数、同一协议应逐行还原目标电压（1e-10 V）。
  const std::vector<double> reproduced =
      dfn::fit::simulate_rows(params, rows, nr, 5.0);
  double reproduce_error = 0.0;
  for (std::size_t i = 0; i < rows.size(); ++i) {
    reproduce_error =
        std::max(reproduce_error, std::abs(reproduced[i] - rows[i].voltage));
  }
  if (reproduce_error > 1e-10) {
    std::cerr << "reproduction error " << reproduce_error << " V\n";
    return fail("simulate_rows does not reproduce the generated rows");
  }
  // 把无量纲 theta 映射回物理参数：theta[0] 同时缩放两侧 k，
  // theta[1]/theta[2] 分别是指数化的固相扩散系数。
  auto apply = [&](const std::vector<double>& theta) {
    dfn::Params p = params;
    const double scale = std::exp(theta[0]);
    p.neg.k *= scale;
    p.pos.k *= scale;
    p.neg.d_s = std::exp(theta[1]);
    p.pos.d_s = std::exp(theta[2]);
    return p;
  };
  // LM 残差：预测电压减目标电压，单位 V；返回 true 表示正演成功。
  const auto residual = [&](const std::vector<double>& theta,
                            std::vector<double>& values) {
    const std::vector<double> predicted =
        dfn::fit::simulate_rows(apply(theta), rows, nr, 5.0);
    values.resize(predicted.size());
    for (std::size_t i = 0; i < predicted.size(); ++i) {
      values[i] = predicted[i] - rows[i].voltage;
    }
    return true;
  };

  // 从真值附近扰动出发（log k 偏 0.3，两个 log D_s 偏 +0.25/-0.25），
  // 检验 LM 能否收敛回全局真值。
  std::vector<double> start = truth;
  start[0] += 0.3;
  start[1] += 0.25;
  start[2] -= 0.25;
  dfn::LmOptions lm;
  lm.max_iterations = 40;
  const dfn::LmResult fit = dfn::levenberg_marquardt(start, residual, lm);
  if (!fit.ok) {
    return fail("LM did not converge");
  }
  // 无噪声问题的最优代价应为 0；0.02 的对数参数误差约合 2% 物理参数偏差。
  // cost > 1e-16（0.5*||r||^2）说明拟合停在非零残差的局部极小。
  double worst = 0.0;
  for (std::size_t j = 0; j < truth.size(); ++j) {
    worst = std::max(worst, std::abs(fit.theta[j] - truth[j]));
  }
  if (worst > 0.02 || fit.cost > 1e-16) {
    std::cerr << "worst parameter error " << worst << ", cost " << fit.cost
              << "\n";
    return fail("LM did not recover the noise-free parameters");
  }

  std::cout << "fit ok (samples=" << rows.size() << ", parameter error "
            << worst << ", cost " << fit.cost << ")\n";
  return 0;
}
