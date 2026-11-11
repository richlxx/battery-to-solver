// 测试文件：DFN 模型的 IDA+KLU 时间积分（对应 models/dfn、solvers/idaklu 模块）。
// 覆盖恒流放电、锂守恒、电压截断事件与恒压保持三种典型工况。

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/solvers/idaklu.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

}  // namespace

int main() {
  dfn::Params params = dfn::make_demo_params();
  // 4/2/4 的 x 网格与 4 层径向网格：精度与运行时间的折中。
  const dfn::StateLayout layout = dfn::make_layout(4, 2, 4, 4, 4);
  dfn::DfnOptions control;
  control.applied_current_density = 8.0;  // 约 1C 放电
  dfn::DfnModel model(params, layout, control);

  const std::vector<double> y0 = model.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(model.size()), 0.0);
  const double lithium0 = model.total_lithium(y0.data());
  const double voltage0 = model.terminal_voltage(y0.data());

  // rtol = 1e-6 与建议 atol：在 60 s 放电中兼顾速度与守恒精度。
  // calculate_consistent_initial 让 IDA 修正 y0/yp0 中的代数不平衡。
  dfn::IdaOptions options;
  options.rtol = 1e-6;
  options.atol = model.suggested_atol();
  options.calculate_consistent_initial = true;
  options.output_times = {15.0, 30.0, 45.0};

  const dfn::SolveReport report =
      dfn::ida_integrate(model, y0, yp0, 0.0, 60.0, options);
  if (!report.ok) {
    std::cerr << report.message << "\n";
    return fail("DFN integration failed");
  }

  // 锂总量在封闭电池中应严格守恒；1e-6 的相对漂移对应积分器容差量级，
  // 超出说明质量守恒型离散被破坏或输出状态不一致。
  const double lithium1 = model.total_lithium(report.y.data());
  const double lithium_error = std::abs(lithium1 - lithium0) / lithium0;
  if (lithium_error > 1e-6) {
    std::cerr << "lithium " << lithium0 << " -> " << lithium1 << "\n";
    return fail("lithium is not conserved");
  }

  // 反应电流积分应平衡外加电流：负极 +I、正极 -I。
  // 2% 容差覆盖反应分布数值积分与局部过电位近似。
  const double expected_current = params.area * control.applied_current_density;
  const double i_neg = model.electrode_current(report.y.data(), false);
  const double i_pos = model.electrode_current(report.y.data(), true);
  if (std::abs(i_neg - expected_current) > 0.02 * expected_current ||
      std::abs(i_pos + expected_current) > 0.02 * expected_current) {
    std::cerr << "expected " << expected_current << ", I_neg " << i_neg
              << ", I_pos " << i_pos << "\n";
    return fail("electrode currents do not balance the applied current");
  }

  // 端电压必须落在物理区间 (0, 4.2] V；越界通常意味着 OCP 或电势参考错误。
  const double voltage = model.terminal_voltage(report.y.data());
  if (!std::isfinite(voltage) || voltage <= 0.0 || voltage > 4.2) {
    std::cerr << "voltage " << voltage << "\n";
    return fail("terminal voltage is implausible");
  }
  // 输出应包含 t0、三个请求时刻和终点 t_end，共至少 5 个样本。
  if (report.output_t.size() < 5) {
    std::cerr << "output samples " << report.output_t.size() << "\n";
    return fail("output sampling");
  }
  // 时间戳必须严格递增，避免插值或输出逻辑重复/倒退。
  for (std::size_t k = 1; k < report.output_t.size(); ++k) {
    if (report.output_t[k] <= report.output_t[k - 1]) {
      return fail("output times are not increasing");
    }
  }

  std::cout << "dfn ida ok: steps " << report.steps << ", jac "
            << report.jacobian_evals << ", V0 " << voltage0 << " -> " << voltage
            << ", Li drift " << std::scientific << std::setprecision(3)
            << lithium_error << std::defaultfloat << ", I_neg " << i_neg
            << " A, samples " << report.output_t.size() << "\n";

  // 工况二：放电到 3.5 V 截止。事件函数 V(y) - 3.5 过零时，
  // IDA 的 rootfinding 应终止积分并报告事件下标 0，且锂仍守恒。
  const double cutoff = 3.5;
  dfn::DfnModel discharge(params, layout, control);
  const std::vector<double> yd0 = discharge.initial_state();
  dfn::IdaOptions cutoff_options;
  cutoff_options.rtol = 1e-6;
  cutoff_options.atol = discharge.suggested_atol();
  cutoff_options.calculate_consistent_initial = true;
  // 事件函数在端电压降到 cutoff 时由正变负，触发停止。
  cutoff_options.events.push_back(
      [&discharge, cutoff](double, const double* y, const double*) {
        return discharge.terminal_voltage(y) - cutoff;
      });
  // 每 1000 s 输出一次，覆盖 0 到 20000 s 的上限时间窗。
  for (int k = 0; k <= 20; ++k) {
    cutoff_options.output_times.push_back(1000.0 * static_cast<double>(k));
  }
  const double lithium_start = discharge.total_lithium(yd0.data());
  const dfn::SolveReport cutoff_report =
      dfn::ida_integrate(discharge, yd0, yp0, 0.0, 20000.0, cutoff_options);
  // 只注册了一个事件，因此成功停止时 event_index 必须为 0。
  if (!cutoff_report.ok || cutoff_report.event_index != 0) {
    std::cerr << cutoff_report.message << " event "
              << cutoff_report.event_index << "\n";
    return fail("voltage cutoff event was not reported");
  }
  // 事件定位精度：1e-3 V 足以确认事件落在截止电压附近，
  // 同时容纳根求解的插值误差；超出说明事件点未被正确精化。
  const double v_event = discharge.terminal_voltage(cutoff_report.y.data());
  if (std::abs(v_event - cutoff) > 1e-3) {
    std::cerr << "event voltage " << v_event << "\n";
    return fail("event voltage is off");
  }
  // 事件停止不得破坏守恒：相对漂移仍需低于 1e-6。
  const double lithium_event = discharge.total_lithium(cutoff_report.y.data());
  if (std::abs(lithium_event - lithium_start) / lithium_start > 1e-6) {
    return fail("lithium drifted during discharge to cutoff");
  }
  // 事件时间必须严格落在积分区间内部，否则说明事件未真正触发。
  if (cutoff_report.t_event <= 0.0 || cutoff_report.t_event >= 20000.0) {
    std::cerr << "t_event " << cutoff_report.t_event << "\n";
    return fail("event time is outside the interval");
  }

  std::cout << "cutoff ok: t_event " << cutoff_report.t_event << " s, steps "
            << cutoff_report.steps << ", V " << v_event << "\n";

  // 工况三：在截止点切换为恒压（CV）控制。端电压被钳在 3.5 V，
  // 电流幅值应随浓度弛豫逐渐衰减，锂总量仍守恒。
  discharge.set_applied_voltage(cutoff);
  const double current0 = discharge.terminal_current(cutoff_report.y.data());
  dfn::IdaOptions cv_options;
  cv_options.rtol = 1e-6;
  cv_options.atol = discharge.suggested_atol();
  cv_options.calculate_consistent_initial = true;
  const double cv_end = cutoff_report.t_event + 300.0;
  const double lithium_cv0 = discharge.total_lithium(cutoff_report.y.data());
  const dfn::SolveReport cv_report =
      dfn::ida_integrate(discharge, cutoff_report.y, cutoff_report.yp,
                         cutoff_report.t_event, cv_end, cv_options);
  if (!cv_report.ok) {
    std::cerr << cv_report.message << "\n";
    return fail("CV hold failed");
  }
  const double current1 = discharge.terminal_current(cv_report.y.data());
  // 恒压阶段颗粒表面浓度趋于均匀，反应电流幅值必须严格下降。
  if (std::abs(current1) >= std::abs(current0)) {
    std::cerr << "current " << current0 << " -> " << current1 << "\n";
    return fail("CV current did not decay");
  }
  // 恒压是代数约束，端电压应以机器精度等于设定值 3.5 V。
  if (std::abs(discharge.terminal_voltage(cv_report.y.data()) - cutoff) > 1e-12) {
    return fail("CV voltage is not held");
  }
  const double lithium_cv1 = discharge.total_lithium(cv_report.y.data());
  if (std::abs(lithium_cv1 - lithium_cv0) / lithium_cv0 > 1e-6) {
    return fail("lithium drifted during CV hold");
  }

  std::cout << "cv ok: I " << current0 << " -> " << current1 << " A, steps "
            << cv_report.steps << "\n";
  return 0;
}
