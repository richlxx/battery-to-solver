#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/solvers/idaklu.hpp"

// 第 9 章配套实验：用完整 DFN 模型跑恒流放电，演示最小端到端流程。
// 第一部分是网格收敛表（V(60 s)、未知量数 n、非零元 nnz、步数、Jacobian 次数、
// 耗时）；第二部分是 CC 放电到 3.5 V 截止后转 CV 保持的电压与电流输出。

namespace {

using Clock = std::chrono::steady_clock;

// 统一的 IDA 调用封装：导数初值取零向量，由 IDA 自行计算一致初值；
// rtol 由调用者给定，atol 用模型建议值，*ms 回传墙钟耗时（毫秒）。
dfn::SolveReport run_cc(dfn::DfnModel& model, const std::vector<double>& y0,
                        double t_end, double rtol, double* ms,
                        const std::vector<double>& output_times = {}) {
  const std::vector<double> yp0(static_cast<std::size_t>(model.size()), 0.0);
  dfn::IdaOptions options;
  options.rtol = rtol;
  options.atol = model.suggested_atol();
  options.calculate_consistent_initial = true;
  options.output_times = output_times;
  const auto start = Clock::now();
  dfn::SolveReport report =
      dfn::ida_integrate(model, y0, yp0, 0.0, t_end, options);
  *ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  return report;
}

// 网格收敛性研究：固定电流密度 8 A/m^2，同时改变 x 方向（负极/隔膜/正极）
// 与径向（nr）网格数，观察 60 s 端电压、步数、Jacobian 次数与耗时的变化。
void grid_convergence(const dfn::Params& params) {
  struct LayoutSpec {
    int nx_neg;
    int nx_sep;
    int nx_pos;
    int nr_neg;
    int nr_pos;
  };

  // 对一组网格规格依次求解并打印表格，以最细网格为参考解，
  // 计算各级电压误差与相邻两级的观测收敛阶 log2(e_k / e_{k+1})。
  const auto study = [&params](const char* name,
                               const std::vector<LayoutSpec>& specs) {
    std::cout << "== " << name << ": V at t = 60 s, 8 A/m^2 ==\n";
    std::cout << "n_x(n/s/p)\tn_r(n/p)\tn\tnnz\tsteps\tjac\tV(60)\tms\n";
    std::vector<double> voltages;
    for (const auto& spec : specs) {
      dfn::DfnOptions control;
      control.applied_current_density = 8.0;
      dfn::DfnModel model(
          params,
          dfn::make_layout(spec.nx_neg, spec.nx_sep, spec.nx_pos, spec.nr_neg,
                           spec.nr_pos),
          control);
      double ms = 0.0;
      const dfn::SolveReport report =
          run_cc(model, model.initial_state(), 60.0, 1e-9, &ms);
      if (!report.ok) {
        std::cout << "failed: " << report.message << "\n";
        return;
      }
      voltages.push_back(model.terminal_voltage(report.y.data()));
      std::cout << spec.nx_neg << "/" << spec.nx_sep << "/" << spec.nx_pos
                << "\t" << spec.nr_neg << "/" << spec.nr_pos << "\t"
                << model.size() << "\t" << model.jacobian_pattern().nnz() << "\t"
                << report.steps << "\t" << report.jacobian_evals << "\t"
                << model.terminal_voltage(report.y.data()) << "\t" << ms << "\n";
    }
    // 最细网格的结果作为参考解；收敛良好时观测阶应为正且趋于格式阶数。
    const double reference = voltages.back();
    std::cout << "differences vs finest:\n";
    for (std::size_t k = 0; k + 1 < voltages.size(); ++k) {
      std::cout << "  level " << k << ": " << std::abs(voltages[k] - reference)
                << "\n";
    }
    for (std::size_t k = 0; k + 2 < voltages.size(); ++k) {
      const double e1 = std::abs(voltages[k] - reference);
      const double e2 = std::abs(voltages[k + 1] - reference);
      if (e2 > 0.0) {
        std::cout << "  observed order " << k << "->" << k + 1 << ": "
                  << std::log2(e1 / e2) << "\n";
      }
    }
    std::cout << "\n";
  };

  // 三级研究分别只加密 x、只加密径向、以及两者同时加密，
  // 用来分离空间离散误差与径向（颗粒内扩散）离散误差。
  study("x refinement (n_r fixed at 8)",
        {{4, 2, 4, 8, 8}, {8, 4, 8, 8, 8}, {16, 8, 16, 8, 8}, {32, 16, 32, 8, 8}});
  study("radial refinement (n_x fixed at 8/4/8)",
        {{8, 4, 8, 4, 4}, {8, 4, 8, 6, 6}, {8, 4, 8, 8, 8}, {8, 4, 8, 12, 12}});
  study("both refined",
        {{4, 2, 4, 4, 4}, {8, 4, 8, 6, 6}, {16, 8, 16, 8, 8}, {32, 16, 32, 12, 12}});
}

// 先以 8 A/m^2 恒流放电，用 IDA 事件定位端电压跌破 3.5 V 的时刻；
// 再从事件状态出发把边界条件切换为 3.5 V 恒压，观察电流如何衰减。
void cc_cv_curve(const dfn::Params& params) {
  dfn::DfnOptions control;
  control.applied_current_density = 8.0;
  dfn::DfnModel model(params, dfn::make_layout(8, 4, 8, 6, 6), control);
  const std::vector<double> y0 = model.initial_state();
  const std::vector<double> yp0(static_cast<std::size_t>(model.size()), 0.0);
  const double cutoff = 3.5;

  dfn::IdaOptions cc;
  cc.rtol = 1e-6;
  cc.atol = model.suggested_atol();
  cc.calculate_consistent_initial = true;
  // 事件函数：端电压与截止电压之差过零时触发，IDA 记录 t_event 并停止。
  cc.events.push_back([&model, cutoff](double, const double* y, const double*) {
    return model.terminal_voltage(y) - cutoff;
  });
  for (int k = 0; k <= 8; ++k) {
    cc.output_times.push_back(2000.0 * static_cast<double>(k));
  }
  const dfn::SolveReport cc_report =
      dfn::ida_integrate(model, y0, yp0, 0.0, 20000.0, cc);
  if (!cc_report.ok) {
    std::cout << "CC failed: " << cc_report.message << "\n";
    return;
  }

  std::cout << "== CC discharge, cutoff " << cutoff << " V ==\n";
  std::cout << "t\tV\tsteps(cumulative)\n";
  for (std::size_t k = 0; k < cc_report.output_t.size(); ++k) {
    std::cout << cc_report.output_t[k] << "\t"
              << model.terminal_voltage(cc_report.output_y[k].data()) << "\n";
  }
  std::cout << "event at " << cc_report.t_event << " s, steps "
            << cc_report.steps << ", jac " << cc_report.jacobian_evals << "\n\n";

  // 切换到恒压模式：初值沿用 CC 结束时的 y 与 yp，保证状态连续。
  model.set_applied_voltage(cutoff);
  dfn::IdaOptions cv;
  cv.rtol = 1e-6;
  cv.atol = model.suggested_atol();
  cv.calculate_consistent_initial = true;
  for (int k = 1; k <= 10; ++k) {
    cv.output_times.push_back(cc_report.t_event + 60.0 * static_cast<double>(k));
  }
  const dfn::SolveReport cv_report = dfn::ida_integrate(
      model, cc_report.y, cc_report.yp, cc_report.t_event,
      cc_report.t_event + 600.0, cv);
  if (!cv_report.ok) {
    std::cout << "CV failed: " << cv_report.message << "\n";
    return;
  }
  std::cout << "== CV hold at " << cutoff << " V ==\n";
  std::cout << "t since switch\tI (A/m^2)\n";
  for (std::size_t k = 0; k < cv_report.output_t.size(); ++k) {
    std::cout << cv_report.output_t[k] - cc_report.t_event << "\t"
              << model.terminal_current(cv_report.output_y[k].data()) << "\n";
  }
  std::cout << "CV steps " << cv_report.steps << "\n\n";
}

}  // namespace

// 实验入口：先做网格收敛分析，再跑 CC→CV 放电曲线。
int main() {
  std::cout << std::setprecision(6);
  const dfn::Params params = dfn::make_demo_params();
  grid_convergence(params);
  cc_cv_curve(params);
  return 0;
}
