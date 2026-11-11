// 仿真驱动：把 CellModel（SPM/SPMe/DFN）交给 IDA 做时间推进，
// 实现 CC / CV / CCCV / pulse 工况的阶段切换、采样输出与截止判定，
// 并统计库仑容量与锂守恒漂移。
#include "dfn/sim/simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <utility>

#include "dfn/models/cell_model.hpp"

#ifdef DFN_USE_KLU
#include "dfn/solvers/idaklu.hpp"
#endif

namespace dfn {

namespace {

using Clock = std::chrono::steady_clock;

#ifdef DFN_USE_KLU

struct PhaseTotals {
  int steps = 0;
  int residual_evals = 0;
  int jacobian_evals = 0;
  int linear_setups = 0;
};

// 生成 (t0, t1) 内形如 k*dt 的采样时刻。ida_integrate 会自行追加 t1，
// 所以这里严格只取 t < t1 的点，避免端点重复；dt <= 0 或区间无效返回空表。
std::vector<double> sample_times(double t0, double t1, double dt) {
  std::vector<double> times;
  if (dt <= 0.0 || t1 <= t0) {
    return times;
  }
  const double first = std::floor(t0 / dt) + 1.0;
  for (double k = first;; k += 1.0) {
    const double t = k * dt;
    if (t >= t1 - 1e-9) {
      break;
    }
    times.push_back(t);
  }
  return times;
}

// 恒流阶段的电流密度（A/m²）：优先使用显式 current_density，
// 否则按 c_rate * i_1c 换算；充电时取负号（项目约定正值 = 放电）。
double controlled_current(const SimulationConfig& config) {
  const double magnitude =
      config.current_density != 0.0
          ? std::abs(config.current_density)
          : std::abs(config.c_rate) * config.i_1c;
  return config.charge ? -magnitude : magnitude;
}

// 把一个阶段的求解输出并入总报告：逐时刻记录端电压、端电流、
// 正负极平均化学计量比与阶段名。不大于上一采样点的时刻会被跳过，
// 保证跨阶段拼接后 samples 的时间仍然单调递增。
void append_samples(const CellModel& model, const SolveReport& report,
                    const std::string& phase, SimulationReport& out) {
  for (std::size_t k = 0; k < report.output_t.size(); ++k) {
    const double t = report.output_t[k];
    if (!out.samples.empty() && t <= out.samples.back().t + 1e-9) {
      continue;
    }
    const double* y = report.output_y[k].data();
    SimulationSample sample;
    sample.t = t;
    sample.voltage = model.terminal_voltage(y);
    sample.current = model.terminal_current(y);
    sample.stoich_neg = model.average_stoichiometry(y, false);
    sample.stoich_pos = model.average_stoichiometry(y, true);
    sample.phase = phase;
    out.samples.push_back(std::move(sample));
  }
}

// 推进一个工况阶段：[t0, t1] 内调用 ida_integrate，统计量累加到 totals；
// 若事件触发则记录停止标志与事件时刻，否则停在 t1。
// phase 是阶段名（"cc" / "cv" / "pulse-on" / "pulse-rest"），写入采样标签。
// 成功后把末态 y、yp 写回，供下一阶段热启动；返回 false 表示积分失败。
bool run_phase(CellModel& model, std::vector<double>& y, std::vector<double>& yp,
               double t0, double t1, const IdaOptions& options, double dt,
               const std::string& phase, SimulationReport& out,
               PhaseTotals& totals, bool* stopped_by_event, double* t_event) {
  IdaOptions phase_options = options;
  phase_options.rtol = options.rtol;
  phase_options.output_times = sample_times(t0, t1, dt);
  const SolveReport report =
      ida_integrate(model, y, yp, t0, t1, phase_options);
  *stopped_by_event = report.event_index >= 0;
  *t_event = report.t_event;

  totals.steps += report.steps;
  totals.residual_evals += report.residual_evals;
  totals.jacobian_evals += report.jacobian_evals;
  totals.linear_setups += report.linear_setups;

  if (!report.ok) {
    out.message = phase + ": " + report.message;
    return false;
  }
  append_samples(model, report, phase, out);
  y = report.y;
  yp = report.yp;
  return true;
}

// 单个电极单位隔膜面积的储锂量，mol/m²：θ * c_max * (1-ε) * L。
double electrode_lithium_per_area(const CellModel& model, const double* y,
                                  bool positive) {
  const Params& params = model.params();
  const ElectrodeParams& electrode = positive ? params.pos : params.neg;
  return model.average_stoichiometry(y, positive) * electrode.c_max *
         (1.0 - electrode.porosity) * electrode.thickness;
}

// 按配置创建具体模型：spm / spme / dfn 对应模型阶梯上的三个层级，
// grid 参数与径向节点数一并传入；DFN 还需用 make_layout 建立状态分块。
std::unique_ptr<CellModel> make_model(const SimulationConfig& config) {
  switch (config.model) {
    case ModelKind::spm:
      return std::make_unique<SpmModel>(config.params, config.n_r_neg,
                                        config.n_r_pos);
    case ModelKind::spme:
      return std::make_unique<SpmeModel>(
          config.params, config.n_x_neg, config.n_x_sep, config.n_x_pos,
          config.n_r_neg, config.n_r_pos);
    case ModelKind::dfn:
      return std::make_unique<DfnModel>(
          config.params,
          make_layout(config.n_x_neg, config.n_x_sep, config.n_x_pos,
                      config.n_r_neg, config.n_r_pos),
          DfnOptions{});
  }
  return nullptr;
}

#endif  // DFN_USE_KLU

}  // namespace

// 顶层驱动：参数校验 → 建立模型与一致初值 → 按工况分阶段推进 →
// 汇总守恒/容量统计 → 记录耗时。返回的 SimulationReport 同时承载
// 采样序列、失败信息（message）与求解器统计量。
SimulationReport run_simulation(const SimulationConfig& config) {
  SimulationReport report;
  const auto start = Clock::now();

  // 先做参数与网格检查，让配置错误以可读信息返回而不是抛异常。
  try {
    validate(config.params);
  } catch (const std::exception& error) {
    report.message = error.what();
    return report;
  }
  if (config.n_x_neg <= 0 || config.n_x_sep <= 0 || config.n_x_pos <= 0 ||
      config.n_r_neg <= 0 || config.n_r_pos <= 0) {
    report.message = "grid sizes must be positive";
    return report;
  }
  if (config.t_max <= 0.0) {
    report.message = "t_max must be positive";
    return report;
  }
  if (config.mode == SimulationMode::pulse &&
      (config.t_pulse <= 0.0 || config.pulse_duty <= 0.0 ||
       config.pulse_duty >= 1.0)) {
    report.message = "pulse period must be positive and duty in (0, 1)";
    return report;
  }

#ifdef DFN_USE_KLU
  std::unique_ptr<CellModel> model_ptr = make_model(config);
  if (!model_ptr) {
    report.message = "unknown model";
    return report;
  }
  CellModel& model = *model_ptr;
  // CV/CCCV 需要模型提供电压控制接口；SPM/SPMe 只有电流控制。
  if ((config.mode == SimulationMode::cv ||
       config.mode == SimulationMode::cccv) &&
      !model.supports_voltage_control()) {
    report.message = std::string(model.name()) +
                     " does not support cv/cccv; use --model dfn";
    return report;
  }
  std::vector<double> y = model.initial_state();
  std::vector<double> yp(static_cast<std::size_t>(model.size()), 0.0);
  report.model = model.name();
  report.unknowns = model.size();
  report.nnz = model.jacobian_pattern().nnz();
  // 初始锂总量与负极储锂量是后面守恒检查和容量统计的基准。
  const double lithium0 = model.total_lithium(y.data());
  const double neg_lithium0 = electrode_lithium_per_area(model, y.data(), false);

  IdaOptions options;
  options.rtol = config.rtol;
  options.atol = model.suggested_atol();
  // 让 IDA 先修正初值，使其满足代数约束，避免起步瞬间出现不连续。
  options.calculate_consistent_initial = true;

  PhaseTotals totals;
  // i_cc：恒流电流；v_cut：截止电压（充电用 v_max、放电用 v_min）；
  // i_cut：CV 段的终止阈值，取恒流幅值的 5%（|I| 降到该值即视为到限）。
  const double i_cc = controlled_current(config);
  const double v_cut = config.charge ? config.v_max : config.v_min;
  const double i_cut = 0.05 * std::abs(i_cc);

  // 截止电压事件：返回值过零即触发。充电检测 V 升到 v_cut，放电检测 V 降到 v_cut。
  auto voltage_cutoff = [&model, v_cut, charge = config.charge](
                            double, const double* state, const double*) {
    const double voltage = model.terminal_voltage(state);
    return charge ? v_cut - voltage : voltage - v_cut;
  };

  double t = 0.0;
  bool stopped = false;
  double t_event = 0.0;

  switch (config.mode) {
    case SimulationMode::cc:
    case SimulationMode::cccv: {
      // CC 阶段：从初值出发恒流推进，事件在端电压触及截止值时停机。
      const double v0 = model.terminal_voltage(y.data());
      if (config.charge ? v0 >= v_cut : v0 <= v_cut) {
        report.message = config.charge
                             ? "initial voltage is already above the cutoff"
                             : "initial voltage is already below the cutoff";
        return report;
      }
      model.set_current_density(i_cc);
      options.events = {voltage_cutoff};
      if (!run_phase(model, y, yp, 0.0, config.t_max, options,
                     config.output_interval, "cc", report, totals, &stopped,
                     &t_event)) {
        return report;
      }
      t = stopped ? t_event : config.t_max;
      report.cutoff_reached = stopped;

      // CCCV 的第二段：保持截止电压，事件换成 |I| 降到 i_cut。
      // 恒流段未触发事件（电压还没到限值）或 t_hold <= 0 时跳过 CV 段。
      if (config.mode == SimulationMode::cccv && stopped &&
          config.t_hold > 0.0) {
        model.set_applied_voltage(v_cut);
        IdaOptions cv_options = options;
        cv_options.events.clear();
        if (i_cut > 0.0) {
          cv_options.events.push_back(
              [&model, i_cut](double, const double* state, const double*) {
                return std::abs(model.terminal_current(state)) - i_cut;
              });
        }
        bool cv_stopped = false;
        double t_cv_event = 0.0;
        if (!run_phase(model, y, yp, t, t + config.t_hold, cv_options,
                       config.output_interval, "cv", report, totals,
                       &cv_stopped, &t_cv_event)) {
          return report;
        }
        t = cv_stopped ? t_cv_event : t + config.t_hold;
      }
      break;
    }
    case SimulationMode::cv: {
      // CV 模式要求初值已经处在保持电压上；若在满电/空电的新电池上
      // 直接跳变电压，应改用 cccv 让电流自然衰减。
      const double v_hold = config.hold_voltage > 0.0 ? config.hold_voltage
                                                      : v_cut;
      if (std::abs(model.terminal_voltage(y.data()) - v_hold) > 1e-4) {
        report.message =
            "cv mode must start at the hold voltage; run cccv instead of "
            "stepping the voltage on a fresh cell";
        return report;
      }
      model.set_applied_voltage(v_hold);
      if (i_cut > 0.0) {
        options.events.push_back(
            [&model, i_cut](double, const double* state, const double*) {
              return std::abs(model.terminal_current(state)) - i_cut;
            });
      }
      bool cv_stopped = false;
      double t_cv_event = 0.0;
      if (!run_phase(model, y, yp, 0.0, config.t_max, options,
                     config.output_interval, "cv", report, totals, &cv_stopped,
                     &t_cv_event)) {
        return report;
      }
      t = cv_stopped ? t_cv_event : config.t_max;
      break;
    }
    case SimulationMode::pulse: {
      // pulse 模式：每个周期拆成“通电 on_time”与“静置 off_time”两段，
      // 两段都带截止电压事件；t_max 截断最后一段。
      const double v0 = model.terminal_voltage(y.data());
      if (config.charge ? v0 >= v_cut : v0 <= v_cut) {
        report.message = config.charge
                             ? "initial voltage is already above the cutoff"
                             : "initial voltage is already below the cutoff";
        return report;
      }
      options.events = {voltage_cutoff};
      // 占空比 duty 决定每个周期内通电时间占比，其余时间电流置零静置。
      const double on_time = config.t_pulse * config.pulse_duty;
      const double off_time = config.t_pulse - on_time;
      while (t < config.t_max - 1e-9) {
        const double t_on_end = std::min(t + on_time, config.t_max);
        model.set_current_density(i_cc);
        if (!run_phase(model, y, yp, t, t_on_end, options,
                       config.output_interval, "pulse-on", report, totals,
                       &stopped, &t_event)) {
          return report;
        }
        t = stopped ? t_event : t_on_end;
        if (stopped) {
          report.cutoff_reached = true;
          break;
        }
        const double t_off_end = std::min(t + off_time, config.t_max);
        model.set_current_density(0.0);
        if (!run_phase(model, y, yp, t, t_off_end, options,
                       config.output_interval, "pulse-rest", report, totals,
                       &stopped, &t_event)) {
          return report;
        }
        t = stopped ? t_event : t_off_end;
        if (stopped) {
          report.cutoff_reached = true;
          break;
        }
      }
      break;
    }
  }

  // 守恒与容量统计：锂只在固相/液相之间搬运，总量相对漂移应停在数值误差量级；
  // 容量取负极储锂变化量乘以法拉第常数，单位为 A s/m²。
  const double lithium1 = model.total_lithium(y.data());
  report.lithium_drift = std::abs(lithium1 - lithium0) / lithium0;
  const double moved =
      neg_lithium0 - electrode_lithium_per_area(model, y.data(), false);
  report.capacity = kFaraday * (config.charge ? -moved : moved);
  report.t_end = report.samples.empty() ? t : report.samples.back().t;
  report.steps = totals.steps;
  report.residual_evals = totals.residual_evals;
  report.jacobian_evals = totals.jacobian_evals;
  report.linear_setups = totals.linear_setups;
  report.ok = true;
#else
  // 未启用 SUNDIALS 时不做积分，直接返回构建提示。
  (void)config;
  report.message =
      "dfn-sim was built without SUNDIALS; configure with -DDFN_USE_KLU=ON";
#endif

  report.wall_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  return report;
}

// 写 CSV：列依次为 t_s, voltage_v, current_a_m2, stoich_neg, stoich_pos, phase。
void write_csv(const SimulationReport& report, std::ostream& out) {
  out << "t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,phase\n";
  out << std::setprecision(10);
  for (const SimulationSample& sample : report.samples) {
    out << sample.t << ',' << sample.voltage << ',' << sample.current << ','
        << sample.stoich_neg << ',' << sample.stoich_pos << ','
        << sample.phase << '\n';
  }
}

void write_csv(const SimulationReport& report, const std::string& path) {
  // --out results/cc.csv 这类路径的父目录可能尚不存在，先建出来。
  const std::filesystem::path parent =
      std::filesystem::path(path).parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("cannot open " + path + " for writing");
  }
  write_csv(report, file);
}

}  // namespace dfn
