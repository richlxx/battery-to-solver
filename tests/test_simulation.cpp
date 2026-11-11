// 测试文件：参数加载与仿真驱动（对应 sim/config、sim/simulation 模块）。
// 验证 YAML 参数、CC-CV/脉冲/充电工况、SPM/SPMe 模型选择、CSV 输出与守恒。

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "dfn/sim/config.hpp"
#include "dfn/sim/simulation.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

// 绝对容差比较：调用方显式给出电压/时间阈值，单位直观。
bool nearly(double a, double b, double tolerance) {
  return std::abs(a - b) <= tolerance;
}

}  // namespace

int main() {
  // YAML 中的 demo 参数必须与内置 make_demo_params 逐项一致；
  // 容差为 0 表示要求解析结果与源码常数完全相同（同一十进制字面量）。
  const dfn::Params builtin = dfn::make_demo_params();
  const dfn::Params loaded =
      dfn::load_params_yaml(std::string(DFN_SOURCE_DIR) + "/configs/demo.yaml");
  if (!nearly(loaded.area, builtin.area, 0.0) ||
      !nearly(loaded.neg.thickness, builtin.neg.thickness, 0.0) ||
      !nearly(loaded.pos.c0, builtin.pos.c0, 0.0) ||
      !nearly(loaded.sep.porosity, builtin.sep.porosity, 0.0) ||
      !nearly(loaded.elyte.t_plus, builtin.elyte.t_plus, 0.0)) {
    return fail("YAML loader does not reproduce the demo parameters");
  }

#ifndef DFN_USE_KLU
  // 未编译 SUNDIALS 时只能检查参数加载，仿真相关分支跳过。
  std::cout << "dfn sim ok (parameter loader only, no SUNDIALS)\n";
  return 0;
#else
  // 工况一：CC-CV 放电。小网格（4/2/4，径向 4）缩短测试时间。
  dfn::SimulationConfig config;
  config.mode = dfn::SimulationMode::cccv;
  config.n_x_neg = 4;
  config.n_x_sep = 2;
  config.n_x_pos = 4;
  config.n_r_neg = 4;
  config.n_r_pos = 4;
  config.t_hold = 400.0;          // CV 段最长 400 s
  config.output_interval = 300.0;  // 每 300 s 采样一次

  const dfn::SimulationReport report = dfn::run_simulation(config);
  if (!report.ok) {
    std::cerr << report.message << "\n";
    return fail("CC-CV simulation failed");
  }
  // 1C 放电应能在 t_max 内降到 3.5 V 截止，否则说明容量或动力学异常。
  if (!report.cutoff_reached) {
    return fail("CC phase did not reach the voltage cutoff");
  }
  // 至少要有 10 个采样点，才能覆盖 CC→CV 过渡并支撑后续检查。
  if (report.samples.size() < 10) {
    return fail("too few output samples");
  }
  // 采样时间必须严格递增，防止重复或乱序输出。
  for (std::size_t k = 1; k < report.samples.size(); ++k) {
    if (report.samples[k].t <= report.samples[k - 1].t) {
      return fail("sample times are not increasing");
    }
  }
  // 首个采样应在 t=0、电压等于初始 OCP ≈ 3.98572 V；2 mV 容差覆盖
  // 初始一致化带来的微小偏移。
  if (!nearly(report.samples.front().t, 0.0, 1e-12) ||
      !nearly(report.samples.front().voltage, 3.98572, 2e-3)) {
    return fail("initial sample is wrong");
  }

  // CV 段电压是代数约束，应以 1e-6 V 的精度保持在 3.5 V。
  bool saw_cv = false;
  for (const dfn::SimulationSample& sample : report.samples) {
    if (sample.phase == "cv") {
      saw_cv = true;
      if (!nearly(sample.voltage, 3.5, 1e-6)) {
        return fail("CV phase does not hold the cutoff voltage");
      }
    }
  }
  if (!saw_cv) {
    return fail("CC-CV did not enter the CV phase");
  }
  // CV 结束时电流应衰减到 2.5 A/m^2 以下（远小于 1C 的 36 A/m^2），
  // 且电压仍为 3.5 V；否则说明恒压控制未生效。
  const dfn::SimulationSample& last = report.samples.back();
  if (std::abs(last.current) >= 2.5 || std::abs(last.voltage - 3.5) > 1e-6) {
    return fail("CV current did not decay at the cutoff voltage");
  }
  // 放电容量应为正且落在 1.2e5-1.45e5 C/m^2（约 33-40 Ah/m^2），
  // 对应 demo 电芯的可用容量；越界说明活性物质利用或计量错误。
  if (report.capacity < 1.2e5 || report.capacity > 1.45e5) {
    std::cerr << "capacity " << report.capacity << "\n";
    return fail("discharged capacity is outside the plausible range");
  }
  // 整个 CC-CV 过程中锂相对漂移应低于 1e-6，验证质量守恒型离散。
  if (report.lithium_drift > 1e-6) {
    return fail("lithium drifted during the CC-CV run");
  }

  // CSV 输出的表头必须固定，且行数 = 采样数 + 1（表头行）。
  std::ostringstream csv;
  dfn::write_csv(report, csv);
  std::size_t lines = 0;
  std::string first_line;
  std::istringstream reader(csv.str());
  std::string line;
  while (std::getline(reader, line)) {
    if (lines == 0) {
      first_line = line;
    }
    ++lines;
  }
  if (first_line !=
          "t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,phase" ||
      lines != report.samples.size() + 1) {
    return fail("CSV output does not match the sample list");
  }

  // 工况二：600 s 恒流，用于检查 SPMe/SPM 的未知量计数与模型选择。
  dfn::SimulationConfig reduced;
  reduced.mode = dfn::SimulationMode::cc;
  reduced.n_x_neg = 4;
  reduced.n_x_sep = 2;
  reduced.n_x_pos = 4;
  reduced.n_r_neg = 4;
  reduced.n_r_pos = 4;
  reduced.t_max = 600.0;
  reduced.output_interval = 300.0;

  // SPMe 未知量 = 两颗粒径向浓度 2*4 + 两个表面浓度 2
  //               + 电解液浓度 N_x_total=10 + 电解液电位 10 = 30。
  reduced.model = dfn::ModelKind::spme;
  const dfn::SimulationReport spme_report = dfn::run_simulation(reduced);
  if (!spme_report.ok || spme_report.model != "spme") {
    std::cerr << spme_report.message << "\n";
    return fail("SPMe simulation failed");
  }
  if (spme_report.unknowns != 4 + 4 + 2 + 10 + 10) {
    return fail("SPMe unknown count is wrong");
  }

  // SPM 没有 x 方向未知量：只有两颗粒浓度 2*4 与两个表面浓度 = 10。
  reduced.model = dfn::ModelKind::spm;
  const dfn::SimulationReport spm_report = dfn::run_simulation(reduced);
  if (!spm_report.ok || spm_report.model != "spm") {
    std::cerr << spm_report.message << "\n";
    return fail("SPM simulation failed");
  }
  if (spm_report.unknowns != 4 + 4 + 2) {
    return fail("SPM unknown count is wrong");
  }

  // SPMe 只实现了恒流控制，切换到 CCCV 必须被明确拒绝（错误信息含
  // "does not support"），而不是静默忽略电压约束。
  reduced.model = dfn::ModelKind::spme;
  reduced.mode = dfn::SimulationMode::cccv;
  const dfn::SimulationReport unsupported = dfn::run_simulation(reduced);
  if (unsupported.ok ||
      unsupported.message.find("does not support") == std::string::npos) {
    return fail("SPMe should reject voltage control");
  }

  // 工况三：脉冲放电。t_pulse = 600 s、占空比 0.5，总时长 1200 s，
  // 每 100 s 采样，应看到 pulse-on 与 pulse-rest 交替出现。
  dfn::SimulationConfig pulse;
  pulse.mode = dfn::SimulationMode::pulse;
  pulse.n_x_neg = 4;
  pulse.n_x_sep = 2;
  pulse.n_x_pos = 4;
  pulse.n_r_neg = 3;
  pulse.n_r_pos = 3;
  pulse.t_max = 1200.0;
  pulse.t_pulse = 600.0;
  pulse.pulse_duty = 0.5;
  pulse.output_interval = 100.0;
  const dfn::SimulationReport pulse_report = dfn::run_simulation(pulse);
  if (!pulse_report.ok) {
    std::cerr << pulse_report.message << "\n";
    return fail("pulse simulation failed");
  }
  // 通电段电流应精确为 1C = 36 A/m^2，静置段电流应为 0。
  bool saw_on = false;
  bool saw_rest = false;
  for (const dfn::SimulationSample& sample : pulse_report.samples) {
    if (sample.phase == "pulse-on" && nearly(sample.current, 36.0, 1e-9)) {
      saw_on = true;
    }
    if (sample.phase == "pulse-rest" && nearly(sample.current, 0.0, 1e-12)) {
      saw_rest = true;
    }
  }
  if (!saw_on || !saw_rest) {
    return fail("pulse phases do not alternate current and rest");
  }
  // 脉冲工况没有电压截止，应在 t_max 准时结束。
  if (!nearly(pulse_report.t_end, 1200.0, 1e-6)) {
    return fail("pulse run did not reach t_max");
  }

  // 工况四：充电 CC-CV，截止电压 4.0 V。充电电流按约定为负值。
  dfn::SimulationConfig charge;
  charge.mode = dfn::SimulationMode::cccv;
  charge.charge = true;
  charge.v_max = 4.0;
  charge.t_hold = 400.0;
  charge.output_interval = 100.0;
  const dfn::SimulationReport charge_report = dfn::run_simulation(charge);
  if (!charge_report.ok) {
    std::cerr << charge_report.message << "\n";
    return fail("charge CC-CV simulation failed");
  }
  if (!charge_report.cutoff_reached || charge_report.capacity <= 0.0) {
    return fail("charge phase did not reach the upper cutoff");
  }
  bool saw_charge_cv = false;
  for (const dfn::SimulationSample& sample : charge_report.samples) {
    // 充电方向电流不得为正，否则符号约定被破坏。
    if (sample.current > 0.0) {
      return fail("charge current must be negative");
    }
    if (sample.phase == "cv") {
      saw_charge_cv = true;
      if (!nearly(sample.voltage, 4.0, 1e-6)) {
        return fail("charge CV phase does not hold v_max");
      }
    }
  }
  if (!saw_charge_cv) {
    return fail("charge CC-CV did not enter the CV phase");
  }
  // 以第 3 个采样（仍处 CC 段）为基准，CV 末端的电流幅值必须更小。
  if (std::abs(charge_report.samples.back().current) >=
      std::abs(charge_report.samples[2].current)) {
    return fail("charge CV current did not decay");
  }

  std::cout << "dfn sim ok: t_event " << report.t_end << " s, capacity "
            << report.capacity / 3600.0 << " Ah/m^2, CV end " << last.current
            << " A/m^2, pulse t_end " << pulse_report.t_end
            << " s, charge capacity " << charge_report.capacity / 3600.0
            << " Ah/m^2, charge CV end "
            << charge_report.samples.back().current << " A/m^2\n";
  return 0;
#endif
}
