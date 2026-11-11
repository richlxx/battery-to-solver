#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/models/spm.hpp"
#include "dfn/models/spme.hpp"
#include "dfn/params.hpp"

namespace dfn {

// 仿真驱动层：把模型选择、工况（CC/CV/CCCV/脉冲）、网格与容差组合成
// 一次完整仿真，输出电压/电流时间序列与求解器统计量。

// 可选模型层级：单颗粒、带电解液的单颗粒、完整 DFN。
enum class ModelKind { spm, spme, dfn };

// 工况类型：恒流、恒压、先恒流后恒压、脉冲（充/放与静置交替）。
enum class SimulationMode { cc, cv, cccv, pulse };

// 一次仿真的全部输入参数。
struct SimulationConfig {
  Params params = make_demo_params();
  ModelKind model = ModelKind::dfn;
  SimulationMode mode = SimulationMode::cccv;

  // 示例电芯约储存 36 Ah/m^2，因此 36 A/m^2 大致相当于 1C。
  // 电流密度单位为 A/m^2，正值表示放电；充电时施加电流取相反符号，
  // 并把 v_max 用作截止条件。
  double i_1c = 36.0;             // 1C 对应电流密度，A/m^2
  double c_rate = 1.0;            // 倍率，实际电流 = c_rate * i_1c
  double current_density = 0.0;   // 非零时覆盖 c_rate 给出的电流
  bool charge = false;            // true 表示充电

  double v_min = 3.5;             // 放电截止电压，V
  double v_max = 4.0;             // 充电截止电压，V
  double hold_voltage = 0.0;      // CV 保持电压；0 表示沿用截止电压

  double t_max = 20000.0;         // 仿真时间上限，s
  double t_hold = 600.0;          // CV 保持段最长时间，s
  double t_pulse = 600.0;         // 单个脉冲段时长，s
  double pulse_duty = 0.5;        // 脉冲占空比（通电时间占比）
  double output_interval = 60.0;  // 采样与输出间隔，s
  double rtol = 1e-6;             // 时间积分相对容差

  int n_x_neg = 8;  // 负极 x 方向控制体数
  int n_x_sep = 4;  // 隔膜 x 方向控制体数
  int n_x_pos = 8;  // 正极 x 方向控制体数
  int n_r_neg = 6;  // 负极颗粒径向控制体数
  int n_r_pos = 6;  // 正极颗粒径向控制体数
};

// 一个输出时刻的采样点。
struct SimulationSample {
  double t = 0.0;           // 时间，s
  double voltage = 0.0;     // 端电压，V
  double current = 0.0;     // 电流密度 A/m^2，正值表示放电
  double stoich_neg = 0.0;  // 负极平均化学计量数 c_s / c_max
  double stoich_pos = 0.0;  // 正极平均化学计量数
  std::string phase;        // 所处阶段标签，如 "cc"、"cv"、"pulse-on"
};

// 仿真结果汇总：采样时间序列加上结束原因与求解器统计量。
struct SimulationReport {
  bool ok = false;             // 是否正常完成
  std::string message;         // 结束/失败原因
  std::string model;           // 实际使用的模型名
  std::vector<SimulationSample> samples;  // 输出时间序列

  double t_end = 0.0;  // 仿真结束时刻，s
  // 容量，A s/m^2，由负极储锂量变化算出；与输出间隔无关。
  double capacity = 0.0;
  bool cutoff_reached = false;  // 是否触发电压截止
  double lithium_drift = 0.0;   // 锂守恒的相对漂移，用于数值检验

  int steps = 0;            // 时间步数
  int residual_evals = 0;   // 残差函数调用次数
  int jacobian_evals = 0;   // Jacobian 计算次数
  int linear_setups = 0;    // 线性求解器分解/重构次数
  int unknowns = 0;         // 未知量个数
  int nnz = 0;              // Jacobian 非零元个数
  double wall_ms = 0.0;     // 墙钟耗时，ms
};

// 按配置运行一次完整仿真（构造模型、一致化初值、推进到截止或 t_max）。
SimulationReport run_simulation(const SimulationConfig& config);

// 输出 CSV，列为 t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,phase。
void write_csv(const SimulationReport& report, std::ostream& out);
void write_csv(const SimulationReport& report, const std::string& path);

}  // namespace dfn
