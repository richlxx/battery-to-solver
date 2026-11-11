#pragma once

#include <stdexcept>
#include <vector>

#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"

namespace dfn {

// 模型层公共接口：SPM、SPMe、DFN 都实现 CellModel，仿真驱动只依赖
// 该接口，因此可以整体替换模型而不改动驱动代码。

// SPM -> SPMe -> DFN 阶梯模型的公共接口。每个模型自行管理状态布局，
// 对外只暴露端电压、端电流与守恒量等“观测量”。
class CellModel : public DaeModel {
 public:
  // t = 0 的状态向量，要求已满足代数约束（一致初值）。
  virtual std::vector<double> initial_state() const = 0;
  // 建议的绝对容差，逐分量给出，量级应与对应未知量匹配。
  virtual std::vector<double> suggested_atol() const = 0;

  // 端电压（V）与端电流密度（A/m^2，正值表示放电）。
  virtual double terminal_voltage(const double* y) const = 0;
  virtual double terminal_current(const double* y) const = 0;

  // 某电极的体积平均化学计量数 c_s / c_max，positive 为 true 表示正极。
  virtual double average_stoichiometry(const double* y, bool positive) const = 0;
  // 全电池锂总量（固相 + 电解液），用于守恒性检查。
  virtual double total_lithium(const double* y) const = 0;

  // 设定外部控制量：电流密度或电压。
  virtual void set_current_density(double value) = 0;
  virtual void set_applied_voltage(double) {
    throw std::logic_error("this model does not support voltage control");
  }
  virtual bool supports_voltage_control() const { return false; }

  virtual const Params& params() const = 0;
  // 模型名（"spm"、"spme"、"dfn"），用于报告与命令行输出。
  virtual const char* name() const = 0;
};

}  // namespace dfn
