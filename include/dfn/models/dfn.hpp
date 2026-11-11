#pragma once

#include <vector>

#include "dfn/grid.hpp"
#include "dfn/models/cell_model.hpp"
#include "dfn/models/reaction.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"
#include "dfn/state.hpp"

namespace dfn {

// 完整 P2D/DFN 模型：x 方向为负极-隔膜-正极，r 方向为球形活性颗粒，
// 耦合固相扩散、电解液扩散、两相电势与 Butler-Volmer 反应。

// x 方向的控制方式：给定电流密度或给定端电压。
enum class ControlMode { current, voltage };

// DFN 的外部控制选项。
struct DfnOptions {
  ControlMode mode = ControlMode::current;
  double applied_current_density = 0.0;  // A/m^2，正值表示放电
  double applied_voltage = 4.0;          // V，仅电压模式使用
};

// 等温 DFN 模型（第 3 章的方程组）。浓溶液修正项 kappa_D 取零，
// 理由见第 9 章。
class DfnModel final : public CellModel {
 public:
  // 由参数、状态布局与控制方式构造模型，内部生成网格、装配算子并建立
  // Jacobian 稀疏模式。
  DfnModel(Params params, StateLayout layout, DfnOptions options = {});

  int size() const override { return size_; }
  // 残差 F(t, y, y')，写入 out；各方程已乘以 eq_scale_* 做量纲缩放。
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  // 变量类型：浓度为 1（微分），电势为 0（代数）。
  std::vector<int> variable_types() const override;
  // 组装 dF/dy + cj * dF/dy'（cj 由 IDA/BDF 给出）。
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

  const Params& params() const override { return params_; }
  const char* name() const override { return "dfn"; }
  bool supports_voltage_control() const override { return true; }
  const StateLayout& layout() const { return layout_; }
  // x 方向分段网格（负极、隔膜、正极连续编号）。
  const Grid1D& x_grid() const { return x_grid_; }

  // 辅助代数未知量（颗粒表面浓度）的全局下标：每个 x 位置一个表面节点，
  // 排在 StateLayout::size() 之后的附加段中。
  int cs_surf_neg(int ix) const { return cs_surf_neg_offset_ + ix; }
  int cs_surf_pos(int ix) const { return cs_surf_pos_offset_ + ix; }

  ControlMode mode() const { return options_.mode; }
  double applied_current_density() const {
    return options_.applied_current_density;
  }
  double applied_voltage() const { return options_.applied_voltage; }
  void set_current_density(double value) override;
  void set_applied_voltage(double value) override;
  // 同时设置控制方式与目标值（电流 A/m^2 或电压 V）。
  void set_control(ControlMode mode, double value);

  std::vector<double> initial_state() const override;
  std::vector<double> suggested_atol() const override;

  // 端电压，单位 V：电压模式下直接返回给定电压，否则由正极集流体处的
  // 固相电势外推得到。
  double terminal_voltage(const double* y) const override;
  // 流过电池的总电流密度，A/m^2，正值表示放电：电流模式直接返回给定值，
  // 电压模式下由集流体处的电势梯度反算。
  double terminal_current(const double* y) const override;

  // 全电池锂总量（固相 + 电解液），mol/m^2。
  double total_lithium(const double* y) const override;
  double solid_lithium(const double* y) const;
  double electrolyte_lithium(const double* y) const;
  // 某电极的体积平均化学计量数 c_s / c_max。
  double average_stoichiometry(const double* y, bool positive) const override;
  // 某电极上反应电流的积分，单位 A；理论上负极应为 +I、正极为 -I。
  double electrode_current(const double* y, bool positive) const;

 private:
  // 把区域内局部 x 下标映射为全局 x 下标（隔膜计入电解液段）。
  int global_index(int local_ix, bool positive) const;
  void build_pattern();

  Params params_;
  StateLayout layout_;
  DfnOptions options_;

  Grid1D x_grid_;      // 整条 x 方向分段网格
  Grid1D sphere_neg_;  // 负极颗粒径向网格
  Grid1D sphere_pos_;  // 正极颗粒径向网格
  Grid1D neg_grid_;    // 负极单独的 x 网格（用于固相电势方程）
  Grid1D pos_grid_;    // 正极单独的 x 网格

  // 装配好的离散算子：径向固相扩散、固相电势、电解液扩散与电势。
  DiffusionOperator radial_neg_;
  DiffusionOperator radial_pos_;
  TridiagonalMatrix phis_neg_op_;
  TridiagonalMatrix phis_pos_op_;
  DiffusionOperator elyte_op_;
  TridiagonalMatrix phie_op_;

  std::vector<double> eps_e_;      // 每个 x 节点的电解液孔隙率
  std::vector<double> d_eff_;      // 每个 x 节点的有效扩散系数
  std::vector<double> kappa_eff_;  // 每个 x 节点的有效离子电导率
  std::vector<double> a_s_;        // 每个 x 节点的比表面积，隔膜内为 0

  double sigma_neg_eff_ = 0.0;   // 负极固相有效电导率，S/m
  double sigma_pos_eff_ = 0.0;   // 正极固相有效电导率，S/m
  double dx_neg_ = 0.0;          // 负极单元宽度，m
  double dx_pos_ = 0.0;          // 正极单元宽度，m
  double neg_dirichlet_ = 0.0;   // 负极集流体边界系数 -2 sigma / dx^2
  double pos_dirichlet_ = 0.0;   // 正极集流体边界系数 -2 sigma / dx^2
  double pos_flux_coeff_ = 0.0;  // 正极集流体处通量系数 -1 / dx
  double surf_coeff_neg_ = 0.0;  // 表面通量源项系数 -A_surf / V_cell（负极）
  double surf_coeff_pos_ = 0.0;  // 同上（正极）
  double alpha_neg_ = 0.0;       // dr / (2 D_s)，由表面浓度重构表面通量
  double alpha_pos_ = 0.0;       // 同上（正极）

  // 方程缩放因子：残差乘以特征时间/浓度，使各方程量级接近、改善条件数。
  double eq_scale_solid_neg_ = 1.0;
  double eq_scale_solid_pos_ = 1.0;
  double eq_scale_electrolyte_ = 1.0;
  double eq_scale_potential_ = 1.0;
  double eq_scale_surface_ = 1.0;

  int size_ = 0;                // 未知量总数（含表面浓度段）
  int cs_surf_neg_offset_ = 0;  // 负极表面浓度段起始下标
  int cs_surf_pos_offset_ = 0;  // 正极表面浓度段起始下标

  // 每个 x 节点的 Butler-Volmer 结果缓存；需在 const 的 residual() 中写入，
  // 故声明为 mutable。
  mutable std::vector<ReactionValue> neg_reaction_;
  mutable std::vector<ReactionValue> pos_reaction_;

  CscPattern pattern_;  // Jacobian 稀疏模式
};

}  // namespace dfn
