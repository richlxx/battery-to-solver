#pragma once

#include <vector>

#include "dfn/grid.hpp"
#include "dfn/models/cell_model.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"

namespace dfn {

// 单颗粒模型（SPM）：每个电极只有一个球形颗粒，反应电流沿 x 方向均匀，
// 忽略电解液。没有 x 网格，规模最小，是 SPM -> SPMe -> DFN 阶梯的起点。

// SPM 模型：未知量为两个颗粒的径向浓度和两个表面浓度（代数约束）。
// 端电压由表面浓度算出：U_pos - U_neg + eta_pos - eta_neg。
class SpmModel final : public CellModel {
 public:
  // nr_neg、nr_pos 分别为负、正极颗粒的径向控制体数。
  SpmModel(Params params, int n_r_neg, int n_r_pos);

  int size() const override { return size_; }
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  // 变量类型：径向浓度为 1，表面浓度为 0（代数）。
  std::vector<int> variable_types() const override;
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

  const Params& params() const override { return params_; }
  const char* name() const override { return "spm"; }

  // 设置外加电流密度，A/m^2（正值表示放电）。
  void set_current_density(double value) override { current_ = value; }
  double current_density() const { return current_; }
  double terminal_current(const double* /*y*/) const override {
    return current_;
  }

  std::vector<double> initial_state() const override;
  std::vector<double> suggested_atol() const override;

  // 端电压（V），对应 set_current_density 设定的电流。
  double terminal_voltage(const double* y) const override;
  // 单颗粒中的固相锂量，mol（每 m^2 电极面积）。
  double particle_lithium(const double* y, bool positive) const;
  double average_stoichiometry(const double* y, bool positive) const override;
  double total_lithium(const double* y) const override {
    return particle_lithium(y, false) + particle_lithium(y, true);
  }

  // 全局未知量下标：负/正极径向浓度连续存放，表面浓度排在最后。
  int cs_neg(int ir) const { return ir; }
  int cs_pos(int ir) const { return n_r_neg_ + ir; }
  int cs_surf_neg() const { return surf_neg_; }
  int cs_surf_pos() const { return surf_pos_; }

  // 某电极的均匀反应通量，mol/(m^2 s)，正值表示锂离开颗粒。
  double reaction_flux(bool positive) const;

 private:
  // 表面过电位 eta 与开路电位 U，均由表面化学计量数算出，单位 V。
  double surface_overpotential(const ElectrodeParams& electrode, bool positive,
                               double cs_surf) const;
  double open_circuit(const ElectrodeParams& electrode, bool positive,
                      double cs_surf) const;
  void build_pattern();

  Params params_;
  int n_r_neg_ = 0;   // 负极径向控制体数
  int n_r_pos_ = 0;   // 正极径向控制体数
  int surf_neg_ = 0;  // 负极表面浓度下标
  int surf_pos_ = 0;  // 正极表面浓度下标
  int size_ = 0;      // 未知量总数
  double current_ = 0.0;  // 外加电流密度，A/m^2

  Grid1D sphere_neg_;  // 负极颗粒径向网格
  Grid1D sphere_pos_;  // 正极颗粒径向网格
  DiffusionOperator radial_neg_;  // 负极径向扩散算子
  DiffusionOperator radial_pos_;  // 正极径向扩散算子

  double surf_coeff_neg_ = 0.0;  // 表面通量源项系数 -A_surf / V_cell
  double surf_coeff_pos_ = 0.0;
  double alpha_neg_ = 0.0;  // dr / (2 D_s)，由表面浓度重构表面通量
  double alpha_pos_ = 0.0;
  double eq_scale_solid_neg_ = 1.0;    // 固相浓度方程缩放
  double eq_scale_solid_pos_ = 1.0;
  double eq_scale_surface_neg_ = 1.0;  // 表面浓度代数方程缩放
  double eq_scale_surface_pos_ = 1.0;

  CscPattern pattern_;  // Jacobian 稀疏模式
};

}  // namespace dfn
