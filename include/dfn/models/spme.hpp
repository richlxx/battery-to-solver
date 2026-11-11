#pragma once

#include <vector>

#include "dfn/grid.hpp"
#include "dfn/models/cell_model.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"

namespace dfn {

// 带电解液的单颗粒模型（SPMe）：每个电极一个颗粒且反应电流沿 x 均匀，
// 另加一维电解液浓度与电势，是 SPM 与 DFN 之间的中间层级。

// SPMe 模型。状态依次为 [c_s 负极, c_s 正极, c_surf 负极, c_surf 正极,
// c_e, phi_e]；phi_e 由负极集流体处的参考电位（gauge）固定。
class SpmeModel final : public CellModel {
 public:
  // n_x_* 为三个区域的 x 控制体数，n_r_* 为颗粒径向控制体数。
  SpmeModel(Params params, int n_x_neg, int n_x_sep, int n_x_pos,
            int n_r_neg, int n_r_pos);

  int size() const override { return size_; }
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  // 变量类型：浓度为 1，电解液电势为 0（代数）。
  std::vector<int> variable_types() const override;
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

  const Params& params() const override { return params_; }
  const char* name() const override { return "spme"; }

  // 设置外加电流密度，A/m^2（正值表示放电）。
  void set_current_density(double value) override;
  double terminal_current(const double* /*y*/) const override {
    return current_;
  }

  std::vector<double> initial_state() const override;
  std::vector<double> suggested_atol() const override;

  // 端电压（V）。
  double terminal_voltage(const double* y) const override;
  double average_stoichiometry(const double* y, bool positive) const override;
  // 全电池锂总量、固相锂量与电解液锂量，单位 mol/m^2。
  double total_lithium(const double* y) const override;
  double solid_lithium(const double* y) const;
  double electrolyte_lithium(const double* y) const;
  double particle_lithium(const double* y, bool positive) const;

  // 全局未知量下标：先放两颗粒浓度，再放两个表面浓度，最后是 c_e 与 phi_e。
  int cs_neg(int ir) const { return ir; }
  int cs_pos(int ir) const { return n_r_neg_ + ir; }
  int cs_surf_neg() const { return cs_pos_offset_ + n_r_pos_; }
  int cs_surf_pos() const { return cs_surf_neg() + 1; }
  int ce(int ix) const { return ce_offset_ + ix; }
  int phie(int ix) const { return phie_offset_ + ix; }

  // 某电极的均匀反应通量，mol/(m^2 s)，正值表示锂离开颗粒。
  double reaction_flux(bool positive) const {
    return positive ? j_pos_ : j_neg_;
  }

 private:
  // 电极的开路电位（V），由表面化学计量数算出。
  double open_circuit(const ElectrodeParams& electrode, bool positive,
                      double cs_surf) const;
  // 第 ix 个 x 节点上电解液浓度方程的源项（反应消耗）。
  double electrolyte_source(int ix) const;
  void build_pattern();

  Params params_;
  int n_x_neg_ = 0;    // 负极 x 控制体数
  int n_x_sep_ = 0;    // 隔膜 x 控制体数
  int n_x_pos_ = 0;    // 正极 x 控制体数
  int n_x_total_ = 0;  // x 方向节点总数
  int n_r_neg_ = 0;    // 负极径向控制体数
  int n_r_pos_ = 0;    // 正极径向控制体数

  int cs_pos_offset_ = 0;  // 正极浓度段起始下标
  int ce_offset_ = 0;      // 电解液浓度段起始下标
  int phie_offset_ = 0;    // 电解液电势段起始下标
  int size_ = 0;           // 未知量总数

  double current_ = 0.0;  // 外加电流密度，A/m^2
  double j_neg_ = 0.0;    // 负极均匀反应通量，mol/(m^2 s)
  double j_pos_ = 0.0;    // 正极均匀反应通量，mol/(m^2 s)
  double a_s_neg_ = 0.0;  // 负极比表面积，1/m
  double a_s_pos_ = 0.0;  // 正极比表面积，1/m

  Grid1D x_grid_;      // x 方向分段网格
  Grid1D sphere_neg_;  // 负极颗粒径向网格
  Grid1D sphere_pos_;  // 正极颗粒径向网格

  DiffusionOperator radial_neg_;  // 负极径向扩散算子
  DiffusionOperator radial_pos_;  // 正极径向扩散算子
  DiffusionOperator elyte_op_;    // 电解液扩散算子
  TridiagonalMatrix phie_op_;     // 电解液电势方程算子

  std::vector<double> eps_e_;  // 每个 x 节点的孔隙率
  std::vector<double> a_s_;    // 每个 x 节点的比表面积

  double surf_coeff_neg_ = 0.0;  // 表面通量源项系数 -A_surf / V_cell
  double surf_coeff_pos_ = 0.0;
  double alpha_neg_ = 0.0;  // dr / (2 D_s)，由表面浓度重构表面通量
  double alpha_pos_ = 0.0;

  // 方程缩放因子，使各类残差量级接近。
  double eq_scale_solid_neg_ = 1.0;
  double eq_scale_solid_pos_ = 1.0;
  double eq_scale_electrolyte_ = 1.0;
  double eq_scale_potential_ = 1.0;
  double eq_scale_surface_ = 1.0;

  CscPattern pattern_;  // Jacobian 稀疏模式
};

}  // namespace dfn
