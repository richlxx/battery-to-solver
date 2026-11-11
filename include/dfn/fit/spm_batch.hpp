#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "dfn/grid.hpp"
#include "dfn/linalg/tridiagonal.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/thermo.hpp"

namespace dfn {
namespace fit {

// 定步长 SPM 批量正演：供参数拟合实验与 CUDA 参考实现对照使用。
// 每个电极一个颗粒，无 x 网格，时间推进为固定步长的后向欧拉。

inline constexpr double kFaraday = 96485.33212;      // 法拉第常数，C/mol
inline constexpr double kGasConstant = 8.314462618;  // 气体常数，J/(mol K)

// 一段恒流工况。
struct BatchSegment {
  double duration = 0.0;  // 持续时间，s
  double current = 0.0;   // A/m^2，正值表示放电
};

// 一行采样：时刻、该区间电流与端电压，用作拟合的目标数据。
struct BatchRow {
  double t = 0.0;        // 时间，s
  double current = 0.0;  // 电流密度，A/m^2
  double voltage = 0.0;  // 端电压，V
};

// 把浓度限制在 (0, c_max) 开区间内，避免幂运算出现 0 或负底数。
inline double clamp_stoich(double cs, double cmax) {
  return std::clamp(cs, 1e-6 * cmax, cmax * (1.0 - 1e-6));
}

// 由 Butler-Volmer 关系用牛顿迭代反解过电位 eta（V）。j 为反应通量
// mol/(m^2 s)，cs_surf 为表面固相浓度，c_e 为电解液浓度。
inline double surface_overpotential(const ElectrodeParams& electrode,
                                    double j, double cs_surf, double c_e,
                                    double temperature) {
  if (j == 0.0) {
    return 0.0;
  }
  const double cs = clamp_stoich(cs_surf, electrode.c_max);
  // 交换电流密度 i0，单位 A/m^2。
  const double i0 =
      kFaraday * electrode.k * std::pow(c_e, electrode.alpha_a) *
      std::pow(electrode.c_max - cs, electrode.alpha_a) *
      std::pow(cs, electrode.alpha_c);
  if (i0 <= 0.0) {
    return 0.0;
  }
  const double a = electrode.alpha_a * kFaraday / (kGasConstant * temperature);
  const double c = electrode.alpha_c * kFaraday / (kGasConstant * temperature);
  double eta = 0.0;
  // 牛顿迭代；指数项限幅在 +-60，防止溢出。
  for (int iteration = 0; iteration < 60; ++iteration) {
    const double ea = std::exp(std::clamp(a * eta, -60.0, 60.0));
    const double ec = std::exp(std::clamp(-c * eta, -60.0, 60.0));
    const double residual = i0 / kFaraday * (ea - ec) - j;
    const double derivative = i0 / kFaraday * (a * ea + c * ec);
    const double step = residual / derivative;
    eta -= step;
    if (std::abs(step) < 1e-12) {
      break;
    }
  }
  return eta;
}

// 单个电极预装配的离散系数：径向扩散矩阵、后向欧拉隐式矩阵，
// 以及表面通量与电流之间的换算系数。
struct ElectrodeGeometry {
  TridiagonalMatrix diffusion;  // 球坐标径向扩散算子（无 dt）
  TridiagonalMatrix implicit;   // I - dt * diffusion，后向欧拉左端
  double surf_coeff = 0.0;      // 表面通量源项系数，-A_surf / V_cell
  double alpha = 0.0;           // dr / (2 D_s)，由表面浓度重构表面通量
  double j_coeff = 0.0;         // 由电流密度换算为反应通量的系数
};

// 为单个电极装配几何系数；nr 为径向单元数，dt 为时间步长。
inline ElectrodeGeometry make_geometry(const ElectrodeParams& electrode, int nr,
                                       double dt) {
  const Grid1D grid = make_uniform_sphere_grid(nr, electrode.radius);
  const DiffusionOperator op = assemble_diffusion(
      grid, std::vector<double>(static_cast<std::size_t>(nr), electrode.d_s));
  ElectrodeGeometry geometry;
  geometry.diffusion = op.matrix;
  geometry.surf_coeff =
      -grid.area[static_cast<std::size_t>(nr)] /
      grid.volume[static_cast<std::size_t>(nr - 1)];
  const double dr = electrode.radius / static_cast<double>(nr);
  geometry.alpha = dr / (2.0 * electrode.d_s);
  const double a_s = 3.0 * (1.0 - electrode.porosity) / electrode.radius;
  geometry.j_coeff = 1.0 / (a_s * kFaraday * electrode.thickness);

  geometry.implicit.n = nr;
  geometry.implicit.lower.resize(static_cast<std::size_t>(nr));
  geometry.implicit.diag.resize(static_cast<std::size_t>(nr));
  geometry.implicit.upper.resize(static_cast<std::size_t>(nr));
  for (int i = 0; i < nr; ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    geometry.implicit.lower[k] = -dt * geometry.diffusion.lower[k];
    geometry.implicit.diag[k] = 1.0 - dt * geometry.diffusion.diag[k];
    geometry.implicit.upper[k] = -dt * geometry.diffusion.upper[k];
  }
  return geometry;
}

// 定步长后向欧拉 SPM，供参数拟合实验与 CUDA 参考实现使用。
// SPM 没有 x 网格，每个电极只有一个球形颗粒。
class SpmBatch {
 public:
  // params 为电池参数，nr 为径向单元数，dt 为时间步长（s）。
  SpmBatch(Params params, int nr, double dt)
      : params_(std::move(params)),
        nr_(nr),
        dt_(dt),
        neg_(make_geometry(params_.neg, nr, dt)),
        pos_(make_geometry(params_.pos, nr, dt)),
        c_neg_(static_cast<std::size_t>(nr), params_.neg.c0),
        c_pos_(static_cast<std::size_t>(nr), params_.pos.c0) {}

  // 把两个颗粒的浓度恢复到初始值 c0。
  void reset() {
    std::fill(c_neg_.begin(), c_neg_.end(), params_.neg.c0);
    std::fill(c_pos_.begin(), c_pos_.end(), params_.pos.c0);
  }

  // 以电流密度 current（A/m^2，正值为放电）推进一个时间步，返回步末电压。
  double step(double current) {
    const double j_neg = current * neg_.j_coeff;
    const double j_pos = -current * pos_.j_coeff;

    std::vector<double> rhs = c_neg_;
    rhs[static_cast<std::size_t>(nr_ - 1)] += dt_ * neg_.surf_coeff * j_neg;
    c_neg_ = solve_tridiagonal(neg_.implicit, rhs);

    rhs = c_pos_;
    rhs[static_cast<std::size_t>(nr_ - 1)] += dt_ * pos_.surf_coeff * j_pos;
    c_pos_ = solve_tridiagonal(pos_.implicit, rhs);
    return voltage(current);
  }

  // 给定电流密度时的端电压：U_pos - U_neg + eta_pos - eta_neg。
  double voltage(double current) const {
    const double j_neg = current * neg_.j_coeff;
    const double j_pos = -current * pos_.j_coeff;
    const double cs_neg = clamp_stoich(
        c_neg_[static_cast<std::size_t>(nr_ - 1)] - neg_.alpha * j_neg,
        params_.neg.c_max);
    const double cs_pos = clamp_stoich(
        c_pos_[static_cast<std::size_t>(nr_ - 1)] - pos_.alpha * j_pos,
        params_.pos.c_max);
    const double u_neg = negative_ocp(cs_neg / params_.neg.c_max).value;
    const double u_pos = positive_ocp(cs_pos / params_.pos.c_max).value;
    const double eta_neg =
        surface_overpotential(params_.neg, j_neg, cs_neg, params_.elyte.c0,
                              params_.temperature);
    const double eta_pos =
        surface_overpotential(params_.pos, j_pos, cs_pos, params_.elyte.c0,
                              params_.temperature);
    return u_pos - u_neg + eta_pos - eta_neg;
  }

  // 由表面通量重构颗粒表面浓度（一阶近似），并限制在物理区间内。
  double surface(bool positive, double current) const {
    const ElectrodeGeometry& geometry = positive ? pos_ : neg_;
    const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
    const std::vector<double>& c = positive ? c_pos_ : c_neg_;
    const double sign = positive ? -1.0 : 1.0;
    const double j = sign * current * geometry.j_coeff;
    return clamp_stoich(
        c[static_cast<std::size_t>(nr_ - 1)] - geometry.alpha * j,
        electrode.c_max);
  }

  // 访问器：参数、径向单元数与时间步长。
  const Params& params() const { return params_; }
  int radial_cells() const { return nr_; }
  double time_step() const { return dt_; }

 private:
  Params params_;
  int nr_ = 0;
  double dt_ = 0.0;
  ElectrodeGeometry neg_;
  ElectrodeGeometry pos_;
  std::vector<double> c_neg_;
  std::vector<double> c_pos_;
};

// 按分段工况推进模型，生成带时间、电流与电压的采样行，作为拟合数据。
// 段长按 dt 取整；段末那一行已经改用下一段的电流计算电压。
inline std::vector<BatchRow> rows_from_protocol(
    const Params& params, const std::vector<BatchSegment>& segments, int nr,
    double dt) {
  if (segments.empty()) {
    throw std::invalid_argument("protocol needs at least one segment");
  }
  SpmBatch model(params, nr, dt);
  std::vector<BatchRow> rows;
  rows.push_back(BatchRow{0.0, segments.front().current,
                          model.voltage(segments.front().current)});
  double t = 0.0;
  for (std::size_t s = 0; s < segments.size(); ++s) {
    const BatchSegment& segment = segments[s];
    const int steps = static_cast<int>(std::lround(segment.duration / dt));
    for (int k = 0; k < steps; ++k) {
      model.step(segment.current);
      t += dt;
      const bool boundary = k + 1 == steps && s + 1 < segments.size();
      const double next =
          boundary ? segments[s + 1].current : segment.current;
      rows.push_back(BatchRow{t, next, model.voltage(next)});
    }
  }
  return rows;
}

// 在给定的采样行（要求时间等距）上推进模型，返回各时刻的电压。
// 常用于检查拟合结果与参考数据的偏差。
inline std::vector<double> simulate_rows(const Params& params,
                                         const std::vector<BatchRow>& rows,
                                         int nr, double dt) {
  SpmBatch model(params, nr, dt);
  std::vector<double> voltages(rows.size(), 0.0);
  voltages[0] = model.voltage(rows[0].current);
  for (std::size_t i = 1; i < rows.size(); ++i) {
    const double expected = rows[0].t + dt * static_cast<double>(i);
    if (std::abs(rows[i].t - expected) > 1e-6) {
      throw std::runtime_error("samples are not on a uniform time grid");
    }
    // 区间 [t_{i-1}, t_i] 使用 t_{i-1} 处记录的电流推进；t_i 时刻的电压
    // 已经反映下一区间的电流。
    model.step(rows[i - 1].current);
    voltages[i] = model.voltage(rows[i].current);
  }
  return voltages;
}

}  // namespace fit
}  // namespace dfn
