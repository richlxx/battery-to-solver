// 单颗粒模型（SPM）实现：每个电极只用一个代表性球形颗粒，
// 反应电流沿电极厚度均匀分布，忽略电解液与固相电位的空间分布。
// 状态量只有两电极球网格上的固相锂浓度，外加两个表面浓度代数变量。
#include "dfn/models/spm.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "dfn/models/reaction.hpp"
#include "dfn/thermo.hpp"

namespace dfn {

namespace {

// 把固相浓度限制在 (1e-6 c_max, (1-1e-6) c_max)，
// 以免 OCV 拟合函数在化学计量比端点处斜率发散。
double clamp_stoich(double cs, double cmax) {
  return std::clamp(cs, 1e-6 * cmax, cmax * (1.0 - 1e-6));
}

}  // namespace

// 构造模型：校验参数与径向节点数，建立两电极的球网格和径向扩散算子，
// 并预计算表面通量系数、表面浓度外推系数与残差缩放因子。
SpmModel::SpmModel(Params params, int n_r_neg, int n_r_pos)
    : params_(std::move(params)), n_r_neg_(n_r_neg), n_r_pos_(n_r_pos) {
  validate(params_);
  if (n_r_neg <= 0 || n_r_pos <= 0) {
    throw std::invalid_argument("SPM radial counts must be positive");
  }

  // 球网格：face 存放半径，volume 为球壳体积，area[i] = 4 pi r_i^2。
  sphere_neg_ = make_uniform_sphere_grid(n_r_neg_, params_.neg.radius);
  sphere_pos_ = make_uniform_sphere_grid(n_r_pos_, params_.pos.radius);
  // 径向固相扩散算子 div(D_s grad c_s)，整个颗粒取同一 D_s。
  radial_neg_ = assemble_diffusion(
      sphere_neg_, std::vector<double>(static_cast<std::size_t>(n_r_neg_),
                                       params_.neg.d_s));
  radial_pos_ = assemble_diffusion(
      sphere_pos_, std::vector<double>(static_cast<std::size_t>(n_r_pos_),
                                       params_.pos.d_s));

  // 表面通量系数 = -A_surf / V_outer：把颗粒表面的摩尔通量
  // [mol/(m^2 s)] 折算成最外层球壳上的体积源项；负号表示通量向外流走。
  surf_coeff_neg_ = -sphere_neg_.area[static_cast<std::size_t>(n_r_neg_)] /
                    sphere_neg_.volume[static_cast<std::size_t>(n_r_neg_ - 1)];
  surf_coeff_pos_ = -sphere_pos_.area[static_cast<std::size_t>(n_r_pos_)] /
                    sphere_pos_.volume[static_cast<std::size_t>(n_r_pos_ - 1)];

  // dr/(2 D_s)：由最外层节点浓度线性外推到颗粒表面的斜率系数，
  // 因为节点位于球壳中心，距表面还有半个网格。
  const double dr_neg = params_.neg.radius / static_cast<double>(n_r_neg_);
  const double dr_pos = params_.pos.radius / static_cast<double>(n_r_pos_);
  alpha_neg_ = dr_neg / (2.0 * params_.neg.d_s);
  alpha_pos_ = dr_pos / (2.0 * params_.pos.d_s);

  // 残差缩放取固相扩散时间尺度 tau_s = R^2 / D_s 再除以 c_max，
  // 使残差的量级与状态量一致，改善求解器的容差判据。
  eq_scale_solid_neg_ = params_.neg.radius * params_.neg.radius /
                        params_.neg.d_s / params_.neg.c_max;
  eq_scale_solid_pos_ = params_.pos.radius * params_.pos.radius /
                        params_.pos.d_s / params_.pos.c_max;
  // 表面浓度行只用 1/c_max 缩放。
  eq_scale_surface_neg_ = 1.0 / params_.neg.c_max;
  eq_scale_surface_pos_ = 1.0 / params_.pos.c_max;

  // 状态排列：[c_s,neg (n_r_neg 个) | c_s,pos (n_r_pos 个) |
  //            c_surf,neg | c_surf,pos]，总长 n_r_neg + n_r_pos + 2。
  surf_neg_ = n_r_neg_ + n_r_pos_;
  surf_pos_ = surf_neg_ + 1;
  size_ = surf_pos_ + 1;

  build_pattern();
}

// 均匀反应通量 j [mol/(m^2 s)]：外加电流密度 I [A/m^2] 除以
// 比表面积 a_s、法拉第常数 F 与电极厚度 L，即界面摩尔通量。
// 符号：负极为 +I/(a_s F L)、正极为 -I/(a_s F L)；j > 0 表示锂离开颗粒，
// 所以放电（I > 0）时负极脱锂、正极嵌锂。
double SpmModel::reaction_flux(bool positive) const {
  const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
  // 比表面积 a_s = 3 (1 - eps) / R，单位 1/m。
  const double a_s = 3.0 * (1.0 - electrode.porosity) / electrode.radius;
  // 同一电流在两个电极上的界面通量方向相反。
  const double sign = positive ? -1.0 : 1.0;
  return sign * current_ / (a_s * kFaraday * electrode.thickness);
}

// 由表面化学计量比查半电池开路电位 U [V]，自变量先夹紧到物理区间。
double SpmModel::open_circuit(const ElectrodeParams& electrode, bool positive,
                              double cs_surf) const {
  const double cs = clamp_stoich(cs_surf, electrode.c_max);
  return positive ? positive_ocp(cs / electrode.c_max).value
                  : negative_ocp(cs / electrode.c_max).value;
}

// 在均匀通量下反解 Butler-Volmer 得到表面过电位 eta [V]。
// SPM 不解析电解液浓度梯度，故 ce 取参考浓度 c0。
double SpmModel::surface_overpotential(const ElectrodeParams& electrode,
                                       bool positive, double cs_surf) const {
  return overpotential_from_flux(params_, electrode, positive,
                                 params_.elyte.c0, cs_surf,
                                 reaction_flux(positive));
}

// 组装 SPM 残差 F(t, y, y') = 0，out 长度等于 size()。
// 固相扩散方程：dc_s/dt = div(D_s grad c_s)，反应通量只进入最外层球壳；
// 表面浓度代数式：c_surf = c_s(最外层) - alpha * j（一阶外推到 r = R）。
void SpmModel::residual(double /*t*/, const double* y, const double* yp,
                        double* out) const {
  std::fill(out, out + size_, 0.0);

  // 两电极的均匀反应通量，由当前外加电流决定。
  const double j_neg = reaction_flux(false);
  const double j_pos = reaction_flux(true);

  for (int ir = 0; ir < n_r_neg_; ++ir) {
    // 径向扩散算子乘以状态，得到 div(D_s grad c_s) 的第 ir 个分量。
    double lap = radial_neg_.matrix.diag[static_cast<std::size_t>(ir)] * y[ir];
    if (ir > 0) {
      lap += radial_neg_.matrix.lower[static_cast<std::size_t>(ir)] *
             y[cs_neg(ir - 1)];
    }
    if (ir + 1 < n_r_neg_) {
      lap += radial_neg_.matrix.upper[static_cast<std::size_t>(ir)] *
             y[cs_neg(ir + 1)];
    }
    // 只有最外层节点接收反应通量；负号已包含在 surf_coeff 中。
    const double source =
        (ir + 1 == n_r_neg_) ? surf_coeff_neg_ * j_neg : 0.0;
    out[cs_neg(ir)] = (yp[cs_neg(ir)] - lap - source) * eq_scale_solid_neg_;
  }

  for (int ir = 0; ir < n_r_pos_; ++ir) {
    double lap = radial_pos_.matrix.diag[static_cast<std::size_t>(ir)] *
                 y[cs_pos(ir)];
    if (ir > 0) {
      lap += radial_pos_.matrix.lower[static_cast<std::size_t>(ir)] *
             y[cs_pos(ir - 1)];
    }
    if (ir + 1 < n_r_pos_) {
      lap += radial_pos_.matrix.upper[static_cast<std::size_t>(ir)] *
             y[cs_pos(ir + 1)];
    }
    const double source =
        (ir + 1 == n_r_pos_) ? surf_coeff_pos_ * j_pos : 0.0;
    out[cs_pos(ir)] = (yp[cs_pos(ir)] - lap - source) * eq_scale_solid_pos_;
  }

  // 代数约束：表面浓度取最外层节点值并向 r = R 线性外推半个网格。
  out[surf_neg_] =
      (y[surf_neg_] - y[cs_neg(n_r_neg_ - 1)] + alpha_neg_ * j_neg) *
      eq_scale_surface_neg_;
  out[surf_pos_] =
      (y[surf_pos_] - y[cs_pos(n_r_pos_ - 1)] + alpha_pos_ * j_pos) *
      eq_scale_surface_pos_;
}

// 变量类型：1 表示微分变量（球网格上的固相浓度），
// 0 表示代数变量（两个表面浓度）。
std::vector<int> SpmModel::variable_types() const {
  std::vector<int> types(static_cast<std::size_t>(size_), 1);
  types[static_cast<std::size_t>(surf_neg_)] = 0;
  types[static_cast<std::size_t>(surf_pos_)] = 0;
  return types;
}

// 解析雅可比 dF/dy + cj * dF/dy'，cj 为求解器给出的时间步系数。
// 扩散部分对角元为 cj - lap 对角系数，邻居节点为 -lap 非对角元；
// 表面行分别是 +1 与 -1 乘以缩放因子。反应通量只依赖外加电流，故不含 y。
void SpmModel::jacobian(double /*t*/, const double* /*y*/, const double* /*yp*/,
                        double cj, CscBuilder& out) const {
  for (int ir = 0; ir < n_r_neg_; ++ir) {
    const int row = cs_neg(ir);
    out.add(row, row,
            (-radial_neg_.matrix.diag[static_cast<std::size_t>(ir)] + cj) *
                eq_scale_solid_neg_);
    if (ir > 0) {
      out.add(row, cs_neg(ir - 1),
              -radial_neg_.matrix.lower[static_cast<std::size_t>(ir)] *
                  eq_scale_solid_neg_);
    }
    if (ir + 1 < n_r_neg_) {
      out.add(row, cs_neg(ir + 1),
              -radial_neg_.matrix.upper[static_cast<std::size_t>(ir)] *
                  eq_scale_solid_neg_);
    }
  }
  for (int ir = 0; ir < n_r_pos_; ++ir) {
    const int row = cs_pos(ir);
    out.add(row, row,
            (-radial_pos_.matrix.diag[static_cast<std::size_t>(ir)] + cj) *
                eq_scale_solid_pos_);
    if (ir > 0) {
      out.add(row, cs_pos(ir - 1),
              -radial_pos_.matrix.lower[static_cast<std::size_t>(ir)] *
                  eq_scale_solid_pos_);
    }
    if (ir + 1 < n_r_pos_) {
      out.add(row, cs_pos(ir + 1),
              -radial_pos_.matrix.upper[static_cast<std::size_t>(ir)] *
                  eq_scale_solid_pos_);
    }
  }

  // 表面浓度行对自身取 +1、对最外层节点取 -1，均带缩放。
  out.add(surf_neg_, surf_neg_, eq_scale_surface_neg_);
  out.add(surf_neg_, cs_neg(n_r_neg_ - 1), -eq_scale_surface_neg_);
  out.add(surf_pos_, surf_pos_, eq_scale_surface_pos_);
  out.add(surf_pos_, cs_pos(n_r_pos_ - 1), -eq_scale_surface_pos_);
}

// 初值：颗粒内部浓度均匀，取设计浓度 c0；表面浓度与内部一致，
// 因此初始过电位为零，外推到表面的修正项也为零。
std::vector<double> SpmModel::initial_state() const {
  std::vector<double> y(static_cast<std::size_t>(size_), 0.0);
  for (int ir = 0; ir < n_r_neg_; ++ir) {
    y[static_cast<std::size_t>(cs_neg(ir))] = params_.neg.c0;
  }
  for (int ir = 0; ir < n_r_pos_; ++ir) {
    y[static_cast<std::size_t>(cs_pos(ir))] = params_.pos.c0;
  }
  y[static_cast<std::size_t>(surf_neg_)] = params_.neg.c0;
  y[static_cast<std::size_t>(surf_pos_)] = params_.pos.c0;
  return y;
}

// 绝对容差：固相浓度量级在 1e4~1e5 mol/m^3，统一取相对宽松的 1e-2。
std::vector<double> SpmModel::suggested_atol() const {
  std::vector<double> atol(static_cast<std::size_t>(size_), 1e-2);
  atol[static_cast<std::size_t>(surf_neg_)] = 1e-2;
  atol[static_cast<std::size_t>(surf_pos_)] = 1e-2;
  return atol;
}

// 端电压 V = U_pos - U_neg + eta_pos - eta_neg。
// 放电时负极脱锂（eta_neg > 0）、正极嵌锂（eta_pos < 0），故 V 低于 OCV 差。
double SpmModel::terminal_voltage(const double* y) const {
  const double cs_neg_surf = y[surf_neg_];
  const double cs_pos_surf = y[surf_pos_];
  const double u_neg = open_circuit(params_.neg, false, cs_neg_surf);
  const double u_pos = open_circuit(params_.pos, true, cs_pos_surf);
  const double eta_neg =
      surface_overpotential(params_.neg, false, cs_neg_surf);
  const double eta_pos = surface_overpotential(params_.pos, true, cs_pos_surf);
  return u_pos - u_neg + eta_pos - eta_neg;
}

// 单个颗粒的锂含量 [mol/m^2 电极面积]：按球壳体积加权求平均浓度，
// 再乘以固相体积分数 (1 - eps) 与电极厚度，得到面容量意义下的锂量。
double SpmModel::particle_lithium(const double* y, bool positive) const {
  const int nr = positive ? n_r_pos_ : n_r_neg_;
  const Grid1D& sphere = positive ? sphere_pos_ : sphere_neg_;
  const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;

  double particle_volume = 0.0;
  double weighted = 0.0;
  for (int ir = 0; ir < nr; ++ir) {
    const double volume = sphere.volume[static_cast<std::size_t>(ir)];
    particle_volume += volume;
    weighted += y[positive ? cs_pos(ir) : cs_neg(ir)] * volume;
  }
  // 体平均浓度为 sum(c_i V_i) / sum(V_i)，各球壳体积自行归一。
  return (weighted / particle_volume) * (1.0 - electrode.porosity) *
         electrode.thickness;
}

// 电极平均化学计量比 = 平均固相浓度 / c_max，用于检查锂守恒。
double SpmModel::average_stoichiometry(const double* y, bool positive) const {
  const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
  return particle_lithium(y, positive) /
         ((1.0 - electrode.porosity) * electrode.thickness * electrode.c_max);
}

// 登记雅可比稀疏结构：两段径向三对角链（含对角），
// 表面行只与自身和最外层节点耦合。
void SpmModel::build_pattern() {
  std::vector<std::pair<int, int>> entries;
  const auto add = [&entries](int row, int col) {
    entries.emplace_back(row, col);
  };

  for (int ir = 0; ir < n_r_neg_; ++ir) {
    add(cs_neg(ir), cs_neg(ir));
    if (ir > 0) {
      add(cs_neg(ir), cs_neg(ir - 1));
    }
    if (ir + 1 < n_r_neg_) {
      add(cs_neg(ir), cs_neg(ir + 1));
    }
  }
  for (int ir = 0; ir < n_r_pos_; ++ir) {
    add(cs_pos(ir), cs_pos(ir));
    if (ir > 0) {
      add(cs_pos(ir), cs_pos(ir - 1));
    }
    if (ir + 1 < n_r_pos_) {
      add(cs_pos(ir), cs_pos(ir + 1));
    }
  }
  add(surf_neg_, surf_neg_);
  add(surf_neg_, cs_neg(n_r_neg_ - 1));
  add(surf_pos_, surf_pos_);
  add(surf_pos_, cs_pos(n_r_pos_ - 1));

  pattern_ = ::dfn::build_pattern(size_, std::move(entries));
}

}  // namespace dfn
