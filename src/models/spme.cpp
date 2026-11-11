// 单颗粒模型加电解液（SPMe）实现：两电极各用一个代表性颗粒，反应电流
// 沿厚度均匀，另解析电解液浓度 c_e 与电位 phi_e 的一维分布。
// phi_e 在负极集流体处取零电位以固定规范；状态排列见构造函数说明。
#include "dfn/models/spme.hpp"

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

// 构造模型：校验参数与各方向节点数，建立 x 方向分段网格与两电极球网格，
// 逐节点计算孔隙率、有效扩散系数与有效电导率，并确定状态排列与缩放因子。
SpmeModel::SpmeModel(Params params, int n_x_neg, int n_x_sep, int n_x_pos,
                     int n_r_neg, int n_r_pos)
    : params_(std::move(params)),
      n_x_neg_(n_x_neg),
      n_x_sep_(n_x_sep),
      n_x_pos_(n_x_pos),
      n_r_neg_(n_r_neg),
      n_r_pos_(n_r_pos) {
  validate(params_);
  if (n_x_neg <= 0 || n_x_sep <= 0 || n_x_pos <= 0 || n_r_neg <= 0 ||
      n_r_pos <= 0) {
    throw std::invalid_argument("SPMe grid sizes must be positive");
  }
  n_x_total_ = n_x_neg_ + n_x_sep_ + n_x_pos_;

  const double l_neg = params_.neg.thickness;
  const double l_sep = params_.sep.thickness;
  const double l_pos = params_.pos.thickness;

  // x 方向网格：负极 | 隔膜 | 正极三段各自等距，块界面精确落在实测厚度上。
  x_grid_ = make_segmented_line_grid({n_x_neg_, n_x_sep_, n_x_pos_},
                                     {l_neg, l_sep, l_pos}, params_.area);
  // 球网格：face 为半径，volume 为球壳体积，area[i] = 4 pi r_i^2。
  sphere_neg_ = make_uniform_sphere_grid(n_r_neg_, params_.neg.radius);
  sphere_pos_ = make_uniform_sphere_grid(n_r_pos_, params_.pos.radius);

  // 固相径向扩散算子 div(D_s grad c_s)，两电极各一套。
  radial_neg_ = assemble_diffusion(
      sphere_neg_, std::vector<double>(static_cast<std::size_t>(n_r_neg_),
                                       params_.neg.d_s));
  radial_pos_ = assemble_diffusion(
      sphere_pos_, std::vector<double>(static_cast<std::size_t>(n_r_pos_),
                                       params_.pos.d_s));

  // 比表面积 a_s = 3 (1 - eps) / R，单位 1/m。
  a_s_neg_ = 3.0 * (1.0 - params_.neg.porosity) / params_.neg.radius;
  a_s_pos_ = 3.0 * (1.0 - params_.pos.porosity) / params_.pos.radius;

  eps_e_.assign(static_cast<std::size_t>(n_x_total_), 0.0);
  std::vector<double> d_eff(static_cast<std::size_t>(n_x_total_), 0.0);
  std::vector<double> kappa_eff(static_cast<std::size_t>(n_x_total_), 0.0);
  a_s_.assign(static_cast<std::size_t>(n_x_total_), 0.0);
  for (int i = 0; i < n_x_total_; ++i) {
    double epsy = 0.0;
    double bruggeman = 1.5;
    // 按节点所属区域（负极、隔膜、正极）取孔隙率与 Bruggeman 指数。
    if (i < n_x_neg_) {
      epsy = params_.neg.porosity;
      bruggeman = params_.neg.bruggeman;
      a_s_[static_cast<std::size_t>(i)] = a_s_neg_;
    } else if (i < n_x_neg_ + n_x_sep_) {
      epsy = params_.sep.porosity;
      bruggeman = params_.sep.bruggeman;
    } else {
      epsy = params_.pos.porosity;
      bruggeman = params_.pos.bruggeman;
      a_s_[static_cast<std::size_t>(i)] = a_s_pos_;
    }
    eps_e_[static_cast<std::size_t>(i)] = epsy;
    // Bruggeman 关系：D_eff = D eps^b，kappa_eff = kappa eps^b。
    d_eff[static_cast<std::size_t>(i)] =
        params_.elyte.d_e * std::pow(epsy, bruggeman);
    kappa_eff[static_cast<std::size_t>(i)] =
        params_.elyte.kappa * std::pow(epsy, bruggeman);
  }
  // 电解液扩散算子 div(D_eff grad c_e)；电位只用其三对角矩阵部分。
  elyte_op_ = assemble_diffusion(x_grid_, d_eff);
  phie_op_ = assemble_diffusion(x_grid_, kappa_eff).matrix;

  // 表面通量系数 = -A_surf / V_outer，把表面摩尔通量折算为最外层球壳的源项。
  surf_coeff_neg_ = -sphere_neg_.area[static_cast<std::size_t>(n_r_neg_)] /
                    sphere_neg_.volume[static_cast<std::size_t>(n_r_neg_ - 1)];
  surf_coeff_pos_ = -sphere_pos_.area[static_cast<std::size_t>(n_r_pos_)] /
                    sphere_pos_.volume[static_cast<std::size_t>(n_r_pos_ - 1)];
  // dr/(2 D_s)：由最外层节点浓度线性外推到颗粒表面的斜率系数。
  const double dr_neg = params_.neg.radius / static_cast<double>(n_r_neg_);
  const double dr_pos = params_.pos.radius / static_cast<double>(n_r_pos_);
  alpha_neg_ = dr_neg / (2.0 * params_.neg.d_s);
  alpha_pos_ = dr_pos / (2.0 * params_.pos.d_s);

  // 固相残差按 tau_s / c_max 缩放，tau_s = R^2 / D_s 为颗粒扩散时间尺度。
  const double tau_s_neg =
      params_.neg.radius * params_.neg.radius / params_.neg.d_s;
  const double tau_s_pos =
      params_.pos.radius * params_.pos.radius / params_.pos.d_s;
  eq_scale_solid_neg_ = tau_s_neg / params_.neg.c_max;
  eq_scale_solid_pos_ = tau_s_pos / params_.pos.c_max;
  // 电解液时间尺度 tau_e = L^2 / D_e，按参考浓度 c0 无量纲化。
  const double length = params_.cell_length();
  const double tau_e = length * length / params_.elyte.d_e;
  eq_scale_electrolyte_ =
      params_.elyte.c0 > 0.0 ? tau_e / params_.elyte.c0 : 1.0;
  // 电位残差按电池厚度 L 缩放，用于平衡各方程的残差量级。
  eq_scale_potential_ = length;
  eq_scale_surface_ =
      1.0 / std::max(params_.neg.c_max, params_.pos.c_max);

  // 状态排列：
  //   [c_s,neg | c_s,pos | c_surf,neg | c_surf,pos | c_e | phi_e]
  // 固相节点分别 n_r_neg、n_r_pos 个，c_e 与 phi_e 各 n_x_total 个，
  // 故总长度 = n_r_neg + n_r_pos + 2 + 2 n_x_total。
  cs_pos_offset_ = n_r_neg_;
  ce_offset_ = cs_surf_pos() + 1;
  phie_offset_ = ce_offset_ + n_x_total_;
  size_ = phie_offset_ + n_x_total_;

  // 先把电流置零，同时确定两侧的均匀反应通量。
  set_current_density(0.0);
  build_pattern();
}

// 设定外加电流密度 I [A/m^2]，并换算成两电极的均匀反应通量：
//   j_neg =  I / (a_s F L)，j_pos = -I / (a_s F L)，
// j > 0 表示锂离开颗粒；放电（I > 0）时负极脱锂、正极嵌锂。
void SpmeModel::set_current_density(double value) {
  current_ = value;
  j_neg_ = current_ / (a_s_neg_ * kFaraday * params_.neg.thickness);
  j_pos_ = -current_ / (a_s_pos_ * kFaraday * params_.pos.thickness);
}

// 由表面化学计量比查半电池开路电位 U [V]，自变量先夹紧到物理区间。
double SpmeModel::open_circuit(const ElectrodeParams& electrode, bool positive,
                               double cs_surf) const {
  const double cs = clamp_stoich(cs_surf, electrode.c_max);
  return positive ? positive_ocp(cs / electrode.c_max).value
                  : negative_ocp(cs / electrode.c_max).value;
}

// 组装 SPMe 残差 F(t, y, y') = 0，out 长度等于 size()。
// 方程组成：
//   1) 固相径向扩散 dc_s/dt = div(D_s grad c_s)，反应通量只作用在最外层球壳；
//   2) 表面浓度代数式 c_surf = c_s(最外层) - alpha * j；
//   3) 电解液质量守恒 eps dc_e/dt = div(D_eff grad c_e) + (1 - t_+) a_s j；
//   4) 电解液电荷守恒 div(kappa_eff grad phi_e) + a_s F j = 0，
//      并在 x = 0（负极集流体）处令 phi_e = 0 作为电位规范。
void SpmeModel::residual(double /*t*/, const double* y, const double* yp,
                         double* out) const {
  std::fill(out, out + size_, 0.0);

  for (int ir = 0; ir < n_r_neg_; ++ir) {
    // 径向扩散算子乘状态，得到 div(D_s grad c_s) 的第 ir 个分量。
    double lap = radial_neg_.matrix.diag[static_cast<std::size_t>(ir)] * y[ir];
    if (ir > 0) {
      lap += radial_neg_.matrix.lower[static_cast<std::size_t>(ir)] *
             y[cs_neg(ir - 1)];
    }
    if (ir + 1 < n_r_neg_) {
      lap += radial_neg_.matrix.upper[static_cast<std::size_t>(ir)] *
             y[cs_neg(ir + 1)];
    }
    // 只有最外层球壳接收反应通量；负号已包含在 surf_coeff 中。
    const double source =
        (ir + 1 == n_r_neg_) ? surf_coeff_neg_ * j_neg_ : 0.0;
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
        (ir + 1 == n_r_pos_) ? surf_coeff_pos_ * j_pos_ : 0.0;
    out[cs_pos(ir)] = (yp[cs_pos(ir)] - lap - source) * eq_scale_solid_pos_;
  }

  // 代数约束：表面浓度取最外层节点值，并向 r = R 线性外推半个网格。
  out[cs_surf_neg()] =
      (y[cs_surf_neg()] - y[cs_neg(n_r_neg_ - 1)] + alpha_neg_ * j_neg_) *
      eq_scale_surface_;
  out[cs_surf_pos()] =
      (y[cs_surf_pos()] - y[cs_pos(n_r_pos_ - 1)] + alpha_pos_ * j_pos_) *
      eq_scale_surface_;

  for (int i = 0; i < n_x_total_; ++i) {
    const int row = ce(i);
    // 电解液扩散：div(D_eff grad c_e)，两端为零通量（Neumann）边界。
    double lap = elyte_op_.matrix.diag[static_cast<std::size_t>(i)] * y[row];
    if (i > 0) {
      lap += elyte_op_.matrix.lower[static_cast<std::size_t>(i)] *
             y[ce(i - 1)];
    }
    if (i + 1 < n_x_total_) {
      lap += elyte_op_.matrix.upper[static_cast<std::size_t>(i)] *
             y[ce(i + 1)];
    }
    // 孔隙率加权的时间导数项，源项来自界面通量中的迁移部分。
    out[row] =
        (eps_e_[static_cast<std::size_t>(i)] * yp[row] - lap -
         electrolyte_source(i)) *
        eq_scale_electrolyte_;
  }

  // 电位规范：负极集流体处 phi_e = 0，只有这一个纯代数行。
  out[phie(0)] = y[phie(0)] * eq_scale_potential_;
  for (int i = 1; i < n_x_total_; ++i) {
    const int row = phie(i);
    // 欧姆定律代入电荷守恒：div(kappa_eff grad phi_e) + a_s F j = 0。
    double lap = phie_op_.diag[static_cast<std::size_t>(i)] * y[row];
    if (i > 0) {
      lap += phie_op_.lower[static_cast<std::size_t>(i)] * y[phie(i - 1)];
    }
    if (i + 1 < n_x_total_) {
      lap += phie_op_.upper[static_cast<std::size_t>(i)] * y[phie(i + 1)];
    }
    // 隔膜内 a_s = 0 因而无源项；负极节点用 j_neg，其余用 j_pos。
    const double source =
        a_s_[static_cast<std::size_t>(i)] * kFaraday *
        (i < n_x_neg_ ? j_neg_ : j_pos_);
    out[row] = (lap + source) * eq_scale_potential_;
  }
}

// 电解液浓度源项 (1 - t_+) a_s j [mol/(m^3 s)]。
// 阳离子迁移数 t_+ < 1，界面反应产生的锂只有 t_+ 部分由迁移带走，
// 剩余 (1 - t_+) 部分留在液相中改变盐浓度；隔膜内 a_s = 0，源项为零。
double SpmeModel::electrolyte_source(int ix) const {
  const double j = ix < n_x_neg_ ? j_neg_ : j_pos_;
  return (1.0 - params_.elyte.t_plus) * a_s_[static_cast<std::size_t>(ix)] * j;
}

// 变量类型：1 表示微分变量（固相节点与电解液浓度），
// 0 表示代数变量（两个表面浓度与整条 phi_e 分布）。
std::vector<int> SpmeModel::variable_types() const {
  std::vector<int> types(static_cast<std::size_t>(size_), 1);
  types[static_cast<std::size_t>(cs_surf_neg())] = 0;
  types[static_cast<std::size_t>(cs_surf_pos())] = 0;
  for (int i = phie_offset_; i < size_; ++i) {
    types[static_cast<std::size_t>(i)] = 0;
  }
  return types;
}

// 解析雅可比 dF/dy + cj * dF/dy'，cj 为求解器给出的时间步系数。
// 固相与电解液浓度行为扩散三对角结构，时间导数项在对角上再加 cj
// （电解液行还要乘孔隙率 eps_e）；电位行为欧姆三对角结构。
// 反应通量由外加电流固定，故矩阵与 y、yp 无关。
void SpmeModel::jacobian(double /*t*/, const double* /*y*/, const double* /*yp*/,
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

  // 表面浓度行对自身取 +1、对最外层节点取 -1，均带统一缩放。
  out.add(cs_surf_neg(), cs_surf_neg(), eq_scale_surface_);
  out.add(cs_surf_neg(), cs_neg(n_r_neg_ - 1), -eq_scale_surface_);
  out.add(cs_surf_pos(), cs_surf_pos(), eq_scale_surface_);
  out.add(cs_surf_pos(), cs_pos(n_r_pos_ - 1), -eq_scale_surface_);

  for (int i = 0; i < n_x_total_; ++i) {
    const int row = ce(i);
    // 对角：cj 乘孔隙率（来自 eps_e * dc_e/dt）再减去扩散对角系数。
    out.add(row, row,
            (eps_e_[static_cast<std::size_t>(i)] * cj -
             elyte_op_.matrix.diag[static_cast<std::size_t>(i)]) *
                eq_scale_electrolyte_);
    if (i > 0) {
      out.add(row, ce(i - 1),
              -elyte_op_.matrix.lower[static_cast<std::size_t>(i)] *
                  eq_scale_electrolyte_);
    }
    if (i + 1 < n_x_total_) {
      out.add(row, ce(i + 1),
              -elyte_op_.matrix.upper[static_cast<std::size_t>(i)] *
                  eq_scale_electrolyte_);
    }
  }

  out.add(phie(0), phie(0), eq_scale_potential_);
  for (int i = 1; i < n_x_total_; ++i) {
    const int row = phie(i);
    out.add(row, row,
            phie_op_.diag[static_cast<std::size_t>(i)] * eq_scale_potential_);
    out.add(row, phie(i - 1),
            phie_op_.lower[static_cast<std::size_t>(i)] * eq_scale_potential_);
    if (i + 1 < n_x_total_) {
      out.add(row, phie(i + 1),
              phie_op_.upper[static_cast<std::size_t>(i)] *
                  eq_scale_potential_);
    }
  }
}

// 初值：颗粒内固相浓度与电解液浓度都取设计值 c0，phi_e 全部取 0，
// 与规范条件一致；此时无浓度梯度，也无欧姆压降。
std::vector<double> SpmeModel::initial_state() const {
  std::vector<double> y(static_cast<std::size_t>(size_), 0.0);
  for (int ir = 0; ir < n_r_neg_; ++ir) {
    y[static_cast<std::size_t>(cs_neg(ir))] = params_.neg.c0;
  }
  for (int ir = 0; ir < n_r_pos_; ++ir) {
    y[static_cast<std::size_t>(cs_pos(ir))] = params_.pos.c0;
  }
  y[static_cast<std::size_t>(cs_surf_neg())] = params_.neg.c0;
  y[static_cast<std::size_t>(cs_surf_pos())] = params_.pos.c0;
  for (int i = 0; i < n_x_total_; ++i) {
    y[static_cast<std::size_t>(ce(i))] = params_.elyte.c0;
  }
  return y;
}

// 绝对容差：内部默认 1e-8；固相浓度按 1e-2、电解液浓度按 1e-3 放宽，
// 电位分量保持 1e-8（初值为零，量级本身很小）。
std::vector<double> SpmeModel::suggested_atol() const {
  std::vector<double> atol(static_cast<std::size_t>(size_), 1e-8);
  for (int ir = 0; ir < n_r_neg_; ++ir) {
    atol[static_cast<std::size_t>(cs_neg(ir))] = 1e-2;
  }
  for (int ir = 0; ir < n_r_pos_; ++ir) {
    atol[static_cast<std::size_t>(cs_pos(ir))] = 1e-2;
  }
  atol[static_cast<std::size_t>(cs_surf_neg())] = 1e-2;
  atol[static_cast<std::size_t>(cs_surf_pos())] = 1e-2;
  for (int i = 0; i < n_x_total_; ++i) {
    atol[static_cast<std::size_t>(ce(i))] = 1e-3;
  }
  return atol;
}

// 端电压 = U_pos - U_neg + eta_pos - eta_neg + (phi_e(正极端) - phi_e(负极端))。
// 电解液电位差即液相欧姆压降与浓度极化带来的电位偏移；放电时 V 低于 OCV 差。
double SpmeModel::terminal_voltage(const double* y) const {
  const double csn = y[cs_surf_neg()];
  const double csp = y[cs_surf_pos()];
  const double u_neg = open_circuit(params_.neg, false, csn);
  const double u_pos = open_circuit(params_.pos, true, csp);

  // 取各电极中点节点作为代表性位置，用该处的电解液浓度求过电位。
  const int ix_neg = n_x_neg_ / 2;
  const int ix_pos = n_x_neg_ + n_x_sep_ + n_x_pos_ / 2;
  const double eta_neg = overpotential_from_flux(
      params_, params_.neg, false, y[ce(ix_neg)], csn, j_neg_);
  const double eta_pos = overpotential_from_flux(
      params_, params_.pos, true, y[ce(ix_pos)], csp, j_pos_);

  // phi_e 沿 x 的增量：正极集流体处减负极集流体处。
  const double phie_drop =
      y[phie(n_x_total_ - 1)] - y[phie(0)];
  return u_pos - u_neg + eta_pos - eta_neg + phie_drop;
}

// 单个颗粒的锂含量 [mol/m^2 电极面积]：按球壳体积加权平均浓度，
// 再乘固相体积分数 (1 - eps) 与电极厚度。
double SpmeModel::particle_lithium(const double* y, bool positive) const {
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
  return (weighted / particle_volume) * (1.0 - electrode.porosity) *
         electrode.thickness;
}

// 固相锂总量 = 负极颗粒锂量 + 正极颗粒锂量 [mol/m^2]。
double SpmeModel::solid_lithium(const double* y) const {
  return particle_lithium(y, false) + particle_lithium(y, true);
}

// 电解液中的锂总量 [mol/m^2 电极面积]：sum(eps_e c_e dx)，
// 体积除以横截面积后即为单元厚度 dx。
double SpmeModel::electrolyte_lithium(const double* y) const {
  double total = 0.0;
  for (int i = 0; i < n_x_total_; ++i) {
    const double dx =
        x_grid_.volume[static_cast<std::size_t>(i)] / params_.area;
    total += eps_e_[static_cast<std::size_t>(i)] *
             y[ce(i)] * dx;
  }
  return total;
}

// 全电池锂总量（固相 + 液相），用于守恒性检查。
double SpmeModel::total_lithium(const double* y) const {
  return solid_lithium(y) + electrolyte_lithium(y);
}

// 电极平均化学计量比 = 平均固相浓度 / c_max。
double SpmeModel::average_stoichiometry(const double* y,
                                        bool positive) const {
  const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
  return particle_lithium(y, positive) /
         ((1.0 - electrode.porosity) * electrode.thickness * electrode.c_max);
}

// 登记雅可比稀疏结构：两段径向三对角、两个表面行、
// 电解液浓度三对角与电位三对角（phi_e(0) 只有对角元）。
void SpmeModel::build_pattern() {
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
  add(cs_surf_neg(), cs_surf_neg());
  add(cs_surf_neg(), cs_neg(n_r_neg_ - 1));
  add(cs_surf_pos(), cs_surf_pos());
  add(cs_surf_pos(), cs_pos(n_r_pos_ - 1));

  for (int i = 0; i < n_x_total_; ++i) {
    add(ce(i), ce(i));
    if (i > 0) {
      add(ce(i), ce(i - 1));
    }
    if (i + 1 < n_x_total_) {
      add(ce(i), ce(i + 1));
    }
  }

  add(phie(0), phie(0));
  for (int i = 1; i < n_x_total_; ++i) {
    add(phie(i), phie(i));
    add(phie(i), phie(i - 1));
    if (i + 1 < n_x_total_) {
      add(phie(i), phie(i + 1));
    }
  }

  pattern_ = ::dfn::build_pattern(size_, std::move(entries));
}

}  // namespace dfn
