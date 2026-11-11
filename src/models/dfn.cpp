// DFN（P2D）模型的残差、Jacobian 与物理量统计实现。
// 状态向量分块顺序：[c_s 负 | c_s 正 | c_e | φ_s 负 | φ_s 正 | φ_e |
//                    c_s,surf 负 | c_s,surf 正]（详细偏移见 StateLayout）。
#include "dfn/models/dfn.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "dfn/models/reaction.hpp"
#include "dfn/thermo.hpp"

namespace dfn {

// 构造函数：一次性建立网格与全部线性算子，缓存与时间无关的系数。
// params：电化学与几何参数；layout：状态向量分块（各区域节点数、径向节点数）；
// options：控制模式（恒流/恒压）及对应的电流密度或端电压。
// 网格分三类：x 方向分段网格（负极|隔膜|正极，用于 c_e、φ_e）、
// 两个电极各自的轴向网格（用于 φ_s）、球形颗粒的径向网格（用于 c_s(r)）。
DfnModel::DfnModel(Params params, StateLayout layout, DfnOptions options)
    : params_(std::move(params)),
      layout_(layout),
      options_(options),
      neg_reaction_(static_cast<std::size_t>(layout.n_x_neg)),
      pos_reaction_(static_cast<std::size_t>(layout.n_x_pos)) {
  validate(params_);

  const double l_neg = params_.neg.thickness;
  const double l_sep = params_.sep.thickness;
  const double l_pos = params_.pos.thickness;

  x_grid_ = make_segmented_line_grid({layout_.n_x_neg, layout_.n_x_sep,
                                      layout_.n_x_pos},
                                     {l_neg, l_sep, l_pos}, params_.area);
  neg_grid_ = make_uniform_line_grid(layout_.n_x_neg, l_neg, params_.area);
  pos_grid_ = make_uniform_line_grid(layout_.n_x_pos, l_pos, params_.area);
  sphere_neg_ = make_uniform_sphere_grid(layout_.n_r_neg, params_.neg.radius);
  sphere_pos_ = make_uniform_sphere_grid(layout_.n_r_pos, params_.pos.radius);

  const std::vector<double> d_s_neg(static_cast<std::size_t>(layout_.n_r_neg),
                                    params_.neg.d_s);
  const std::vector<double> d_s_pos(static_cast<std::size_t>(layout_.n_r_pos),
                                    params_.pos.d_s);
  // 球形颗粒扩散算子：径向节点上取同一 D_s，两端面通量暂置零，
  // 颗粒表面的反应通量在残差里作为最外层源项单独加入。
  radial_neg_ = assemble_diffusion(sphere_neg_, d_s_neg);
  radial_pos_ = assemble_diffusion(sphere_pos_, d_s_pos);

  // Bruggeman 修正：σ_eff = σ (1-ε)^b，固相体积分数为 1-ε。
  sigma_neg_eff_ = params_.neg.sigma *
                   std::pow(1.0 - params_.neg.porosity,
                            params_.neg.bruggeman);
  sigma_pos_eff_ = params_.pos.sigma *
                   std::pow(1.0 - params_.pos.porosity,
                            params_.pos.bruggeman);

  const std::vector<double> sigma_neg(static_cast<std::size_t>(layout_.n_x_neg),
                                      sigma_neg_eff_);
  const std::vector<double> sigma_pos(static_cast<std::size_t>(layout_.n_x_pos),
                                      sigma_pos_eff_);
  // 固相电位算子：均匀的电导率下与扩散算子同构，只取矩阵部分。
  phis_neg_op_ = assemble_diffusion(neg_grid_, sigma_neg).matrix;
  phis_pos_op_ = assemble_diffusion(pos_grid_, sigma_pos).matrix;

  // 逐节点缓存区域属性：负极在前、隔膜居中、正极在后（与全局 x 编号一致）。
  const int nx = layout_.n_x_total();
  eps_e_.assign(static_cast<std::size_t>(nx), 0.0);
  d_eff_.assign(static_cast<std::size_t>(nx), 0.0);
  kappa_eff_.assign(static_cast<std::size_t>(nx), 0.0);
  a_s_.assign(static_cast<std::size_t>(nx), 0.0);
  for (int i = 0; i < nx; ++i) {
    double epsy = 0.0;
    double bruggeman = 1.5;
    if (i < layout_.n_x_neg) {
      epsy = params_.neg.porosity;
      bruggeman = params_.neg.bruggeman;
      // 比表面积 a_s = 3(1-ε)/R，球形颗粒单位体积的界面面积。
      a_s_[static_cast<std::size_t>(i)] =
          3.0 * (1.0 - epsy) / params_.neg.radius;
    } else if (i < layout_.n_x_neg + layout_.n_x_sep) {
      epsy = params_.sep.porosity;
      bruggeman = params_.sep.bruggeman;
    } else {
      epsy = params_.pos.porosity;
      bruggeman = params_.pos.bruggeman;
      a_s_[static_cast<std::size_t>(i)] =
          3.0 * (1.0 - epsy) / params_.pos.radius;
    }
    eps_e_[static_cast<std::size_t>(i)] = epsy;
    // 有效扩散系数与有效离子电导率均按 Bruggeman 关系 ε^b 缩放。
    d_eff_[static_cast<std::size_t>(i)] =
        params_.elyte.d_e * std::pow(epsy, bruggeman);
    kappa_eff_[static_cast<std::size_t>(i)] =
        params_.elyte.kappa * std::pow(epsy, bruggeman);
  }
  // assemble_diffusion 给出有限体积形式的 div(D ∇c) 三对角算子（对角元为负），
  // 残差里直接以“算子乘状态向量”出现。
  elyte_op_ = assemble_diffusion(x_grid_, d_eff_);
  phie_op_ = assemble_diffusion(x_grid_, kappa_eff_).matrix;

  // 以下系数用于集流体边界条件、颗粒表面通量与表面浓度外推。
  dx_neg_ = l_neg / static_cast<double>(layout_.n_x_neg);
  dx_pos_ = l_pos / static_cast<double>(layout_.n_x_pos);
  // 半格虚拟单元实现 Dirichlet 条件：面电导 σ/dx/2 折算后为 -2σ/dx²。
  neg_dirichlet_ = -2.0 * sigma_neg_eff_ / (dx_neg_ * dx_neg_);
  pos_dirichlet_ = -2.0 * sigma_pos_eff_ / (dx_pos_ * dx_pos_);
  // 恒流边界通量项：面通量除以单元体积（面积/体积 = 1/dx）。
  pos_flux_coeff_ = -1.0 / dx_pos_;
  // 颗粒表面反应通量按面积/壳层体积折算到最外一圈球壳。
  surf_coeff_neg_ = -sphere_neg_.area[static_cast<std::size_t>(layout_.n_r_neg)] /
                    sphere_neg_.volume[static_cast<std::size_t>(layout_.n_r_neg - 1)];
  surf_coeff_pos_ = -sphere_pos_.area[static_cast<std::size_t>(layout_.n_r_pos)] /
                    sphere_pos_.volume[static_cast<std::size_t>(layout_.n_r_pos - 1)];
  const double dr_neg = params_.neg.radius / static_cast<double>(layout_.n_r_neg);
  const double dr_pos = params_.pos.radius / static_cast<double>(layout_.n_r_pos);
  // 表面浓度外推系数 α = dr/(2 D_s)，由表面通量修正最外层浓度。
  alpha_neg_ = dr_neg / (2.0 * params_.neg.d_s);
  alpha_pos_ = dr_pos / (2.0 * params_.pos.d_s);

  // 方程缩放因子：用特征时间与特征浓度把各方程残差化为无量纲量级，
  // 避免量纲悬殊导致 Newton 收敛判据失衡。
  const double tau_s_neg = params_.neg.radius * params_.neg.radius /
                           params_.neg.d_s;
  const double tau_s_pos = params_.pos.radius * params_.pos.radius /
                           params_.pos.d_s;
  eq_scale_solid_neg_ = tau_s_neg / params_.neg.c_max;
  eq_scale_solid_pos_ = tau_s_pos / params_.pos.c_max;
  const double length = params_.cell_length();
  const double tau_e = length * length / params_.elyte.d_e;
  eq_scale_electrolyte_ = params_.elyte.c0 > 0.0
                              ? tau_e / params_.elyte.c0
                              : 1.0;
  // 参考电流密度 1 A/m^2；其具体取值只影响缩放，不影响解。
  eq_scale_potential_ = length;

  // 表面浓度作为附加代数未知量，紧接在 StateLayout 的各块之后追加。
  cs_surf_neg_offset_ = layout_.size();
  cs_surf_pos_offset_ = layout_.size() + layout_.n_x_neg;
  size_ = layout_.size() + layout_.n_x_neg + layout_.n_x_pos;
  eq_scale_surface_ =
      1.0 / std::max(params_.neg.c_max, params_.pos.c_max);

  // 稀疏模式只需建立一次，之后 Jacobian 只负责填值。
  build_pattern();
}

// 局部轴向节点 ix 映射到全局 x 编号：负极从 0 起，正极位于隔膜之后。
int DfnModel::global_index(int local_ix, bool positive) const {
  return positive ? layout_.n_x_neg + layout_.n_x_sep + local_ix : local_ix;
}

// 组装 DAE 残差 F(t, y, y') = 0，是 IDA 每步 Newton 迭代与误差估计的入口。
// t：当前时间（残差本身不显含时间，故不使用）；y / yp：状态向量及其时间导数；
// out：长度 size() 的残差，函数内部先清零再按分块顺序逐块写入。
// 微分变量为 c_s 与 c_e（出现 yp），代数变量为 c_s,surf 与全部电位。
// 每条方程都乘以对应的 eq_scale_*，把量纲差异压缩到相近量级。
void DfnModel::residual(double /*t*/, const double* y, const double* yp,
                        double* out) const {
  const int n = size();
  std::fill(out, out + n, 0.0);

  const double f = kFaraday;
  const double one_minus_t = 1.0 - params_.elyte.t_plus;
  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;

  // 先按当前状态算出每个涂布节点上的 Butler-Volmer 反应通量 j（mol/(m² s)，
  // j > 0 表示锂离开固相进入电解液），后面的扩散与守恒方程都复用这份结果。
  // 注意 reaction_at 内部已完成 i0、过电位与 OCP 斜率的计算。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int gx = global_index(ix, false);
    neg_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_, params_.neg, false, y[layout_.ce(gx)],
        y[cs_surf_neg(ix)], y[layout_.phis_neg(ix)],
        y[layout_.phie(gx)]);
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int gx = global_index(ix, true);
    pos_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_, params_.pos, true, y[layout_.ce(gx)],
        y[cs_surf_pos(ix)], y[layout_.phis_pos(ix)],
        y[layout_.phie(gx)]);
  }

  // 固相球形颗粒扩散：∂c_s/∂t = ∇·(D_s ∇c_s)。径向用有限体积离散，
  // 最外圈球壳额外带上反应通量源项；j > 0 时锂从颗粒表面被抽走，故为负源。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const double j_surf =
        neg_reaction_[static_cast<std::size_t>(ix)].j;
    for (int ir = 0; ir < nr_neg; ++ir) {
      const int row = layout_.cs_neg(ix, ir);
      double lap = radial_neg_.matrix.diag[static_cast<std::size_t>(ir)] * y[row];
      if (ir > 0) {
        lap += radial_neg_.matrix.lower[static_cast<std::size_t>(ir)] *
               y[layout_.cs_neg(ix, ir - 1)];
      }
      if (ir + 1 < nr_neg) {
        lap += radial_neg_.matrix.upper[static_cast<std::size_t>(ir)] *
               y[layout_.cs_neg(ix, ir + 1)];
      }
      // 只有最外层球壳接收反应通量，内部节点为纯扩散。
      const double source =
          (ir + 1 == nr_neg) ? surf_coeff_neg_ * j_surf : 0.0;
      out[row] = (yp[row] - lap - source) * eq_scale_solid_neg_;
    }
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const double j_surf =
        pos_reaction_[static_cast<std::size_t>(ix)].j;
    for (int ir = 0; ir < nr_pos; ++ir) {
      const int row = layout_.cs_pos(ix, ir);
      double lap = radial_pos_.matrix.diag[static_cast<std::size_t>(ir)] * y[row];
      if (ir > 0) {
        lap += radial_pos_.matrix.lower[static_cast<std::size_t>(ir)] *
               y[layout_.cs_pos(ix, ir - 1)];
      }
      if (ir + 1 < nr_pos) {
        lap += radial_pos_.matrix.upper[static_cast<std::size_t>(ir)] *
               y[layout_.cs_pos(ix, ir + 1)];
      }
      const double source =
          (ir + 1 == nr_pos) ? surf_coeff_pos_ * j_surf : 0.0;
      out[row] = (yp[row] - lap - source) * eq_scale_solid_pos_;
    }
  }

  // 表面浓度由最外层节点浓度与表面通量线性外推：c_surf = c_last - α j。
  // 该式把“颗粒表面”从体积平均中还原出来，供 Butler-Volmer 使用。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int row = cs_surf_neg(ix);
    out[row] = (y[row] - y[layout_.cs_neg(ix, nr_neg - 1)] +
                alpha_neg_ * neg_reaction_[static_cast<std::size_t>(ix)].j) *
               eq_scale_surface_;
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int row = cs_surf_pos(ix);
    out[row] = (y[row] - y[layout_.cs_pos(ix, nr_pos - 1)] +
                alpha_pos_ * pos_reaction_[static_cast<std::size_t>(ix)].j) *
               eq_scale_surface_;
  }

  // 电解液质量守恒：ε ∂c_e/∂t = ∇·(D_eff ∇c_e) + (1-t_+) a_s j。
  // 源项只存在于两个电极区域，隔膜内 a_s = 0，只剩扩散。
  const int nx = layout_.n_x_total();
  for (int i = 0; i < nx; ++i) {
    const int row = layout_.ce(i);
    double lap = elyte_op_.matrix.diag[static_cast<std::size_t>(i)] *
                 y[row];
    if (i > 0) {
      lap += elyte_op_.matrix.lower[static_cast<std::size_t>(i)] *
             y[layout_.ce(i - 1)];
    }
    if (i + 1 < nx) {
      lap += elyte_op_.matrix.upper[static_cast<std::size_t>(i)] *
             y[layout_.ce(i + 1)];
    }
    double source = 0.0;
    if (i < layout_.n_x_neg) {
      source = one_minus_t * a_s_[static_cast<std::size_t>(i)] *
               neg_reaction_[static_cast<std::size_t>(i)].j;
    } else if (i >= layout_.n_x_neg + layout_.n_x_sep) {
      source = one_minus_t * a_s_[static_cast<std::size_t>(i)] *
               pos_reaction_[static_cast<std::size_t>(i - layout_.n_x_neg -
                                                      layout_.n_x_sep)]
                   .j;
    }
    out[row] = (eps_e_[static_cast<std::size_t>(i)] * yp[row] - lap - source) *
               eq_scale_electrolyte_;
  }

  // 负极固相电荷守恒：∇·(σ_eff ∇φ_s) = a_s F j。
  // 负极为参考电极：ix = 0 处（负极集流体）施加 φ_s = 0 的 Dirichlet 条件。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int row = layout_.phis_neg(ix);
    double lap = phis_neg_op_.diag[static_cast<std::size_t>(ix)] * y[row];
    if (ix > 0) {
      lap += phis_neg_op_.lower[static_cast<std::size_t>(ix)] *
             y[layout_.phis_neg(ix - 1)];
    }
    if (ix + 1 < layout_.n_x_neg) {
      lap += phis_neg_op_.upper[static_cast<std::size_t>(ix)] *
             y[layout_.phis_neg(ix + 1)];
    }
    if (ix == 0) {
      lap += neg_dirichlet_ * y[row];
    }
    const double source = a_s_[static_cast<std::size_t>(ix)] * f *
                          neg_reaction_[static_cast<std::size_t>(ix)].j;
    out[row] = (lap - source) * eq_scale_potential_;
  }

  // 正极固相电荷守恒：∇·(σ_eff ∇φ_s) = a_s F j。
  // 最末节点贴着正极集流体：恒流模式施加 Neumann 通量（外加电流密度），
  // 恒压模式施加 Dirichlet 电位 φ_s(L) = applied_voltage。
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int row = layout_.phis_pos(ix);
    double lap = phis_pos_op_.diag[static_cast<std::size_t>(ix)] * y[row];
    if (ix > 0) {
      lap += phis_pos_op_.lower[static_cast<std::size_t>(ix)] *
             y[layout_.phis_pos(ix - 1)];
    }
    if (ix + 1 < layout_.n_x_pos) {
      lap += phis_pos_op_.upper[static_cast<std::size_t>(ix)] *
             y[layout_.phis_pos(ix + 1)];
    }
    // 集流体边界项：恒流时为 Neumann 通量，恒压时为半格虚拟单元的 Dirichlet。
    const bool last = (ix + 1 == layout_.n_x_pos);
    if (last) {
      if (options_.mode == ControlMode::current) {
        lap += pos_flux_coeff_ * options_.applied_current_density;
      } else {
        lap += pos_dirichlet_ * y[row] - pos_dirichlet_ * options_.applied_voltage;
      }
    }
    const int gx = global_index(ix, true);
    const double source = a_s_[static_cast<std::size_t>(gx)] * f *
                          pos_reaction_[static_cast<std::size_t>(ix)].j;
    out[row] = (lap - source) * eq_scale_potential_;
  }

  // 液相电荷守恒：∇·(κ_eff ∇φ_e) + a_s F j = 0。
  // 两端集流体处离子电流为零，故 φ_e 无额外边界项（自然边界）。
  for (int i = 0; i < nx; ++i) {
    const int row = layout_.phie(i);
    double lap = phie_op_.diag[static_cast<std::size_t>(i)] * y[row];
    if (i > 0) {
      lap += phie_op_.lower[static_cast<std::size_t>(i)] * y[layout_.phie(i - 1)];
    }
    if (i + 1 < nx) {
      lap += phie_op_.upper[static_cast<std::size_t>(i)] * y[layout_.phie(i + 1)];
    }
    double source = 0.0;
    if (i < layout_.n_x_neg) {
      source = a_s_[static_cast<std::size_t>(i)] * f *
               neg_reaction_[static_cast<std::size_t>(i)].j;
    } else if (i >= layout_.n_x_neg + layout_.n_x_sep) {
      source = a_s_[static_cast<std::size_t>(i)] * f *
               pos_reaction_[static_cast<std::size_t>(i - layout_.n_x_neg -
                                                      layout_.n_x_sep)]
                   .j;
    }
    out[row] = (lap + source) * eq_scale_potential_;
  }
}

// 标记变量的微分/代数属性：1 = 微分变量，0 = 代数变量。
// 电位块（φ_s 负、φ_s 正、φ_e）来自电荷守恒的椭圆方程，不含 ∂/∂t；
// 表面浓度由外推式直接定义，同样从 layout_.size() 起整块置 0。
std::vector<int> DfnModel::variable_types() const {
  std::vector<int> types(static_cast<std::size_t>(size()), 1);
  for (int i = layout_.phis_neg_offset(); i < layout_.size(); ++i) {
    types[static_cast<std::size_t>(i)] = 0;
  }
  for (int i = cs_surf_neg_offset_; i < size(); ++i) {
    types[static_cast<std::size_t>(i)] = 0;
  }
  return types;
}

// 组装 Jacobian J = dF/dy + cj * dF/dy'，写入预分配的 CSC 结构。
// cj 是 BDF 方法的标量系数（IDA 传入）；y 在当前迭代点上取值；
// yp 不参与（残差对 y' 的依赖是线性的，只有 cj 与 eq_scale 的乘积）。
// 结构分三类：各块沿轴向/径向的三对角耦合、同一 x 节点上
// c_s,surf / c_e / φ_s / φ_e 之间的 Butler-Volmer 局部耦合、集流体边界项。
// 耦合系数由反应通量的偏导 dj_dcs、dj_dce、dj_deta 经链式法则得到。
void DfnModel::jacobian(double /*t*/, const double* y, const double* /*yp*/,
                        double cj, CscBuilder& out) const {
  const double f = kFaraday;
  const double one_minus_t = 1.0 - params_.elyte.t_plus;
  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;

  // 与 residual 同步重算反应通量及其偏导，保证 Jacobian 与残差在同一线性化点。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int gx = global_index(ix, false);
    neg_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_, params_.neg, false, y[layout_.ce(gx)],
        y[cs_surf_neg(ix)], y[layout_.phis_neg(ix)],
        y[layout_.phie(gx)]);
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int gx = global_index(ix, true);
    pos_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_, params_.pos, true, y[layout_.ce(gx)],
        y[cs_surf_pos(ix)], y[layout_.phis_pos(ix)],
        y[layout_.phie(gx)]);
  }

  // 固相扩散行：对角元含 cj（来自 dF/dy'），径向相邻节点为纯扩散耦合。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const auto& reaction = neg_reaction_[static_cast<std::size_t>(ix)];
    for (int ir = 0; ir < nr_neg; ++ir) {
      const int row = layout_.cs_neg(ix, ir);
      out.add(row, row,
              (-radial_neg_.matrix.diag[static_cast<std::size_t>(ir)] + cj) *
                  eq_scale_solid_neg_);
      if (ir > 0) {
        out.add(row, layout_.cs_neg(ix, ir - 1),
                -radial_neg_.matrix.lower[static_cast<std::size_t>(ir)] *
                    eq_scale_solid_neg_);
      }
      if (ir + 1 < nr_neg) {
        out.add(row, layout_.cs_neg(ix, ir + 1),
                -radial_neg_.matrix.upper[static_cast<std::size_t>(ir)] *
                    eq_scale_solid_neg_);
      }
    }
    // 最外层球壳的源项含 j，因此对 c_surf、c_e、φ_s、φ_e 都有非零偏导；
    // φ_s 与 φ_e 的符号相反，因为过电位 η = φ_s - φ_e - U。
    const int last_cell = layout_.cs_neg(ix, nr_neg - 1);
    const int gx = global_index(ix, false);
    out.add(last_cell, cs_surf_neg(ix),
            -surf_coeff_neg_ * reaction.dj_dcs * eq_scale_solid_neg_);
    out.add(last_cell, layout_.ce(gx),
            -surf_coeff_neg_ * reaction.dj_dce * eq_scale_solid_neg_);
    out.add(last_cell, layout_.phis_neg(ix),
            -surf_coeff_neg_ * reaction.dj_deta * eq_scale_solid_neg_);
    out.add(last_cell, layout_.phie(gx),
            surf_coeff_neg_ * reaction.dj_deta * eq_scale_solid_neg_);
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const auto& reaction = pos_reaction_[static_cast<std::size_t>(ix)];
    for (int ir = 0; ir < nr_pos; ++ir) {
      const int row = layout_.cs_pos(ix, ir);
      out.add(row, row,
              (-radial_pos_.matrix.diag[static_cast<std::size_t>(ir)] + cj) *
                  eq_scale_solid_pos_);
      if (ir > 0) {
        out.add(row, layout_.cs_pos(ix, ir - 1),
                -radial_pos_.matrix.lower[static_cast<std::size_t>(ir)] *
                    eq_scale_solid_pos_);
      }
      if (ir + 1 < nr_pos) {
        out.add(row, layout_.cs_pos(ix, ir + 1),
                -radial_pos_.matrix.upper[static_cast<std::size_t>(ir)] *
                    eq_scale_solid_pos_);
      }
    }
    const int last_cell = layout_.cs_pos(ix, nr_pos - 1);
    const int gx = global_index(ix, true);
    out.add(last_cell, cs_surf_pos(ix),
            -surf_coeff_pos_ * reaction.dj_dcs * eq_scale_solid_pos_);
    out.add(last_cell, layout_.ce(gx),
            -surf_coeff_pos_ * reaction.dj_dce * eq_scale_solid_pos_);
    out.add(last_cell, layout_.phis_pos(ix),
            -surf_coeff_pos_ * reaction.dj_deta * eq_scale_solid_pos_);
    out.add(last_cell, layout_.phie(gx),
            surf_coeff_pos_ * reaction.dj_deta * eq_scale_solid_pos_);
  }

  // 表面浓度代数方程（无 cj 项）：c_surf = c_last - α j，
  // 对角元是 1 + α dj/dc_surf，其余为 α 乘以 j 对各电位的偏导。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int row = cs_surf_neg(ix);
    const int gx = global_index(ix, false);
    const auto& reaction = neg_reaction_[static_cast<std::size_t>(ix)];
    out.add(row, row,
            (1.0 + alpha_neg_ * reaction.dj_dcs) * eq_scale_surface_);
    out.add(row, layout_.cs_neg(ix, nr_neg - 1), -eq_scale_surface_);
    out.add(row, layout_.ce(gx),
            alpha_neg_ * reaction.dj_dce * eq_scale_surface_);
    out.add(row, layout_.phis_neg(ix),
            alpha_neg_ * reaction.dj_deta * eq_scale_surface_);
    out.add(row, layout_.phie(gx),
            -alpha_neg_ * reaction.dj_deta * eq_scale_surface_);
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int row = cs_surf_pos(ix);
    const int gx = global_index(ix, true);
    const auto& reaction = pos_reaction_[static_cast<std::size_t>(ix)];
    out.add(row, row,
            (1.0 + alpha_pos_ * reaction.dj_dcs) * eq_scale_surface_);
    out.add(row, layout_.cs_pos(ix, nr_pos - 1), -eq_scale_surface_);
    out.add(row, layout_.ce(gx),
            alpha_pos_ * reaction.dj_dce * eq_scale_surface_);
    out.add(row, layout_.phis_pos(ix),
            alpha_pos_ * reaction.dj_deta * eq_scale_surface_);
    out.add(row, layout_.phie(gx),
            -alpha_pos_ * reaction.dj_deta * eq_scale_surface_);
  }

  // 电解液行：对角含 ε cj（时间项），相邻节点来自有效扩散系数；
  // 电极区域再叠加反应源对 c_surf、c_e、φ_s、φ_e 的偏导，隔膜内无此耦合。
  const int nx = layout_.n_x_total();
  for (int i = 0; i < nx; ++i) {
    const int row = layout_.ce(i);
    out.add(row, row,
            (eps_e_[static_cast<std::size_t>(i)] * cj -
             elyte_op_.matrix.diag[static_cast<std::size_t>(i)]) *
                eq_scale_electrolyte_);
    if (i > 0) {
      out.add(row, layout_.ce(i - 1),
              -elyte_op_.matrix.lower[static_cast<std::size_t>(i)] *
                  eq_scale_electrolyte_);
    }
    if (i + 1 < nx) {
      out.add(row, layout_.ce(i + 1),
              -elyte_op_.matrix.upper[static_cast<std::size_t>(i)] *
                  eq_scale_electrolyte_);
    }
    double coupling = 0.0;
    const ReactionValue* reaction = nullptr;
    int surface = -1;
    int phis_index = -1;
    if (i < layout_.n_x_neg) {
      reaction = &neg_reaction_[static_cast<std::size_t>(i)];
      surface = cs_surf_neg(i);
      phis_index = layout_.phis_neg(i);
      coupling = -(one_minus_t) * a_s_[static_cast<std::size_t>(i)] *
                 eq_scale_electrolyte_;
    } else if (i >= layout_.n_x_neg + layout_.n_x_sep) {
      const int ix = i - layout_.n_x_neg - layout_.n_x_sep;
      reaction = &pos_reaction_[static_cast<std::size_t>(ix)];
      surface = cs_surf_pos(ix);
      phis_index = layout_.phis_pos(ix);
      coupling = -(one_minus_t) * a_s_[static_cast<std::size_t>(i)] *
                 eq_scale_electrolyte_;
    }
    if (reaction != nullptr) {
      out.add(row, surface, coupling * reaction->dj_dcs);
      out.add(row, row, coupling * reaction->dj_dce);
      out.add(row, phis_index, coupling * reaction->dj_deta);
      out.add(row, layout_.phie(i), -coupling * reaction->dj_deta);
    }
  }

  // 固相电位行：Laplacian 与集流体边界项，另经反应通量耦合到
  // c_surf、c_e、φ_s 自身与 φ_e（coupling = -a_s F * eq_scale）。
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int row = layout_.phis_neg(ix);
    double diagonal = phis_neg_op_.diag[static_cast<std::size_t>(ix)];
    if (ix == 0) {
      diagonal += neg_dirichlet_;
    }
    out.add(row, row, diagonal * eq_scale_potential_);
    if (ix > 0) {
      out.add(row, layout_.phis_neg(ix - 1),
              phis_neg_op_.lower[static_cast<std::size_t>(ix)] *
                  eq_scale_potential_);
    }
    if (ix + 1 < layout_.n_x_neg) {
      out.add(row, layout_.phis_neg(ix + 1),
              phis_neg_op_.upper[static_cast<std::size_t>(ix)] *
                  eq_scale_potential_);
    }
    const auto& reaction = neg_reaction_[static_cast<std::size_t>(ix)];
    const int gx = global_index(ix, false);
    const double coupling =
        -a_s_[static_cast<std::size_t>(gx)] * f * eq_scale_potential_;
    out.add(row, cs_surf_neg(ix), coupling * reaction.dj_dcs);
    out.add(row, layout_.ce(gx), coupling * reaction.dj_dce);
    out.add(row, row, coupling * reaction.dj_deta);
    out.add(row, layout_.phie(gx), -coupling * reaction.dj_deta);
  }
  // 正极固相电位行：结构与负极对称；恒压模式下末节点附加 Dirichlet 对角项。
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int row = layout_.phis_pos(ix);
    double diagonal = phis_pos_op_.diag[static_cast<std::size_t>(ix)];
    if (ix + 1 == layout_.n_x_pos && options_.mode == ControlMode::voltage) {
      diagonal += pos_dirichlet_;
    }
    out.add(row, row, diagonal * eq_scale_potential_);
    if (ix > 0) {
      out.add(row, layout_.phis_pos(ix - 1),
              phis_pos_op_.lower[static_cast<std::size_t>(ix)] *
                  eq_scale_potential_);
    }
    if (ix + 1 < layout_.n_x_pos) {
      out.add(row, layout_.phis_pos(ix + 1),
              phis_pos_op_.upper[static_cast<std::size_t>(ix)] *
                  eq_scale_potential_);
    }
    const auto& reaction = pos_reaction_[static_cast<std::size_t>(ix)];
    const int gx = global_index(ix, true);
    const double coupling =
        -a_s_[static_cast<std::size_t>(gx)] * f * eq_scale_potential_;
    out.add(row, cs_surf_pos(ix), coupling * reaction.dj_dcs);
    out.add(row, layout_.ce(gx), coupling * reaction.dj_dce);
    out.add(row, row, coupling * reaction.dj_deta);
    out.add(row, layout_.phie(gx), -coupling * reaction.dj_deta);
  }

  // 液相电位行：纯 Laplacian 加反应耦合。φ_e 出现在 η 的负号侧，
  // 故对角项为 -coupling * dj_deta，且与 φ_s 的耦合符号相反。
  for (int i = 0; i < nx; ++i) {
    const int row = layout_.phie(i);
    out.add(row, row,
            phie_op_.diag[static_cast<std::size_t>(i)] * eq_scale_potential_);
    if (i > 0) {
      out.add(row, layout_.phie(i - 1),
              phie_op_.lower[static_cast<std::size_t>(i)] *
                  eq_scale_potential_);
    }
    if (i + 1 < nx) {
      out.add(row, layout_.phie(i + 1),
              phie_op_.upper[static_cast<std::size_t>(i)] *
                  eq_scale_potential_);
    }
    const ReactionValue* reaction = nullptr;
    if (i < layout_.n_x_neg) {
      reaction = &neg_reaction_[static_cast<std::size_t>(i)];
    } else if (i >= layout_.n_x_neg + layout_.n_x_sep) {
      reaction = &pos_reaction_[static_cast<std::size_t>(i - layout_.n_x_neg -
                                                        layout_.n_x_sep)];
    }
    if (reaction != nullptr) {
      const double coupling =
          a_s_[static_cast<std::size_t>(i)] * f * eq_scale_potential_;
      if (i < layout_.n_x_neg) {
        out.add(row, cs_surf_neg(i), coupling * reaction->dj_dcs);
        out.add(row, layout_.phis_neg(i), coupling * reaction->dj_deta);
      } else {
        const int ix = i - layout_.n_x_neg - layout_.n_x_sep;
        out.add(row, cs_surf_pos(ix), coupling * reaction->dj_dcs);
        out.add(row, layout_.phis_pos(ix), coupling * reaction->dj_deta);
      }
      out.add(row, layout_.ce(i), coupling * reaction->dj_dce);
      out.add(row, row, -coupling * reaction->dj_deta);
    }
  }
}

// 切换到恒流控制：value 为外加电流密度（A/m²，正值表示放电、负值表示充电）。
void DfnModel::set_current_density(double value) {
  options_.mode = ControlMode::current;
  options_.applied_current_density = value;
}

// 切换到恒压控制：value 为集流体处的端电压（V）。
void DfnModel::set_applied_voltage(double value) {
  options_.mode = ControlMode::voltage;
  options_.applied_voltage = value;
}

// 按模式分发到恒流/恒压设置，供仿真驱动在 CC→CV 阶段切换时调用。
void DfnModel::set_control(ControlMode mode, double value) {
  if (mode == ControlMode::current) {
    set_current_density(value);
  } else {
    set_applied_voltage(value);
  }
}

// 构造一致的初始状态。基础部分（c_s、c_e）由 make_initial_state 填成均匀值；
// 表面浓度取各颗粒最外层径向节点；φ_e 取 -U_neg、φ_s 正极取 U_pos - U_neg，
// 使初始 η = φ_s - φ_e - U(θ) ≈ 0，端电压等于开路电压。
std::vector<double> DfnModel::initial_state() const {
  std::vector<double> base = make_initial_state(params_, layout_);
  std::vector<double> y(static_cast<std::size_t>(size()), 0.0);
  std::copy(base.begin(), base.end(), y.begin());

  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    y[static_cast<std::size_t>(cs_surf_neg(ix))] =
        base[static_cast<std::size_t>(layout_.cs_neg(ix, nr_neg - 1))];
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    y[static_cast<std::size_t>(cs_surf_pos(ix))] =
        base[static_cast<std::size_t>(layout_.cs_pos(ix, nr_pos - 1))];
  }

  // 按初始化学计量比查 OCP 曲线，得到与开路状态一致的电位分布。
  const double x0 = params_.neg.c0 / params_.neg.c_max;
  const double y0 = params_.pos.c0 / params_.pos.c_max;
  const double u_neg = negative_ocp(x0).value;
  const double u_pos = positive_ocp(y0).value;
  for (int i = 0; i < layout_.n_x_total(); ++i) {
    y[static_cast<std::size_t>(layout_.phie(i))] = -u_neg;
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    y[static_cast<std::size_t>(layout_.phis_pos(ix))] = u_pos - u_neg;
  }
  return y;
}

// 逐分量的建议绝对容差：浓度量级约 1e4 mol/m³，故 c_s 与 c_surf 取 1e-2、
// c_e 取 1e-3；电位（V）保持默认 1e-8。
std::vector<double> DfnModel::suggested_atol() const {
  std::vector<double> atol(static_cast<std::size_t>(size()), 1e-8);
  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    for (int ir = 0; ir < nr_neg; ++ir) {
      atol[static_cast<std::size_t>(layout_.cs_neg(ix, ir))] = 1e-2;
    }
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    for (int ir = 0; ir < nr_pos; ++ir) {
      atol[static_cast<std::size_t>(layout_.cs_pos(ix, ir))] = 1e-2;
    }
  }
  for (int i = 0; i < layout_.n_x_total(); ++i) {
    atol[static_cast<std::size_t>(layout_.ce(i))] = 1e-3;
  }
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    atol[static_cast<std::size_t>(cs_surf_neg(ix))] = 1e-2;
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    atol[static_cast<std::size_t>(cs_surf_pos(ix))] = 1e-2;
  }
  return atol;
}

// 端电压（V）：恒压模式下直接返回设定值；恒流模式下由正极集流体的
// Neumann 边界条件，把最末单元的 φ_s 外推半个单元到边界。
double DfnModel::terminal_voltage(const double* y) const {
  if (options_.mode == ControlMode::voltage) {
    return options_.applied_voltage;
  }
  const int last = layout_.n_x_pos - 1;
  const double phi_last = y[layout_.phis_pos(last)];
  return phi_last -
         options_.applied_current_density * dx_pos_ / (2.0 * sigma_pos_eff_);
}

// 端电流密度（A/m²，正值表示放电）：恒流模式返回设定值；
// 恒压模式由边界电位差与半格电导反算，与 terminal_voltage 互为逆运算。
double DfnModel::terminal_current(const double* y) const {
  if (options_.mode == ControlMode::current) {
    return options_.applied_current_density;
  }
  const int last = layout_.n_x_pos - 1;
  const double phi_last = y[layout_.phis_pos(last)];
  return 2.0 * sigma_pos_eff_ *
         (phi_last - options_.applied_voltage) / dx_pos_;
}

// 固相储锂总量（mol）：每个 x 单元先按球壳体积加权求颗粒平均浓度，
// 再乘固相体积分数 (1-ε) 与单元体积，最后对整片电极求和。
double DfnModel::solid_lithium(const double* y) const {
  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;
  double total = 0.0;
  double particle_volume = 0.0;
  for (double v : sphere_neg_.volume) {
    particle_volume += v;
  }
  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    double weighted = 0.0;
    for (int ir = 0; ir < nr_neg; ++ir) {
      weighted += y[layout_.cs_neg(ix, ir)] *
                  sphere_neg_.volume[static_cast<std::size_t>(ir)];
    }
    total += (weighted / particle_volume) * (1.0 - params_.neg.porosity) *
             neg_grid_.volume[static_cast<std::size_t>(ix)];
  }
  particle_volume = 0.0;
  for (double v : sphere_pos_.volume) {
    particle_volume += v;
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    double weighted = 0.0;
    for (int ir = 0; ir < nr_pos; ++ir) {
      weighted += y[layout_.cs_pos(ix, ir)] *
                  sphere_pos_.volume[static_cast<std::size_t>(ir)];
    }
    total += (weighted / particle_volume) * (1.0 - params_.pos.porosity) *
             pos_grid_.volume[static_cast<std::size_t>(ix)];
  }
  return total;
}

// 电解液储锂总量（mol）：Σ ε c_e ΔV，与 solid_lithium 一起用于守恒检查。
double DfnModel::electrolyte_lithium(const double* y) const {
  double total = 0.0;
  for (int i = 0; i < layout_.n_x_total(); ++i) {
    total += eps_e_[static_cast<std::size_t>(i)] * y[layout_.ce(i)] *
             x_grid_.volume[static_cast<std::size_t>(i)];
  }
  return total;
}

// 全电池储锂总量（mol）= 固相 + 液相。锂只在两相之间转移，总量应守恒。
double DfnModel::total_lithium(const double* y) const {
  return solid_lithium(y) + electrolyte_lithium(y);
}

// 电极的体积平均化学计量比 θ = <c_s> / c_max；positive 选择正极或负极。
// 输出曲线与容量核算都用它表示电极的荷电状态。
double DfnModel::average_stoichiometry(const double* y, bool positive) const {
  const int nx = positive ? layout_.n_x_pos : layout_.n_x_neg;
  const int nr = positive ? layout_.n_r_pos : layout_.n_r_neg;
  const Grid1D& sphere = positive ? sphere_pos_ : sphere_neg_;
  const Grid1D& grid = positive ? pos_grid_ : neg_grid_;
  const double cmax = positive ? params_.pos.c_max : params_.neg.c_max;
  const double porosity =
      positive ? params_.pos.porosity : params_.neg.porosity;

  double particle_volume = 0.0;
  for (double v : sphere.volume) {
    particle_volume += v;
  }

  double lithium = 0.0;
  double volume = 0.0;
  for (int ix = 0; ix < nx; ++ix) {
    double weighted = 0.0;
    for (int ir = 0; ir < nr; ++ir) {
      const int row = positive ? layout_.cs_pos(ix, ir)
                               : layout_.cs_neg(ix, ir);
      weighted += y[row] * sphere.volume[static_cast<std::size_t>(ir)];
    }
    const double solid_volume =
        (1.0 - porosity) * grid.volume[static_cast<std::size_t>(ix)];
    lithium += (weighted / particle_volume) * solid_volume;
    volume += solid_volume;
  }
  return lithium / (volume * cmax);
}

// 对整个电极积分反应电流（A）：Σ a_s F j ΔV。
// 按项目约定应为负极 +I、正极 -I（I 为外加电流），可用于自检。
double DfnModel::electrode_current(const double* y, bool positive) const {
  const int nx = positive ? layout_.n_x_pos : layout_.n_x_neg;
  const int gx0 = positive ? layout_.n_x_neg + layout_.n_x_sep : 0;
  const double volume = params_.area * (positive ? dx_pos_ : dx_neg_);
  double total = 0.0;
  for (int ix = 0; ix < nx; ++ix) {
    const int gx = gx0 + ix;
    const ReactionValue reaction = reaction_at(
        params_, positive ? params_.pos : params_.neg, positive,
        y[layout_.ce(gx)],
        y[positive ? cs_surf_pos(ix) : cs_surf_neg(ix)],
        y[positive ? layout_.phis_pos(ix) : layout_.phis_neg(ix)],
        y[layout_.phie(gx)]);
    total += a_s_[static_cast<std::size_t>(gx)] * kFaraday * reaction.j *
             volume;
  }
  return total;
}

// 建立 Jacobian 的 CSC 稀疏模式。逐块列出所有可能非零的位置：
// 径向三对角（c_s）、轴向三对角（c_e、φ_s、φ_e），以及同一 x 节点上
// c_s 最外层 / c_surf / c_e / φ_s / φ_e 之间的 Butler-Volmer 局部耦合。
// 隔膜只有 c_e 与 φ_e 两条轴向三对角链；模式固定后 CscBuilder 只写值。
void DfnModel::build_pattern() {
  std::vector<std::pair<int, int>> entries;
  entries.reserve(static_cast<std::size_t>(size()) * 7);

  const auto add = [&entries](int row, int col) {
    entries.emplace_back(row, col);
  };

  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;
  const int nx = layout_.n_x_total();

  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int gx = global_index(ix, false);
    for (int ir = 0; ir < nr_neg; ++ir) {
      const int row = layout_.cs_neg(ix, ir);
      add(row, row);
      if (ir > 0) {
        add(row, layout_.cs_neg(ix, ir - 1));
      }
      if (ir + 1 < nr_neg) {
        add(row, layout_.cs_neg(ix, ir + 1));
      }
    }
    const int last_cell = layout_.cs_neg(ix, nr_neg - 1);
    const int surf = cs_surf_neg(ix);
    add(last_cell, surf);
    add(last_cell, layout_.ce(gx));
    add(last_cell, layout_.phis_neg(ix));
    add(last_cell, layout_.phie(gx));

    add(surf, surf);
    add(surf, last_cell);
    add(surf, layout_.ce(gx));
    add(surf, layout_.phis_neg(ix));
    add(surf, layout_.phie(gx));

    add(layout_.ce(gx), layout_.ce(gx));
    if (gx > 0) {
      add(layout_.ce(gx), layout_.ce(gx - 1));
    }
    if (gx + 1 < nx) {
      add(layout_.ce(gx), layout_.ce(gx + 1));
    }
    add(layout_.ce(gx), surf);
    add(layout_.ce(gx), layout_.phis_neg(ix));
    add(layout_.ce(gx), layout_.phie(gx));

    add(layout_.phis_neg(ix), layout_.phis_neg(ix));
    if (ix > 0) {
      add(layout_.phis_neg(ix), layout_.phis_neg(ix - 1));
    }
    if (ix + 1 < layout_.n_x_neg) {
      add(layout_.phis_neg(ix), layout_.phis_neg(ix + 1));
    }
    add(layout_.phis_neg(ix), surf);
    add(layout_.phis_neg(ix), layout_.ce(gx));
    add(layout_.phis_neg(ix), layout_.phie(gx));

    add(layout_.phie(gx), layout_.phie(gx));
    if (gx > 0) {
      add(layout_.phie(gx), layout_.phie(gx - 1));
    }
    if (gx + 1 < nx) {
      add(layout_.phie(gx), layout_.phie(gx + 1));
    }
    add(layout_.phie(gx), surf);
    add(layout_.phie(gx), layout_.ce(gx));
    add(layout_.phie(gx), layout_.phis_neg(ix));
  }

  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int gx = global_index(ix, true);
    for (int ir = 0; ir < nr_pos; ++ir) {
      const int row = layout_.cs_pos(ix, ir);
      add(row, row);
      if (ir > 0) {
        add(row, layout_.cs_pos(ix, ir - 1));
      }
      if (ir + 1 < nr_pos) {
        add(row, layout_.cs_pos(ix, ir + 1));
      }
    }
    const int last_cell = layout_.cs_pos(ix, nr_pos - 1);
    const int surf = cs_surf_pos(ix);
    add(last_cell, surf);
    add(last_cell, layout_.ce(gx));
    add(last_cell, layout_.phis_pos(ix));
    add(last_cell, layout_.phie(gx));

    add(surf, surf);
    add(surf, last_cell);
    add(surf, layout_.ce(gx));
    add(surf, layout_.phis_pos(ix));
    add(surf, layout_.phie(gx));

    add(layout_.ce(gx), layout_.ce(gx));
    if (gx > 0) {
      add(layout_.ce(gx), layout_.ce(gx - 1));
    }
    if (gx + 1 < nx) {
      add(layout_.ce(gx), layout_.ce(gx + 1));
    }
    add(layout_.ce(gx), surf);
    add(layout_.ce(gx), layout_.phis_pos(ix));
    add(layout_.ce(gx), layout_.phie(gx));

    add(layout_.phis_pos(ix), layout_.phis_pos(ix));
    if (ix > 0) {
      add(layout_.phis_pos(ix), layout_.phis_pos(ix - 1));
    }
    if (ix + 1 < layout_.n_x_pos) {
      add(layout_.phis_pos(ix), layout_.phis_pos(ix + 1));
    }
    add(layout_.phis_pos(ix), surf);
    add(layout_.phis_pos(ix), layout_.ce(gx));
    add(layout_.phis_pos(ix), layout_.phie(gx));

    add(layout_.phie(gx), layout_.phie(gx));
    if (gx > 0) {
      add(layout_.phie(gx), layout_.phie(gx - 1));
    }
    if (gx + 1 < nx) {
      add(layout_.phie(gx), layout_.phie(gx + 1));
    }
    add(layout_.phie(gx), surf);
    add(layout_.phie(gx), layout_.ce(gx));
    add(layout_.phie(gx), layout_.phis_pos(ix));
  }

  for (int i = layout_.n_x_neg; i < layout_.n_x_neg + layout_.n_x_sep; ++i) {
    add(layout_.ce(i), layout_.ce(i));
    if (i > 0) {
      add(layout_.ce(i), layout_.ce(i - 1));
    }
    if (i + 1 < nx) {
      add(layout_.ce(i), layout_.ce(i + 1));
    }
    add(layout_.phie(i), layout_.phie(i));
    if (i > 0) {
      add(layout_.phie(i), layout_.phie(i - 1));
    }
    if (i + 1 < nx) {
      add(layout_.phie(i), layout_.phie(i + 1));
    }
  }

  pattern_ = ::dfn::build_pattern(size(), std::move(entries));
}

}  // namespace dfn
