#include "dfn/models/dfn.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "dfn/thermo.hpp"

namespace dfn {

namespace {

constexpr double kPi = 3.14159265358979323846;

}  // namespace

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
  radial_neg_ = assemble_diffusion(sphere_neg_, d_s_neg);
  radial_pos_ = assemble_diffusion(sphere_pos_, d_s_pos);

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
  phis_neg_op_ = assemble_diffusion(neg_grid_, sigma_neg).matrix;
  phis_pos_op_ = assemble_diffusion(pos_grid_, sigma_pos).matrix;

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
    d_eff_[static_cast<std::size_t>(i)] =
        params_.elyte.d_e * std::pow(epsy, bruggeman);
    kappa_eff_[static_cast<std::size_t>(i)] =
        params_.elyte.kappa * std::pow(epsy, bruggeman);
  }
  elyte_op_ = assemble_diffusion(x_grid_, d_eff_);
  phie_op_ = assemble_diffusion(x_grid_, kappa_eff_).matrix;

  dx_neg_ = l_neg / static_cast<double>(layout_.n_x_neg);
  dx_pos_ = l_pos / static_cast<double>(layout_.n_x_pos);
  neg_dirichlet_ = -2.0 * sigma_neg_eff_ / (dx_neg_ * dx_neg_);
  pos_dirichlet_ = -2.0 * sigma_pos_eff_ / (dx_pos_ * dx_pos_);
  pos_flux_coeff_ = -1.0 / dx_pos_;
  surf_coeff_neg_ = -sphere_neg_.area[static_cast<std::size_t>(layout_.n_r_neg)] /
                    sphere_neg_.volume[static_cast<std::size_t>(layout_.n_r_neg - 1)];
  surf_coeff_pos_ = -sphere_pos_.area[static_cast<std::size_t>(layout_.n_r_pos)] /
                    sphere_pos_.volume[static_cast<std::size_t>(layout_.n_r_pos - 1)];
  const double dr_neg = params_.neg.radius / static_cast<double>(layout_.n_r_neg);
  const double dr_pos = params_.pos.radius / static_cast<double>(layout_.n_r_pos);
  alpha_neg_ = dr_neg / (2.0 * params_.neg.d_s);
  alpha_pos_ = dr_pos / (2.0 * params_.pos.d_s);

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
  // Reference current density 1 A/m^2; the exact value only affects scaling.
  eq_scale_potential_ = length;

  cs_surf_neg_offset_ = layout_.size();
  cs_surf_pos_offset_ = layout_.size() + layout_.n_x_neg;
  size_ = layout_.size() + layout_.n_x_neg + layout_.n_x_pos;
  eq_scale_surface_ =
      1.0 / std::max(params_.neg.c_max, params_.pos.c_max);

  build_pattern();
}

int DfnModel::global_index(int local_ix, bool positive) const {
  return positive ? layout_.n_x_neg + layout_.n_x_sep + local_ix : local_ix;
}

DfnModel::ReactionValue DfnModel::reaction_at(const ElectrodeParams& electrode,
                                              bool positive, double ce,
                                              double cs_surf, double phis,
                                              double phie) const {
  const double cmax = electrode.c_max;
  const double cs = std::clamp(cs_surf, 1e-6 * cmax, cmax * (1.0 - 1e-6));
  const double ce_clamped = std::max(ce, 1e-6 * params_.elyte.c0);
  const double stoich = cs / cmax;
  const OcpValue ocp =
      positive ? positive_ocp(stoich) : negative_ocp(stoich);

  const double i0 = kFaraday * electrode.k *
                    std::pow(ce_clamped, electrode.alpha_a) *
                    std::pow(cmax - cs, electrode.alpha_a) *
                    std::pow(cs, electrode.alpha_c);
  const double eta = phis - phie - ocp.value;
  const double a =
      electrode.alpha_a * kFaraday / (kGasConstant * params_.temperature);
  const double c =
      electrode.alpha_c * kFaraday / (kGasConstant * params_.temperature);
  const double ea = std::exp(std::clamp(a * eta, -60.0, 60.0));
  const double ec = std::exp(std::clamp(-c * eta, -60.0, 60.0));

  ReactionValue result;
  result.j = i0 / kFaraday * (ea - ec);
  result.dj_deta = i0 / kFaraday * (a * ea + c * ec);
  result.dj_dce = result.j * electrode.alpha_a / ce_clamped;
  result.dj_dcs =
      result.j * (-electrode.alpha_a / (cmax - cs) +
                  electrode.alpha_c / cs) -
      result.dj_deta * ocp.slope / cmax;
  return result;
}

void DfnModel::residual(double /*t*/, const double* y, const double* yp,
                        double* out) const {
  const int n = size();
  std::fill(out, out + n, 0.0);

  const double f = kFaraday;
  const double one_minus_t = 1.0 - params_.elyte.t_plus;
  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;

  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int gx = global_index(ix, false);
    neg_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_.neg, false, y[layout_.ce(gx)],
        y[cs_surf_neg(ix)], y[layout_.phis_neg(ix)],
        y[layout_.phie(gx)]);
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int gx = global_index(ix, true);
    pos_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_.pos, true, y[layout_.ce(gx)],
        y[cs_surf_pos(ix)], y[layout_.phis_pos(ix)],
        y[layout_.phie(gx)]);
  }

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

  // Surface concentration: c_surf = c_last - alpha * j.
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

void DfnModel::jacobian(double /*t*/, const double* y, const double* /*yp*/,
                        double cj, CscBuilder& out) const {
  const double f = kFaraday;
  const double one_minus_t = 1.0 - params_.elyte.t_plus;
  const int nr_neg = layout_.n_r_neg;
  const int nr_pos = layout_.n_r_pos;

  for (int ix = 0; ix < layout_.n_x_neg; ++ix) {
    const int gx = global_index(ix, false);
    neg_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_.neg, false, y[layout_.ce(gx)],
        y[cs_surf_neg(ix)], y[layout_.phis_neg(ix)],
        y[layout_.phie(gx)]);
  }
  for (int ix = 0; ix < layout_.n_x_pos; ++ix) {
    const int gx = global_index(ix, true);
    pos_reaction_[static_cast<std::size_t>(ix)] = reaction_at(
        params_.pos, true, y[layout_.ce(gx)],
        y[cs_surf_pos(ix)], y[layout_.phis_pos(ix)],
        y[layout_.phie(gx)]);
  }

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

  // Surface-concentration equations.
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

void DfnModel::set_current_density(double value) {
  options_.mode = ControlMode::current;
  options_.applied_current_density = value;
}

void DfnModel::set_applied_voltage(double value) {
  options_.mode = ControlMode::voltage;
  options_.applied_voltage = value;
}

void DfnModel::set_control(ControlMode mode, double value) {
  if (mode == ControlMode::current) {
    set_current_density(value);
  } else {
    set_applied_voltage(value);
  }
}

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

double DfnModel::terminal_voltage(const double* y) const {
  if (options_.mode == ControlMode::voltage) {
    return options_.applied_voltage;
  }
  const int last = layout_.n_x_pos - 1;
  const double phi_last = y[layout_.phis_pos(last)];
  return phi_last -
         options_.applied_current_density * dx_pos_ / (2.0 * sigma_pos_eff_);
}

double DfnModel::terminal_current(const double* y) const {
  if (options_.mode == ControlMode::current) {
    return options_.applied_current_density;
  }
  const int last = layout_.n_x_pos - 1;
  const double phi_last = y[layout_.phis_pos(last)];
  return 2.0 * sigma_pos_eff_ *
         (phi_last - options_.applied_voltage) / dx_pos_;
}

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

double DfnModel::electrolyte_lithium(const double* y) const {
  double total = 0.0;
  for (int i = 0; i < layout_.n_x_total(); ++i) {
    total += eps_e_[static_cast<std::size_t>(i)] * y[layout_.ce(i)] *
             x_grid_.volume[static_cast<std::size_t>(i)];
  }
  return total;
}

double DfnModel::total_lithium(const double* y) const {
  return solid_lithium(y) + electrolyte_lithium(y);
}

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

double DfnModel::electrode_current(const double* y, bool positive) const {
  const int nx = positive ? layout_.n_x_pos : layout_.n_x_neg;
  const int gx0 = positive ? layout_.n_x_neg + layout_.n_x_sep : 0;
  const double volume = params_.area * (positive ? dx_pos_ : dx_neg_);
  double total = 0.0;
  for (int ix = 0; ix < nx; ++ix) {
    const int gx = gx0 + ix;
    const ReactionValue reaction = reaction_at(
        positive ? params_.pos : params_.neg, positive, y[layout_.ce(gx)],
        y[positive ? cs_surf_pos(ix) : cs_surf_neg(ix)],
        y[positive ? layout_.phis_pos(ix) : layout_.phis_neg(ix)],
        y[layout_.phie(gx)]);
    total += a_s_[static_cast<std::size_t>(gx)] * kFaraday * reaction.j *
             volume;
  }
  return total;
}

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
