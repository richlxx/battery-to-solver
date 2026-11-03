#include "dfn/models/spme.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "dfn/models/reaction.hpp"
#include "dfn/thermo.hpp"

namespace dfn {

namespace {

double clamp_stoich(double cs, double cmax) {
  return std::clamp(cs, 1e-6 * cmax, cmax * (1.0 - 1e-6));
}

}  // namespace

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

  x_grid_ = make_segmented_line_grid({n_x_neg_, n_x_sep_, n_x_pos_},
                                     {l_neg, l_sep, l_pos}, params_.area);
  sphere_neg_ = make_uniform_sphere_grid(n_r_neg_, params_.neg.radius);
  sphere_pos_ = make_uniform_sphere_grid(n_r_pos_, params_.pos.radius);

  radial_neg_ = assemble_diffusion(
      sphere_neg_, std::vector<double>(static_cast<std::size_t>(n_r_neg_),
                                       params_.neg.d_s));
  radial_pos_ = assemble_diffusion(
      sphere_pos_, std::vector<double>(static_cast<std::size_t>(n_r_pos_),
                                       params_.pos.d_s));

  a_s_neg_ = 3.0 * (1.0 - params_.neg.porosity) / params_.neg.radius;
  a_s_pos_ = 3.0 * (1.0 - params_.pos.porosity) / params_.pos.radius;

  eps_e_.assign(static_cast<std::size_t>(n_x_total_), 0.0);
  std::vector<double> d_eff(static_cast<std::size_t>(n_x_total_), 0.0);
  std::vector<double> kappa_eff(static_cast<std::size_t>(n_x_total_), 0.0);
  a_s_.assign(static_cast<std::size_t>(n_x_total_), 0.0);
  for (int i = 0; i < n_x_total_; ++i) {
    double epsy = 0.0;
    double bruggeman = 1.5;
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
    d_eff[static_cast<std::size_t>(i)] =
        params_.elyte.d_e * std::pow(epsy, bruggeman);
    kappa_eff[static_cast<std::size_t>(i)] =
        params_.elyte.kappa * std::pow(epsy, bruggeman);
  }
  elyte_op_ = assemble_diffusion(x_grid_, d_eff);
  phie_op_ = assemble_diffusion(x_grid_, kappa_eff).matrix;

  surf_coeff_neg_ = -sphere_neg_.area[static_cast<std::size_t>(n_r_neg_)] /
                    sphere_neg_.volume[static_cast<std::size_t>(n_r_neg_ - 1)];
  surf_coeff_pos_ = -sphere_pos_.area[static_cast<std::size_t>(n_r_pos_)] /
                    sphere_pos_.volume[static_cast<std::size_t>(n_r_pos_ - 1)];
  const double dr_neg = params_.neg.radius / static_cast<double>(n_r_neg_);
  const double dr_pos = params_.pos.radius / static_cast<double>(n_r_pos_);
  alpha_neg_ = dr_neg / (2.0 * params_.neg.d_s);
  alpha_pos_ = dr_pos / (2.0 * params_.pos.d_s);

  const double tau_s_neg =
      params_.neg.radius * params_.neg.radius / params_.neg.d_s;
  const double tau_s_pos =
      params_.pos.radius * params_.pos.radius / params_.pos.d_s;
  eq_scale_solid_neg_ = tau_s_neg / params_.neg.c_max;
  eq_scale_solid_pos_ = tau_s_pos / params_.pos.c_max;
  const double length = params_.cell_length();
  const double tau_e = length * length / params_.elyte.d_e;
  eq_scale_electrolyte_ =
      params_.elyte.c0 > 0.0 ? tau_e / params_.elyte.c0 : 1.0;
  eq_scale_potential_ = length;
  eq_scale_surface_ =
      1.0 / std::max(params_.neg.c_max, params_.pos.c_max);

  cs_pos_offset_ = n_r_neg_;
  ce_offset_ = cs_surf_pos() + 1;
  phie_offset_ = ce_offset_ + n_x_total_;
  size_ = phie_offset_ + n_x_total_;

  set_current_density(0.0);
  build_pattern();
}

void SpmeModel::set_current_density(double value) {
  current_ = value;
  j_neg_ = current_ / (a_s_neg_ * kFaraday * params_.neg.thickness);
  j_pos_ = -current_ / (a_s_pos_ * kFaraday * params_.pos.thickness);
}

double SpmeModel::open_circuit(const ElectrodeParams& electrode, bool positive,
                               double cs_surf) const {
  const double cs = clamp_stoich(cs_surf, electrode.c_max);
  return positive ? positive_ocp(cs / electrode.c_max).value
                  : negative_ocp(cs / electrode.c_max).value;
}

void SpmeModel::residual(double /*t*/, const double* y, const double* yp,
                         double* out) const {
  std::fill(out, out + size_, 0.0);

  for (int ir = 0; ir < n_r_neg_; ++ir) {
    double lap = radial_neg_.matrix.diag[static_cast<std::size_t>(ir)] * y[ir];
    if (ir > 0) {
      lap += radial_neg_.matrix.lower[static_cast<std::size_t>(ir)] *
             y[cs_neg(ir - 1)];
    }
    if (ir + 1 < n_r_neg_) {
      lap += radial_neg_.matrix.upper[static_cast<std::size_t>(ir)] *
             y[cs_neg(ir + 1)];
    }
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

  out[cs_surf_neg()] =
      (y[cs_surf_neg()] - y[cs_neg(n_r_neg_ - 1)] + alpha_neg_ * j_neg_) *
      eq_scale_surface_;
  out[cs_surf_pos()] =
      (y[cs_surf_pos()] - y[cs_pos(n_r_pos_ - 1)] + alpha_pos_ * j_pos_) *
      eq_scale_surface_;

  for (int i = 0; i < n_x_total_; ++i) {
    const int row = ce(i);
    double lap = elyte_op_.matrix.diag[static_cast<std::size_t>(i)] * y[row];
    if (i > 0) {
      lap += elyte_op_.matrix.lower[static_cast<std::size_t>(i)] *
             y[ce(i - 1)];
    }
    if (i + 1 < n_x_total_) {
      lap += elyte_op_.matrix.upper[static_cast<std::size_t>(i)] *
             y[ce(i + 1)];
    }
    out[row] =
        (eps_e_[static_cast<std::size_t>(i)] * yp[row] - lap -
         electrolyte_source(i)) *
        eq_scale_electrolyte_;
  }

  out[phie(0)] = y[phie(0)] * eq_scale_potential_;
  for (int i = 1; i < n_x_total_; ++i) {
    const int row = phie(i);
    double lap = phie_op_.diag[static_cast<std::size_t>(i)] * y[row];
    if (i > 0) {
      lap += phie_op_.lower[static_cast<std::size_t>(i)] * y[phie(i - 1)];
    }
    if (i + 1 < n_x_total_) {
      lap += phie_op_.upper[static_cast<std::size_t>(i)] * y[phie(i + 1)];
    }
    const double source =
        a_s_[static_cast<std::size_t>(i)] * kFaraday *
        (i < n_x_neg_ ? j_neg_ : j_pos_);
    out[row] = (lap + source) * eq_scale_potential_;
  }
}

double SpmeModel::electrolyte_source(int ix) const {
  const double j = ix < n_x_neg_ ? j_neg_ : j_pos_;
  return (1.0 - params_.elyte.t_plus) * a_s_[static_cast<std::size_t>(ix)] * j;
}

std::vector<int> SpmeModel::variable_types() const {
  std::vector<int> types(static_cast<std::size_t>(size_), 1);
  types[static_cast<std::size_t>(cs_surf_neg())] = 0;
  types[static_cast<std::size_t>(cs_surf_pos())] = 0;
  for (int i = phie_offset_; i < size_; ++i) {
    types[static_cast<std::size_t>(i)] = 0;
  }
  return types;
}

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

  out.add(cs_surf_neg(), cs_surf_neg(), eq_scale_surface_);
  out.add(cs_surf_neg(), cs_neg(n_r_neg_ - 1), -eq_scale_surface_);
  out.add(cs_surf_pos(), cs_surf_pos(), eq_scale_surface_);
  out.add(cs_surf_pos(), cs_pos(n_r_pos_ - 1), -eq_scale_surface_);

  for (int i = 0; i < n_x_total_; ++i) {
    const int row = ce(i);
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

double SpmeModel::terminal_voltage(const double* y) const {
  const double csn = y[cs_surf_neg()];
  const double csp = y[cs_surf_pos()];
  const double u_neg = open_circuit(params_.neg, false, csn);
  const double u_pos = open_circuit(params_.pos, true, csp);

  const int ix_neg = n_x_neg_ / 2;
  const int ix_pos = n_x_neg_ + n_x_sep_ + n_x_pos_ / 2;
  const double eta_neg = overpotential_from_flux(
      params_, params_.neg, false, y[ce(ix_neg)], csn, j_neg_);
  const double eta_pos = overpotential_from_flux(
      params_, params_.pos, true, y[ce(ix_pos)], csp, j_pos_);

  const double phie_drop =
      y[phie(n_x_total_ - 1)] - y[phie(0)];
  return u_pos - u_neg + eta_pos - eta_neg + phie_drop;
}

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

double SpmeModel::solid_lithium(const double* y) const {
  return particle_lithium(y, false) + particle_lithium(y, true);
}

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

double SpmeModel::total_lithium(const double* y) const {
  return solid_lithium(y) + electrolyte_lithium(y);
}

double SpmeModel::average_stoichiometry(const double* y,
                                        bool positive) const {
  const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
  return particle_lithium(y, positive) /
         ((1.0 - electrode.porosity) * electrode.thickness * electrode.c_max);
}

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
