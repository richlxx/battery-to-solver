#include "dfn/reduction/spm.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "dfn/thermo.hpp"

namespace dfn {

namespace {

double clamp_stoich(double cs, double cmax) {
  return std::clamp(cs, 1e-6 * cmax, cmax * (1.0 - 1e-6));
}

}  // namespace

SpmModel::SpmModel(Params params, int n_r_neg, int n_r_pos)
    : params_(std::move(params)), n_r_neg_(n_r_neg), n_r_pos_(n_r_pos) {
  validate(params_);
  if (n_r_neg <= 0 || n_r_pos <= 0) {
    throw std::invalid_argument("SPM radial counts must be positive");
  }

  sphere_neg_ = make_uniform_sphere_grid(n_r_neg_, params_.neg.radius);
  sphere_pos_ = make_uniform_sphere_grid(n_r_pos_, params_.pos.radius);
  radial_neg_ = assemble_diffusion(
      sphere_neg_, std::vector<double>(static_cast<std::size_t>(n_r_neg_),
                                       params_.neg.d_s));
  radial_pos_ = assemble_diffusion(
      sphere_pos_, std::vector<double>(static_cast<std::size_t>(n_r_pos_),
                                       params_.pos.d_s));

  surf_coeff_neg_ = -sphere_neg_.area[static_cast<std::size_t>(n_r_neg_)] /
                    sphere_neg_.volume[static_cast<std::size_t>(n_r_neg_ - 1)];
  surf_coeff_pos_ = -sphere_pos_.area[static_cast<std::size_t>(n_r_pos_)] /
                    sphere_pos_.volume[static_cast<std::size_t>(n_r_pos_ - 1)];

  const double dr_neg = params_.neg.radius / static_cast<double>(n_r_neg_);
  const double dr_pos = params_.pos.radius / static_cast<double>(n_r_pos_);
  alpha_neg_ = dr_neg / (2.0 * params_.neg.d_s);
  alpha_pos_ = dr_pos / (2.0 * params_.pos.d_s);

  eq_scale_solid_neg_ = params_.neg.radius * params_.neg.radius /
                        params_.neg.d_s / params_.neg.c_max;
  eq_scale_solid_pos_ = params_.pos.radius * params_.pos.radius /
                        params_.pos.d_s / params_.pos.c_max;
  eq_scale_surface_neg_ = 1.0 / params_.neg.c_max;
  eq_scale_surface_pos_ = 1.0 / params_.pos.c_max;

  surf_neg_ = n_r_neg_ + n_r_pos_;
  surf_pos_ = surf_neg_ + 1;
  size_ = surf_pos_ + 1;

  build_pattern();
}

double SpmModel::reaction_flux(bool positive) const {
  const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
  const double a_s = 3.0 * (1.0 - electrode.porosity) / electrode.radius;
  const double sign = positive ? -1.0 : 1.0;
  return sign * current_ / (a_s * kFaraday * electrode.thickness);
}

double SpmModel::open_circuit(const ElectrodeParams& electrode, bool positive,
                              double cs_surf) const {
  const double cs = clamp_stoich(cs_surf, electrode.c_max);
  return positive ? positive_ocp(cs / electrode.c_max).value
                  : negative_ocp(cs / electrode.c_max).value;
}

double SpmModel::surface_overpotential(const ElectrodeParams& electrode,
                                       bool positive, double cs_surf) const {
  const double j = reaction_flux(positive);
  if (j == 0.0) {
    return 0.0;
  }
  const double cs = clamp_stoich(cs_surf, electrode.c_max);
  const double i0 = kFaraday * electrode.k *
                    std::pow(params_.elyte.c0, electrode.alpha_a) *
                    std::pow(electrode.c_max - cs, electrode.alpha_a) *
                    std::pow(cs, electrode.alpha_c);
  if (i0 <= 0.0) {
    return 0.0;
  }
  const double a =
      electrode.alpha_a * kFaraday / (kGasConstant * params_.temperature);
  const double c =
      electrode.alpha_c * kFaraday / (kGasConstant * params_.temperature);

  double eta = 0.0;
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

void SpmModel::residual(double /*t*/, const double* y, const double* yp,
                        double* out) const {
  std::fill(out, out + size_, 0.0);

  const double j_neg = reaction_flux(false);
  const double j_pos = reaction_flux(true);

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

  out[surf_neg_] =
      (y[surf_neg_] - y[cs_neg(n_r_neg_ - 1)] + alpha_neg_ * j_neg) *
      eq_scale_surface_neg_;
  out[surf_pos_] =
      (y[surf_pos_] - y[cs_pos(n_r_pos_ - 1)] + alpha_pos_ * j_pos) *
      eq_scale_surface_pos_;
}

std::vector<int> SpmModel::variable_types() const {
  std::vector<int> types(static_cast<std::size_t>(size_), 1);
  types[static_cast<std::size_t>(surf_neg_)] = 0;
  types[static_cast<std::size_t>(surf_pos_)] = 0;
  return types;
}

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

  out.add(surf_neg_, surf_neg_, eq_scale_surface_neg_);
  out.add(surf_neg_, cs_neg(n_r_neg_ - 1), -eq_scale_surface_neg_);
  out.add(surf_pos_, surf_pos_, eq_scale_surface_pos_);
  out.add(surf_pos_, cs_pos(n_r_pos_ - 1), -eq_scale_surface_pos_);
}

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

std::vector<double> SpmModel::suggested_atol() const {
  std::vector<double> atol(static_cast<std::size_t>(size_), 1e-2);
  atol[static_cast<std::size_t>(surf_neg_)] = 1e-2;
  atol[static_cast<std::size_t>(surf_pos_)] = 1e-2;
  return atol;
}

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
  return (weighted / particle_volume) * (1.0 - electrode.porosity) *
         electrode.thickness;
}

double SpmModel::average_stoichiometry(const double* y, bool positive) const {
  const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
  return particle_lithium(y, positive) /
         ((1.0 - electrode.porosity) * electrode.thickness * electrode.c_max);
}

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
