#pragma once

#include <vector>

#include "dfn/grid.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"

namespace dfn {

// Single particle model: one spherical particle per electrode, uniform
// reaction current and no electrolyte. The terminal voltage is
// U_pos - U_neg + eta_pos - eta_neg, evaluated from the surface
// concentrations; there is no spatial x grid.
class SpmModel final : public DaeModel {
 public:
  SpmModel(Params params, int n_r_neg, int n_r_pos);

  int size() const override { return size_; }
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  std::vector<int> variable_types() const override;
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

  const Params& params() const { return params_; }

  void set_current_density(double value) { current_ = value; }
  double current_density() const { return current_; }

  std::vector<double> initial_state() const;
  std::vector<double> suggested_atol() const;

  // Terminal voltage for the current set by set_current_density.
  double terminal_voltage(const double* y) const;
  // Solid-phase lithium of one particle, mol per m^2 of electrode area.
  double particle_lithium(const double* y, bool positive) const;
  double average_stoichiometry(const double* y, bool positive) const;

  int cs_neg(int ir) const { return ir; }
  int cs_pos(int ir) const { return n_r_neg_ + ir; }
  int cs_surf_neg() const { return surf_neg_; }
  int cs_surf_pos() const { return surf_pos_; }

  // Uniform reaction current of one electrode, mol/(m^2 s), positive
  // means lithium leaves the particle.
  double reaction_flux(bool positive) const;

 private:
  double surface_overpotential(const ElectrodeParams& electrode, bool positive,
                               double cs_surf) const;
  double open_circuit(const ElectrodeParams& electrode, bool positive,
                      double cs_surf) const;
  void build_pattern();

  Params params_;
  int n_r_neg_ = 0;
  int n_r_pos_ = 0;
  int surf_neg_ = 0;
  int surf_pos_ = 0;
  int size_ = 0;
  double current_ = 0.0;

  Grid1D sphere_neg_;
  Grid1D sphere_pos_;
  DiffusionOperator radial_neg_;
  DiffusionOperator radial_pos_;

  double surf_coeff_neg_ = 0.0;
  double surf_coeff_pos_ = 0.0;
  double alpha_neg_ = 0.0;
  double alpha_pos_ = 0.0;
  double eq_scale_solid_neg_ = 1.0;
  double eq_scale_solid_pos_ = 1.0;
  double eq_scale_surface_neg_ = 1.0;
  double eq_scale_surface_pos_ = 1.0;

  CscPattern pattern_;
};

}  // namespace dfn
