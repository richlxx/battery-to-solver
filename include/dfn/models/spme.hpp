#pragma once

#include <vector>

#include "dfn/grid.hpp"
#include "dfn/models/cell_model.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"

namespace dfn {

// Single particle model with electrolyte: one spherical particle per
// electrode with uniform reaction current, plus 1D electrolyte concentration
// and potential. The state is [c_s neg, c_s pos, c_surf neg, c_surf pos,
// c_e, phi_e]; phi_e is pinned by a gauge at the negative current collector.
class SpmeModel final : public CellModel {
 public:
  SpmeModel(Params params, int n_x_neg, int n_x_sep, int n_x_pos,
            int n_r_neg, int n_r_pos);

  int size() const override { return size_; }
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  std::vector<int> variable_types() const override;
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

  const Params& params() const override { return params_; }
  const char* name() const override { return "spme"; }

  void set_current_density(double value) override;
  double terminal_current(const double* /*y*/) const override {
    return current_;
  }

  std::vector<double> initial_state() const override;
  std::vector<double> suggested_atol() const override;

  double terminal_voltage(const double* y) const override;
  double average_stoichiometry(const double* y, bool positive) const override;
  double total_lithium(const double* y) const override;
  double solid_lithium(const double* y) const;
  double electrolyte_lithium(const double* y) const;
  double particle_lithium(const double* y, bool positive) const;

  int cs_neg(int ir) const { return ir; }
  int cs_pos(int ir) const { return n_r_neg_ + ir; }
  int cs_surf_neg() const { return cs_pos_offset_ + n_r_pos_; }
  int cs_surf_pos() const { return cs_surf_neg() + 1; }
  int ce(int ix) const { return ce_offset_ + ix; }
  int phie(int ix) const { return phie_offset_ + ix; }

  // Uniform reaction flux of one electrode, mol/(m^2 s), positive means
  // lithium leaves the particle.
  double reaction_flux(bool positive) const {
    return positive ? j_pos_ : j_neg_;
  }

 private:
  double open_circuit(const ElectrodeParams& electrode, bool positive,
                      double cs_surf) const;
  double electrolyte_source(int ix) const;
  void build_pattern();

  Params params_;
  int n_x_neg_ = 0;
  int n_x_sep_ = 0;
  int n_x_pos_ = 0;
  int n_x_total_ = 0;
  int n_r_neg_ = 0;
  int n_r_pos_ = 0;

  int cs_pos_offset_ = 0;
  int ce_offset_ = 0;
  int phie_offset_ = 0;
  int size_ = 0;

  double current_ = 0.0;
  double j_neg_ = 0.0;
  double j_pos_ = 0.0;
  double a_s_neg_ = 0.0;
  double a_s_pos_ = 0.0;

  Grid1D x_grid_;
  Grid1D sphere_neg_;
  Grid1D sphere_pos_;

  DiffusionOperator radial_neg_;
  DiffusionOperator radial_pos_;
  DiffusionOperator elyte_op_;
  TridiagonalMatrix phie_op_;

  std::vector<double> eps_e_;
  std::vector<double> a_s_;

  double surf_coeff_neg_ = 0.0;
  double surf_coeff_pos_ = 0.0;
  double alpha_neg_ = 0.0;
  double alpha_pos_ = 0.0;

  double eq_scale_solid_neg_ = 1.0;
  double eq_scale_solid_pos_ = 1.0;
  double eq_scale_electrolyte_ = 1.0;
  double eq_scale_potential_ = 1.0;
  double eq_scale_surface_ = 1.0;

  CscPattern pattern_;
};

}  // namespace dfn
