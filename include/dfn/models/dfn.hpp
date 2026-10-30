#pragma once

#include <vector>

#include "dfn/grid.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"
#include "dfn/state.hpp"

namespace dfn {

enum class ControlMode { current, voltage };

struct DfnOptions {
  ControlMode mode = ControlMode::current;
  double applied_current_density = 0.0;  // A/m^2, positive = discharge
  double applied_voltage = 4.0;          // V, used in voltage mode
};

// Isothermal DFN with the equations of chapter 3. The concentrated-solution
// correction kappa_D is set to zero; see chapter 9 for the rationale.
class DfnModel final : public DaeModel {
 public:
  DfnModel(Params params, StateLayout layout, DfnOptions options = {});

  int size() const override { return size_; }
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  std::vector<int> variable_types() const override;
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

  const Params& params() const { return params_; }
  const StateLayout& layout() const { return layout_; }
  const Grid1D& x_grid() const { return x_grid_; }

  // Auxiliary algebraic unknowns: surface concentration of each particle.
  int cs_surf_neg(int ix) const { return cs_surf_neg_offset_ + ix; }
  int cs_surf_pos(int ix) const { return cs_surf_pos_offset_ + ix; }

  ControlMode mode() const { return options_.mode; }
  double applied_current_density() const {
    return options_.applied_current_density;
  }
  double applied_voltage() const { return options_.applied_voltage; }
  void set_current_density(double value);
  void set_applied_voltage(double value);
  void set_control(ControlMode mode, double value);

  std::vector<double> initial_state() const;
  std::vector<double> suggested_atol() const;

  // Terminal voltage at the positive current-collector face.
  double terminal_voltage(const double* y) const;
  // Total current density through the cell, A/m^2 (positive = discharge).
  double terminal_current(const double* y) const;

  double total_lithium(const double* y) const;
  double solid_lithium(const double* y) const;
  double electrolyte_lithium(const double* y) const;
  // Volume-averaged stoichiometry c_s / c_max of one electrode.
  double average_stoichiometry(const double* y, bool positive) const;
  // Integrated reaction current over one electrode, in amperes; should be
  // +I for the negative electrode and -I for the positive one.
  double electrode_current(const double* y, bool positive) const;

 private:
  struct ReactionValue {
    double j = 0.0;       // mol/(m^2 s)
    double dj_deta = 0.0;
    double dj_dce = 0.0;
    double dj_dcs = 0.0;
  };

  ReactionValue reaction_at(const ElectrodeParams& electrode, bool positive,
                            double ce, double cs_surf, double phis,
                            double phie) const;
  int global_index(int local_ix, bool positive) const;
  void build_pattern();

  Params params_;
  StateLayout layout_;
  DfnOptions options_;

  Grid1D x_grid_;
  Grid1D sphere_neg_;
  Grid1D sphere_pos_;
  Grid1D neg_grid_;
  Grid1D pos_grid_;

  DiffusionOperator radial_neg_;
  DiffusionOperator radial_pos_;
  TridiagonalMatrix phis_neg_op_;
  TridiagonalMatrix phis_pos_op_;
  DiffusionOperator elyte_op_;
  TridiagonalMatrix phie_op_;

  std::vector<double> eps_e_;      // per x node
  std::vector<double> d_eff_;      // per x node
  std::vector<double> kappa_eff_;  // per x node
  std::vector<double> a_s_;        // per x node, zero in the separator

  double sigma_neg_eff_ = 0.0;
  double sigma_pos_eff_ = 0.0;
  double dx_neg_ = 0.0;
  double dx_pos_ = 0.0;
  double neg_dirichlet_ = 0.0;  // -2 sigma / dx^2
  double pos_dirichlet_ = 0.0;  // -2 sigma / dx^2
  double pos_flux_coeff_ = 0.0; // -1 / dx
  double surf_coeff_neg_ = 0.0;
  double surf_coeff_pos_ = 0.0;
  double alpha_neg_ = 0.0;  // dr / (2 D_s), surface reconstruction
  double alpha_pos_ = 0.0;

  double eq_scale_solid_neg_ = 1.0;
  double eq_scale_solid_pos_ = 1.0;
  double eq_scale_electrolyte_ = 1.0;
  double eq_scale_potential_ = 1.0;
  double eq_scale_surface_ = 1.0;

  int size_ = 0;
  int cs_surf_neg_offset_ = 0;
  int cs_surf_pos_offset_ = 0;

  mutable std::vector<ReactionValue> neg_reaction_;
  mutable std::vector<ReactionValue> pos_reaction_;

  CscPattern pattern_;
};

}  // namespace dfn
