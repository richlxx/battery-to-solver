#pragma once

#include <vector>

#include "dfn/solvers/dae.hpp"

namespace dfn {

struct ChainDaeOptions {
  int cells = 1;
  double lambda = 100.0;
  double diffusion = 0.0;
};

// Teaching DAE: one stiff relaxation equation per cell plus one algebraic
// variable per cell. Unknowns are interleaved as (c_0, w_0, c_1, w_1, ...).
//
//   c_i' + lambda (c_i - w_i) - diffusion * lap(c)_i = 0
//   w_i - sin(t) = 0
//
// With one cell and no diffusion the solution is known in closed form.
class ChainDae final : public DaeModel {
 public:
  explicit ChainDae(ChainDaeOptions options = {});

  int cells() const { return options_.cells; }
  double lambda() const { return options_.lambda; }
  double diffusion() const { return options_.diffusion; }

  int size() const override { return 2 * options_.cells; }
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  std::vector<int> variable_types() const override;
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

 private:
  ChainDaeOptions options_;
  CscPattern pattern_;
};

// Uniform profile that satisfies both equations at t0 when paired with
// chain_dae_consistent_derivative(): every cell takes the single-cell
// closed-form concentration, every algebraic variable takes sin(t0).
std::vector<double> chain_dae_consistent_state(const ChainDae& model, double t0);

// Matching derivative of chain_dae_consistent_state, so that F(t0, y, y') = 0.
std::vector<double> chain_dae_consistent_derivative(const ChainDae& model,
                                                    double t0);

// Sine-shaped concentration perturbation on top of the uniform profile; pair
// it with chain_dae_derivative_from_state() to get F(t0, y, y') = 0:
//   c_i = exact(t0) + amplitude * sin((i + 1) * pi / (cells + 1))
std::vector<double> chain_dae_mode_state(const ChainDae& model, double t0,
                                         double amplitude);

// Given y, returns the y' that makes F(t, y, y') = 0.
std::vector<double> chain_dae_derivative_from_state(
    const ChainDae& model, double t, const std::vector<double>& y);

// Exact concentration for a single cell with diffusion = 0.
double chain_dae_exact_concentration(double lambda, double t);

// Time derivative of the exact concentration.
double chain_dae_exact_concentration_rate(double lambda, double t);

}  // namespace dfn
