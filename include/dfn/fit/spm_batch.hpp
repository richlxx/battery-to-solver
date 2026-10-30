#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "dfn/grid.hpp"
#include "dfn/linalg/tridiagonal.hpp"
#include "dfn/operators.hpp"
#include "dfn/params.hpp"
#include "dfn/thermo.hpp"

namespace dfn {
namespace fit {

inline constexpr double kFaraday = 96485.33212;
inline constexpr double kGasConstant = 8.314462618;

struct BatchSegment {
  double duration = 0.0;
  double current = 0.0;  // A/m^2, positive = discharge
};

struct BatchRow {
  double t = 0.0;
  double current = 0.0;
  double voltage = 0.0;
};

inline double clamp_stoich(double cs, double cmax) {
  return std::clamp(cs, 1e-6 * cmax, cmax * (1.0 - 1e-6));
}

inline double surface_overpotential(const ElectrodeParams& electrode,
                                    double j, double cs_surf, double c_e,
                                    double temperature) {
  if (j == 0.0) {
    return 0.0;
  }
  const double cs = clamp_stoich(cs_surf, electrode.c_max);
  const double i0 =
      kFaraday * electrode.k * std::pow(c_e, electrode.alpha_a) *
      std::pow(electrode.c_max - cs, electrode.alpha_a) *
      std::pow(cs, electrode.alpha_c);
  if (i0 <= 0.0) {
    return 0.0;
  }
  const double a = electrode.alpha_a * kFaraday / (kGasConstant * temperature);
  const double c = electrode.alpha_c * kFaraday / (kGasConstant * temperature);
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

struct ElectrodeGeometry {
  TridiagonalMatrix diffusion;
  TridiagonalMatrix implicit;
  double surf_coeff = 0.0;
  double alpha = 0.0;
  double j_coeff = 0.0;
};

inline ElectrodeGeometry make_geometry(const ElectrodeParams& electrode, int nr,
                                       double dt) {
  const Grid1D grid = make_uniform_sphere_grid(nr, electrode.radius);
  const DiffusionOperator op = assemble_diffusion(
      grid, std::vector<double>(static_cast<std::size_t>(nr), electrode.d_s));
  ElectrodeGeometry geometry;
  geometry.diffusion = op.matrix;
  geometry.surf_coeff =
      -grid.area[static_cast<std::size_t>(nr)] /
      grid.volume[static_cast<std::size_t>(nr - 1)];
  const double dr = electrode.radius / static_cast<double>(nr);
  geometry.alpha = dr / (2.0 * electrode.d_s);
  const double a_s = 3.0 * (1.0 - electrode.porosity) / electrode.radius;
  geometry.j_coeff = 1.0 / (a_s * kFaraday * electrode.thickness);

  geometry.implicit.n = nr;
  geometry.implicit.lower.resize(static_cast<std::size_t>(nr));
  geometry.implicit.diag.resize(static_cast<std::size_t>(nr));
  geometry.implicit.upper.resize(static_cast<std::size_t>(nr));
  for (int i = 0; i < nr; ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    geometry.implicit.lower[k] = -dt * geometry.diffusion.lower[k];
    geometry.implicit.diag[k] = 1.0 - dt * geometry.diffusion.diag[k];
    geometry.implicit.upper[k] = -dt * geometry.diffusion.upper[k];
  }
  return geometry;
}

// Fixed-step backward-Euler SPM used by the fitting lab and the CUDA
// reference. SPM has no x grid; each electrode is one particle.
class SpmBatch {
 public:
  SpmBatch(Params params, int nr, double dt)
      : params_(std::move(params)),
        nr_(nr),
        dt_(dt),
        neg_(make_geometry(params_.neg, nr, dt)),
        pos_(make_geometry(params_.pos, nr, dt)),
        c_neg_(static_cast<std::size_t>(nr), params_.neg.c0),
        c_pos_(static_cast<std::size_t>(nr), params_.pos.c0) {}

  void reset() {
    std::fill(c_neg_.begin(), c_neg_.end(), params_.neg.c0);
    std::fill(c_pos_.begin(), c_pos_.end(), params_.pos.c0);
  }

  double step(double current) {
    const double j_neg = current * neg_.j_coeff;
    const double j_pos = -current * pos_.j_coeff;

    std::vector<double> rhs = c_neg_;
    rhs[static_cast<std::size_t>(nr_ - 1)] += dt_ * neg_.surf_coeff * j_neg;
    c_neg_ = solve_tridiagonal(neg_.implicit, rhs);

    rhs = c_pos_;
    rhs[static_cast<std::size_t>(nr_ - 1)] += dt_ * pos_.surf_coeff * j_pos;
    c_pos_ = solve_tridiagonal(pos_.implicit, rhs);
    return voltage(current);
  }

  double voltage(double current) const {
    const double j_neg = current * neg_.j_coeff;
    const double j_pos = -current * pos_.j_coeff;
    const double cs_neg = clamp_stoich(
        c_neg_[static_cast<std::size_t>(nr_ - 1)] - neg_.alpha * j_neg,
        params_.neg.c_max);
    const double cs_pos = clamp_stoich(
        c_pos_[static_cast<std::size_t>(nr_ - 1)] - pos_.alpha * j_pos,
        params_.pos.c_max);
    const double u_neg = negative_ocp(cs_neg / params_.neg.c_max).value;
    const double u_pos = positive_ocp(cs_pos / params_.pos.c_max).value;
    const double eta_neg =
        surface_overpotential(params_.neg, j_neg, cs_neg, params_.elyte.c0,
                              params_.temperature);
    const double eta_pos =
        surface_overpotential(params_.pos, j_pos, cs_pos, params_.elyte.c0,
                              params_.temperature);
    return u_pos - u_neg + eta_pos - eta_neg;
  }

  double surface(bool positive, double current) const {
    const ElectrodeGeometry& geometry = positive ? pos_ : neg_;
    const ElectrodeParams& electrode = positive ? params_.pos : params_.neg;
    const std::vector<double>& c = positive ? c_pos_ : c_neg_;
    const double sign = positive ? -1.0 : 1.0;
    const double j = sign * current * geometry.j_coeff;
    return clamp_stoich(
        c[static_cast<std::size_t>(nr_ - 1)] - geometry.alpha * j,
        electrode.c_max);
  }

  const Params& params() const { return params_; }
  int radial_cells() const { return nr_; }
  double time_step() const { return dt_; }

 private:
  Params params_;
  int nr_ = 0;
  double dt_ = 0.0;
  ElectrodeGeometry neg_;
  ElectrodeGeometry pos_;
  std::vector<double> c_neg_;
  std::vector<double> c_pos_;
};

inline std::vector<BatchRow> rows_from_protocol(
    const Params& params, const std::vector<BatchSegment>& segments, int nr,
    double dt) {
  if (segments.empty()) {
    throw std::invalid_argument("protocol needs at least one segment");
  }
  SpmBatch model(params, nr, dt);
  std::vector<BatchRow> rows;
  rows.push_back(BatchRow{0.0, segments.front().current,
                          model.voltage(segments.front().current)});
  double t = 0.0;
  for (std::size_t s = 0; s < segments.size(); ++s) {
    const BatchSegment& segment = segments[s];
    const int steps = static_cast<int>(std::lround(segment.duration / dt));
    for (int k = 0; k < steps; ++k) {
      model.step(segment.current);
      t += dt;
      const bool boundary = k + 1 == steps && s + 1 < segments.size();
      const double next =
          boundary ? segments[s + 1].current : segment.current;
      rows.push_back(BatchRow{t, next, model.voltage(next)});
    }
  }
  return rows;
}

inline std::vector<double> simulate_rows(const Params& params,
                                         const std::vector<BatchRow>& rows,
                                         int nr, double dt) {
  SpmBatch model(params, nr, dt);
  std::vector<double> voltages(rows.size(), 0.0);
  voltages[0] = model.voltage(rows[0].current);
  for (std::size_t i = 1; i < rows.size(); ++i) {
    const double expected = rows[0].t + dt * static_cast<double>(i);
    if (std::abs(rows[i].t - expected) > 1e-6) {
      throw std::runtime_error("samples are not on a uniform time grid");
    }
    // The interval [t_{i-1}, t_i] uses the current stored at t_{i-1}; the
    // voltage at t_i already reflects the current of the next interval.
    model.step(rows[i - 1].current);
    voltages[i] = model.voltage(rows[i].current);
  }
  return voltages;
}

}  // namespace fit
}  // namespace dfn
