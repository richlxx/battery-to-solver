#include "dfn/models/reaction.hpp"

#include <algorithm>
#include <cmath>

#include "dfn/thermo.hpp"

namespace dfn {

namespace {

double exchange_current(const Params& params, const ElectrodeParams& electrode,
                        double ce, double cs_surf) {
  const double cmax = electrode.c_max;
  const double cs = std::clamp(cs_surf, 1e-6 * cmax, cmax * (1.0 - 1e-6));
  const double ce_clamped = std::max(ce, 1e-6 * params.elyte.c0);
  return kFaraday * electrode.k * std::pow(ce_clamped, electrode.alpha_a) *
         std::pow(cmax - cs, electrode.alpha_a) *
         std::pow(cs, electrode.alpha_c);
}

}  // namespace

ReactionValue reaction_at(const Params& params, const ElectrodeParams& electrode,
                          bool positive, double ce, double cs_surf, double phis,
                          double phie) {
  const double cmax = electrode.c_max;
  const double cs = std::clamp(cs_surf, 1e-6 * cmax, cmax * (1.0 - 1e-6));
  const double ce_clamped = std::max(ce, 1e-6 * params.elyte.c0);
  const double stoich = cs / cmax;
  const OcpValue ocp =
      positive ? positive_ocp(stoich) : negative_ocp(stoich);

  const double i0 = exchange_current(params, electrode, ce_clamped, cs);
  const double eta = phis - phie - ocp.value;
  const double a =
      electrode.alpha_a * kFaraday / (kGasConstant * params.temperature);
  const double c =
      electrode.alpha_c * kFaraday / (kGasConstant * params.temperature);
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

double overpotential_from_flux(const Params& params,
                               const ElectrodeParams& electrode,
                               bool /*positive*/, double ce, double cs_surf,
                               double j) {
  if (j == 0.0) {
    return 0.0;
  }
  const double i0 = exchange_current(params, electrode, ce, cs_surf);
  if (i0 <= 0.0) {
    return 0.0;
  }
  const double a =
      electrode.alpha_a * kFaraday / (kGasConstant * params.temperature);
  const double c =
      electrode.alpha_c * kFaraday / (kGasConstant * params.temperature);

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

}  // namespace dfn
