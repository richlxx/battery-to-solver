#pragma once

#include "dfn/params.hpp"

namespace dfn {

// Butler-Volmer reaction at one interface, with the derivatives needed by
// the DFN Jacobian. j is positive when lithium leaves the solid.
struct ReactionValue {
  double j = 0.0;       // mol/(m^2 s)
  double dj_deta = 0.0;
  double dj_dce = 0.0;
  double dj_dcs = 0.0;
};

ReactionValue reaction_at(const Params& params, const ElectrodeParams& electrode,
                          bool positive, double ce, double cs_surf, double phis,
                          double phie);

// Inverts Butler-Volmer for a prescribed flux j. Used by the reduced models,
// where j is fixed by the applied current instead of the local potential.
double overpotential_from_flux(const Params& params,
                               const ElectrodeParams& electrode, bool positive,
                               double ce, double cs_surf, double j);

}  // namespace dfn
