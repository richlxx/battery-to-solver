#pragma once

namespace dfn {

// Open-circuit potential and its derivative with respect to stoichiometry.
struct OcpValue {
  double value = 0.0;  // V
  double slope = 0.0;  // V per unit stoichiometry
};

// Illustrative half-cell fits shared by the labs and the DFN model. The
// Chen2020 parameter set will replace them once the configuration loader
// lands; the shape (slope signs, range) is the part that matters here.
OcpValue negative_ocp(double stoich);
OcpValue positive_ocp(double stoich);

}  // namespace dfn
