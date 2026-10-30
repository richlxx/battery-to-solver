#include "dfn/thermo.hpp"

#include <algorithm>
#include <cmath>

namespace dfn {

OcpValue negative_ocp(double stoich) {
  const double x = std::clamp(stoich, 0.0, 1.0);
  const double one_minus_x = 1.0 - x;
  OcpValue result;
  result.value =
      0.05 + 0.8 * std::exp(-10.0 * x) + 0.1 * one_minus_x * one_minus_x;
  result.slope = -8.0 * std::exp(-10.0 * x) - 0.2 * one_minus_x;
  return result;
}

OcpValue positive_ocp(double stoich) {
  const double y = std::clamp(stoich, 0.0, 1.0);
  OcpValue result;
  result.value = 3.9 + 0.2 * (1.0 - y);
  result.slope = -0.2;
  return result;
}

}  // namespace dfn
