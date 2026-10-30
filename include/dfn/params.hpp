#pragma once

#include <string>
#include <vector>

namespace dfn {

inline constexpr double kFaraday = 96485.33212;  // C/mol
inline constexpr double kGasConstant = 8.314462618;  // J/(mol K)

struct ElectrodeParams {
  double thickness = 0.0;    // m
  double radius = 0.0;       // m
  double porosity = 0.0;     // -
  double bruggeman = 1.5;    // -
  double c_max = 0.0;        // mol/m^3
  double c0 = 0.0;           // mol/m^3
  double d_s = 0.0;          // m^2/s
  double k = 0.0;            // reaction rate constant
  double alpha_a = 0.5;      // anodic transfer coefficient
  double alpha_c = 0.5;      // cathodic transfer coefficient
  double sigma = 0.0;        // S/m
};

struct SeparatorParams {
  double thickness = 0.0;    // m
  double porosity = 0.0;     // -
  double bruggeman = 1.5;    // -
};

struct ElectrolyteParams {
  double c0 = 0.0;           // mol/m^3
  double d_e = 0.0;          // m^2/s
  double t_plus = 0.0;       // -
  double kappa = 0.0;        // S/m
};

struct Params {
  double area = 0.0;         // m^2
  double temperature = 0.0;  // K
  ElectrodeParams neg;
  ElectrodeParams pos;
  SeparatorParams sep;
  ElectrolyteParams elyte;

  double cell_length() const {
    return neg.thickness + sep.thickness + pos.thickness;
  }
};

// Illustrative values for wiring tests; replace with the Chen2020 table.
Params make_demo_params();

void validate(const Params& params);

}  // namespace dfn
