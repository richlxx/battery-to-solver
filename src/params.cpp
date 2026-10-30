#include "dfn/params.hpp"

#include <stdexcept>

namespace dfn {

Params make_demo_params() {
  Params p;
  p.area = 0.0165;
  p.temperature = 298.15;

  p.neg.thickness = 85.0e-6;
  p.neg.radius = 5.0e-6;
  p.neg.porosity = 0.25;
  p.neg.bruggeman = 1.5;
  p.neg.c_max = 30000.0;
  p.neg.c0 = 0.8 * p.neg.c_max;
  p.neg.d_s = 3.0e-14;
  p.neg.k = 1.0e-6;
  p.neg.sigma = 100.0;

  p.pos.thickness = 75.0e-6;
  p.pos.radius = 5.0e-6;
  p.pos.porosity = 0.33;
  p.pos.bruggeman = 1.5;
  p.pos.c_max = 60000.0;
  p.pos.c0 = 0.3 * p.pos.c_max;
  p.pos.d_s = 4.0e-15;
  p.pos.k = 1.0e-6;
  p.pos.sigma = 10.0;

  p.sep.thickness = 12.0e-6;
  p.sep.porosity = 0.47;
  p.sep.bruggeman = 1.5;

  p.elyte.c0 = 1000.0;
  p.elyte.d_e = 5.0e-10;
  p.elyte.t_plus = 0.26;
  p.elyte.kappa = 1.0;

  return p;
}

void validate(const Params& params) {
  const auto require = [](bool condition, const std::string& message) {
    if (!condition) {
      throw std::invalid_argument(message);
    }
  };

  const auto check_electrode = [&require](const ElectrodeParams& e,
                                          const std::string& name) {
    require(e.thickness > 0.0, name + ": thickness must be positive");
    require(e.radius > 0.0, name + ": particle radius must be positive");
    require(e.porosity > 0.0 && e.porosity < 1.0,
            name + ": porosity must be in (0, 1)");
    require(e.c_max > 0.0, name + ": c_max must be positive");
    require(e.c0 >= 0.0 && e.c0 <= e.c_max,
            name + ": c0 must be in [0, c_max]");
    require(e.d_s > 0.0, name + ": d_s must be positive");
    require(e.k > 0.0, name + ": reaction rate must be positive");
    require(e.alpha_a > 0.0 && e.alpha_c > 0.0,
            name + ": transfer coefficients must be positive");
    require(e.sigma > 0.0, name + ": sigma must be positive");
  };

  require(params.area > 0.0, "cell area must be positive");
  require(params.temperature > 0.0, "temperature must be positive");
  check_electrode(params.neg, "negative electrode");
  check_electrode(params.pos, "positive electrode");

  require(params.sep.thickness > 0.0, "separator: thickness must be positive");
  require(params.sep.porosity > 0.0 && params.sep.porosity < 1.0,
          "separator: porosity must be in (0, 1)");

  require(params.elyte.c0 > 0.0, "electrolyte: c0 must be positive");
  require(params.elyte.d_e > 0.0, "electrolyte: d_e must be positive");
  require(params.elyte.t_plus >= 0.0 && params.elyte.t_plus < 1.0,
          "electrolyte: t_plus must be in [0, 1)");
  require(params.elyte.kappa > 0.0, "electrolyte: kappa must be positive");
}

}  // namespace dfn
