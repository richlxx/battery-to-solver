#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/solvers/bdf2.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

}  // namespace

int main() {
  const auto rhs = [](double, const std::vector<double>& y) {
    return std::vector<double>{-y[0]};
  };
  const auto jacobian = [](double, const std::vector<double>&) {
    return std::vector<std::vector<double>>{{-1.0}};
  };

  const auto run = [&](double h) {
    return dfn::bdf2_integrate(rhs, jacobian, {1.0}, 0.0, 1.0, h);
  };

  const double exact = std::exp(-1.0);
  const dfn::Bdf2Result coarse = run(0.05);
  const dfn::Bdf2Result medium = run(0.025);
  const dfn::Bdf2Result fine = run(0.0125);

  if (!coarse.ok || !medium.ok || !fine.ok) {
    return fail("BDF2 integration failed");
  }
  if (coarse.steps != 20 || medium.steps != 40 || fine.steps != 80) {
    return fail("unexpected step counts");
  }
  if (std::abs(coarse.t.back() - 1.0) > 1e-12) {
    return fail("final time");
  }

  const double e1 = std::abs(coarse.y.back()[0] - exact);
  const double e2 = std::abs(medium.y.back()[0] - exact);
  const double e3 = std::abs(fine.y.back()[0] - exact);
  const double ratio1 = e1 / e2;
  const double ratio2 = e2 / e3;
  if (ratio1 < 3.0 || ratio1 > 5.5 || ratio2 < 3.0 || ratio2 > 5.5) {
    return fail("BDF2 observed order");
  }

  const auto stiff_rhs = [](double, const std::vector<double>& y) {
    return std::vector<double>{-100.0 * y[0]};
  };
  const auto stiff_jacobian = [](double, const std::vector<double>&) {
    return std::vector<std::vector<double>>{{-100.0}};
  };
  const dfn::Bdf2Result stiff = dfn::bdf2_integrate(
      stiff_rhs, stiff_jacobian, {1.0}, 0.0, 0.2, 0.05);
  if (!stiff.ok || !std::isfinite(stiff.y.back()[0]) ||
      std::abs(stiff.y.back()[0]) > 1.0) {
    return fail("BDF2 stiff stability");
  }

  std::cout << "bdf2 ok (order ratios " << ratio1 << ", " << ratio2 << ")\n";
  return 0;
}
