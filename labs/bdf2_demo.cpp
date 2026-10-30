#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/solvers/bdf2.hpp"

int main() {
  const auto rhs = [](double, const std::vector<double>& y) {
    return std::vector<double>{-y[0]};
  };
  const auto jacobian = [](double, const std::vector<double>&) {
    return std::vector<std::vector<double>>{{-1.0}};
  };

  const double exact = std::exp(-1.0);
  std::cout << "h,steps,error,ratio\n" << std::scientific << std::setprecision(4);

  double previous_error = 0.0;
  for (double h : {0.2, 0.1, 0.05, 0.025}) {
    const dfn::Bdf2Result result =
        dfn::bdf2_integrate(rhs, jacobian, {1.0}, 0.0, 1.0, h);
    const double error = std::abs(result.y.back()[0] - exact);
    const double ratio = previous_error > 0.0 ? previous_error / error : 0.0;
    std::cout << h << ',' << result.steps << ',' << error << ',' << ratio << '\n';
    previous_error = error;
  }

  return 0;
}
