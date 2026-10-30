#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

int main() {
  const double lambda = -100.0;
  const double t_end = 0.2;
  const double exact = std::exp(lambda * t_end);

  std::cout << "h,steps,explicit_final,implicit_final,exact_final,explicit_max\n"
            << std::scientific << std::setprecision(6);

  for (double h : {0.1, 0.05, 0.02, 0.01, 0.005}) {
    const int steps = static_cast<int>(
        std::lround(t_end / h));

    double explicit_value = 1.0;
    double explicit_max = 1.0;
    double implicit_value = 1.0;
    for (int i = 0; i < steps; ++i) {
      explicit_value *= 1.0 + lambda * h;
      explicit_max = std::max(explicit_max, std::abs(explicit_value));
      implicit_value /= 1.0 - lambda * h;
    }

    std::cout << h << ',' << steps << ',' << explicit_value << ','
              << implicit_value << ',' << exact << ',' << explicit_max << '\n';
  }

  return 0;
}
