#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
  const double x = 1.0;
  const double exact = std::cos(x);

  std::cout << "h,forward,central,forward_err,central_err\n"
            << std::scientific << std::setprecision(6);

  for (int k = 0; k < 48; ++k) {
    const double h = 0.1 * std::pow(0.5, k);
    const double forward = (std::sin(x + h) - std::sin(x)) / h;
    const double central = (std::sin(x + h) - std::sin(x - h)) / (2.0 * h);
    std::cout << h << ',' << forward << ',' << central << ','
              << std::abs(forward - exact) << ','
              << std::abs(central - exact) << '\n';
  }

  return 0;
}
