#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

std::vector<double> euler(double lambda, double h, int steps) {
  std::vector<double> y(steps + 1, 1.0);
  for (int i = 0; i < steps; ++i) {
    y[i + 1] = y[i] + h * lambda * y[i];
  }
  return y;
}

double max_abs(const std::vector<double>& y) {
  double m = 0.0;
  for (double v : y) {
    m = std::max(m, std::abs(v));
  }
  return m;
}

}  // namespace

int main() {
  std::cout << std::scientific << std::setprecision(6);

  std::cout << "decay: y'=-1*y, y(0)=1, h=0.1\n";
  std::cout << "t,euler,exact\n";
  const auto decay = euler(-1.0, 0.1, 5);
  for (int i = 0; i < static_cast<int>(decay.size()); ++i) {
    const double t = 0.1 * i;
    std::cout << t << ',' << decay[i] << ',' << std::exp(-t) << '\n';
  }

  std::cout << "stability: y'=-100*y, t_end=0.2, exact="
            << std::exp(-100.0 * 0.2) << '\n';
  std::cout << "h,steps,max_abs,end\n";
  for (double h : {0.05, 0.02, 0.01, 0.005}) {
    const int steps = static_cast<int>(std::lround(0.2 / h));
    const auto y = euler(-100.0, h, steps);
    std::cout << h << ',' << steps << ',' << max_abs(y) << ',' << y.back()
              << '\n';
  }

  std::cout << "two modes: lambda=(-1,-100), y(0)=(1,1), t_end=0.2\n";
  std::cout << "h,end_slow,end_fast\n";
  for (double h : {0.05, 0.02, 0.01}) {
    const int steps = static_cast<int>(std::lround(0.2 / h));
    const auto slow = euler(-1.0, h, steps);
    const auto fast = euler(-100.0, h, steps);
    std::cout << h << ',' << slow.back() << ',' << fast.back() << '\n';
  }

  return 0;
}
