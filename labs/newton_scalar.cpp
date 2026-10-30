#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  double x = 1.0;
  if (argc > 1) {
    x = std::atof(argv[1]);
  }

  std::cout << "iter,x,f(x)\n" << std::setprecision(10);
  for (int iteration = 0; iteration < 10; ++iteration) {
    const double f = x * x - 2.0;
    std::cout << iteration << ',' << x << ',' << f << '\n';
    if (std::abs(f) < 1e-14) {
      break;
    }
    x -= f / (2.0 * x);
  }
  return 0;
}
