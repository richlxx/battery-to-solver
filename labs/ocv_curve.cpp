#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string_view>

namespace {

// Illustrative half-cell OCP fits; replaced by the Chen2020 data in chapter 3.
double negative_ocp(double x) {
  const double one_minus_x = 1.0 - x;
  return 0.05 + 0.8 * std::exp(-10.0 * x) + 0.1 * one_minus_x * one_minus_x;
}

double positive_ocp(double y) {
  return 3.9 + 0.2 * (1.0 - y);
}

double open_circuit_voltage(double soc) {
  const double x = 0.10 + 0.75 * soc;
  const double y = 0.90 - 0.80 * soc;
  return positive_ocp(y) - negative_ocp(x);
}

void print_usage() {
  std::cout << "Usage: ocv_curve [--n POINTS]\n"
               "Print soc, electrode stoichiometries and OCV as CSV.\n";
}

}  // namespace

int main(int argc, char** argv) {
  int n = 101;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      print_usage();
      return 0;
    }
    if (arg == "--n") {
      if (i + 1 >= argc) {
        std::cerr << "error: --n requires a value\n";
        return 2;
      }
      n = std::atoi(argv[++i]);
      if (n < 2) {
        std::cerr << "error: --n must be at least 2\n";
        return 2;
      }
      continue;
    }
    std::cerr << "error: unknown option " << arg << "\n";
    return 2;
  }

  std::cout << "soc,sto_neg,sto_pos,ocv_v\n"
            << std::fixed << std::setprecision(6);

  for (int k = 0; k < n; ++k) {
    const double soc = static_cast<double>(k) / static_cast<double>(n - 1);
    const double x = 0.10 + 0.75 * soc;
    const double y = 0.90 - 0.80 * soc;
    std::cout << soc << ',' << x << ',' << y << ','
              << open_circuit_voltage(soc) << '\n';
  }

  return 0;
}
