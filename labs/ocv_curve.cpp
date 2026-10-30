#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string_view>

#include "dfn/thermo.hpp"

namespace {

double open_circuit_voltage(double soc) {
  const double x = 0.10 + 0.75 * soc;
  const double y = 0.90 - 0.80 * soc;
  return dfn::positive_ocp(y).value - dfn::negative_ocp(x).value;
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
