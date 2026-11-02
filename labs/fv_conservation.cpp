#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/grid.hpp"
#include "dfn/operators.hpp"

namespace {

double total(const dfn::Grid1D& grid, const std::vector<double>& c) {
  double sum = 0.0;
  for (int i = 0; i < grid.n; ++i) {
    sum += grid.volume[static_cast<std::size_t>(i)] *
           c[static_cast<std::size_t>(i)];
  }
  return sum;
}

std::vector<double> advance(const dfn::DiffusionOperator& op,
                            std::vector<double> c, double dt) {
  const std::vector<double> rhs = op.matrix.multiply(c);
  for (std::size_t i = 0; i < c.size(); ++i) {
    c[i] += dt * (rhs[i] + op.source[i]);
  }
  return c;
}

void run(const char* name, const dfn::Grid1D& grid,
         const dfn::DiffusionOperator& op, double surface_flux, double dt,
         int steps) {
  std::vector<double> c = {1200.0, 1000.0, 800.0};
  const double initial = total(grid, c);
  const double area = grid.area[static_cast<std::size_t>(grid.n)];
  std::cout << "case,step,total,drift,expected\n";
  for (int step = 0; step <= steps; ++step) {
    const double now = total(grid, c);
    const double expected = -area * surface_flux * dt * step;
    std::cout << name << ',' << step << ',' << now << ',' << (now - initial)
              << ',' << expected << '\n';
    if (step < steps) {
      c = advance(op, c, dt);
    }
  }
}

}  // namespace

int main() {
  constexpr double radius = 5.0e-6;
  constexpr double diffusivity = 3.0e-14;
  constexpr double dt = 10.0;
  constexpr int steps = 5;
  constexpr double surface_flux = 1.0e-5;
  const int n = 3;

  const dfn::Grid1D grid = dfn::make_uniform_sphere_grid(n, radius);
  const std::vector<double> D(static_cast<std::size_t>(n), diffusivity);
  const dfn::DiffusionOperator sealed = dfn::assemble_diffusion(grid, D);
  const dfn::DiffusionOperator open =
      dfn::assemble_diffusion(grid, D, 0.0, surface_flux);

  std::cout << std::scientific << std::setprecision(6);
  std::cout << "n=" << n << " R=" << radius << " D=" << diffusivity
            << " dt=" << dt << " j=" << surface_flux << '\n';
  run("sealed", grid, sealed, 0.0, dt, steps);
  run("open", grid, open, surface_flux, dt, steps);

  std::cout << "\ncell,volume,dc_dt,ledger\n";
  const std::vector<double> c = {1200.0, 1000.0, 800.0};
  const std::vector<double> rhs = sealed.matrix.multiply(c);
  double ledger_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const double ledger = grid.volume[static_cast<std::size_t>(i)] *
                          rhs[static_cast<std::size_t>(i)];
    ledger_sum += ledger;
    std::cout << i << ',' << grid.volume[static_cast<std::size_t>(i)] << ','
              << rhs[static_cast<std::size_t>(i)] << ',' << ledger << '\n';
  }
  std::cout << "ledger_sum," << ledger_sum << '\n';
  return 0;
}
