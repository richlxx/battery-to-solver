#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "dfn/grid.hpp"
#include "dfn/operators.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

bool close_relative(double a, double b, double tolerance) {
  const double scale = std::max(std::abs(a), std::abs(b));
  if (scale == 0.0) {
    return true;
  }
  return std::abs(a - b) <= tolerance * scale;
}

}  // namespace

int main() {
  constexpr double pi = 3.14159265358979323846;
  const double radius = 5.0e-6;
  const int n = 8;

  const dfn::Grid1D sphere = dfn::make_uniform_sphere_grid(n, radius);
  double volume = 0.0;
  for (double v : sphere.volume) {
    volume += v;
  }
  const double exact_volume = 4.0 / 3.0 * pi * radius * radius * radius;
  if (!close_relative(volume, exact_volume, 1e-12)) {
    return fail("sphere volume");
  }
  if (std::abs(sphere.area[0]) > 0.0) {
    return fail("center face area must be zero");
  }
  if (!close_relative(sphere.area[static_cast<std::size_t>(n)],
                      4.0 * pi * radius * radius, 1e-12)) {
    return fail("sphere surface area");
  }

  const dfn::Grid1D line = dfn::make_uniform_line_grid(10, 100.0e-6, 1.0e-4);
  double line_volume = 0.0;
  for (double v : line.volume) {
    line_volume += v;
  }
  if (!close_relative(line_volume, 1.0e-8, 1e-12)) {
    return fail("line volume");
  }

  const std::vector<double> D_sphere(static_cast<std::size_t>(sphere.n), 3.0e-14);
  const dfn::DiffusionOperator op = dfn::assemble_diffusion(sphere, D_sphere);

  const std::vector<double> uniform(static_cast<std::size_t>(sphere.n), 1.0);
  const std::vector<double> uniform_rhs = op.matrix.multiply(uniform);
  for (double value : uniform_rhs) {
    if (std::abs(value) > 1e-12) {
      return fail("uniform state must be steady with zero flux");
    }
  }

  std::mt19937 generator(42);
  std::uniform_real_distribution<double> distribution(0.5, 1.5);
  std::vector<double> c(static_cast<std::size_t>(sphere.n));
  for (double& value : c) {
    value = distribution(generator);
  }

  const std::vector<double> rhs = op.matrix.multiply(c);
  double balance = 0.0;
  double scale = 0.0;
  for (int i = 0; i < sphere.n; ++i) {
    balance += sphere.volume[static_cast<std::size_t>(i)] *
               (rhs[static_cast<std::size_t>(i)] +
                op.source[static_cast<std::size_t>(i)]);
    scale += std::abs(sphere.volume[static_cast<std::size_t>(i)] *
                      rhs[static_cast<std::size_t>(i)]);
  }
  if (std::abs(balance) > 1e-10 * scale) {
    return fail("sphere conservation");
  }

  const double surface_flux = 1.0e-6;
  const dfn::DiffusionOperator with_flux =
      dfn::assemble_diffusion(sphere, D_sphere, 0.0, surface_flux);
  const std::vector<double> flux_rhs = with_flux.matrix.multiply(c);
  double flux_balance = 0.0;
  for (int i = 0; i < sphere.n; ++i) {
    flux_balance += sphere.volume[static_cast<std::size_t>(i)] *
                    (flux_rhs[static_cast<std::size_t>(i)] +
                     with_flux.source[static_cast<std::size_t>(i)]);
  }
  const double expected = -sphere.area[static_cast<std::size_t>(n)] * surface_flux;
  if (!close_relative(flux_balance, expected, 1e-10)) {
    return fail("surface flux balance");
  }

  const std::vector<double> D_line(static_cast<std::size_t>(line.n), 1.0e-9);
  const dfn::DiffusionOperator line_op = dfn::assemble_diffusion(line, D_line);
  std::vector<double> line_c(static_cast<std::size_t>(line.n));
  for (double& value : line_c) {
    value = distribution(generator);
  }
  const std::vector<double> line_rhs = line_op.matrix.multiply(line_c);
  double line_balance = 0.0;
  double line_scale = 0.0;
  for (int i = 0; i < line.n; ++i) {
    line_balance += line.volume[static_cast<std::size_t>(i)] *
                    line_rhs[static_cast<std::size_t>(i)];
    line_scale += std::abs(line.volume[static_cast<std::size_t>(i)] *
                           line_rhs[static_cast<std::size_t>(i)]);
  }
  if (std::abs(line_balance) > 1e-10 * line_scale) {
    return fail("line conservation");
  }

  std::cout << "discretization ok (sphere n=" << sphere.n
            << ", line n=" << line.n << ")\n";
  return 0;
}
