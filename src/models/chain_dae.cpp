#include "dfn/models/chain_dae.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace dfn {

namespace {

int cell_index(int cell) { return 2 * cell; }
int aux_index(int cell) { return 2 * cell + 1; }

}  // namespace

ChainDae::ChainDae(ChainDaeOptions options) : options_(options) {
  if (options_.cells < 1) {
    throw std::invalid_argument("ChainDae: cells must be >= 1");
  }
  if (options_.lambda <= 0.0) {
    throw std::invalid_argument("ChainDae: lambda must be positive");
  }
  if (options_.diffusion < 0.0) {
    throw std::invalid_argument("ChainDae: diffusion must be non-negative");
  }

  std::vector<std::pair<int, int>> entries;
  entries.reserve(static_cast<std::size_t>(8 * options_.cells));
  for (int i = 0; i < options_.cells; ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    entries.emplace_back(c, c);
    entries.emplace_back(c, w);
    if (i > 0) {
      entries.emplace_back(c, cell_index(i - 1));
    }
    if (i + 1 < options_.cells) {
      entries.emplace_back(c, cell_index(i + 1));
    }
    entries.emplace_back(w, w);
  }
  pattern_ = build_pattern(size(), std::move(entries));
}

void ChainDae::residual(double t, const double* y, const double* yp,
                        double* out) const {
  const double d = options_.diffusion;
  for (int i = 0; i < options_.cells; ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    const double left = (i > 0) ? y[cell_index(i - 1)] : y[c];
    const double right = (i + 1 < options_.cells) ? y[cell_index(i + 1)] : y[c];
    const double laplacian = left + right - 2.0 * y[c];

    out[c] = yp[c] + options_.lambda * (y[c] - y[w]) - d * laplacian;
    out[w] = y[w] - std::sin(t);
  }
}

std::vector<int> ChainDae::variable_types() const {
  std::vector<int> types(static_cast<std::size_t>(size()), 0);
  for (int i = 0; i < options_.cells; ++i) {
    types[static_cast<std::size_t>(cell_index(i))] = 1;
  }
  return types;
}

void ChainDae::jacobian(double t, const double* y, const double* yp, double cj,
                        CscBuilder& out) const {
  (void)t;
  (void)y;
  (void)yp;
  const double d = options_.diffusion;
  for (int i = 0; i < options_.cells; ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    const bool interior = (i > 0) && (i + 1 < options_.cells);
    const double laplacian_diagonal =
        interior ? 2.0 * d : ((options_.cells > 1) ? d : 0.0);
    out.add(c, c, cj + options_.lambda + laplacian_diagonal);
    if (i > 0) {
      out.add(c, cell_index(i - 1), -d);
    }
    if (i + 1 < options_.cells) {
      out.add(c, cell_index(i + 1), -d);
    }
    out.add(c, w, -options_.lambda);
    out.add(w, w, 1.0);
  }
}

std::vector<double> chain_dae_consistent_state(const ChainDae& model,
                                               double t0) {
  return chain_dae_mode_state(model, t0, 0.0);
}

std::vector<double> chain_dae_mode_state(const ChainDae& model, double t0,
                                         double amplitude) {
  std::vector<double> y(static_cast<std::size_t>(model.size()), 0.0);
  const double concentration =
      chain_dae_exact_concentration(model.lambda(), t0);
  const double spacing =
      3.14159265358979323846 / static_cast<double>(model.cells() + 1);
  for (int i = 0; i < model.cells(); ++i) {
    const double shape = std::sin(static_cast<double>(i + 1) * spacing);
    y[static_cast<std::size_t>(cell_index(i))] =
        concentration + amplitude * shape;
    y[static_cast<std::size_t>(aux_index(i))] = std::sin(t0);
  }
  return y;
}

std::vector<double> chain_dae_derivative_from_state(
    const ChainDae& model, double t, const std::vector<double>& y) {
  std::vector<double> yp(static_cast<std::size_t>(model.size()), 0.0);
  const double d = model.diffusion();
  for (int i = 0; i < model.cells(); ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    const double left = (i > 0) ? y[static_cast<std::size_t>(cell_index(i - 1))]
                                : y[static_cast<std::size_t>(c)];
    const double right =
        (i + 1 < model.cells())
            ? y[static_cast<std::size_t>(cell_index(i + 1))]
            : y[static_cast<std::size_t>(c)];
    const double laplacian =
        left + right - 2.0 * y[static_cast<std::size_t>(c)];
    yp[static_cast<std::size_t>(c)] =
        -model.lambda() * (y[static_cast<std::size_t>(c)] -
                           y[static_cast<std::size_t>(w)]) +
        d * laplacian;
    yp[static_cast<std::size_t>(w)] = std::cos(t);
  }
  return yp;
}

std::vector<double> chain_dae_consistent_derivative(const ChainDae& model,
                                                    double t0) {
  return chain_dae_derivative_from_state(
      model, t0, chain_dae_mode_state(model, t0, 0.0));
}

double chain_dae_exact_concentration(double lambda, double t) {
  if (lambda <= 0.0) {
    throw std::invalid_argument("chain_dae_exact_concentration: lambda > 0");
  }
  const double denominator = lambda * lambda + 1.0;
  return lambda *
         (lambda * std::sin(t) - std::cos(t) + std::exp(-lambda * t)) /
         denominator;
}

double chain_dae_exact_concentration_rate(double lambda, double t) {
  if (lambda <= 0.0) {
    throw std::invalid_argument(
        "chain_dae_exact_concentration_rate: lambda > 0");
  }
  const double denominator = lambda * lambda + 1.0;
  return lambda *
         (lambda * std::cos(t) + std::sin(t) - lambda * std::exp(-lambda * t)) /
         denominator;
}

}  // namespace dfn
