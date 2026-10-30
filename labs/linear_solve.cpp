#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/linalg/dense.hpp"
#include "dfn/linalg/tridiagonal.hpp"

namespace {

dfn::TridiagonalMatrix make_matrix(int n) {
  dfn::TridiagonalMatrix matrix;
  matrix.n = n;
  matrix.lower.assign(static_cast<std::size_t>(n), -1.0);
  matrix.diag.assign(static_cast<std::size_t>(n), 4.0);
  matrix.upper.assign(static_cast<std::size_t>(n), -1.0);
  matrix.lower[0] = 0.0;
  matrix.upper[static_cast<std::size_t>(n - 1)] = 0.0;
  return matrix;
}

std::vector<double> make_exact(int n) {
  std::vector<double> exact(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    exact[static_cast<std::size_t>(i)] = 1.0 + 0.01 * static_cast<double>(i);
  }
  return exact;
}

std::vector<double> make_rhs(const dfn::TridiagonalMatrix& matrix,
                             const std::vector<double>& exact) {
  const int n = matrix.n;
  std::vector<double> rhs(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    rhs[static_cast<std::size_t>(i)] =
        matrix.diag[static_cast<std::size_t>(i)] * exact[static_cast<std::size_t>(i)];
    if (i > 0) {
      rhs[static_cast<std::size_t>(i)] +=
          matrix.lower[static_cast<std::size_t>(i)] *
          exact[static_cast<std::size_t>(i - 1)];
    }
    if (i + 1 < n) {
      rhs[static_cast<std::size_t>(i)] +=
          matrix.upper[static_cast<std::size_t>(i)] *
          exact[static_cast<std::size_t>(i + 1)];
    }
  }
  return rhs;
}

std::vector<std::vector<double>> to_dense(const dfn::TridiagonalMatrix& matrix) {
  const int n = matrix.n;
  std::vector<std::vector<double>> dense(
      static_cast<std::size_t>(n),
      std::vector<double>(static_cast<std::size_t>(n), 0.0));
  for (int i = 0; i < n; ++i) {
    dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] =
        matrix.diag[static_cast<std::size_t>(i)];
    if (i > 0) {
      dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(i - 1)] =
          matrix.lower[static_cast<std::size_t>(i)];
    }
    if (i + 1 < n) {
      dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(i + 1)] =
          matrix.upper[static_cast<std::size_t>(i)];
    }
  }
  return dense;
}

double max_error(const std::vector<double>& x, const std::vector<double>& exact) {
  double error = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    error = std::max(error, std::abs(x[i] - exact[i]));
  }
  return error;
}

}  // namespace

int main() {
  std::cout << "n,dense_ms,thomas_ms,dense_max_error,thomas_max_error\n"
            << std::scientific << std::setprecision(3);

  for (int n : {50, 100, 200, 400}) {
    const dfn::TridiagonalMatrix matrix = make_matrix(n);
    const std::vector<double> exact = make_exact(n);
    const std::vector<double> rhs = make_rhs(matrix, exact);
    const std::vector<std::vector<double>> dense = to_dense(matrix);

    const auto dense_start = std::chrono::steady_clock::now();
    const std::vector<double> x_dense = dfn::solve_dense(dense, rhs);
    const auto dense_end = std::chrono::steady_clock::now();

    const auto thomas_start = std::chrono::steady_clock::now();
    const std::vector<double> x_thomas = dfn::solve_tridiagonal(matrix, rhs);
    const auto thomas_end = std::chrono::steady_clock::now();

    const double dense_ms =
        std::chrono::duration<double, std::milli>(dense_end - dense_start).count();
    const double thomas_ms =
        std::chrono::duration<double, std::milli>(thomas_end - thomas_start).count();

    std::cout << n << ',' << dense_ms << ',' << thomas_ms << ','
              << max_error(x_dense, exact) << ','
              << max_error(x_thomas, exact) << '\n';
  }

  return 0;
}
