#include "dfn/linalg/tridiagonal.hpp"

#include <cmath>
#include <stdexcept>

namespace dfn {

namespace {
constexpr double kPivotTolerance = 1e-14;
}

std::vector<double> solve_tridiagonal(const TridiagonalMatrix& A,
                                      const std::vector<double>& b) {
  const int n = A.n;
  if (n <= 0 ||
      static_cast<int>(A.lower.size()) != n ||
      static_cast<int>(A.diag.size()) != n ||
      static_cast<int>(A.upper.size()) != n ||
      static_cast<int>(b.size()) != n) {
    throw std::invalid_argument("solve_tridiagonal: size mismatch");
  }
  if (std::abs(A.diag[0]) < kPivotTolerance) {
    throw std::runtime_error("solve_tridiagonal: zero pivot");
  }

  std::vector<double> c(static_cast<std::size_t>(n), 0.0);
  std::vector<double> d(static_cast<std::size_t>(n), 0.0);
  std::vector<double> x(static_cast<std::size_t>(n), 0.0);

  c[0] = (n > 1) ? A.upper[0] / A.diag[0] : 0.0;
  d[0] = b[0] / A.diag[0];

  for (int i = 1; i < n; ++i) {
    const double denominator = A.diag[static_cast<std::size_t>(i)] -
                               A.lower[static_cast<std::size_t>(i)] *
                                   c[static_cast<std::size_t>(i - 1)];
    if (std::abs(denominator) < kPivotTolerance) {
      throw std::runtime_error("solve_tridiagonal: zero pivot");
    }
    c[static_cast<std::size_t>(i)] =
        (i + 1 < n) ? A.upper[static_cast<std::size_t>(i)] / denominator : 0.0;
    d[static_cast<std::size_t>(i)] =
        (b[static_cast<std::size_t>(i)] -
         A.lower[static_cast<std::size_t>(i)] *
             d[static_cast<std::size_t>(i - 1)]) /
        denominator;
  }

  x[static_cast<std::size_t>(n - 1)] = d[static_cast<std::size_t>(n - 1)];
  for (int i = n - 2; i >= 0; --i) {
    x[static_cast<std::size_t>(i)] =
        d[static_cast<std::size_t>(i)] -
        c[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i + 1)];
  }

  return x;
}

}  // namespace dfn
