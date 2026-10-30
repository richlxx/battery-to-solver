#include "dfn/linalg/dense.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace dfn {

namespace {
constexpr double kPivotTolerance = 1e-14;
}

LuFactorization lu_factor(std::vector<std::vector<double>> A) {
  const int n = static_cast<int>(A.size());
  if (n == 0) {
    throw std::invalid_argument("lu_factor: empty matrix");
  }
  for (const auto& row : A) {
    if (static_cast<int>(row.size()) != n) {
      throw std::invalid_argument("lu_factor: matrix must be square");
    }
  }

  LuFactorization result;
  result.n = n;
  result.perm.resize(static_cast<std::size_t>(n));
  std::iota(result.perm.begin(), result.perm.end(), 0);

  for (int k = 0; k < n; ++k) {
    int pivot = k;
    double best = std::abs(A[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)]);
    for (int i = k + 1; i < n; ++i) {
      const double value = std::abs(A[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)]);
      if (value > best) {
        best = value;
        pivot = i;
      }
    }
    if (best < kPivotTolerance) {
      throw std::runtime_error("lu_factor: matrix is singular");
    }
    if (pivot != k) {
      std::swap(A[static_cast<std::size_t>(pivot)], A[static_cast<std::size_t>(k)]);
      std::swap(result.perm[static_cast<std::size_t>(pivot)],
                result.perm[static_cast<std::size_t>(k)]);
    }

    for (int i = k + 1; i < n; ++i) {
      A[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] /=
          A[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)];
      const double multiplier =
          A[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
      for (int j = k + 1; j < n; ++j) {
        A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] -=
            multiplier * A[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
      }
    }
  }

  result.lu = std::move(A);
  return result;
}

std::vector<double> lu_solve(const LuFactorization& factorization,
                             const std::vector<double>& b) {
  if (static_cast<int>(b.size()) != factorization.n) {
    throw std::invalid_argument("lu_solve: right-hand side size mismatch");
  }

  const int n = factorization.n;
  std::vector<double> x(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    x[static_cast<std::size_t>(i)] =
        b[static_cast<std::size_t>(factorization.perm[static_cast<std::size_t>(i)])];
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      x[static_cast<std::size_t>(i)] -=
          factorization.lu[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
          x[static_cast<std::size_t>(j)];
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j < n; ++j) {
      x[static_cast<std::size_t>(i)] -=
          factorization.lu[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
          x[static_cast<std::size_t>(j)];
    }
    x[static_cast<std::size_t>(i)] /=
        factorization.lu[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
  }

  return x;
}

std::vector<double> solve_dense(std::vector<std::vector<double>> A,
                                std::vector<double> b) {
  const LuFactorization factorization = lu_factor(std::move(A));
  return lu_solve(factorization, b);
}

}  // namespace dfn
