#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dfn/linalg/dense.hpp"
#include "dfn/linalg/tridiagonal.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

bool close(double a, double b, double tolerance = 1e-10) {
  return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

}  // namespace

int main() {
  const std::vector<double> x =
      dfn::solve_dense({{2.0, 1.0}, {1.0, 3.0}}, {5.0, 10.0});
  if (!close(x[0], 1.0) || !close(x[1], 3.0)) {
    return fail("dense 2x2 solve");
  }

  const std::vector<double> x_pivot =
      dfn::solve_dense({{0.0, 1.0}, {1.0, 0.0}}, {2.0, 3.0});
  if (!close(x_pivot[0], 3.0) || !close(x_pivot[1], 2.0)) {
    return fail("partial pivoting");
  }

  const dfn::LuFactorization factorization =
      dfn::lu_factor({{2.0, 1.0}, {1.0, 3.0}});
  const std::vector<double> x_reuse =
      dfn::lu_solve(factorization, {1.0, 0.0});
  if (!close(x_reuse[0], 0.6) || !close(x_reuse[1], -0.2)) {
    return fail("LU reuse");
  }

  const int n = 20;
  dfn::TridiagonalMatrix matrix;
  matrix.n = n;
  matrix.lower.assign(static_cast<std::size_t>(n), -1.0);
  matrix.diag.assign(static_cast<std::size_t>(n), 4.0);
  matrix.upper.assign(static_cast<std::size_t>(n), -1.0);
  matrix.lower[0] = 0.0;
  matrix.upper[static_cast<std::size_t>(n - 1)] = 0.0;

  std::vector<double> exact(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    exact[static_cast<std::size_t>(i)] = 1.0 + 0.1 * static_cast<double>(i);
  }

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

  const std::vector<double> x_tridiagonal = dfn::solve_tridiagonal(matrix, rhs);
  std::vector<std::vector<double>> dense(static_cast<std::size_t>(n),
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
  const std::vector<double> x_dense = dfn::solve_dense(dense, rhs);

  for (int i = 0; i < n; ++i) {
    if (!close(x_tridiagonal[static_cast<std::size_t>(i)],
               exact[static_cast<std::size_t>(i)])) {
      return fail("Thomas solve");
    }
    if (!close(x_dense[static_cast<std::size_t>(i)],
               exact[static_cast<std::size_t>(i)])) {
      return fail("dense tridiagonal solve");
    }
  }

  bool singular_threw = false;
  try {
    (void)dfn::solve_dense({{1.0, 1.0}, {1.0, 1.0}}, {1.0, 2.0});
  } catch (const std::exception&) {
    singular_threw = true;
  }
  if (!singular_threw) {
    return fail("singular matrix must throw");
  }

  bool shape_threw = false;
  try {
    (void)dfn::solve_dense({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}, {1.0, 2.0});
  } catch (const std::exception&) {
    shape_threw = true;
  }
  if (!shape_threw) {
    return fail("non-square matrix must throw");
  }

  std::cout << "linear solvers ok (n=" << n << ")\n";
  return 0;
}
