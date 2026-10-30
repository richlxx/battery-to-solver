#include "dfn/operators.hpp"

#include <stdexcept>

namespace dfn {

std::vector<double> TridiagonalMatrix::multiply(
    const std::vector<double>& x) const {
  if (static_cast<int>(x.size()) != n) {
    throw std::invalid_argument("tridiagonal multiply: size mismatch");
  }

  std::vector<double> y(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    double value = diag[i] * x[i];
    if (i > 0) {
      value += lower[i] * x[i - 1];
    }
    if (i + 1 < n) {
      value += upper[i] * x[i + 1];
    }
    y[i] = value;
  }
  return y;
}

DiffusionOperator assemble_diffusion(const Grid1D& grid,
                                     const std::vector<double>& D_cell,
                                     double left_flux,
                                     double right_flux) {
  if (static_cast<int>(D_cell.size()) != grid.n) {
    throw std::invalid_argument("diffusion operator: D_cell size mismatch");
  }

  const int n = grid.n;
  DiffusionOperator op;
  op.matrix.n = n;
  op.matrix.lower.assign(static_cast<std::size_t>(n), 0.0);
  op.matrix.diag.assign(static_cast<std::size_t>(n), 0.0);
  op.matrix.upper.assign(static_cast<std::size_t>(n), 0.0);
  op.source.assign(static_cast<std::size_t>(n), 0.0);

  std::vector<double> D_face(static_cast<std::size_t>(n + 1), 0.0);
  D_face[0] = D_cell[0];
  D_face[static_cast<std::size_t>(n)] = D_cell[static_cast<std::size_t>(n - 1)];
  for (int i = 1; i < n; ++i) {
    const double a = D_cell[static_cast<std::size_t>(i - 1)];
    const double b = D_cell[static_cast<std::size_t>(i)];
    D_face[static_cast<std::size_t>(i)] = 2.0 * a * b / (a + b);
  }

  for (int i = 0; i < n; ++i) {
    double west = 0.0;
    if (i > 0) {
      west = D_face[static_cast<std::size_t>(i)] * grid.area[static_cast<std::size_t>(i)] /
             grid.spacing[static_cast<std::size_t>(i)] /
             grid.volume[static_cast<std::size_t>(i)];
    }
    double east = 0.0;
    if (i + 1 < n) {
      east = D_face[static_cast<std::size_t>(i + 1)] *
             grid.area[static_cast<std::size_t>(i + 1)] /
             grid.spacing[static_cast<std::size_t>(i + 1)] /
             grid.volume[static_cast<std::size_t>(i)];
    }
    op.matrix.diag[static_cast<std::size_t>(i)] = -(west + east);
    op.matrix.lower[static_cast<std::size_t>(i)] = west;
    op.matrix.upper[static_cast<std::size_t>(i)] = east;
  }

  op.source[0] -=
      grid.area[0] * left_flux / grid.volume[0];
  op.source[static_cast<std::size_t>(n - 1)] -=
      grid.area[static_cast<std::size_t>(n)] * right_flux /
      grid.volume[static_cast<std::size_t>(n - 1)];

  return op;
}

}  // namespace dfn
