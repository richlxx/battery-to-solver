#pragma once

#include <vector>

#include "dfn/grid.hpp"

namespace dfn {

struct TridiagonalMatrix {
  int n = 0;
  std::vector<double> lower;  // lower[0] is unused
  std::vector<double> diag;
  std::vector<double> upper;  // upper[n - 1] is unused

  std::vector<double> multiply(const std::vector<double>& x) const;
};

struct DiffusionOperator {
  TridiagonalMatrix matrix;
  std::vector<double> source;  // per unit volume
};

DiffusionOperator assemble_diffusion(const Grid1D& grid,
                                     const std::vector<double>& D_cell,
                                     double left_flux = 0.0,
                                     double right_flux = 0.0);

}  // namespace dfn
