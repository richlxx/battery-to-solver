#pragma once

#include <vector>

#include "dfn/operators.hpp"

namespace dfn {

std::vector<double> solve_tridiagonal(const TridiagonalMatrix& A,
                                      const std::vector<double>& b);

}  // namespace dfn
