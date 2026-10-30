#pragma once

#include <vector>

namespace dfn {

struct LuFactorization {
  int n = 0;
  std::vector<std::vector<double>> lu;
  std::vector<int> perm;
};

LuFactorization lu_factor(std::vector<std::vector<double>> A);

std::vector<double> lu_solve(const LuFactorization& factorization,
                             const std::vector<double>& b);

std::vector<double> solve_dense(std::vector<std::vector<double>> A,
                                std::vector<double> b);

}  // namespace dfn
