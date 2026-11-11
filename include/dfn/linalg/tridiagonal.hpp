#pragma once

#include <vector>

#include "dfn/operators.hpp"

namespace dfn {

// 三对角线性系统求解（Thomas 算法），左端为 operators.hpp 中的
// TridiagonalMatrix，是扩散隐式步的默认求解器。

// 求解 A x = b：A 为 n 阶三对角矩阵，b 长度 n，返回 x 长度 n。
std::vector<double> solve_tridiagonal(const TridiagonalMatrix& A,
                                      const std::vector<double>& b);

}  // namespace dfn
