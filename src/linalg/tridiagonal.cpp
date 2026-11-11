// 本文件实现三对角线性方程组的追赶法（Thomas 算法），即不做主元选取的
// 高斯消元：先向前消元得到上对角系数 c 与变换后的右端 d，再向后回代。

#include "dfn/linalg/tridiagonal.hpp"

#include <cmath>
#include <stdexcept>

namespace dfn {

namespace {
// 主元绝对值下限，与稠密 LU 保持一致；低于它视为零主元，直接报错。
constexpr double kPivotTolerance = 1e-14;
}

// 求解三对角系统 A x = b。
// 参数：A 由 lower/diag/upper 三条对角线描述（lower[0]、upper[n-1] 不使用）；
//       b 为右端向量，长度必须为 n。
// 返回：解向量 x（长度 n）。
// 前提：各步主元非零（扩散类矩阵对角占优、负定，通常满足）；算法不做行交换，
//       因此对角元接近 0 时会抛出异常而不是继续计算。
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

  // 前代所需的工作数组：c 存归一化后的上对角，d 存消元后的右端，
  // 二者相当于 L U x = b 中的 U 与中间向量 y。
  std::vector<double> c(static_cast<std::size_t>(n), 0.0);
  std::vector<double> d(static_cast<std::size_t>(n), 0.0);
  std::vector<double> x(static_cast<std::size_t>(n), 0.0);

  // 首行归一化：除以 diag[0]；n = 1 时系统只有一条方程，没有上对角可用。
  c[0] = (n > 1) ? A.upper[0] / A.diag[0] : 0.0;
  d[0] = b[0] / A.diag[0];

  for (int i = 1; i < n; ++i) {
    // denominator 为消去 x[i-1] 后的主元 diag[i] - lower[i]*c[i-1]；
    // 它接近 0 说明矩阵（数值上）奇异。
    const double denominator = A.diag[static_cast<std::size_t>(i)] -
                               A.lower[static_cast<std::size_t>(i)] *
                                   c[static_cast<std::size_t>(i - 1)];
    if (std::abs(denominator) < kPivotTolerance) {
      throw std::runtime_error("solve_tridiagonal: zero pivot");
    }
    // 末行的上对角不参与运算，c[n-1] 置 0 以统一回代公式。
    c[static_cast<std::size_t>(i)] =
        (i + 1 < n) ? A.upper[static_cast<std::size_t>(i)] / denominator : 0.0;
    d[static_cast<std::size_t>(i)] =
        (b[static_cast<std::size_t>(i)] -
         A.lower[static_cast<std::size_t>(i)] *
             d[static_cast<std::size_t>(i - 1)]) /
        denominator;
  }

  // 向后回代：末行的解就是 d[n-1]，上一行逐次减去上对角的贡献。
  x[static_cast<std::size_t>(n - 1)] = d[static_cast<std::size_t>(n - 1)];
  for (int i = n - 2; i >= 0; --i) {
    x[static_cast<std::size_t>(i)] =
        d[static_cast<std::size_t>(i)] -
        c[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i + 1)];
  }

  return x;
}

}  // namespace dfn
