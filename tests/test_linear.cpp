// 测试文件：稠密与三对角线性求解器（对应 linalg 模块）。
// 验证部分主元、LU 复用、Thomas 算法以及奇异/非方阵的错误处理。

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

// 混合绝对/相对容差：既允许接近 0 的分量有 1e-10 的绝对误差，
// 也允许大分量按相对比例缩放，避免量纲差异导致误判。
bool close(double a, double b, double tolerance = 1e-10) {
  return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

}  // namespace

int main() {
  // 2x2 方程组：2x + y = 5，x + 3y = 10，解析解为 (1, 3)。
  const std::vector<double> x =
      dfn::solve_dense({{2.0, 1.0}, {1.0, 3.0}}, {5.0, 10.0});
  if (!close(x[0], 1.0) || !close(x[1], 3.0)) {
    return fail("dense 2x2 solve");
  }

  // 首行主元为 0，必须交换行才能消元；解为 (3, 2)。
  const std::vector<double> x_pivot =
      dfn::solve_dense({{0.0, 1.0}, {1.0, 0.0}}, {2.0, 3.0});
  if (!close(x_pivot[0], 3.0) || !close(x_pivot[1], 2.0)) {
    return fail("partial pivoting");
  }

  // 同一系数矩阵只分解一次、多次回代：A x = (1, 0) 的解析解为 (0.6, -0.2)。
  const dfn::LuFactorization factorization =
      dfn::lu_factor({{2.0, 1.0}, {1.0, 3.0}});
  const std::vector<double> x_reuse =
      dfn::lu_solve(factorization, {1.0, 0.0});
  if (!close(x_reuse[0], 0.6) || !close(x_reuse[1], -0.2)) {
    return fail("LU reuse");
  }

  // 构造对角占优的三对角系统，精确解取线性分布 x_i = 1 + 0.1 i，
  // 再由矩阵-向量乘生成右端项，用于校验 Thomas 算法与稠密求解一致。
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

  // 两种求解路径应给出同一精确解：Thomas 为 O(n)，稠密为 O(n^3) 参考实现。
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

  // 奇异矩阵（两行相同、行列式为 0）必须抛异常，避免静默返回 NaN/Inf。
  bool singular_threw = false;
  try {
    (void)dfn::solve_dense({{1.0, 1.0}, {1.0, 1.0}}, {1.0, 2.0});
  } catch (const std::exception&) {
    singular_threw = true;
  }
  if (!singular_threw) {
    return fail("singular matrix must throw");
  }

  // 非方阵无法求解，必须在入口处拒绝，而不是越界读取。
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
