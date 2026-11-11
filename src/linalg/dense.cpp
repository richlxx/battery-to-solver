// 本文件实现稠密线性方程组的直接求解：带部分主元选取的 LU 分解，
// 以及基于该分解的前代/回代。用于小规模稠密系统与教学对照。

#include "dfn/linalg/dense.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace dfn {

namespace {
// 主元绝对值下限：小于它即视为矩阵奇异。1e-14 取在双精度舍入误差
// （约 2.2e-16）之上，只用于检测结构性奇异，不代表可接受的病态程度。
constexpr double kPivotTolerance = 1e-14;
}

// 对按行存储的方阵 A 做带部分主元的 LU 分解。
// 参数 A：按值传入，函数内就地分解，调用方的矩阵不受影响。
// 返回：LuFactorization。lu 的严格下三角存放消元乘子（L 的对角为 1）、
//       上三角为 U，即 P*A = L*U；perm[i] 表示第 i 行方程来自原始第 perm[i] 行。
// 若某列主元绝对值小于 kPivotTolerance，抛出 std::runtime_error。
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
  // perm 初始为单位排列；之后每次行交换都要同步交换，保证右端项重排一致。
  std::iota(result.perm.begin(), result.perm.end(), 0);

  for (int k = 0; k < n; ++k) {
    // 部分主元选取：在第 k 列对角线及以下找绝对值最大的元素并换到第 k 行，
    // 以抑制消元过程中的舍入误差放大。
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

    // 消元：A[i][k] 先除以主元归一化为乘子（随即写回原位），
    // 再用该乘子更新第 k 行右侧的子矩阵。
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

// 利用 lu_factor 的结果求解 A x = b。
// 参数：factorization 为已完成的 LU 分解；b 为右端向量，长度必须为 n。
// 返回：解向量 x（长度 n）。
// 步骤：先按 perm 重排 b（即 P*b），再前代解 L y = P*b，最后回代解 U x = y。
std::vector<double> lu_solve(const LuFactorization& factorization,
                             const std::vector<double>& b) {
  if (static_cast<int>(b.size()) != factorization.n) {
    throw std::invalid_argument("lu_solve: right-hand side size mismatch");
  }

  const int n = factorization.n;
  std::vector<double> x(static_cast<std::size_t>(n), 0.0);
  // 载入右端项时直接按 perm 重排，x 随后先存 y、再被覆盖为最终解。
  for (int i = 0; i < n; ++i) {
    x[static_cast<std::size_t>(i)] =
        b[static_cast<std::size_t>(factorization.perm[static_cast<std::size_t>(i)])];
  }

  // 前代：L 为单位下三角（对角为 1），乘子存放在 lu 的严格下三角，
  // 因此只需减去已知项，不做除法。
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      x[static_cast<std::size_t>(i)] -=
          factorization.lu[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
          x[static_cast<std::size_t>(j)];
    }
  }

  // 回代：从末行向上，先减去已知项，再除以对角元 U[i][i]。
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

// 便捷包装：对 A 的一份拷贝做 LU 分解并求解；A 与 b 均按值传入，调用方数据不变。
std::vector<double> solve_dense(std::vector<std::vector<double>> A,
                                std::vector<double> b) {
  const LuFactorization factorization = lu_factor(std::move(A));
  return lu_solve(factorization, b);
}

}  // namespace dfn
