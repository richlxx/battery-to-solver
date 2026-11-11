#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/linalg/dense.hpp"
#include "dfn/linalg/tridiagonal.hpp"

// 实验（第 5 章 线性求解）：对 n 阶对角占优三对角系统（对角 4、次/超对角 -1）
// 比较稠密 LU（O(n^3)）与 Thomas 追赶法（O(n)）的耗时与最大误差。
// 输出 CSV 列：n, dense_ms, thomas_ms, dense_max_error, thomas_max_error。
namespace {

// 构造三对角矩阵：对角 4、次对角与超对角 -1，首行无下对角、末行无上对角。
// 对角占优（4 > 1 + 1）保证追赶法无需选主元即可稳定求解。
dfn::TridiagonalMatrix make_matrix(int n) {
  dfn::TridiagonalMatrix matrix;
  matrix.n = n;
  matrix.lower.assign(static_cast<std::size_t>(n), -1.0);
  matrix.diag.assign(static_cast<std::size_t>(n), 4.0);
  matrix.upper.assign(static_cast<std::size_t>(n), -1.0);
  matrix.lower[0] = 0.0;
  matrix.upper[static_cast<std::size_t>(n - 1)] = 0.0;
  return matrix;
}

// 预设精确解 x_i = 1 + 0.01 i，用来生成右端项并度量两种解法的误差。
std::vector<double> make_exact(int n) {
  std::vector<double> exact(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    exact[static_cast<std::size_t>(i)] = 1.0 + 0.01 * static_cast<double>(i);
  }
  return exact;
}

// 右端项 b = A * x_exact：按三对角行展开，只累加实际存在的三个非零元。
std::vector<double> make_rhs(const dfn::TridiagonalMatrix& matrix,
                             const std::vector<double>& exact) {
  const int n = matrix.n;
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
  return rhs;
}

// 把三对角存储展开成稠密矩阵，作为 solve_dense 的输入（仅用于性能对比）。
std::vector<std::vector<double>> to_dense(const dfn::TridiagonalMatrix& matrix) {
  const int n = matrix.n;
  std::vector<std::vector<double>> dense(
      static_cast<std::size_t>(n),
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
  return dense;
}

// 分量最大绝对误差（无穷范数），两种解法都应只有舍入误差量级。
double max_error(const std::vector<double>& x, const std::vector<double>& exact) {
  double error = 0.0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    error = std::max(error, std::abs(x[i] - exact[i]));
  }
  return error;
}

}  // namespace

int main() {
  std::cout << "n,dense_ms,thomas_ms,dense_max_error,thomas_max_error\n"
            << std::scientific << std::setprecision(3);

  // 规模依次翻倍：稠密 LU 耗时约按 8 倍增长（O(n^3)），
  // Thomas 约按 2 倍增长（O(n)），两者差距随 n 迅速拉开。
  for (int n : {50, 100, 200, 400}) {
    const dfn::TridiagonalMatrix matrix = make_matrix(n);
    const std::vector<double> exact = make_exact(n);
    const std::vector<double> rhs = make_rhs(matrix, exact);
    const std::vector<std::vector<double>> dense = to_dense(matrix);

    // 稠密 LU 分解 + 前代回代（O(n^3)），用 steady_clock 计时（毫秒）。
    const auto dense_start = std::chrono::steady_clock::now();
    const std::vector<double> x_dense = dfn::solve_dense(dense, rhs);
    const auto dense_end = std::chrono::steady_clock::now();

    // Thomas 追赶法（O(n)）：只消元两条对角线，无需存储完整矩阵。
    const auto thomas_start = std::chrono::steady_clock::now();
    const std::vector<double> x_thomas = dfn::solve_tridiagonal(matrix, rhs);
    const auto thomas_end = std::chrono::steady_clock::now();

    const double dense_ms =
        std::chrono::duration<double, std::milli>(dense_end - dense_start).count();
    const double thomas_ms =
        std::chrono::duration<double, std::milli>(thomas_end - thomas_start).count();

    std::cout << n << ',' << dense_ms << ',' << thomas_ms << ','
              << max_error(x_dense, exact) << ','
              << max_error(x_thomas, exact) << '\n';
  }

  return 0;
}
