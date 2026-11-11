// DAE 基础设施：CSC 稀疏模式的构建与查询、按固定模式装配数值的 CscBuilder，
// 以及有限差分 Jacobian 与稠密 Jacobian 的组装（dF/dy + cj·dF/dy'）。

#include "dfn/solvers/dae.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace dfn {

// 在 CSC 模式里查 (row, col) 的存储位置，返回 values 数组下标。
// 每列内的行号已升序排列，所以用二分查找；越界或不在模式中时返回 -1。
int CscPattern::index_of(int row, int col) const {
  if (col < 0 || col >= n || row < 0 || row >= n) {
    return -1;
  }
  const auto begin = row_index.begin() + col_ptr[static_cast<std::size_t>(col)];
  const auto end = row_index.begin() + col_ptr[static_cast<std::size_t>(col + 1)];
  const auto it = std::lower_bound(begin, end, row);
  if (it == end || *it != row) {
    return -1;
  }
  return static_cast<int>(it - row_index.begin());
}

// 由 (row, col) 列表构建 CSC 模式（列优先存储）：越界直接报错，
// 然后按列、列内按行排序并去掉重复项，最后用前缀和得到列起始偏移。
// entries 按值传入，排序去重不会影响调用方数据。
CscPattern build_pattern(int n, std::vector<std::pair<int, int>> entries) {
  if (n <= 0) {
    throw std::invalid_argument("build_pattern: n must be positive");
  }
  for (const auto& [row, col] : entries) {
    if (row < 0 || row >= n || col < 0 || col >= n) {
      throw std::out_of_range("build_pattern: entry (" + std::to_string(row) +
                              "," + std::to_string(col) + ") outside " +
                              std::to_string(n) + "x" + std::to_string(n));
    }
  }
  std::sort(entries.begin(), entries.end(),
            [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
              return a.second != b.second ? a.second < b.second
                                          : a.first < b.first;
            });
  entries.erase(std::unique(entries.begin(), entries.end()), entries.end());

  CscPattern pattern;
  pattern.n = n;
  pattern.col_ptr.assign(static_cast<std::size_t>(n) + 1, 0);
  pattern.row_index.reserve(entries.size());
  for (const auto& [row, col] : entries) {
    ++pattern.col_ptr[static_cast<std::size_t>(col) + 1];
    pattern.row_index.push_back(row);
  }
  // 前缀和：把“每列非零个数”转换成 CSC 约定的列起始偏移 col_ptr。
  for (int col = 0; col < n; ++col) {
    pattern.col_ptr[static_cast<std::size_t>(col) + 1] +=
        pattern.col_ptr[static_cast<std::size_t>(col)];
  }
  return pattern;
}

// 生成 n×n 全稠密 CSC 模式：按列优先顺序列出所有 (row, col) 位置。
CscPattern dense_pattern(int n) {
  std::vector<std::pair<int, int>> entries;
  entries.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      entries.emplace_back(row, col);
    }
  }
  return build_pattern(n, std::move(entries));
}

// 按给定模式预分配数值数组（初值全 0）；模式在构造时固定，之后只改数值。
CscBuilder::CscBuilder(const CscPattern& pattern)
    : pattern_(pattern),
      values_(static_cast<std::size_t>(pattern.nnz()), 0.0) {}

// 累加一个元素：有限差分 Jacobian 会把 dF/dy 与 cj·dF/dy' 两块加在同一位置，
// 所以必须用 add 而不是 set。模式中没有该位置说明装配与模式不一致，直接抛错。
void CscBuilder::add(int row, int col, double value) {
  const int k = pattern_.index_of(row, col);
  if (k < 0) {
    throw std::out_of_range(
        "CscBuilder: entry (" + std::to_string(row) + "," + std::to_string(col) +
        ") is not in the Jacobian pattern; the pattern and the assembly "
        "disagree");
  }
  values_[static_cast<std::size_t>(k)] += value;
}

// 覆盖一个元素；位置不在模式中同样视为装配错误并抛错。
void CscBuilder::set(int row, int col, double value) {
  const int k = pattern_.index_of(row, col);
  if (k < 0) {
    throw std::out_of_range(
        "CscBuilder: entry (" + std::to_string(row) + "," + std::to_string(col) +
        ") is not in the Jacobian pattern; the pattern and the assembly "
        "disagree");
  }
  values_[static_cast<std::size_t>(k)] = value;
}

// 稠密有限差分 Jacobian：向 out 累加 Jacobian = dF/dy + cj·dF/dy'。
// 对每一列 j 依次扰动 y_j（得 ∂F/∂y_j）和 y'_j（得 ∂F/∂y'_j，再乘 cj），
// 差分步长取 1e-7·max(1, |分量|)，兼顾舍入误差与截断误差。
// cj 是 IDA 中的标量：y' = cj·y + 常数，故 dF/dy 实际是 dF/dy + cj·dF/dy'。
void finite_difference_jacobian(const DaeModel& model, double t, const double* y,
                                const double* yp, double cj,
                                CscBuilder& out) {
  const int n = model.size();
  std::vector<double> y_base(y, y + n);
  std::vector<double> yp_base(yp, yp + n);
  std::vector<double> y_pert = y_base;
  std::vector<double> yp_pert = yp_base;
  std::vector<double> r0(static_cast<std::size_t>(n));
  std::vector<double> r1(static_cast<std::size_t>(n));

  // 基准残差 R0 = F(t, y, y')，后续所有差分都与它相减。
  model.residual(t, y_base.data(), yp_base.data(), r0.data());

  constexpr double kStep = 1e-7;
  for (int j = 0; j < n; ++j) {
    // 第 j 列：先扰动 y_j 得到 ∂F/∂y_j 的整列，再扰动 y'_j 并乘 cj 叠加。
    const double dj = kStep * std::max(1.0, std::abs(y_base[static_cast<std::size_t>(j)]));
    y_pert[static_cast<std::size_t>(j)] += dj;
    model.residual(t, y_pert.data(), yp_base.data(), r1.data());
    for (int i = 0; i < n; ++i) {
      out.add(i, j,
              (r1[static_cast<std::size_t>(i)] - r0[static_cast<std::size_t>(i)]) / dj);
    }
    y_pert[static_cast<std::size_t>(j)] = y_base[static_cast<std::size_t>(j)];

    const double djp =
        kStep * std::max(1.0, std::abs(yp_base[static_cast<std::size_t>(j)]));
    yp_pert[static_cast<std::size_t>(j)] += djp;
    model.residual(t, y_base.data(), yp_pert.data(), r1.data());
    for (int i = 0; i < n; ++i) {
      out.add(i, j,
              cj * (r1[static_cast<std::size_t>(i)] - r0[static_cast<std::size_t>(i)]) /
                  djp);
    }
    yp_pert[static_cast<std::size_t>(j)] = yp_base[static_cast<std::size_t>(j)];
  }
}

// 稠密版 Jacobian：借用稠密 CSC 模式让 model.jacobian 装配数值，
// 再转成行优先的嵌套 vector（J[i][j] = ∂F_i/∂y_j），供 newton_solve 直接使用。
std::vector<std::vector<double>> dense_jacobian(const DaeModel& model, double t,
                                                const double* y, const double* yp,
                                                double cj) {
  const int n = model.size();
  CscBuilder builder(dense_pattern(n));
  model.jacobian(t, y, yp, cj, builder);
  const std::vector<double>& values = builder.values();

  std::vector<std::vector<double>> matrix(
      static_cast<std::size_t>(n),
      std::vector<double>(static_cast<std::size_t>(n), 0.0));
  // CSC 数值按列连续存放：第 col 列的 n 个值下标为 col*n + row（行已升序）。
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] =
          values[static_cast<std::size_t>(col * n + row)];
    }
  }
  return matrix;
}

}  // namespace dfn
