#include "dfn/solvers/dae.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace dfn {

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
  for (int col = 0; col < n; ++col) {
    pattern.col_ptr[static_cast<std::size_t>(col) + 1] +=
        pattern.col_ptr[static_cast<std::size_t>(col)];
  }
  return pattern;
}

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

CscBuilder::CscBuilder(const CscPattern& pattern)
    : pattern_(pattern),
      values_(static_cast<std::size_t>(pattern.nnz()), 0.0) {}

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

  model.residual(t, y_base.data(), yp_base.data(), r0.data());

  constexpr double kStep = 1e-7;
  for (int j = 0; j < n; ++j) {
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
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] =
          values[static_cast<std::size_t>(col * n + row)];
    }
  }
  return matrix;
}

}  // namespace dfn
