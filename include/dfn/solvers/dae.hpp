#pragma once

#include <utility>
#include <vector>

namespace dfn {

// CSC pattern: column j owns rows row_index[col_ptr[j] .. col_ptr[j+1]).
// Rows are sorted inside each column so that lookup is a binary search.
struct CscPattern {
  int n = 0;
  std::vector<int> col_ptr;
  std::vector<int> row_index;

  int nnz() const { return static_cast<int>(row_index.size()); }
  int index_of(int row, int col) const;
};

// Builds a CSC pattern from (row, col) pairs; duplicates are merged.
CscPattern build_pattern(int n, std::vector<std::pair<int, int>> entries);
CscPattern dense_pattern(int n);

// Writes values into a fixed pattern. An entry that is missing from the
// pattern is an assembly bug, so it throws instead of silently growing.
class CscBuilder {
 public:
  explicit CscBuilder(const CscPattern& pattern);

  void add(int row, int col, double value);
  void set(int row, int col, double value);

  const std::vector<double>& values() const { return values_; }
  std::vector<double> take_values() { return std::move(values_); }

 private:
  CscPattern pattern_;
  std::vector<double> values_;
};

// Residual form F(t, y, y') = 0. Differential rows depend on y', algebraic
// rows do not; the model itself does not care which is which.
class DaeModel {
 public:
  virtual ~DaeModel() = default;

  virtual int size() const = 0;
  virtual void residual(double t, const double* y, const double* yp,
                        double* out) const = 0;

  virtual CscPattern jacobian_pattern() const = 0;

  // 1 = differential variable, 0 = algebraic variable. IDA needs this to
  // compute consistent initial conditions.
  virtual std::vector<int> variable_types() const {
    return std::vector<int>(static_cast<std::size_t>(size()), 1);
  }

  // Must add dF/dy + cj * dF/dy' to `out`.
  virtual void jacobian(double t, const double* y, const double* yp, double cj,
                        CscBuilder& out) const = 0;
};

// Dense finite-difference approximation of dF/dy + cj * dF/dy'.
// The pattern must be dense; used for teaching and as a correctness check.
void finite_difference_jacobian(const DaeModel& model, double t, const double* y,
                                const double* yp, double cj,
                                CscBuilder& out);

std::vector<std::vector<double>> dense_jacobian(const DaeModel& model, double t,
                                                const double* y, const double* yp,
                                                double cj);

}  // namespace dfn
