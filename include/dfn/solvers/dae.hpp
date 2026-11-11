#pragma once

#include <utility>
#include <vector>

namespace dfn {

// DAE（微分-代数方程组）求解基础设施：CSC 稀疏模式、Jacobian 装配器、
// 模型抽象基类 DaeModel，以及有限差分 Jacobian 校验工具。

// 压缩稀疏列（CSC）模式：第 j 列占用 row_index[col_ptr[j] .. col_ptr[j+1])。
// 每列内行号递增排列，便于二分查找。这里只记录结构，不存数值。
struct CscPattern {
  int n = 0;                   // 矩阵阶数（方阵）
  std::vector<int> col_ptr;    // 列指针，长度 n + 1
  std::vector<int> row_index;  // 行号数组，长度 nnz

  // 非零元个数。
  int nnz() const { return static_cast<int>(row_index.size()); }
  // 返回 (row, col) 在 row_index 中的位置；不在模式中时返回 -1。
  int index_of(int row, int col) const;
};

// 由 (row, col) 列表构造 CSC 模式，重复项会被合并；越界项抛异常。
CscPattern build_pattern(int n, std::vector<std::pair<int, int>> entries);
// 构造 n 阶稠密模式（所有位置都视为非零），用于教学与校验。
CscPattern dense_pattern(int n);

// 按固定模式写入数值的装配器。若 (row, col) 不在模式中，说明模式与
// 装配代码不一致，直接抛异常而不是静默扩展稀疏结构。
class CscBuilder {
 public:
  explicit CscBuilder(const CscPattern& pattern);

  // 累加或覆盖写入一个矩阵元素；同一位置有多个贡献时用 add。
  void add(int row, int col, double value);
  void set(int row, int col, double value);

  // 按 CSC 列序排列的数值数组，长度 nnz。
  const std::vector<double>& values() const { return values_; }
  std::vector<double> take_values() { return std::move(values_); }

 private:
  CscPattern pattern_;
  std::vector<double> values_;
};

// 残差形式 F(t, y, y') = 0 的模型抽象基类。微分行依赖 y'，代数行不依赖；
// 模型本身不必区分，由 variable_types() 告知求解器。
class DaeModel {
 public:
  virtual ~DaeModel() = default;

  // 未知量个数（= 残差行数）。
  virtual int size() const = 0;
  // 计算残差 out（长度 size()）：t 为时间，y 为状态，yp 为 dy/dt。
  virtual void residual(double t, const double* y, const double* yp,
                        double* out) const = 0;

  // Jacobian 的稀疏结构，不随 t 与 y 改变。
  virtual CscPattern jacobian_pattern() const = 0;

  // 1 = 微分变量，0 = 代数变量。IDA 需要它来计算一致初值。
  virtual std::vector<int> variable_types() const {
    return std::vector<int>(static_cast<std::size_t>(size()), 1);
  }

  // 必须把 dF/dy + cj * dF/dy' 累加到 out（cj 为 BDF 的差分系数）。
  virtual void jacobian(double t, const double* y, const double* yp, double cj,
                        CscBuilder& out) const = 0;
};

// dF/dy + cj * dF/dy' 的稠密有限差分近似（前向差分，步长按分量缩放）。
// out 必须是稠密模式；用于教学演示与解析 Jacobian 的正确性核对。
void finite_difference_jacobian(const DaeModel& model, double t, const double* y,
                                const double* yp, double cj,
                                CscBuilder& out);

// 返回稠密 Jacobian（按行存放的二维数组），内部调用模型提供的 jacobian()。
std::vector<std::vector<double>> dense_jacobian(const DaeModel& model, double t,
                                                const double* y, const double* yp,
                                                double cj);

}  // namespace dfn
