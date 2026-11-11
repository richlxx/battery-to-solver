// ChainDae：用于教学与求解器自测的链式 DAE 模型。
// 每个单元含一个刚性松弛微分方程和一个代数约束，未知量按
// (c_0, w_0, c_1, w_1, ...) 交错排列；另附单单元解析解作为参照。
#include "dfn/models/chain_dae.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace dfn {

namespace {

// 索引约定：偶数下标存放微分变量 c_i，奇数下标存放代数变量 w_i。
int cell_index(int cell) { return 2 * cell; }
int aux_index(int cell) { return 2 * cell + 1; }

}  // namespace

// 校验参数并预先声明雅可比稀疏结构。
// 第 i 个微分行 c_i 与 c_{i-1}、c_i、c_{i+1}、w_i 耦合；
// 代数行 w_i 只与自身耦合，因此矩阵只有少量非零元。
ChainDae::ChainDae(ChainDaeOptions options) : options_(options) {
  if (options_.cells < 1) {
    throw std::invalid_argument("ChainDae: cells must be >= 1");
  }
  if (options_.lambda <= 0.0) {
    throw std::invalid_argument("ChainDae: lambda must be positive");
  }
  if (options_.diffusion < 0.0) {
    throw std::invalid_argument("ChainDae: diffusion must be non-negative");
  }

  std::vector<std::pair<int, int>> entries;
  entries.reserve(static_cast<std::size_t>(8 * options_.cells));
  // 按行登记非零元：先对角，再左、右邻居，最后登记松弛耦合项 (c, w)。
  for (int i = 0; i < options_.cells; ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    entries.emplace_back(c, c);
    entries.emplace_back(c, w);
    if (i > 0) {
      entries.emplace_back(c, cell_index(i - 1));
    }
    if (i + 1 < options_.cells) {
      entries.emplace_back(c, cell_index(i + 1));
    }
    entries.emplace_back(w, w);
  }
  pattern_ = build_pattern(size(), std::move(entries));
}

// 计算残差 F(t, y, y') = 0，y 为状态、yp 为时间导数、out 为残差。
// 第 i 个单元的两条方程为
//   c_i' + lambda (c_i - w_i) - diffusion * lap(c)_i = 0
//   w_i - sin(t) = 0
// 边界用镜像（零通量）方式处理：域外邻居直接取自身值。
void ChainDae::residual(double t, const double* y, const double* yp,
                        double* out) const {
  const double d = options_.diffusion;
  for (int i = 0; i < options_.cells; ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    const double left = (i > 0) ? y[cell_index(i - 1)] : y[c];
    const double right = (i + 1 < options_.cells) ? y[cell_index(i + 1)] : y[c];
    // 三点差分：lap(c)_i = c_{i-1} - 2 c_i + c_{i+1}。
    const double laplacian = left + right - 2.0 * y[c];

    // 刚性松弛项把 c_i 拉向 w_i（lambda 越大越刚），扩散项带负号。
    out[c] = yp[c] + options_.lambda * (y[c] - y[w]) - d * laplacian;
    // 代数约束：w_i 由时间显式给定，无时间导数。
    out[w] = y[w] - std::sin(t);
  }
}

// 变量类型标记：1 表示微分变量（c_i），0 表示代数变量（w_i）。
std::vector<int> ChainDae::variable_types() const {
  std::vector<int> types(static_cast<std::size_t>(size()), 0);
  for (int i = 0; i < options_.cells; ++i) {
    types[static_cast<std::size_t>(cell_index(i))] = 1;
  }
  return types;
}

// 解析雅可比，形式为 dF/dy + cj * dF/dy'。
// cj 是求解器给出的时间步系数（BDF 主对角贡献，约等于 1/h）；
// 本模型 F 对 y、yp 都是线性的，因此不需要 t、y、yp。
void ChainDae::jacobian(double t, const double* y, const double* yp, double cj,
                        CscBuilder& out) const {
  (void)t;
  (void)y;
  (void)yp;
  const double d = options_.diffusion;
  for (int i = 0; i < options_.cells; ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    const bool interior = (i > 0) && (i + 1 < options_.cells);
    // lap 对 c_i 的导数：内部单元为 2d；只有单侧邻居（镜像边界）时为 d；
    // 单单元时差分恒为零，故取 0。
    const double laplacian_diagonal =
        interior ? 2.0 * d : ((options_.cells > 1) ? d : 0.0);
    out.add(c, c, cj + options_.lambda + laplacian_diagonal);
    // 左右邻居的耦合系数均为 -d。
    if (i > 0) {
      out.add(c, cell_index(i - 1), -d);
    }
    if (i + 1 < options_.cells) {
      out.add(c, cell_index(i + 1), -d);
    }
    // 松弛项 -lambda * w_i 对 w_i 的导数。
    out.add(c, w, -options_.lambda);
    // 代数行 w_i - sin(t) 对 w_i 的导数为 1。
    out.add(w, w, 1.0);
  }
}

// 构造 t0 时刻的一致初始状态：所有单元取单单元解析解，
// 代数变量取 sin(t0)，与 chain_dae_consistent_derivative 配套使用。
std::vector<double> chain_dae_consistent_state(const ChainDae& model,
                                               double t0) {
  return chain_dae_mode_state(model, t0, 0.0);
}

// 在均匀解上叠加正弦模态扰动：
//   c_i = exact(t0) + amplitude * sin((i + 1) pi / (cells + 1))
// 该形状是离散 Laplace 算子的特征向量，便于观察误差随时间的衰减。
std::vector<double> chain_dae_mode_state(const ChainDae& model, double t0,
                                         double amplitude) {
  std::vector<double> y(static_cast<std::size_t>(model.size()), 0.0);
  const double concentration =
      chain_dae_exact_concentration(model.lambda(), t0);
  // 相位间隔：特征向量把半波长均分到 cells + 1 个区间。
  const double spacing =
      3.14159265358979323846 / static_cast<double>(model.cells() + 1);
  for (int i = 0; i < model.cells(); ++i) {
    // 第 i 个单元上的模态形状；amplitude = 0 时退化为均匀解。
    const double shape = std::sin(static_cast<double>(i + 1) * spacing);
    y[static_cast<std::size_t>(cell_index(i))] =
        concentration + amplitude * shape;
    y[static_cast<std::size_t>(aux_index(i))] = std::sin(t0);
  }
  return y;
}

// 由给定状态 y 反解出满足 F(t, y, y') = 0 的导数 y'，
// 即把残差分别对 c_i' 和 w_i' 整理成显式形式，用于生成一致初值。
std::vector<double> chain_dae_derivative_from_state(
    const ChainDae& model, double t, const std::vector<double>& y) {
  std::vector<double> yp(static_cast<std::size_t>(model.size()), 0.0);
  const double d = model.diffusion();
  for (int i = 0; i < model.cells(); ++i) {
    const int c = cell_index(i);
    const int w = aux_index(i);
    const double left = (i > 0) ? y[static_cast<std::size_t>(cell_index(i - 1))]
                                : y[static_cast<std::size_t>(c)];
    const double right =
        (i + 1 < model.cells())
            ? y[static_cast<std::size_t>(cell_index(i + 1))]
            : y[static_cast<std::size_t>(c)];
    const double laplacian =
        left + right - 2.0 * y[static_cast<std::size_t>(c)];
    // 由 c_i' = -lambda (c_i - w_i) + d * lap(c)_i 解出时间导数。
    yp[static_cast<std::size_t>(c)] =
        -model.lambda() * (y[static_cast<std::size_t>(c)] -
                           y[static_cast<std::size_t>(w)]) +
        d * laplacian;
    // 由代数约束 w_i = sin(t) 求导得 w_i' = cos(t)。
    yp[static_cast<std::size_t>(w)] = std::cos(t);
  }
  return yp;
}

// 均匀初始状态对应的一致导数。
std::vector<double> chain_dae_consistent_derivative(const ChainDae& model,
                                                    double t0) {
  return chain_dae_derivative_from_state(
      model, t0, chain_dae_mode_state(model, t0, 0.0));
}

// 单单元且无扩散时 c' = -lambda (c - sin t)、c(0) = 0 的解析解：
//   c(t) = lambda (lambda sin t - cos t + exp(-lambda t)) / (lambda^2 + 1)
// 其中 exp(-lambda t) 是刚性瞬态，随时间迅速衰减到周期解上。
double chain_dae_exact_concentration(double lambda, double t) {
  if (lambda <= 0.0) {
    throw std::invalid_argument("chain_dae_exact_concentration: lambda > 0");
  }
  const double denominator = lambda * lambda + 1.0;
  return lambda *
         (lambda * std::sin(t) - std::cos(t) + std::exp(-lambda * t)) /
         denominator;
}

// 对上式求时间导数得到的解析变化率：
//   c'(t) = lambda (lambda cos t + sin t - lambda exp(-lambda t)) / (lambda^2 + 1)
double chain_dae_exact_concentration_rate(double lambda, double t) {
  if (lambda <= 0.0) {
    throw std::invalid_argument(
        "chain_dae_exact_concentration_rate: lambda > 0");
  }
  const double denominator = lambda * lambda + 1.0;
  return lambda *
         (lambda * std::cos(t) + std::sin(t) - lambda * std::exp(-lambda * t)) /
         denominator;
}

}  // namespace dfn
