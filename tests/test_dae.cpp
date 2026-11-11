// 测试文件：CSC 稀疏结构与链式教学 DAE（对应 solvers/dae、models/chain_dae、
// solvers/bdf2_dae 模块）。验证稀疏模式构建、残差/Jacobian 一致性与 BDF2 精度。

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dfn/models/chain_dae.hpp"
#include "dfn/solvers/bdf2_dae.hpp"
#include "dfn/solvers/dae.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

// 混合绝对/相对容差：链条状态为 O(1)，1 + |a| + |b| 使容差随量级自适应。
bool close(double a, double b, double tolerance) {
  return std::abs(a - b) <= tolerance * (1.0 + std::abs(a) + std::abs(b));
}

// 验证 CSC 模式的三要素：非零元计数、行列查找、越界写入必须报错。
int test_pattern() {
  // 3x3 三对角模式（对角 + 次对角 + 超对角）共 7 个非零元。
  const dfn::CscPattern tridiagonal = dfn::build_pattern(3, {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {2, 1}, {1, 2}, {2, 2}});
  if (tridiagonal.nnz() != 7) {
    return fail("tridiagonal nnz");
  }
  // (1,0) 在模式内应能查到有效下标，(0,2) 不在模式内应返回负值。
  if (tridiagonal.index_of(1, 0) < 0 || tridiagonal.index_of(0, 2) >= 0) {
    return fail("tridiagonal lookup");
  }

  // 向模式外写入是装配错误，CscBuilder 必须抛 out_of_range，
  // 否则会悄悄改变稀疏结构、导致求解器内存越界。
  dfn::CscBuilder builder(tridiagonal);
  bool threw = false;
  try {
    builder.set(0, 2, 1.0);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  if (!threw) {
    return fail("missing pattern entry must throw");
  }

  // 稠密模式 4x4 应包含全部 16 个位置。
  const dfn::CscPattern dense = dfn::dense_pattern(4);
  if (dense.nnz() != 16) {
    return fail("dense nnz");
  }
  return 0;
}

// 验证 3 单元链式模型的规模、稀疏模式、一致初值残差与解析 Jacobian。
int test_chain_model() {
  dfn::ChainDaeOptions options;
  options.cells = 3;
  options.lambda = 50.0;   // 中等刚性，便于有限差分 Jacobian 保持良态
  options.diffusion = 4.0;  // 打开扩散耦合，使相邻单元产生非零通量
  const dfn::ChainDae model(options);

  // 每单元 2 个未知量 (c_i, w_i)，3 单元共 6 个。
  if (model.size() != 6) {
    return fail("chain size");
  }
  // 微分行为对角 + w 耦合 + 左右邻居（3+4+3=10），代数行各 1 个对角元，
  // 合计 10 + 3 = 13 个非零元。
  if (model.jacobian_pattern().nnz() != 13) {
    return fail("chain pattern nnz");
  }

  // 一致初值应使残差 F(t, y, y') 精确为 0；阈值取 1e-12，
  // 因为 y 与 y' 由同一闭式解生成，只剩舍入误差。
  double t = 0.3;
  std::vector<double> y = dfn::chain_dae_consistent_state(model, t);
  std::vector<double> yp = dfn::chain_dae_consistent_derivative(model, t);
  std::vector<double> residual(static_cast<std::size_t>(model.size()), 0.0);
  model.residual(t, y.data(), yp.data(), residual.data());
  for (double value : residual) {
    if (std::abs(value) > 1e-12) {
      return fail("consistent state does not zero the residual");
    }
  }

  // 非均匀模态叠加让扩散项真正参与计算，同样必须零残差；
  // 若此处失败，说明 lap(c) 的离散或导数反解有符号/系数错误。
  y = dfn::chain_dae_mode_state(model, t, 0.1);
  yp = dfn::chain_dae_derivative_from_state(model, t, y);
  model.residual(t, y.data(), yp.data(), residual.data());
  for (double value : residual) {
    if (std::abs(value) > 1e-12) {
      return fail("non-uniform consistent state does not zero the residual");
    }
  }

  // 解析 Jacobian 必须与有限差分一致（含 cj*dF/dy' 项）。
  // 人为扰动 y/yp，使各块导数都不为零，避免漏检。
  y[0] = 0.4;
  y[2] = -0.2;
  yp[0] = 1.5;
  yp[4] = -0.7;
  const double cj = 120.0;  // 相当于 h ~ 1/120，接近典型隐式步长

  dfn::CscBuilder analytic(dfn::dense_pattern(model.size()));
  model.jacobian(t, y.data(), yp.data(), cj, analytic);

  dfn::CscBuilder numeric(dfn::dense_pattern(model.size()));
  dfn::finite_difference_jacobian(model, t, y.data(), yp.data(), cj, numeric);

  // 1e-6 的相对容差覆盖前向差分的截断误差与抵消误差；
  // 超差通常意味着解析导数漏项或 cj 乘错位置。
  for (std::size_t k = 0; k < analytic.values().size(); ++k) {
    if (!close(analytic.values()[k], numeric.values()[k], 1e-6)) {
      std::cerr << "analytic " << analytic.values()[k] << " vs numeric "
                << numeric.values()[k] << " at " << k << "\n";
      return fail("Jacobian mismatch");
    }
  }
  return 0;
}

// 单单元闭式解对照：检验 BDF2-DAE 的精度、代数约束与二阶收敛。
int test_single_cell_solution() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;  // 刚性比约 100，显式方法需极小步长
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);

  // 辅助函数自检：t = 0.4 处的一致导数应等于解析变化率（相对容差 1e-12）。
  const std::vector<double> yp_at_start =
      dfn::chain_dae_consistent_derivative(model, 0.4);
  if (!close(yp_at_start[0],
             dfn::chain_dae_exact_concentration_rate(100.0, 0.4), 1e-12)) {
    return fail("consistent derivative disagrees with the closed form");
  }

  // 步长 h = 0.005，Newton 容差收紧到 1e-12，使时间离散误差占主导。
  dfn::Bdf2DaeOptions bdf;
  bdf.newton_tolerance = 1e-12;
  bdf.newton_max_iterations = 30;
  const dfn::Bdf2DaeResult result =
      dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, 0.005, bdf);
  if (!result.ok) {
    return fail(result.message.c_str());
  }

  // 该步长下时间离散误差应远小于 2e-4；容差同时留出启动步与
  // 刚性瞬态的影响，又能抓住方法退化为一阶的情况。
  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);
  const double computed = result.y.back()[0];
  if (!close(computed, exact, 2e-4)) {
    std::cerr << "computed " << computed << " exact " << exact << "\n";
    return fail("single cell solution");
  }
  // 代数约束 w = sin(t) 是强制的，1e-10 仅允许 Newton 收敛残差。
  const double w = result.y.back()[1];
  if (std::abs(w - std::sin(t_end)) > 1e-10) {
    return fail("algebraic constraint at t_end");
  }

  // 步长减半时误差应缩小约 4 倍；[3.0, 5.5] 区间容纳渐近前行为。
  const auto error = [&](double h) {
    const dfn::Bdf2DaeResult run =
        dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, h, bdf);
    if (!run.ok) {
      throw std::runtime_error(run.message);
    }
    return std::abs(run.y.back()[0] - exact);
  };
  const double e1 = error(0.04);
  const double e2 = error(0.02);
  const double e3 = error(0.01);
  const double ratio1 = e1 / e2;
  const double ratio2 = e2 / e3;
  if (ratio1 < 3.0 || ratio1 > 5.5 || ratio2 < 3.0 || ratio2 > 5.5) {
    std::cerr << "ratios " << ratio1 << ", " << ratio2 << "\n";
    return fail("BDF2-DAE observed order");
  }

  // diffusion = 0 时各单元完全解耦，每个单元都应等于单单元闭式解；
  // 若某个单元偏差超限，说明链式装配错误地引入了跨单元耦合。
  dfn::ChainDaeOptions chain_options;
  chain_options.cells = 4;
  chain_options.lambda = 100.0;
  chain_options.diffusion = 0.0;
  const dfn::ChainDae chain(chain_options);
  const dfn::Bdf2DaeResult chain_result = dfn::bdf2_dae_integrate(
      chain, dfn::chain_dae_consistent_state(chain, 0.0), 0.0, t_end, 0.005, bdf);
  if (!chain_result.ok) {
    return fail(chain_result.message.c_str());
  }
  for (int i = 0; i < chain.cells(); ++i) {
    if (!close(chain_result.y.back()[static_cast<std::size_t>(2 * i)], exact, 2e-4)) {
      return fail("decoupled chain cells disagree with closed form");
    }
  }

  std::cout << "dae ok (BDF2-DAE order ratios " << ratio1 << ", " << ratio2
            << ")\n";
  return 0;
}

}  // namespace

int main() {
  if (const int rc = test_pattern(); rc != 0) {
    return rc;
  }
  if (const int rc = test_chain_model(); rc != 0) {
    return rc;
  }
  return test_single_cell_solution();
}
