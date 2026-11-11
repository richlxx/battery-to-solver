// 测试文件：Newton 迭代求解器（对应 solvers/newton 模块）。
// 验证耦合非线性系统收敛、变量缩放、残差历史与奇异 Jacobian 的错误处理。

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dfn/solvers/newton.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

// 混合绝对/相对容差 1e-8，适应 O(1) 到 O(1e4) 的变量量级。
bool close(double a, double b, double tolerance = 1e-8) {
  return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

}  // namespace

int main() {
  // 非线性方程组：x0^2 + x1 = 3，x0 + x1^2 = 5，精确解为 (1, 2)。
  const dfn::ResidualFunction residual =
      [](const std::vector<double>& x) {
        return std::vector<double>{x[0] * x[0] + x[1] - 3.0,
                                   x[0] + x[1] * x[1] - 5.0};
      };
  // 解析 Jacobian：[[2*x0, 1], [1, 2*x1]]，用于检验牛顿方向是否正确。
  const dfn::JacobianFunction jacobian =
      [](const std::vector<double>& x) {
        return std::vector<std::vector<double>>{{2.0 * x[0], 1.0},
                                                {1.0, 2.0 * x[1]}};
      };

  // 从 (0.5, 1.5) 出发，要求残差降到 1e-12；牛顿法二次收敛，迭代数应远小于 30。
  dfn::NewtonOptions options;
  options.tolerance = 1e-12;
  options.max_iterations = 30;
  const dfn::NewtonResult result =
      dfn::newton_solve(residual, jacobian, {0.5, 1.5}, options);

  if (!result.converged) {
    return fail("coupled system did not converge");
  }
  if (!close(result.x[0], 1.0) || !close(result.x[1], 2.0)) {
    return fail("coupled system solution");
  }
  // 收敛标志还须与残差历史一致：最后一次残差不得超过设定容差。
  if (result.residual_history.empty() ||
      result.residual_history.back() > 1e-12) {
    return fail("residual history");
  }
  // 二次收敛通常 5-8 步完成，超过 15 步说明线搜索或 Jacobian 有问题。
  if (result.iterations > 15) {
    return fail("Newton took too many iterations");
  }

  // 病态缩放系统：x0 的量级为 1e4、x1 为 1，若不做变量缩放，
  // Jacobian 两列范数相差 4 个数量级，消元容易损失精度。
  // 精确解：x1 = 1，x0 = 1e4。
  const dfn::ResidualFunction scaled_residual =
      [](const std::vector<double>& x) {
        return std::vector<double>{x[0] / 1.0e4 + x[1] - 2.0,
                                   x[0] / 1.0e4 - x[1]};
      };
  const dfn::JacobianFunction scaled_jacobian =
      [](const std::vector<double>&) {
        return std::vector<std::vector<double>>{{1.0e-4, 1.0},
                                                {1.0e-4, -1.0}};
      };
  // variable_scale 按变量量级归一化步长，使两个分量在求解中同等重要。
  dfn::NewtonOptions scaled_options;
  scaled_options.tolerance = 1e-12;
  scaled_options.variable_scale = {1.0e4, 1.0};
  const dfn::NewtonResult scaled_result = dfn::newton_solve(
      scaled_residual, scaled_jacobian, {0.5e4, 0.5}, scaled_options);
  if (!scaled_result.converged ||
      !close(scaled_result.x[0], 1.0e4, 1e-8) ||
      !close(scaled_result.x[1], 1.0, 1e-8)) {
    return fail("scaled system");
  }

  // Jacobian 恒为 0（秩亏），线性子问题无解，必须抛出异常而非继续迭代。
  bool singular_threw = false;
  try {
    (void)dfn::newton_solve(
        [](const std::vector<double>& x) {
          return std::vector<double>{x[0] - 1.0};
        },
        [](const std::vector<double>&) {
          return std::vector<std::vector<double>>{{0.0}};
        },
        {0.0});
  } catch (const std::exception&) {
    singular_threw = true;
  }
  if (!singular_threw) {
    return fail("singular Jacobian must throw");
  }

  std::cout << "newton ok (" << result.iterations << " iterations)\n";
  return 0;
}
