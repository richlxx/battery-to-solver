// 测试文件：定步长 BDF2 常微分方程积分器（对应 solvers/bdf2 模块）。
// 验证二阶收敛阶、步数/终点时间，以及刚性衰减问题的数值稳定性。

#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/solvers/bdf2.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

}  // namespace

int main() {
  // 标量测试问题 y' = -y，y(0) = 1，解析解 y(t) = exp(-t)。
  const auto rhs = [](double, const std::vector<double>& y) {
    return std::vector<double>{-y[0]};
  };
  const auto jacobian = [](double, const std::vector<double>&) {
    return std::vector<std::vector<double>>{{-1.0}};
  };

  const auto run = [&](double h) {
    return dfn::bdf2_integrate(rhs, jacobian, {1.0}, 0.0, 1.0, h);
  };

  const double exact = std::exp(-1.0);  // t = 1 时的精确解，约 0.367879
  // 步长逐次减半，用于估计观测收敛阶。
  const dfn::Bdf2Result coarse = run(0.05);
  const dfn::Bdf2Result medium = run(0.025);
  const dfn::Bdf2Result fine = run(0.0125);

  if (!coarse.ok || !medium.ok || !fine.ok) {
    return fail("BDF2 integration failed");
  }
  // 定步长积分：[0,1] 上 h = 0.05/0.025/0.0125 应分别走 20/40/80 步。
  if (coarse.steps != 20 || medium.steps != 40 || fine.steps != 80) {
    return fail("unexpected step counts");
  }
  // 终点时间由 t0 + steps*h 得到，1e-12 只允许浮点累加误差。
  if (std::abs(coarse.t.back() - 1.0) > 1e-12) {
    return fail("final time");
  }

  const double e1 = std::abs(coarse.y.back()[0] - exact);
  const double e2 = std::abs(medium.y.back()[0] - exact);
  const double e3 = std::abs(fine.y.back()[0] - exact);
  // 二阶方法的误差应随 h 减半缩小约 4 倍；放宽到 [3.0, 5.5] 是为容纳
  // 启动步（第一步只能用一阶方法）和渐近区间的偏差。
  const double ratio1 = e1 / e2;
  const double ratio2 = e2 / e3;
  if (ratio1 < 3.0 || ratio1 > 5.5 || ratio2 < 3.0 || ratio2 > 5.5) {
    return fail("BDF2 observed order");
  }

  // 刚性衰减 y' = -100 y：h*lambda = -5，显式方法会失稳，
  // 而 BDF2 是 A-稳定的，解必须有限且幅值不超过初始值 1。
  const auto stiff_rhs = [](double, const std::vector<double>& y) {
    return std::vector<double>{-100.0 * y[0]};
  };
  const auto stiff_jacobian = [](double, const std::vector<double>&) {
    return std::vector<std::vector<double>>{{-100.0}};
  };
  const dfn::Bdf2Result stiff = dfn::bdf2_integrate(
      stiff_rhs, stiff_jacobian, {1.0}, 0.0, 0.2, 0.05);
  if (!stiff.ok || !std::isfinite(stiff.y.back()[0]) ||
      std::abs(stiff.y.back()[0]) > 1.0) {
    return fail("BDF2 stiff stability");
  }

  std::cout << "bdf2 ok (order ratios " << ratio1 << ", " << ratio2 << ")\n";
  return 0;
}
