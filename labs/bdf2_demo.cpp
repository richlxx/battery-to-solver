#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/solvers/bdf2.hpp"

// 实验（第 7 章 时间推进）：用 BDF2 求解标量衰减方程 y' = -y, y(0) = 1，
// 在 t = 1 处与解析解 exp(-1) 比较，检验二阶方法的收敛阶。
// 输出 CSV 列：h 步长(s), steps 步数, error 终点绝对误差, ratio 误差比(≈4)。
int main() {
  // 右端 f(t, y) = -y：线性问题，解析解为指数衰减，便于精确度量误差。
  const auto rhs = [](double, const std::vector<double>& y) {
    return std::vector<double>{-y[0]};
  };
  // 解析 Jacobian J = df/dy = -1，供每个 BDF2 时间步内的 Newton 迭代使用。
  const auto jacobian = [](double, const std::vector<double>&) {
    return std::vector<std::vector<double>>{{-1.0}};
  };

  // t = 1 处的精确解 e^{-1}（无量纲）。
  const double exact = std::exp(-1.0);
  // 表头；第一个步长没有前一次误差，ratio 记为 0。
  std::cout << "h,steps,error,ratio\n" << std::scientific << std::setprecision(4);

  double previous_error = 0.0;
  // 步长逐次减半：二阶方法在 h 减半时误差应缩小到约 1/4。
  for (double h : {0.2, 0.1, 0.05, 0.025}) {
    // 从 t = 0 积分到 t = 1，初值 y(0) = 1。
    const dfn::Bdf2Result result =
        dfn::bdf2_integrate(rhs, jacobian, {1.0}, 0.0, 1.0, h);
    // 终点误差：数值解 y(1) 与解析解之差的绝对值。
    const double error = std::abs(result.y.back()[0] - exact);
    // 相邻两次误差之比；理想二阶方法取 h 减半时该比值约为 4。
    const double ratio = previous_error > 0.0 ? previous_error / error : 0.0;
    std::cout << h << ',' << result.steps << ',' << error << ',' << ratio << '\n';
    previous_error = error;
  }

  return 0;
}
