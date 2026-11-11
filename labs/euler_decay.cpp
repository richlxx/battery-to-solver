#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// 实验（第 1 章 稳定性 / 第 7 章 时间推进）：显式欧拉法求解 y' = λ y,
// y(0) = 1，观察衰减速率与稳定条件 |1 + hλ| ≤ 1。分三段：单模式衰减表、
// 步长扫描、快慢两个模式对比，均输出 CSV（t,euler,exact 等）。
namespace {

// 显式欧拉推进 steps 步：y_{n+1} = y_n + h λ y_n = (1 + hλ) y_n。
// λ 的单位为 1/s，h 的单位为 s，返回值含 steps + 1 个时间点，初值 y(0) = 1。
std::vector<double> euler(double lambda, double h, int steps) {
  std::vector<double> y(steps + 1, 1.0);
  for (int i = 0; i < steps; ++i) {
    y[i + 1] = y[i] + h * lambda * y[i];
  }
  return y;
}

// 解向量中绝对值最大的分量，用于判断显式欧拉是否把数值解放大（不稳定）。
double max_abs(const std::vector<double>& y) {
  double m = 0.0;
  for (double v : y) {
    m = std::max(m, std::abs(v));
  }
  return m;
}

}  // namespace

int main() {
  std::cout << std::scientific << std::setprecision(6);

  // 第一段：慢模式 λ = -1，h = 0.1 s，放大因子 1 + hλ = 0.9（|·| < 1，稳定），
  // 5 步后与解析解 exp(-t) 逐点对照，体现一阶精度与单调衰减。
  std::cout << "decay: y'=-1*y, y(0)=1, h=0.1\n";
  std::cout << "t,euler,exact\n";
  const auto decay = euler(-1.0, 0.1, 5);
  for (int i = 0; i < static_cast<int>(decay.size()); ++i) {
    const double t = 0.1 * i;
    std::cout << t << ',' << decay[i] << ',' << std::exp(-t) << '\n';
  }

  // 第二段：刚性模式 λ = -100，t_end = 0.2 s，精确终值 exp(-20) ≈ 2e-9。
  // 稳定要求 h ≤ 2/|λ| = 0.02 s：h = 0.05 时 |1 + hλ| = 4 > 1 会发散，
  // max_abs 记录全过程最大模，> 1 即说明数值解被人为放大。
  std::cout << "stability: y'=-100*y, t_end=0.2, exact="
            << std::exp(-100.0 * 0.2) << '\n';
  std::cout << "h,steps,max_abs,end\n";
  for (double h : {0.05, 0.02, 0.01, 0.005}) {
    // 步数取整，保证终点落在 t_end 上。
    const int steps = static_cast<int>(std::lround(0.2 / h));
    const auto y = euler(-100.0, h, steps);
    std::cout << h << ',' << steps << ',' << max_abs(y) << ',' << y.back()
              << '\n';
  }

  // 第三段：同一时间区间的两个模式 λ = (-1, -100)，各自独立积分。
  // 步长必须同时满足两个模式的稳定条件，实际受快模式 λ = -100 支配。
  std::cout << "two modes: lambda=(-1,-100), y(0)=(1,1), t_end=0.2\n";
  std::cout << "h,end_slow,end_fast\n";
  for (double h : {0.05, 0.02, 0.01}) {
    const int steps = static_cast<int>(std::lround(0.2 / h));
    const auto slow = euler(-1.0, h, steps);
    const auto fast = euler(-100.0, h, steps);
    std::cout << h << ',' << slow.back() << ',' << fast.back() << '\n';
  }

  return 0;
}
