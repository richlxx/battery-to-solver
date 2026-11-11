#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// 实验（第 7 章 刚性方程）：对 y' = -100 y, y(0) = 1 在 [0, 0.2] 上对比显式与
// 隐式欧拉。显式放大因子 1 + hλ 在 h > 0.02 s 时模大于 1 而振荡发散，隐式因子
// 1/(1 - hλ) 恒有界。CSV 列：h, steps, 显式/隐式终点值, 精确值, 显式最大模。
int main() {
  // 刚性特征值 λ = -100 1/s（对应时间常数 0.01 s），积分到 t_end = 0.2 s。
  const double lambda = -100.0;
  const double t_end = 0.2;
  // 精确终点值 exp(-20) ≈ 2.06e-9，数值上已接近 0。
  const double exact = std::exp(lambda * t_end);

  std::cout << "h,steps,explicit_final,implicit_final,exact_final,explicit_max\n"
            << std::scientific << std::setprecision(6);

  for (double h : {0.1, 0.05, 0.02, 0.01, 0.005}) {
    // 步数取整，使终点恰好落在 t_end 上。
    const int steps = static_cast<int>(
        std::lround(t_end / h));

    double explicit_value = 1.0;
    // 显式解全过程 |y| 的最大值，用于判断是否发生数值放大。
    double explicit_max = 1.0;
    double implicit_value = 1.0;
    for (int i = 0; i < steps; ++i) {
      // 显式欧拉：y_{n+1} = (1 + λh) y_n，|1 + λh| 即每步放大因子。
      explicit_value *= 1.0 + lambda * h;
      explicit_max = std::max(explicit_max, std::abs(explicit_value));
      // 隐式欧拉：y_{n+1} = y_n / (1 - λh)，对任意 h > 0 都有 1/(1 - λh) < 1。
      implicit_value /= 1.0 - lambda * h;
    }

    std::cout << h << ',' << steps << ',' << explicit_value << ','
              << implicit_value << ',' << exact << ',' << explicit_max << '\n';
  }

  return 0;
}
