#include <cmath>
#include <iomanip>
#include <iostream>

// 实验（第 2 章 数值微分）：用前向差分与中心差分近似 sin'(1) = cos(1)，
// 步长 h 从 0.1 起连续减半 48 次，观察截断误差与舍入误差的权衡。
// 输出 CSV 列：h, forward, central（两种近似值）, forward_err, central_err。
int main() {
  // 求导点 x = 1 与解析导数值 cos(1)，误差以此为基准（无量纲）。
  const double x = 1.0;
  const double exact = std::cos(x);

  std::cout << "h,forward,central,forward_err,central_err\n"
            << std::scientific << std::setprecision(6);

  // h_k = 0.1 * 2^{-k}，k 到 47 时 h ≈ 7e-16，已接近双精度机器精度。
  for (int k = 0; k < 48; ++k) {
    const double h = 0.1 * std::pow(0.5, k);
    // 前向差分：(sin(x+h) - sin(x)) / h，截断误差 O(h)，一阶精度。
    const double forward = (std::sin(x + h) - std::sin(x)) / h;
    // 中心差分：(sin(x+h) - sin(x-h)) / (2h)，截断误差 O(h^2)，二阶精度。
    const double central = (std::sin(x + h) - std::sin(x - h)) / (2.0 * h);
    // 两个绝对误差：h 过大时截断误差占优，h 过小时舍入误差 ~ eps/h 占优。
    std::cout << h << ',' << forward << ',' << central << ','
              << std::abs(forward - exact) << ','
              << std::abs(central - exact) << '\n';
  }

  return 0;
}
