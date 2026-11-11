#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

// 实验（第 6 章 Newton 法）：用 Newton 迭代求 f(x) = x^2 - 2 = 0 的正根 √2。
// 迭代式 x_{k+1} = x_k - f(x_k) / f'(x_k)，其中 f'(x) = 2x。
// 输出 CSV 列：iter 迭代序号, x 当前近似, f(x) 残差；收敛时有效数字每步翻倍。
int main(int argc, char** argv) {
  // 默认初值 x0 = 1；也可用命令行参数实验不同初值（x0 = 0 处导数为零会失败）。
  double x = 1.0;
  if (argc > 1) {
    x = std::atof(argv[1]);
  }

  std::cout << "iter,x,f(x)\n" << std::setprecision(10);
  // 最多迭代 10 次；二次收敛，通常 6 次迭代内达到双精度，先打印再判收敛。
  for (int iteration = 0; iteration < 10; ++iteration) {
    // 残差 f(x) = x^2 - 2，零点即 √2。
    const double f = x * x - 2.0;
    std::cout << iteration << ',' << x << ',' << f << '\n';
    // |f| < 1e-14 视为已达机器精度，提前结束。
    if (std::abs(f) < 1e-14) {
      break;
    }
    // Newton 校正步：x <- x - f / (2x)，其中 2x 是 f 在 x 处的导数。
    x -= f / (2.0 * x);
  }
  return 0;
}
