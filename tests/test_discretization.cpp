// 测试文件：有限体积网格与扩散算子（对应 grid / operators 模块）。
// 验证球/线段网格的几何量、均匀态的零通量、离散守恒性与表面通量符号。

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "dfn/grid.hpp"
#include "dfn/operators.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

// 纯相对容差：几何量都是正数且量级相近，用相对误差可直接对应离散精度；
// 两者同时为 0 时视为一致，避免除零。
bool close_relative(double a, double b, double tolerance) {
  const double scale = std::max(std::abs(a), std::abs(b));
  if (scale == 0.0) {
    return true;
  }
  return std::abs(a - b) <= tolerance * scale;
}

}  // namespace

int main() {
  constexpr double pi = 3.14159265358979323846;
  const double radius = 5.0e-6;  // 5 µm，典型活性颗粒半径
  const int n = 8;

  // 球壳网格的体积之和必须精确等于整球体积；1e-12 只容纳舍入误差。
  const dfn::Grid1D sphere = dfn::make_uniform_sphere_grid(n, radius);
  double volume = 0.0;
  for (double v : sphere.volume) {
    volume += v;
  }
  const double exact_volume = 4.0 / 3.0 * pi * radius * radius * radius;
  if (!close_relative(volume, exact_volume, 1e-12)) {
    return fail("sphere volume");
  }
  // 球心处面积为 0，保证径向通量在原点不自发产生源项。
  if (std::abs(sphere.area[0]) > 0.0) {
    return fail("center face area must be zero");
  }
  // 外表面面积应为 4*pi*r^2，它是颗粒与电解液交换锂离子的界面。
  if (!close_relative(sphere.area[static_cast<std::size_t>(n)],
                      4.0 * pi * radius * radius, 1e-12)) {
    return fail("sphere surface area");
  }

  // 线段网格：长 100 µm、截面积 1e-4 m^2，总体积应为 1e-8 m^3。
  const dfn::Grid1D line = dfn::make_uniform_line_grid(10, 100.0e-6, 1.0e-4);
  double line_volume = 0.0;
  for (double v : line.volume) {
    line_volume += v;
  }
  if (!close_relative(line_volume, 1.0e-8, 1e-12)) {
    return fail("line volume");
  }

  // 固相扩散系数取 3e-14 m^2/s（典型石墨量级），用于装配隐式扩散矩阵。
  const std::vector<double> D_sphere(static_cast<std::size_t>(sphere.n), 3.0e-14);
  const dfn::DiffusionOperator op = dfn::assemble_diffusion(sphere, D_sphere);

  // 均匀浓度场没有梯度，任何界面的扩散通量都应为 0；
  // 1e-12 的绝对阈值对应矩阵行和中的浮点消去误差。
  const std::vector<double> uniform(static_cast<std::size_t>(sphere.n), 1.0);
  const std::vector<double> uniform_rhs = op.matrix.multiply(uniform);
  for (double value : uniform_rhs) {
    if (std::abs(value) > 1e-12) {
      return fail("uniform state must be steady with zero flux");
    }
  }

  // 固定随机种子 42，生成 [0.5, 1.5] 的非均匀浓度场，使所有界面通量非零；
  // 固定种子保证失败可复现。
  std::mt19937 generator(42);
  std::uniform_real_distribution<double> distribution(0.5, 1.5);
  std::vector<double> c(static_cast<std::size_t>(sphere.n));
  for (double& value : c) {
    value = distribution(generator);
  }

  // 守恒性：把每个单元的体积乘上净通量并求和，内部界面通量两两抵消，
  // 只剩边界项（此处为 0）。阈值按通量绝对量级缩放，避免只看绝对误差。
  const std::vector<double> rhs = op.matrix.multiply(c);
  double balance = 0.0;
  double scale = 0.0;
  for (int i = 0; i < sphere.n; ++i) {
    balance += sphere.volume[static_cast<std::size_t>(i)] *
               (rhs[static_cast<std::size_t>(i)] +
                op.source[static_cast<std::size_t>(i)]);
    scale += std::abs(sphere.volume[static_cast<std::size_t>(i)] *
                      rhs[static_cast<std::size_t>(i)]);
  }
  if (std::abs(balance) > 1e-10 * scale) {
    return fail("sphere conservation");
  }

  // 给定外表面通量 1e-6（外法向为正），体积加权总变化应等于
  // -A_surface * flux：正通量为流出，颗粒内锂总量相应减少。
  const double surface_flux = 1.0e-6;
  const dfn::DiffusionOperator with_flux =
      dfn::assemble_diffusion(sphere, D_sphere, 0.0, surface_flux);
  const std::vector<double> flux_rhs = with_flux.matrix.multiply(c);
  double flux_balance = 0.0;
  for (int i = 0; i < sphere.n; ++i) {
    flux_balance += sphere.volume[static_cast<std::size_t>(i)] *
                    (flux_rhs[static_cast<std::size_t>(i)] +
                     with_flux.source[static_cast<std::size_t>(i)]);
  }
  const double expected = -sphere.area[static_cast<std::size_t>(n)] * surface_flux;
  if (!close_relative(flux_balance, expected, 1e-10)) {
    return fail("surface flux balance");  // 符号或面积归一化错误都会在此暴露
  }

  // 线段网格同样检查守恒：两端给定 0 通量，所有内部通量必须抵消。
  const std::vector<double> D_line(static_cast<std::size_t>(line.n), 1.0e-9);
  const dfn::DiffusionOperator line_op = dfn::assemble_diffusion(line, D_line);
  std::vector<double> line_c(static_cast<std::size_t>(line.n));
  for (double& value : line_c) {
    value = distribution(generator);
  }
  const std::vector<double> line_rhs = line_op.matrix.multiply(line_c);
  double line_balance = 0.0;
  double line_scale = 0.0;
  for (int i = 0; i < line.n; ++i) {
    line_balance += line.volume[static_cast<std::size_t>(i)] *
                    line_rhs[static_cast<std::size_t>(i)];
    line_scale += std::abs(line.volume[static_cast<std::size_t>(i)] *
                           line_rhs[static_cast<std::size_t>(i)]);
  }
  if (std::abs(line_balance) > 1e-10 * line_scale) {
    return fail("line conservation");
  }

  std::cout << "discretization ok (sphere n=" << sphere.n
            << ", line n=" << line.n << ")\n";
  return 0;
}
