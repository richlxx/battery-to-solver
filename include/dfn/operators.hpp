#pragma once

#include <functional>
#include <vector>

#include "dfn/grid.hpp"

namespace dfn {

// 三对角矩阵容器、基于 Grid1D 的扩散算子 div(D * grad(c)) 装配，以及
// 配套的线性 BDF2 推进，是空间离散（有限体积与有限元）与三对角线性
// 求解之间的接口层。

// 三对角矩阵，第 i 行表示
//   lower[i] * x[i - 1] + diag[i] * x[i] + upper[i] * x[i + 1]。
// 三个数组长度均为 n，首尾各有一个未使用的元素，便于统一循环边界。
struct TridiagonalMatrix {
  int n = 0;                  // 阶数（= 未知量个数）
  std::vector<double> lower;  // 次对角元，lower[0] 未使用
  std::vector<double> diag;   // 主对角元
  std::vector<double> upper;  // 超对角元，upper[n - 1] 未使用

  // 计算矩阵与向量 x 的乘积，x 长度应为 n，返回长度 n 的 Ax。
  std::vector<double> multiply(const std::vector<double>& x) const;
};

// Grid1D 上 div(D * grad(c)) 的离散结果。matrix 存放单元之间的耦合系数，
// source 存放边界通量除以单元体积；装配后的方程为 matrix * c = source。
struct DiffusionOperator {
  TridiagonalMatrix matrix;
  std::vector<double> source;  // 单位体积的边界源项
};

// D_cell 给出每个单元的扩散系数，界面值在函数内部重构（相邻单元取
// 调和平均，保证 D 跃变处通量连续）。left_flux 与 right_flux 是穿过
// 首、末界面的给定通量，正方向取外法向，它们进入源项。
DiffusionOperator assemble_diffusion(const Grid1D& grid,
                                     const std::vector<double>& D_cell,
                                     double left_flux = 0.0,
                                     double right_flux = 0.0);

// n 阶单位阵。有限体积的未知量已是“单位体积的变化率”，把半离散系统
// 统一写成 X dc/dt = Y c + f 时，它的 X 退化为单位阵。
TridiagonalMatrix identity_matrix(int n);

// 取负副本：把 assemble_diffusion 给出的 dc/dt = A c + source 统一写成
// X dc/dt = Y c + f 时，Y = -A。
TridiagonalMatrix negated_matrix(const TridiagonalMatrix& a);

// 线性半离散系统 X dc/dt = Y c + f 的 BDF2 推进（首步隐式 Euler 起步，
// 压掉初值激发的高频分量）：首步 (X + dt·Y) c¹ = X c⁰ + dt·f，之后
// (1.5X + dt·Y) cⁿ⁺¹ = X(2cⁿ − 0.5cⁿ⁻¹) + dt·f。X = I 时退化为有限
// 体积的隐式推进；X 为（行缩放的）质量阵时即有限元的推进。线性问题、
// 等步长，不需要 Newton 迭代。
// record(step, c) 在初始时刻（step = 0）与每个步末被调用，可缺省。
std::vector<double> march_tridiagonal(
    const TridiagonalMatrix& x_matrix, const TridiagonalMatrix& y_matrix,
    const std::vector<double>& source, std::vector<double> c, double dt,
    int steps,
    const std::function<void(int, const std::vector<double>&)>& record = {});

}  // namespace dfn
