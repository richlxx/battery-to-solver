#pragma once

#include <vector>

#include "dfn/grid.hpp"

namespace dfn {

// 三对角矩阵容器与基于 Grid1D 的扩散算子 div(D * grad(c)) 装配，
// 是有限体积离散与三对角线性求解之间的接口层。

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

}  // namespace dfn
