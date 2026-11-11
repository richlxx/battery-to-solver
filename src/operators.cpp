// 本文件实现一维有限体积离散用到的三对角算子：三对角矩阵与向量相乘，
// 以及扩散算子 div(D grad c) 的装配（矩阵耦合 + 边界通量源项）。

#include "dfn/operators.hpp"

#include <stdexcept>

namespace dfn {

// 计算 y = A x，A 为本对象描述的三对角矩阵。
// 参数 x：长度必须为 n；返回长度 n 的向量 y。
// 约定：lower[0] 与 upper[n-1] 不参与运算（首行无左邻、末行无右邻）。
std::vector<double> TridiagonalMatrix::multiply(
    const std::vector<double>& x) const {
  if (static_cast<int>(x.size()) != n) {
    throw std::invalid_argument("tridiagonal multiply: size mismatch");
  }

  // 第 i 行：lower[i]*x[i-1] + diag[i]*x[i] + upper[i]*x[i+1]，
  // 越界的相邻项自动跳过。
  std::vector<double> y(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    double value = diag[i] * x[i];
    if (i > 0) {
      value += lower[i] * x[i - 1];
    }
    if (i + 1 < n) {
      value += upper[i] * x[i + 1];
    }
    y[i] = value;
  }
  return y;
}

// 在 Grid1D 上装配 div(D grad c) 的有限体积离散。
// 参数：grid 网格；D_cell 每个控制体的扩散系数（m^2/s；电导率 S/m 的方程
//       形式相同，只是物理含义不同）；left_flux / right_flux 为左、右边界面上
//       给定的通量密度，正方向为流出计算域。
// 返回：DiffusionOperator。matrix 为对称负定的三对角矩阵（对角元为负），
//       source 为已除以控制体体积的边界通量项，满足 matrix * c = source。
DiffusionOperator assemble_diffusion(const Grid1D& grid,
                                     const std::vector<double>& D_cell,
                                     double left_flux,
                                     double right_flux) {
  if (static_cast<int>(D_cell.size()) != grid.n) {
    throw std::invalid_argument("diffusion operator: D_cell size mismatch");
  }

  const int n = grid.n;
  DiffusionOperator op;
  op.matrix.n = n;
  op.matrix.lower.assign(static_cast<std::size_t>(n), 0.0);
  op.matrix.diag.assign(static_cast<std::size_t>(n), 0.0);
  op.matrix.upper.assign(static_cast<std::size_t>(n), 0.0);
  op.source.assign(static_cast<std::size_t>(n), 0.0);

  // 界面扩散系数：内部界面取相邻两控制体值的调和平均 2ab/(a+b)，
  // 等价于两个半格串联，可保证 D 跳变处通量连续；
  // 边界界面直接取相邻控制体的值（边界面不参与内部差分）。
  std::vector<double> D_face(static_cast<std::size_t>(n + 1), 0.0);
  D_face[0] = D_cell[0];
  D_face[static_cast<std::size_t>(n)] = D_cell[static_cast<std::size_t>(n - 1)];
  for (int i = 1; i < n; ++i) {
    const double a = D_cell[static_cast<std::size_t>(i - 1)];
    const double b = D_cell[static_cast<std::size_t>(i)];
    D_face[static_cast<std::size_t>(i)] = 2.0 * a * b / (a + b);
  }

  // 对控制体 i 做有限体积积分：每个界面的传导系数为
  // D_face * area / spacing / volume，即（面积 × 扩散系数）/ 中心距 / 体积。
  // west、east 存放非对角元的幅值（正数），对角元取其负和，
  // 使 matrix 对应 -div(D grad c) 这一负定形式。
  for (int i = 0; i < n; ++i) {
    double west = 0.0;
    if (i > 0) {
      west = D_face[static_cast<std::size_t>(i)] * grid.area[static_cast<std::size_t>(i)] /
             grid.spacing[static_cast<std::size_t>(i)] /
             grid.volume[static_cast<std::size_t>(i)];
    }
    double east = 0.0;
    if (i + 1 < n) {
      east = D_face[static_cast<std::size_t>(i + 1)] *
             grid.area[static_cast<std::size_t>(i + 1)] /
             grid.spacing[static_cast<std::size_t>(i + 1)] /
             grid.volume[static_cast<std::size_t>(i)];
    }
    op.matrix.diag[static_cast<std::size_t>(i)] = -(west + east);
    op.matrix.lower[static_cast<std::size_t>(i)] = west;
    op.matrix.upper[static_cast<std::size_t>(i)] = east;
  }

  // 边界面上的通量是已知量（Neumann 型条件），不引入新未知量，
  // 因此移入右端项：通量以流出为正，故按 -area*flux/volume 写入 source。
  op.source[0] -=
      grid.area[0] * left_flux / grid.volume[0];
  op.source[static_cast<std::size_t>(n - 1)] -=
      grid.area[static_cast<std::size_t>(n)] * right_flux /
      grid.volume[static_cast<std::size_t>(n - 1)];

  return op;
}

}  // namespace dfn
