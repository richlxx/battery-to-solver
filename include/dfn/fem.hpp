#pragma once

#include <vector>

#include "dfn/operators.hpp"

namespace dfn {

// 球坐标径向扩散的 P1（线性帽函数）有限元装配，与 operators.hpp 的
// 有限体积装配构成同一问题的两条离散路线。未知量在节点上（含球心与
// 表面），配合 operators.hpp 的 march_tridiagonal 即可推进。

// 装配结果：节点半径与两个对称三对角矩阵。r² 权重来自球坐标体积元，
// 使 4π · (1ᵀ M c) 恰为颗粒内总量。
struct FemSphere {
  int nodes = 0;  // 节点数 = 单元数 + 1，未知量在节点上
  std::vector<double> r;  // 节点半径：r[0] = 0，r[nodes-1] = radius

  // 质量阵 M_ij = ∫ r² φ_i φ_j dr [m³]，时变项的权重。
  TridiagonalMatrix mass;

  // 刚度阵 K_ij = D ∫ r² φ'_i φ'_j dr [m³/s]，扩散的弱形式。
  TridiagonalMatrix stiffness;
};

// 在 [0, radius] 上用 elements 个等距线性单元装配。基函数为分段线性
// 帽函数，单元积分全部用解析公式（φ' 在单元内是常数，∫r² dr 用立方差），
// 不引入数值积分误差。
FemSphere assemble_fem_sphere(int elements, double radius, double diffusivity);

// 颗粒内总量 4π ∫ c r² dr ≈ 4π · (1ᵀ M c)。取 v = 1（全部基函数之和，
// 它仍在有限元空间里）时，弱形式就是这份总量的收支方程，因此推进中
// 它随时间的变化应精确等于边界流出——全局守恒的代数根源。
double fem_sphere_inventory(const FemSphere& fem, const std::vector<double>& c);

// 行缩放后的推进系统 X dc/dt = Y c + f（配合 march_tridiagonal）：
// 把质量阵、刚度阵、右端的每一行除以质量阵行和（节点体积权重
// ∫ r² φ_i dr）。行缩放不改变解，只是把 m³ 量级的矩阵拉回 1 附近，
// 避开 Thomas 算法的主元下限；assemble_diffusion 除以控制体体积
// 干的就是同一个动作，它的质量阵因此退化成单位阵。
struct FemMarchSystem {
  TridiagonalMatrix x_matrix;  // 行缩放质量阵，每行行和为 1
  TridiagonalMatrix y_matrix;  // 行缩放刚度阵，量级 ~ 1/s
  std::vector<double> source;  // 行缩放右端，量级 ~ mol/(m³ s)
};

// 组装行缩放系统；surface_flux 为外表面通量 [mol/(m² s)]，正方向朝外。
// 给定通量（Neumann）条件经弱形式的分部积分自动出现在右端，只落在
// 末节点；球心处 r² = 0，对称边界同样不需要任何处理。
FemMarchSystem make_fem_march_system(const FemSphere& fem, double surface_flux);

}  // namespace dfn
