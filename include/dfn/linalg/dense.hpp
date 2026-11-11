#pragma once

#include <vector>

namespace dfn {

// 稠密线性代数：带部分主元的 LU 分解与直接求解。用于教学规模的小系统，
// 例如 Levenberg-Marquardt 的正规方程和稀疏 Jacobian 的正确性核对。

// LU 分解结果。lu 保存原地分解后的 L 与 U（L 的对角元隐含为 1），
// perm 记录行交换，满足 P * A = L * U。
struct LuFactorization {
  int n = 0;                            // 矩阵阶数
  std::vector<std::vector<double>> lu;  // 分解后的矩阵，按行存放
  std::vector<int> perm;  // 行置换：第 i 行对应原矩阵的第 perm[i] 行
};

// 对 A 做带部分主元的 LU 分解（A 按值传入，内部被覆盖）。
// 空矩阵、非方阵或数值奇异时抛异常。
LuFactorization lu_factor(std::vector<std::vector<double>> A);

// 用已分解结果求解 A x = b，b 长度应为 n，返回长度 n 的解向量。
std::vector<double> lu_solve(const LuFactorization& factorization,
                             const std::vector<double>& b);

// 一次性完成分解与回代，适合小规模稠密问题。
std::vector<double> solve_dense(std::vector<std::vector<double>> A,
                                std::vector<double> b);

}  // namespace dfn
