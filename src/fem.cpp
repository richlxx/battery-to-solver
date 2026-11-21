// 本文件实现球坐标径向扩散的 P1 有限元装配：质量阵与刚度阵（解析单元
// 积分）、总量泛函与行缩放后的推进系统，与 operators.cpp 的有限体积
// 装配平行。

#include "dfn/fem.hpp"

#include <cmath>
#include <stdexcept>

namespace dfn {

namespace {
constexpr double kPi = 3.14159265358979323846;

TridiagonalMatrix zero_matrix(int n) {
  TridiagonalMatrix m;
  m.n = n;
  m.lower.assign(static_cast<std::size_t>(n), 0.0);
  m.diag.assign(static_cast<std::size_t>(n), 0.0);
  m.upper.assign(static_cast<std::size_t>(n), 0.0);
  return m;
}
}  // namespace

// 装配：遍历每个单元，把 2×2 的单元质量阵与单元刚度阵累加进全局三对角
// 矩阵。相邻基函数的支集只覆盖一个公共单元，因此全局矩阵是三对角的，
// 与有限体积的装配结果同构，可复用 Thomas 求解。
FemSphere assemble_fem_sphere(int elements, double radius, double diffusivity) {
  if (elements < 1) {
    throw std::invalid_argument("fem sphere: need at least one element");
  }

  FemSphere fem;
  fem.nodes = elements + 1;
  const double h = radius / elements;
  fem.r.resize(static_cast<std::size_t>(fem.nodes));
  for (int k = 0; k < fem.nodes; ++k) {
    fem.r[static_cast<std::size_t>(k)] = k * h;
  }
  fem.mass = zero_matrix(fem.nodes);
  fem.stiffness = zero_matrix(fem.nodes);

  for (int e = 0; e < elements; ++e) {
    const double rk = fem.r[static_cast<std::size_t>(e)];
    const double rk1 = fem.r[static_cast<std::size_t>(e + 1)];

    // 单元刚度：φ' 在单元内是常数 ±1/h，因此 ∫ r² φ'_i φ'_j dr 逐个
    // 元素都是 ±D(rk1³ - rk³)/(3h²)。局部 2×2 块为 [[w,-w],[-w,w]]。
    const double w =
        diffusivity * (rk1 * rk1 * rk1 - rk * rk * rk) / (3.0 * h * h);
    fem.stiffness.diag[static_cast<std::size_t>(e)] += w;
    fem.stiffness.upper[static_cast<std::size_t>(e)] -= w;
    fem.stiffness.lower[static_cast<std::size_t>(e + 1)] -= w;
    fem.stiffness.diag[static_cast<std::size_t>(e + 1)] += w;

    // 单元质量：局部坐标 s ∈ [0,1]、r = rk + s·h 下把 r² 与 φ_iφ_j 的
    // 乘积逐项积分（r² 是 s 的二次式，结果是上面的解析公式）。
    const double m00 = h * (rk * rk / 3.0 + rk * h / 6.0 + h * h / 30.0);
    const double m01 = h * (rk * rk / 6.0 + rk * h / 6.0 + h * h / 20.0);
    const double m11 = h * (rk * rk / 3.0 + rk * h / 2.0 + h * h / 5.0);
    fem.mass.diag[static_cast<std::size_t>(e)] += m00;
    fem.mass.upper[static_cast<std::size_t>(e)] += m01;
    fem.mass.lower[static_cast<std::size_t>(e + 1)] += m01;
    fem.mass.diag[static_cast<std::size_t>(e + 1)] += m11;
  }
  return fem;
}

double fem_sphere_inventory(const FemSphere& fem,
                            const std::vector<double>& c) {
  if (static_cast<int>(c.size()) != fem.nodes) {
    throw std::invalid_argument("fem inventory: size mismatch");
  }
  const std::vector<double> weighted = fem.mass.multiply(c);
  double sum = 0.0;
  for (double value : weighted) {
    sum += value;
  }
  return 4.0 * kPi * sum;
}

FemMarchSystem make_fem_march_system(const FemSphere& fem, double surface_flux) {
  FemMarchSystem out;
  out.x_matrix = fem.mass;
  out.y_matrix = fem.stiffness;
  out.source.assign(static_cast<std::size_t>(fem.nodes), 0.0);
  // 右端只含末节点：分部积分后表面项为 -R² j v(R)，v = φ_N 只在表面为 1。
  const double radius = fem.r.back();
  out.source[static_cast<std::size_t>(fem.nodes - 1)] =
      -radius * radius * surface_flux;
  for (int i = 0; i < fem.nodes; ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    // 质量阵行和 = ∫ r² φ_i dr（∑_j φ_j ≡ 1 处处成立）。
    const double w = fem.mass.diag[k] + fem.mass.lower[k] + fem.mass.upper[k];
    out.x_matrix.diag[k] /= w;
    out.x_matrix.lower[k] /= w;
    out.x_matrix.upper[k] /= w;
    out.y_matrix.diag[k] /= w;
    out.y_matrix.lower[k] /= w;
    out.y_matrix.upper[k] /= w;
    out.source[k] /= w;
  }
  return out;
}

}  // namespace dfn
