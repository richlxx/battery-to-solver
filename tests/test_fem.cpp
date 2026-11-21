// 测试文件：有限元装配与 FEM/FV 对照（对应 fem 模块，第 12 章）。
// 验证质量阵/刚度阵的对称性与总量恒等式、行缩放系统、推进中的总锂量
// 守恒，以及第 12 章的两个定量结论：两条离散路线的表面浓度都是二阶
// 收敛，且同一网格下 FV 的误差常数明显大于 FEM。

#include <cmath>
#include <iostream>
#include <vector>

#include "dfn/fem.hpp"
#include "dfn/grid.hpp"
#include "dfn/operators.hpp"

namespace {

constexpr double kPi = 3.14159265358979323846;

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

// 纯相对容差：装配量都是正数且量级相近；两者同时为 0 时视为一致。
bool close_relative(double a, double b, double tolerance) {
  const double scale = std::max(std::abs(a), std::abs(b));
  if (scale == 0.0) {
    return true;
  }
  return std::abs(a - b) <= tolerance * scale;
}

// 与 labs/fem_sphere.cpp 相同的物理参数与推进设置。
struct Problem {
  double radius = 5.0e-6;       // 颗粒半径 [m]
  double diffusivity = 3.0e-14; // 固相扩散系数 [m^2/s]
  double flux = 1.0e-6;         // 外表面通量 [mol/(m^2 s)]，正方向朝外
  double c0 = 1000.0;           // 初始均匀浓度 [mol/m^3]
  double dt = 0.5;              // 时间步长 [s]
};

std::vector<double> run_fem(const Problem& p, int elements, int steps) {
  const dfn::FemSphere fem =
      dfn::assemble_fem_sphere(elements, p.radius, p.diffusivity);
  const dfn::FemMarchSystem sys = dfn::make_fem_march_system(fem, p.flux);
  const std::vector<double> c(static_cast<std::size_t>(fem.nodes), p.c0);
  return dfn::march_tridiagonal(sys.x_matrix, sys.y_matrix, sys.source, c,
                                p.dt, steps);
}

std::vector<double> run_fv(const Problem& p, int cells, int steps) {
  const dfn::Grid1D grid = dfn::make_uniform_sphere_grid(cells, p.radius);
  const std::vector<double> d(static_cast<std::size_t>(cells), p.diffusivity);
  const dfn::DiffusionOperator op =
      dfn::assemble_diffusion(grid, d, 0.0, p.flux);
  const std::vector<double> c(static_cast<std::size_t>(cells), p.c0);
  return dfn::march_tridiagonal(dfn::identity_matrix(cells),
                                dfn::negated_matrix(op.matrix), op.source, c,
                                p.dt, steps);
}

double fv_surface(const Problem& p, int cells,
                  const std::vector<double>& c) {
  const double dr = p.radius / cells;
  return c[static_cast<std::size_t>(cells - 1)] -
         p.flux * dr / (2.0 * p.diffusivity);
}

}  // namespace

int main() {
  constexpr Problem p{};
  const int steps = 1200;  // t_end = 600 s，与第 12 章实验一致

  // —— 装配的结构性质（n = 8 个单元、9 个节点）——
  const dfn::FemSphere fem =
      dfn::assemble_fem_sphere(8, p.radius, p.diffusivity);

  // 对称性：弱形式的双线性型对称，M 与 K 都应满足 lower[i] = upper[i-1]。
  for (int i = 1; i < fem.nodes; ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    if (!close_relative(fem.mass.lower[k], fem.mass.upper[k - 1], 1e-12)) {
      return fail("mass matrix symmetry");
    }
    if (!close_relative(fem.stiffness.lower[k], fem.stiffness.upper[k - 1],
                        1e-12)) {
      return fail("stiffness matrix symmetry");
    }
  }

  // 总量恒等式：1ᵀM1 = ∫ r² dr = R³/3，这是全局守恒的代数根源。
  double mass_sum = 0.0;
  for (int i = 0; i < fem.nodes; ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    mass_sum += fem.mass.diag[k] + fem.mass.lower[k] + fem.mass.upper[k];
  }
  if (!close_relative(mass_sum, p.radius * p.radius * p.radius / 3.0, 1e-12)) {
    return fail("mass matrix total must equal R^3/3");
  }

  // 常数是扩散的稳态：K 乘均匀场应得零（舍入量级）。
  const std::vector<double> uniform(static_cast<std::size_t>(fem.nodes), 1.0);
  const std::vector<double> steady = fem.stiffness.multiply(uniform);
  double steady_scale = 0.0;
  for (double value : fem.stiffness.diag) {
    steady_scale += std::abs(value);
  }
  for (double value : steady) {
    if (std::abs(value) > 1e-12 * steady_scale) {
      return fail("uniform field must be steady");
    }
  }

  // —— 行缩放系统：X 每行行和为 1，右端只落在末节点且为负（流出）——
  const dfn::FemMarchSystem sys = dfn::make_fem_march_system(fem, p.flux);
  for (int i = 0; i < fem.nodes; ++i) {
    const std::size_t k = static_cast<std::size_t>(i);
    const double row = sys.x_matrix.diag[k] + sys.x_matrix.lower[k] +
                       sys.x_matrix.upper[k];
    if (std::abs(row - 1.0) > 1e-12) {
      return fail("row-scaled mass rows must sum to one");
    }
    if (i + 1 < fem.nodes && sys.source[k] != 0.0) {
      return fail("surface flux must only enter the last node");
    }
  }
  if (sys.source.back() >= 0.0) {
    return fail("outward flux must give a negative source");
  }

  // —— 推进中的总锂量守恒：变化应精确等于 -4πR²j t ——
  {
    const dfn::FemSphere fem8 =
        dfn::assemble_fem_sphere(8, p.radius, p.diffusivity);
    const dfn::FemMarchSystem sys8 = dfn::make_fem_march_system(fem8, p.flux);
    const std::vector<double> c0(static_cast<std::size_t>(fem8.nodes), p.c0);
    const double initial = dfn::fem_sphere_inventory(fem8, c0);
    const int n_steps = 100;  // 50 s
    const std::vector<double> c =
        dfn::march_tridiagonal(sys8.x_matrix, sys8.y_matrix, sys8.source, c0,
                               p.dt, n_steps);
    const double now = dfn::fem_sphere_inventory(fem8, c);
    const double expected = -4.0 * kPi * p.radius * p.radius * p.flux *
                            p.dt * n_steps;
    if (!close_relative(now - initial, expected, 1e-10)) {
      return fail("fem global conservation");
    }
  }

  // —— 第 12 章的定量结论：二阶收敛与误差常数差 ——
  const std::vector<double> c_ref = run_fem(p, 2048, steps);
  const double surface_ref = c_ref.back();
  // 参考解本身应在物理范围内（初始 1000，600 s 内平均降约 360）。
  if (!(500.0 < surface_ref && surface_ref < 1000.0)) {
    return fail("reference surface concentration out of physical range");
  }

  double err_fem = 0.0;
  double err_fv = 0.0;
  double prev_fem = 0.0;
  double prev_fv = 0.0;
  for (int n : {16, 32, 64, 128}) {
    prev_fem = err_fem;
    prev_fv = err_fv;
    err_fem = std::abs(run_fem(p, n, steps).back() - surface_ref);
    err_fv = std::abs(fv_surface(p, n, run_fv(p, n, steps)) - surface_ref);
    if (n > 16) {
      // 误差必须单调下降，且 64 档的观测阶落在二阶附近（实测约 2.000）。
      if (err_fem >= prev_fem || err_fv >= prev_fv) {
        return fail("surface error must decrease under refinement");
      }
      if (n == 64) {
        const double order_fem = std::log2(prev_fem / err_fem);
        const double order_fv = std::log2(prev_fv / err_fv);
        if (!(1.85 <= order_fem && order_fem <= 2.15)) {
          return fail("fem observed order must be second order");
        }
        if (!(1.85 <= order_fv && order_fv <= 2.15)) {
          return fail("fv observed order must be second order");
        }
      }
    }
  }
  // 同一网格（循环结束时为 128 单元）下 FV 的表面误差约为 FEM 的 12 倍：
  // 外推半格的代价体现在常数上。断言放宽到 8–16 倍，避免过度约束。
  const double ratio = err_fv / err_fem;
  if (!(8.0 <= ratio && ratio <= 16.0)) {
    return fail("fv error constant must dominate fem at the same mesh");
  }

  std::cout << "fem ok (fem/fv order ~2, fv/fem error ratio "
            << ratio << ")\n";
  return 0;
}
