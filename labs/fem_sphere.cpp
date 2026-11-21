#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/fem.hpp"
#include "dfn/grid.hpp"
#include "dfn/operators.hpp"

// 实验（第 12 章 有限元离散）：在球颗粒径向扩散这个第 4 章的老问题上，
// 用 P1 有限元重新离散，与有限体积逐项对照。半离散系统统一写成
// X dc/dt = Y c + f：FEM 的 X 是（行缩放的）一致质量阵、Y 是刚度阵；
// FV 的 X 是单位阵（质量被体积除掉）、Y 是 -A。两条路共用同一条 BDF2
// 推进，收敛研究量到的因此只有空间差。
// 输出三段 CSV：conv（表面浓度收敛研究）、ledger（总锂量账本）、
// profile（末时刻浓度剖面，供配图）。
namespace {

constexpr double kPi = 3.14159265358979323846;

// 物理与时间推进参数，单位见各字段注释。
struct Problem {
  double radius;        // 颗粒半径 [m]
  double diffusivity;   // 固相扩散系数 [m^2/s]
  double flux;          // 外表面通量 [mol/(m^2 s)]，正值表示锂流出颗粒
  double c0;            // 初始均匀浓度 [mol/m^3]
  double dt;            // 时间步长 [s]
  int steps;            // 总步数，t_end = dt * steps
};

// 装配 FEM 并推进到 t_end，返回节点浓度（末元素即表面浓度）。
std::vector<double> run_fem(const Problem& p, int elements) {
  const dfn::FemSphere fem =
      dfn::assemble_fem_sphere(elements, p.radius, p.diffusivity);
  const dfn::FemMarchSystem sys = dfn::make_fem_march_system(fem, p.flux);
  const std::vector<double> c(static_cast<std::size_t>(fem.nodes), p.c0);
  return dfn::march_tridiagonal(sys.x_matrix, sys.y_matrix, sys.source, c,
                                p.dt, p.steps);
}

// 装配 FV 并推进到 t_end，返回单元（格心）浓度。装配与第 4 章完全相同：
// assemble_diffusion 已把每个格子的方程除以体积，所以 X 是单位阵。
std::vector<double> run_fv(const Problem& p, int cells) {
  const dfn::Grid1D grid = dfn::make_uniform_sphere_grid(cells, p.radius);
  const std::vector<double> d(static_cast<std::size_t>(cells), p.diffusivity);
  const dfn::DiffusionOperator op =
      dfn::assemble_diffusion(grid, d, 0.0, p.flux);
  const std::vector<double> c(static_cast<std::size_t>(cells), p.c0);
  return dfn::march_tridiagonal(dfn::identity_matrix(cells),
                                dfn::negated_matrix(op.matrix), op.source, c,
                                p.dt, p.steps);
}

// FV 的表面浓度，即第 4 章式 (4.12) 的单侧线性外推：
// c_surf = c_{n-1} - j Δr / (2D)。表面在格外半格处，不是未知量。
double fv_surface(const Problem& p, int cells,
                  const std::vector<double>& c) {
  const double dr = p.radius / cells;
  return c[static_cast<std::size_t>(cells - 1)] -
         p.flux * dr / (2.0 * p.diffusivity);
}

}  // namespace

int main() {
  // 与第 4 章 fv_conservation 同一颗颗粒：R = 5 µm、D = 3e-14 m^2/s；
  // 通量取得更小（1e-6），保证 600 s 内浓度远未耗尽。dt = 0.5 s 时 BDF2
  // 的时间误差在 1e-6 mol/m^3 量级，低于最细网格的空间误差，收敛研究
  // 量到的才是空间阶数（对照 4.9 节“时间误差没跟上会盖住空间误差”）。
  constexpr Problem p{5.0e-6, 3.0e-14, 1.0e-6, 1000.0, 0.5, 1200};

  std::cout << std::scientific << std::setprecision(6);
  std::cout << "# fem_sphere: R=" << p.radius << " D=" << p.diffusivity
            << " j=" << p.flux << " c0=" << p.c0 << " dt=" << p.dt
            << " steps=" << p.steps << '\n';

  // —— 收敛研究：表面浓度对参考解（FEM、2048 单元）的误差与观测阶 ——
  // FEM 表面就是节点，浓度直接取节点值；FV 表面在格子外半格处，须外推。
  const std::vector<double> c_ref = run_fem(p, 2048);
  const double surface_ref =
      c_ref.back();  // N=2048 的空间误差已低于时间误差，可当参考
  std::cout << "conv,method,n,dof,c_surf,err,order\n";
  for (const bool fem : {true, false}) {
    const char* method = fem ? "fem" : "fv";
    double prev_err = 0.0;
    for (int n : {4, 8, 16, 32, 64, 128, 256}) {
      // FEM：n 个单元、n+1 个自由度；FV：n 个格子、n 个自由度。
      const double surface =
          fem ? run_fem(p, n).back() : fv_surface(p, n, run_fv(p, n));
      const int dof = fem ? n + 1 : n;
      const double err = std::abs(surface - surface_ref);
      // 观测阶 p = log2(e(h) / e(h/2))，与第 2、9 章同一个定义。
      std::cout << "conv," << method << ',' << n << ',' << dof << ',' << surface
                << ',' << err;
      if (n > 4) {
        std::cout << ',' << std::log2(prev_err / err);
      }
      std::cout << '\n';
      prev_err = err;
    }
  }

  // —— 总锂量账本：两种方法都应精确守恒（舍入误差量级）——
  // FEM 的总量来自 v = 1 的弱形式，FV 的总量来自逐格求和；期望值都是
  // -4πR² j t（外表面以恒定速率带走锂）。
  std::cout << "\nledger,method,step,time,total,drift,expected\n";
  {
    const int n = 16;
    const dfn::FemSphere fem =
        dfn::assemble_fem_sphere(n, p.radius, p.diffusivity);
    const dfn::FemMarchSystem sys = dfn::make_fem_march_system(fem, p.flux);
    const std::vector<double> c0(static_cast<std::size_t>(fem.nodes), p.c0);
    const double initial = dfn::fem_sphere_inventory(fem, c0);
    dfn::march_tridiagonal(
        sys.x_matrix, sys.y_matrix, sys.source, c0, p.dt, p.steps,
        [&](int step, const std::vector<double>& c) {
          if (step % 240 != 0) {
            return;
          }
          const double now = dfn::fem_sphere_inventory(fem, c);
          std::cout << "ledger,fem," << step << ',' << step * p.dt << ','
                    << now << ',' << (now - initial) << ','
                    << -4.0 * kPi * p.radius * p.radius * p.flux * step * p.dt
                    << '\n';
        });
  }
  {
    const int n = 16;
    const dfn::Grid1D grid = dfn::make_uniform_sphere_grid(n, p.radius);
    const std::vector<double> d(static_cast<std::size_t>(n), p.diffusivity);
    const dfn::DiffusionOperator op =
        dfn::assemble_diffusion(grid, d, 0.0, p.flux);
    const std::vector<double> c0(static_cast<std::size_t>(n), p.c0);
    double initial = 0.0;
    for (int i = 0; i < n; ++i) {
      initial += grid.volume[static_cast<std::size_t>(i)] *
                 c0[static_cast<std::size_t>(i)];
    }
    dfn::march_tridiagonal(
        dfn::identity_matrix(n), dfn::negated_matrix(op.matrix), op.source,
        c0, p.dt, p.steps, [&](int step, const std::vector<double>& c) {
          if (step % 240 != 0) {
            return;
          }
          double now = 0.0;
          for (int i = 0; i < n; ++i) {
            now += grid.volume[static_cast<std::size_t>(i)] *
                   c[static_cast<std::size_t>(i)];
          }
          std::cout << "ledger,fv," << step << ',' << step * p.dt << ','
                    << now << ',' << (now - initial) << ','
                    << -grid.area[static_cast<std::size_t>(n)] * p.flux *
                           step * p.dt
                    << '\n';
        });
  }

  // —— 末时刻剖面（供配图）：粗网格上 FEM 节点值、FV 格心值与外推表面
  // 点，叠加参考解。FV 的外推点相对偏低，正是 (4.12) 的误差来源。
  std::cout << "\nprofile,method,r,c\n";
  {
    const int n = 16;
    const std::vector<double> cf = run_fem(p, n);
    for (int k = 0; k <= n; ++k) {
      std::cout << "profile,fem," << p.radius * k / n << ','
                << cf[static_cast<std::size_t>(k)] << '\n';
    }
    const dfn::Grid1D grid = dfn::make_uniform_sphere_grid(n, p.radius);
    const std::vector<double> cv = run_fv(p, n);
    for (int i = 0; i < n; ++i) {
      std::cout << "profile,fv," << grid.center[static_cast<std::size_t>(i)]
                << ',' << cv[static_cast<std::size_t>(i)] << '\n';
    }
    std::cout << "profile,fv_surf," << p.radius << ','
              << fv_surface(p, n, cv) << '\n';
    for (int k = 0; k <= 2048; k += 128) {
      std::cout << "profile,ref," << p.radius * k / 2048.0 << ','
                << c_ref[static_cast<std::size_t>(k)] << '\n';
    }
  }
  return 0;
}
