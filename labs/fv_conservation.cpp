#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/grid.hpp"
#include "dfn/operators.hpp"

// 实验（第 4 章 有限体积离散）：3 层洋葱球对称网格上的显式扩散，检验通量守恒。
// sealed 段外表面零通量，总锂量应当不变；open 段给定外表面通量 j，总量应按
// -A*j*dt*step 线性下降。输出 CSV：case,step,total,drift,expected，末尾为逐格账本。
namespace {

// 体积加权总量 Σ V_i c_i：c 为浓度 mol/m^3，V 为壳体积 m^3，结果为 mol。
double total(const dfn::Grid1D& grid, const std::vector<double>& c) {
  double sum = 0.0;
  for (int i = 0; i < grid.n; ++i) {
    sum += grid.volume[static_cast<std::size_t>(i)] *
           c[static_cast<std::size_t>(i)];
  }
  return sum;
}

// 显式欧拉推进一步：c^{n+1} = c^n + dt * (M c^n + source)。
// M c 是体积平均后的 dc/dt（mol/(m^3 s)），source 为边界通量除以体积。
std::vector<double> advance(const dfn::DiffusionOperator& op,
                            std::vector<double> c, double dt) {
  const std::vector<double> rhs = op.matrix.multiply(c);
  for (std::size_t i = 0; i < c.size(); ++i) {
    c[i] += dt * (rhs[i] + op.source[i]);
  }
  return c;
}

// 运行一个算例并逐步打印总量：drift = 当前总量 - 初始总量（mol），
// expected = -A * j * dt * step 是外表面通量累计带走的锂量，二者应吻合。
void run(const char* name, const dfn::Grid1D& grid,
         const dfn::DiffusionOperator& op, double surface_flux, double dt,
         int steps) {
  // 三层壳的初始浓度（mol/m^3），由内向外递减。
  std::vector<double> c = {1200.0, 1000.0, 800.0};
  const double initial = total(grid, c);
  // grid.area[n] 是最外层（颗粒表面）面积，单位 m^2。
  const double area = grid.area[static_cast<std::size_t>(grid.n)];
  std::cout << "case,step,total,drift,expected\n";
  for (int step = 0; step <= steps; ++step) {
    const double now = total(grid, c);
    const double expected = -area * surface_flux * dt * step;
    std::cout << name << ',' << step << ',' << now << ',' << (now - initial)
              << ',' << expected << '\n';
    if (step < steps) {
      c = advance(op, c, dt);
    }
  }
}

}  // namespace

int main() {
  // 球颗粒半径 5 um，固相扩散系数 3e-14 m^2/s，时间步 10 s，共 5 步；
  // 表面通量 1e-5 mol/(m^2 s)，正值表示从颗粒内部向外流出。
  constexpr double radius = 5.0e-6;
  constexpr double diffusivity = 3.0e-14;
  constexpr double dt = 10.0;
  constexpr int steps = 5;
  constexpr double surface_flux = 1.0e-5;
  const int n = 3;

  // 3 个控制体的均匀球网格：face 为半径，area = 4*pi*r^2，volume 为球壳体积。
  const dfn::Grid1D grid = dfn::make_uniform_sphere_grid(n, radius);
  const std::vector<double> D(static_cast<std::size_t>(n), diffusivity);
  // sealed：两端零通量；open：仅外表面施加通量 surface_flux。
  const dfn::DiffusionOperator sealed = dfn::assemble_diffusion(grid, D);
  const dfn::DiffusionOperator open =
      dfn::assemble_diffusion(grid, D, 0.0, surface_flux);

  std::cout << std::scientific << std::setprecision(6);
  std::cout << "n=" << n << " R=" << radius << " D=" << diffusivity
            << " dt=" << dt << " j=" << surface_flux << '\n';
  run("sealed", grid, sealed, 0.0, dt, steps);
  run("open", grid, open, surface_flux, dt, steps);

  // 逐格账本：V_i * (dc/dt)_i，单位 mol/s。内部面通量对相邻两格一进一出，
  // 求和后应相互抵消，故 sealed 情形的 ledger_sum 应为 0（舍入误差量级）。
  std::cout << "\ncell,volume,dc_dt,ledger\n";
  const std::vector<double> c = {1200.0, 1000.0, 800.0};
  const std::vector<double> rhs = sealed.matrix.multiply(c);
  double ledger_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const double ledger = grid.volume[static_cast<std::size_t>(i)] *
                          rhs[static_cast<std::size_t>(i)];
    ledger_sum += ledger;
    std::cout << i << ',' << grid.volume[static_cast<std::size_t>(i)] << ','
              << rhs[static_cast<std::size_t>(i)] << ',' << ledger << '\n';
  }
  std::cout << "ledger_sum," << ledger_sum << '\n';
  return 0;
}
