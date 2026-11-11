#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "dfn/models/chain_dae.hpp"
#include "dfn/solvers/bdf2_dae.hpp"
#include "dfn/solvers/idaklu.hpp"

// 第 8 章配套实验：把自写 BDF2 与 IDA+KLU 放在同一个 DAE 上对比。
// 三组实验依次是容差扫描、BDF2 与 IDA 的精度/步数/耗时对比、以及链式 DAE
// 规模放大时 KLU 稀疏求解与稠密求解的对照。误差均相对解析解计算。

namespace {

using Clock = std::chrono::steady_clock;

// 返回从 start 到当前时刻的墙钟毫秒数。
double elapsed_ms(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// 在单电池 DAE（cells=1、lambda=100、无扩散）上扫描 rtol 从 1e-3 到 1e-9，
// 观察 IDA 的步数、残差/Jacobian 评估次数和终值误差如何随容差收紧而增长。
void tolerance_sweep() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);
  const std::vector<double> yp0 = dfn::chain_dae_consistent_derivative(model, 0.0);
  // 解析解在 t_end 处给出，是唯一可信的误差基准。
  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);

  std::cout << "== single cell DAE: tolerance sweep (IDA + KLU) ==\n";
  std::cout << "rtol\tatol\tsteps\tres_evals\tjac_evals\terror\n";
  for (const double rtol : {1e-3, 1e-5, 1e-7, 1e-9}) {
    dfn::IdaOptions ida;
    ida.rtol = rtol;
    // atol 随 rtol 同比例收紧，避免绝对容差成为瓶颈。
    ida.atol = {1e-3 * rtol};
    ida.calculate_consistent_initial = false;
    const dfn::SolveReport report =
        dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
    if (!report.ok) {
      std::cout << rtol << "\t" << ida.atol[0] << "\tfailed: " << report.message
                << "\n";
      continue;
    }
    const double error = std::abs(report.y[0] - exact);
    std::cout << rtol << "\t" << ida.atol[0] << "\t" << report.steps << "\t"
              << report.residual_evals << "\t\t" << report.jacobian_evals << "\t"
              << error << "\n";
  }
  std::cout << "\n";
}

// 固定步长 BDF2（h = 0.01 / 0.005 / 0.0025）与自适应步长 IDA 对比：
// 同一模型、同一精确解，比较步数、Newton 迭代、Jacobian 评估与墙钟耗时。
void bdf2_versus_ida() {
  dfn::ChainDaeOptions options;
  options.cells = 1;
  options.lambda = 100.0;
  const dfn::ChainDae model(options);

  const double t_end = 1.0;
  const std::vector<double> y0 = dfn::chain_dae_consistent_state(model, 0.0);
  const std::vector<double> yp0 = dfn::chain_dae_consistent_derivative(model, 0.0);
  const double exact = dfn::chain_dae_exact_concentration(100.0, t_end);

  std::cout << "== hand-written BDF2 versus IDA (single cell) ==\n";
  std::cout << "solver\tparams\tsteps\tnewton\tjac\terror\tms\n";

  // BDF2 走 h 减半序列，可以观察时间离散误差随 h 的二阶下降。
  for (const double h : {0.01, 0.005, 0.0025}) {
    dfn::Bdf2DaeOptions bdf;
    bdf.newton_tolerance = 1e-12;
    bdf.newton_max_iterations = 30;
    const auto start = Clock::now();
    const dfn::Bdf2DaeResult result =
        dfn::bdf2_dae_integrate(model, y0, 0.0, t_end, h, bdf);
    const double ms = elapsed_ms(start);
    if (!result.ok) {
      std::cout << "BDF2\th=" << h << "\tfailed: " << result.message << "\n";
      continue;
    }
    std::cout << "BDF2\th=" << h << "\t" << result.steps << "\t"
              << result.newton_iterations << "\t-\t"
              << std::abs(result.y.back()[0] - exact) << "\t" << ms << "\n";
  }

  // IDA 只给容差，步长由误差控制器自适应选择，因此参数列是 rtol。
  for (const double rtol : {1e-5, 1e-7, 1e-9}) {
    dfn::IdaOptions ida;
    ida.rtol = rtol;
    ida.atol = {1e-3 * rtol};
    ida.calculate_consistent_initial = false;
    const auto start = Clock::now();
    const dfn::SolveReport report =
        dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
    const double ms = elapsed_ms(start);
    if (!report.ok) {
      std::cout << "IDA\trtol=" << rtol << "\tfailed: " << report.message << "\n";
      continue;
    }
    std::cout << "IDA+KLU\trtol=" << rtol << "\t" << report.steps << "\t"
              << report.nonlinear_iterations << "\t" << report.jacobian_evals
              << "\t" << std::abs(report.y[0] - exact) << "\t" << ms << "\n";
  }
  std::cout << "\n";
}

// 规模放大实验：cells 从 80 增到 1280，比较 KLU 稀疏 LU 与稠密 LU
// 在步数、线性求解准备次数（setups）和总耗时上的差异。
void chain_scale_up() {
  std::cout << "== chain DAE: KLU versus dense linear solver ==\n";
  std::cout << "cells\tn\tlin_solver\tsteps\tjac\tsetups\tms\n";

  for (const int cells : {80, 320, 1280}) {
    dfn::ChainDaeOptions options;
    options.cells = cells;
    // lambda 取 1e4 让系统足够刚性，再加上较强扩散形成三对角耦合。
    options.lambda = 1e4;
    options.diffusion = 100.0;
    const dfn::ChainDae model(options);

    // 初值取正弦扰动模态，并配套算出满足 F(t0, y, y') = 0 的导数。
    const double t_end = 0.02;
    const std::vector<double> y0 = dfn::chain_dae_mode_state(model, 0.0, 0.1);
    const std::vector<double> yp0 =
        dfn::chain_dae_derivative_from_state(model, 0.0, y0);

    // 同一初值、同一容差下只切换线性求解器，保证对比公平。
    for (const auto kind :
         {dfn::LinearSolverKind::klu, dfn::LinearSolverKind::dense}) {
      dfn::IdaOptions ida;
      ida.rtol = 1e-6;
      ida.atol = {1e-8};
      ida.calculate_consistent_initial = false;
      ida.linear_solver = kind;
      const auto start = Clock::now();
      const dfn::SolveReport report =
          dfn::ida_integrate(model, y0, yp0, 0.0, t_end, ida);
      const double ms = elapsed_ms(start);
      const char* name =
          kind == dfn::LinearSolverKind::klu ? "KLU" : "dense";
      if (!report.ok) {
        std::cout << cells << "\t" << model.size() << "\t" << name
                  << "\tfailed: " << report.message << "\n";
        continue;
      }
      std::cout << cells << "\t" << model.size() << "\t" << name << "\t\t"
                << report.steps << "\t" << report.jacobian_evals << "\t"
                << report.linear_setups << "\t" << ms << "\n";
    }
  }
  std::cout << "\n";
}

}  // namespace

// 实验入口：依次运行三组对比，输出均带表头，方便直接贴进教程。
int main() {
  std::cout << std::setprecision(4);
  tolerance_sweep();
  bdf2_versus_ida();
  chain_scale_up();
  return 0;
}
