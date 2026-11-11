#pragma once

#include <functional>
#include <string>
#include <vector>

#include "dfn/solvers/dae.hpp"

namespace dfn {

// SUNDIALS IDA 封装：求解 F(t, y, y') = 0，支持稀疏 Jacobian（KLU）、
// 事件（根函数）检测与指定时刻输出。本文件不引入 SUNDIALS 头文件。

// 线性求解后端：klu 为稀疏直接法，dense 为稠密直接法。
enum class LinearSolverKind { klu, dense };

// 根函数 g(t, y, y')。g 变号时 IDA 报告一次事件并停止积分。
using EventFunction =
    std::function<double(double, const double*, const double*)>;

struct IdaOptions {
  double rtol = 1e-6;  // 相对容差
  // 绝对容差：空表示每个分量都取 1e-8；长度 1 表示全体统一值；
  // 长度 n 表示逐分量设置。
  std::vector<double> atol;
  double max_step = 0.0;  // 最大步长，<= 0 表示由 IDA 自适应决定
  bool supply_jacobian = true;  // 是否提供解析 Jacobian（否则用内部差分）
  bool calculate_consistent_initial = true;  // 是否先求一致初值
  LinearSolverKind linear_solver = LinearSolverKind::klu;  // 线性求解后端
  std::vector<EventFunction> events;  // 事件（根）函数列表
  std::vector<double> output_times;   // 除 t0 与终点外额外采样的时刻
};

struct SolveReport {
  bool ok = false;              // 是否成功返回
  std::string message;          // 结果说明或错误信息
  int error_code = 0;           // IDA 返回标志
  int steps = 0;                // 内部时间步数
  int residual_evals = 0;       // 残差求值次数
  int jacobian_evals = 0;       // Jacobian 求值次数
  int linear_setups = 0;        // 线性求解器设置次数
  int nonlinear_iterations = 0; // 牛顿迭代累计次数
  double t_final = 0.0;         // 实际到达的结束时刻
  std::vector<double> y;        // 终点状态
  std::vector<double> yp;       // 终点导数
  int event_index = -1;         // 触发事件在 IdaOptions::events 中的下标
  double t_event = 0.0;         // 事件发生时刻，s
  std::vector<double> output_t;               // 采样时刻
  std::vector<std::vector<double>> output_y;  // 各采样时刻的状态
};

// 用 SUNDIALS IDA 积分 F(t, y, y') = 0；稀疏路径需要编译时打开
// DFN_USE_KLU。
SolveReport ida_integrate(const DaeModel& model, const std::vector<double>& y0,
                          const std::vector<double>& yp0, double t0,
                          double t_end, const IdaOptions& options = {});

// 把 IDA 返回标志转换成可读的错误描述。
std::string ida_error_message(int flag);

}  // namespace dfn
