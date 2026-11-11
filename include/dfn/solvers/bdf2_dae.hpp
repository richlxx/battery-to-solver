#pragma once

#include <string>
#include <vector>

#include "dfn/solvers/dae.hpp"

namespace dfn {

// DAE 形式 F(t, y, y') = 0 的定步长 BDF2 时间推进，是 IDA 的可读替代品，
// 便于在教程中逐步跟踪 Newton 与 Jacobian 装配。

struct Bdf2DaeOptions {
  double newton_tolerance = 1e-10;  // 每步牛顿迭代的收敛阈值
  int newton_max_iterations = 20;   // 每步牛顿迭代次数上限
};

struct Bdf2DaeResult {
  std::vector<double> t;               // 采样时间点
  std::vector<std::vector<double>> y;  // 各时间点对应的状态
  bool ok = true;                      // 是否所有步都成功
  int steps = 0;                       // 实际时间步数
  int newton_iterations = 0;           // 牛顿迭代累计次数
  std::string message;                 // 失败原因或结束说明
};

// 用 BDF2 对 F(t, y, y') = 0 做定步长积分：第一步用 BDF1 启动，之后用
// 二阶格式。y0 需与模型一致（代数约束已满足）。
Bdf2DaeResult bdf2_dae_integrate(const DaeModel& model, std::vector<double> y0,
                                 double t0, double t_end, double h,
                                 const Bdf2DaeOptions& options = {});

}  // namespace dfn
