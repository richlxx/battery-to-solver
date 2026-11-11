#pragma once

#include <vector>

#include "dfn/solvers/dae.hpp"

namespace dfn {

// 教学用刚性 DAE：一维 "链式" 松弛-扩散-代数约束系统。规模可调、有闭式解，
// 便于演示 BDF2、Newton 与刚性比等概念。

struct ChainDaeOptions {
  int cells = 1;           // 单元个数
  double lambda = 100.0;   // 松弛速率，越大方程越刚性
  double diffusion = 0.0;  // 相邻单元之间的扩散耦合强度
};

// 每个单元一个刚性松弛方程加一个代数变量，未知量按
// (c_0, w_0, c_1, w_1, ...) 交错排列：
//
//   c_i' + lambda (c_i - w_i) - diffusion * lap(c)_i = 0
//   w_i - sin(t) = 0
//
// 单单元且 diffusion = 0 时解有闭式表达式，可直接验证数值误差。
class ChainDae final : public DaeModel {
 public:
  explicit ChainDae(ChainDaeOptions options = {});

  int cells() const { return options_.cells; }        // 单元个数
  double lambda() const { return options_.lambda; }   // 松弛速率
  double diffusion() const { return options_.diffusion; }  // 扩散耦合强度

  int size() const override { return 2 * options_.cells; }
  // 残差 F(t, y, y')，写入 out（长度 size()）。
  void residual(double t, const double* y, const double* yp,
                double* out) const override;
  CscPattern jacobian_pattern() const override { return pattern_; }
  // 变量类型：c_i 为微分变量（1），w_i 为代数变量（0）。
  std::vector<int> variable_types() const override;
  // 组装 dF/dy + cj * dF/dy'（cj 为 BDF/IDA 的差分系数）。
  void jacobian(double t, const double* y, const double* yp, double cj,
                CscBuilder& out) const override;

 private:
  ChainDaeOptions options_;
  CscPattern pattern_;
};

// 与 chain_dae_consistent_derivative() 配对时在 t0 同时满足两个方程的
// 均匀初值：每个单元取单单元闭式浓度，代数变量取 sin(t0)。
std::vector<double> chain_dae_consistent_state(const ChainDae& model, double t0);

// 与上述初值匹配的时间导数，使 F(t0, y, y') = 0。
std::vector<double> chain_dae_consistent_derivative(const ChainDae& model,
                                                    double t0);

// 在均匀解上叠加正弦扰动，配合 chain_dae_derivative_from_state() 仍满足
// F(t0, y, y') = 0：
//   c_i = exact(t0) + amplitude * sin((i + 1) * pi / (cells + 1))
std::vector<double> chain_dae_mode_state(const ChainDae& model, double t0,
                                         double amplitude);

// 给定 y，反解出使 F(t, y, y') = 0 的 y'（用于构造一致初值）。
std::vector<double> chain_dae_derivative_from_state(
    const ChainDae& model, double t, const std::vector<double>& y);

// diffusion = 0 时单单元的精确浓度。
double chain_dae_exact_concentration(double lambda, double t);

// 上述精确浓度的时间导数。
double chain_dae_exact_concentration_rate(double lambda, double t);

}  // namespace dfn
