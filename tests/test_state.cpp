// 测试文件：参数与状态布局（对应 params / state 模块）。
// 验证 cell_length 计算、make_layout 的索引覆盖、初始状态取值与 validate 的错误处理。

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dfn/params.hpp"
#include "dfn/state.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

}  // namespace

int main() {
  const dfn::Params params = dfn::make_demo_params();
  dfn::validate(params);  // demo 参数必须通过合法性校验（正厚度、正孔隙率等）

  // 电芯总厚度应为负极、隔膜、正极厚度之和；
  // 三者都是同一量级的有限小数，1e-18 的绝对容差只允许舍入误差。
  const double expected_length = 85.0e-6 + 12.0e-6 + 75.0e-6;
  if (std::abs(params.cell_length() - expected_length) > 1e-18) {
    return fail("cell length");
  }

  const dfn::StateLayout layout = dfn::make_layout(3, 2, 4, 5, 6);
  // 未知量数 = 负极固相浓度 N_x_neg*N_r_neg + 正极固相浓度 N_x_pos*N_r_pos
  //          + 电解液浓度与电位 2*N_x_total + 两电极固相电位 N_x_neg+N_x_pos。
  const int expected_size = 3 * 5 + 4 * 6 + 9 + 3 + 4 + 9;
  if (layout.size() != expected_size) {
    return fail("layout size");
  }

  // 遍历全部索引，统计每个偏移被访问的次数；正确布局必须恰好覆盖一次。
  std::vector<int> hits(static_cast<std::size_t>(layout.size()), 0);
  for (int i_x = 0; i_x < layout.n_x_neg; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_neg; ++i_r) {
      ++hits[layout.cs_neg(i_x, i_r)];
    }
  }
  for (int i_x = 0; i_x < layout.n_x_pos; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_pos; ++i_r) {
      ++hits[layout.cs_pos(i_x, i_r)];
    }
  }
  for (int i_x = 0; i_x < layout.n_x_total(); ++i_x) {
    ++hits[layout.ce(i_x)];
    ++hits[layout.phie(i_x)];
  }
  for (int i_x = 0; i_x < layout.n_x_neg; ++i_x) {
    ++hits[layout.phis_neg(i_x)];
  }
  for (int i_x = 0; i_x < layout.n_x_pos; ++i_x) {
    ++hits[layout.phis_pos(i_x)];
  }
  for (int count : hits) {
    if (count != 1) {
      return fail("index coverage");  // 漏访问或重叠访问都说明偏移计算有误
    }
  }

  // 初始状态：固相与电解液浓度取平衡浓度 c0，电位从 0 起步。
  const std::vector<double> y = dfn::make_initial_state(params, layout);
  if (y.size() != static_cast<std::size_t>(layout.size())) {
    return fail("initial state size");
  }
  if (std::abs(y[layout.cs_neg(2, 4)] - params.neg.c0) > 1e-12) {
    return fail("negative concentration");
  }
  if (std::abs(y[layout.cs_pos(3, 5)] - params.pos.c0) > 1e-12) {
    return fail("positive concentration");
  }
  if (std::abs(y[layout.ce(4)] - params.elyte.c0) > 1e-12) {
    return fail("electrolyte concentration");
  }
  // 电位初值精确为 0，可直接用 != 比较（这里是有意检查赋值而非数值近似）。
  if (y[layout.phis_neg(0)] != 0.0 || y[layout.phie(0)] != 0.0) {
    return fail("potentials should start at zero");
  }

  // 负厚度在物理上无意义，validate 必须抛出 std::invalid_argument。
  dfn::Params bad = params;
  bad.neg.thickness = -1.0;
  bool threw = false;
  try {
    dfn::validate(bad);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  if (!threw) {
    return fail("validation should reject negative thickness");
  }

  std::cout << "state layout ok (" << layout.size() << " unknowns)\n";
  return 0;
}
