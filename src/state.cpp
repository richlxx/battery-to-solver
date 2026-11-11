// 本文件实现 DFN 状态向量布局（StateLayout）的构造与初值生成。
// 状态按块依次排列：负极固相浓度、正极固相浓度、电解液浓度、
// 负极固相电位、正极固相电位、电解液电位；各块内按 x（再按 r）展平。

#include "dfn/state.hpp"

#include <stdexcept>

namespace dfn {

// 校验网格划分数并构造状态布局。
// 参数：n_x_neg / n_x_sep / n_x_pos 为负极、隔膜、正极沿厚度方向的控制体数；
//       n_r_neg / n_r_pos 为负极、正极活性颗粒沿半径方向的球壳数。
// 返回：StateLayout，其中各 *_offset() 给出对应块在 y 中的起始下标，
//       总长度 layout.size() 为固相浓度、电解液浓度与三个电位块长度之和。
// 几何约定：x 从负极集流体指向正极集流体；r = 0 为颗粒中心，r = R 为表面。
StateLayout make_layout(int n_x_neg, int n_x_sep, int n_x_pos,
                        int n_r_neg, int n_r_pos) {
  // 任一方向划分数非正都会让状态长度为 0，或使后续索引越界，直接拒绝。
  if (n_x_neg <= 0 || n_x_sep <= 0 || n_x_pos <= 0 ||
      n_r_neg <= 0 || n_r_pos <= 0) {
    throw std::invalid_argument("grid counts must be positive");
  }
  StateLayout layout;
  layout.n_x_neg = n_x_neg;
  layout.n_x_sep = n_x_sep;
  layout.n_x_pos = n_x_pos;
  layout.n_r_neg = n_r_neg;
  layout.n_r_pos = n_r_pos;
  return layout;
}

// 生成初始状态向量 y（长度 layout.size()，先整体置 0）。
// 初始假设：颗粒内部固相浓度均匀、电解液浓度沿厚度方向均匀，即平衡态；
// 电位未知量保持 0，待模型装配方程后由 Newton 迭代确定。
// 参数：params 提供三个初始浓度 c0；layout 决定各块下标。
std::vector<double> make_initial_state(const Params& params,
                                       const StateLayout& layout) {
  std::vector<double> y(static_cast<std::size_t>(layout.size()), 0.0);

  // 负极固相：同一 x 位置的所有半径球壳连续存放，
  // 下标为 cs_neg = i_x*n_r_neg + i_r，初值取均匀的 neg.c0。
  for (int i_x = 0; i_x < layout.n_x_neg; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_neg; ++i_r) {
      y[layout.cs_neg(i_x, i_r)] = params.neg.c0;
    }
  }
  // 正极固相：布局与负极相同，但球壳数 n_r_pos 可以不同。
  for (int i_x = 0; i_x < layout.n_x_pos; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_pos; ++i_r) {
      y[layout.cs_pos(i_x, i_r)] = params.pos.c0;
    }
  }
  // 电解液浓度：覆盖负极、隔膜、正极的全部 x 位置，初值取自 elyte.c0。
  for (int i_x = 0; i_x < layout.n_x_total(); ++i_x) {
    y[layout.ce(i_x)] = params.elyte.c0;
  }

  return y;
}

}  // namespace dfn
