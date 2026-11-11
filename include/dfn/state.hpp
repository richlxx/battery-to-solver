#pragma once

#include <vector>

#include "dfn/params.hpp"

namespace dfn {

// DFN 状态向量（未知量数组）的分段布局与初值构造。布局把各物理量映射
// 为全局下标，残差、Jacobian 与输出代码都依赖同一套索引约定。

// DFN 未知量分段布局，顺序为：
//   c_s(neg) -> c_s(pos) -> c_e -> phi_s(neg) -> phi_s(pos) -> phi_e
// 其中两段固相浓度按 (i_x, i_r) 二维展平，r 方向连续存放（同一 x 位置
// 的径向节点相邻）。浓度单位 mol/m^3，电势单位 V。
struct StateLayout {
  int n_x_neg = 0;  // 负极 x 方向控制体数
  int n_x_sep = 0;  // 隔膜 x 方向控制体数
  int n_x_pos = 0;  // 正极 x 方向控制体数
  int n_r_neg = 0;  // 负极颗粒径向控制体数
  int n_r_pos = 0;  // 正极颗粒径向控制体数

  // 整条 x 方向的节点总数（负极 + 隔膜 + 正极）。
  int n_x_total() const { return n_x_neg + n_x_sep + n_x_pos; }

  // 以下 offset() 给出各段在全局向量中的起始下标；第一段从 0 开始。
  int cs_neg_offset() const { return 0; }
  int cs_pos_offset() const { return n_x_neg * n_r_neg; }
  int ce_offset() const { return cs_pos_offset() + n_x_pos * n_r_pos; }
  int phis_neg_offset() const { return ce_offset() + n_x_total(); }
  int phis_pos_offset() const { return phis_neg_offset() + n_x_neg; }
  int phie_offset() const { return phis_pos_offset() + n_x_pos; }
  // 状态向量总长度。
  int size() const { return phie_offset() + n_x_total(); }

  // 把 (x 节点, 径向节点) 映射为全局下标：i_x * n_r + i_r 加上段偏移。
  int cs_neg(int i_x, int i_r) const {
    return cs_neg_offset() + i_x * n_r_neg + i_r;
  }
  int cs_pos(int i_x, int i_r) const {
    return cs_pos_offset() + i_x * n_r_pos + i_r;
  }
  // 全局 x 节点索引（0 起，跨负极/隔膜/正极连续编号）的标量未知量下标。
  int ce(int i_x) const { return ce_offset() + i_x; }
  int phis_neg(int i_x) const { return phis_neg_offset() + i_x; }
  int phis_pos(int i_x) const { return phis_pos_offset() + i_x; }
  int phie(int i_x) const { return phie_offset() + i_x; }
};

// 依据各方向控制体数建立布局；任一网格数为非正时抛异常。
StateLayout make_layout(int n_x_neg, int n_x_sep, int n_x_pos,
                        int n_r_neg, int n_r_pos);

// 生成 t = 0 的初始状态：固相与电解液浓度取各自 c0，电势全为 0。
std::vector<double> make_initial_state(const Params& params,
                                       const StateLayout& layout);

}  // namespace dfn
