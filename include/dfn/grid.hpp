#pragma once

#include <vector>

namespace dfn {

// 一维有限体积网格及网格生成函数，是扩散离散与模型几何的共同基础。
// 区域被划分为 n 个控制体（cell），相邻控制体之间共 n + 1 个界面（face）。

// 一维有限体积网格。未知量位于单元中心，通量在界面上计算。
struct Grid1D {
  int n = 0;  // 单元（控制体）个数

  // 界面坐标（长度 n + 1）。face[0] 与 face[n] 是区域边界；
  // 第 i 个单元占据 [face[i], face[i + 1]]。球坐标网格中存的是半径。
  std::vector<double> face;

  // 单元中心坐标（长度 n）。
  std::vector<double> center;

  // 单元体积（长度 n），一维时为界面面积乘单元长度。
  std::vector<double> volume;

  // 界面面积（长度 n + 1），area[i] 对应 face[i]。
  std::vector<double> area;

  // 每个界面看到的相邻中心距（长度 n + 1），用于面上梯度差分。内部
  // 元素参与扩散模板；两个边界元素只是半单元距离，不被内部模板使用。
  std::vector<double> spacing;
};

// 在 [0, length] 上生成 n 个单元的均匀网格，每个单元体积均为
// cross_section * length / n（cross_section 为横截面积）。
Grid1D make_uniform_line_grid(int n, double length, double cross_section);

// 生成分段均匀网格：counts[k] 个长度为 lengths[k] 的块首尾相接，
// 块与块的界面被精确保存在 face 中，块内其余界面取相邻单元中心的中点。
Grid1D make_segmented_line_grid(const std::vector<int>& counts,
                                const std::vector<double>& lengths,
                                double cross_section);

// 半径为 radius 的球体径向均匀网格。此处 face 存的是半径，
// area[i] = 4 * pi * face[i]^2，volume[i] 为 face[i] 与 face[i + 1]
// 之间球壳的体积。
Grid1D make_uniform_sphere_grid(int n, double radius);

}  // namespace dfn
