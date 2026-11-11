// 本文件实现一维有限体积网格的构造：均匀直线网格、分段（多区域）网格
// 与球坐标径向网格。未知量存放在控制体（cell）中心，通量在界面上计算。

#include "dfn/grid.hpp"

#include <cmath>
#include <stdexcept>

namespace dfn {

namespace {
// 球面积与球体积公式所需的 π，直接给出双精度字面量，避免依赖 M_PI 宏。
constexpr double kPi = 3.14159265358979323846;
}

// 在 [0, length] 上生成 n 个等宽控制体的直线网格。
// 参数：n 控制体数；length 区域长度（m）；cross_section 截面面积（m^2），
//       平面几何下所有界面与所有控制体共用该面积。
// 返回：Grid1D，其中 face[i] = i*dx，center[i] = (i+0.5)*dx，volume[i] = A*dx。
Grid1D make_uniform_line_grid(int n, double length, double cross_section) {
  if (n <= 0 || length <= 0.0 || cross_section <= 0.0) {
    throw std::invalid_argument("line grid needs positive counts and sizes");
  }

  Grid1D grid;
  grid.n = n;
  grid.face.resize(static_cast<std::size_t>(n + 1));
  grid.center.resize(static_cast<std::size_t>(n));
  grid.volume.resize(static_cast<std::size_t>(n));
  grid.area.assign(static_cast<std::size_t>(n + 1), cross_section);
  grid.spacing.resize(static_cast<std::size_t>(n + 1));

  // 界面位于 i*dx，控制体中心取相邻两界面的中点；
  // 均匀网格下每个界面的中心距都等于 dx。
  const double dx = length / static_cast<double>(n);
  for (int i = 0; i <= n; ++i) {
    grid.face[i] = i * dx;
    grid.spacing[i] = dx;
  }
  for (int i = 0; i < n; ++i) {
    grid.center[i] = (i + 0.5) * dx;
    grid.volume[i] = cross_section * dx;
  }
  return grid;
}

// 生成按物理区域分段的直线网格，用于负极/隔膜/正极串联的一维电芯。
// 参数：counts[k] 为第 k 段的控制体数；lengths[k] 为第 k 段厚度（m）；
//       cross_section 为极片面积（m^2），全部界面共用。
// 返回：共 total = Σ counts[k] 个控制体。段与段之间的物理界面精确落在 face 上，
//       段内其余内部界面取相邻两个控制体中心的中点。
Grid1D make_segmented_line_grid(const std::vector<int>& counts,
                                const std::vector<double>& lengths,
                                double cross_section) {
  if (counts.empty() || counts.size() != lengths.size() ||
      cross_section <= 0.0) {
    throw std::invalid_argument("segmented grid needs matching, non-empty blocks");
  }
  int total = 0;
  for (std::size_t k = 0; k < counts.size(); ++k) {
    if (counts[k] <= 0 || lengths[k] <= 0.0) {
      throw std::invalid_argument("segmented grid needs positive blocks");
    }
    total += counts[k];
  }

  Grid1D grid;
  grid.n = total;
  grid.face.assign(static_cast<std::size_t>(total) + 1, 0.0);
  grid.center.assign(static_cast<std::size_t>(total), 0.0);
  grid.volume.assign(static_cast<std::size_t>(total), 0.0);
  grid.area.assign(static_cast<std::size_t>(total) + 1, cross_section);
  grid.spacing.assign(static_cast<std::size_t>(total) + 1, 0.0);

  // 第一遍：按段填充控制体中心与体积。每段内部等宽，宽度为
  // lengths[r] / counts[r]；boundary[k] 记录第 k 段起点的 x 坐标。
  std::vector<double> boundary(counts.size() + 1, 0.0);
  int offset = 0;
  for (std::size_t r = 0; r < counts.size(); ++r) {
    const double dx = lengths[r] / static_cast<double>(counts[r]);
    boundary[r + 1] = boundary[r] + lengths[r];
    for (int k = 0; k < counts[r]; ++k) {
      const int i = offset + k;
      grid.center[static_cast<std::size_t>(i)] =
          boundary[r] + (static_cast<double>(k) + 0.5) * dx;
      grid.volume[static_cast<std::size_t>(i)] = cross_section * dx;
    }
    offset += counts[r];
  }

  // 第二遍：段边界（负极/隔膜、隔膜/正极）精确保留物理位置，
  // 段内其余内部界面放在相邻两个控制体中心的中点。
  grid.face[0] = 0.0;
  grid.face[static_cast<std::size_t>(total)] = boundary.back();
  offset = counts[0];
  for (std::size_t r = 1; r < counts.size(); ++r) {
    grid.face[static_cast<std::size_t>(offset)] = boundary[r];
    offset += counts[r];
  }
  for (int i = 1; i < total; ++i) {
    if (grid.face[static_cast<std::size_t>(i)] == 0.0) {
      grid.face[static_cast<std::size_t>(i)] =
          0.5 * (grid.center[static_cast<std::size_t>(i - 1)] +
                 grid.center[static_cast<std::size_t>(i)]);
    }
    grid.spacing[static_cast<std::size_t>(i)] =
        grid.center[static_cast<std::size_t>(i)] -
        grid.center[static_cast<std::size_t>(i - 1)];
  }

  // 边界界面只有一侧有控制体中心，其 spacing 取半个控制体宽度；
  // 内部扩散模板不会引用这两个值，保留它们只为数组长度一致。
  grid.spacing[0] = 0.5 * (grid.center[0] - grid.face[0]);
  grid.spacing[static_cast<std::size_t>(total)] =
      grid.face[static_cast<std::size_t>(total)] -
      grid.center[static_cast<std::size_t>(total - 1)];
  return grid;
}

// 生成球坐标下的均匀径向网格，描述活性颗粒内部的锂扩散。
// 参数：n 半径方向上的球壳数；radius 颗粒半径（m）。
// 返回：face[i] = i*dr 为半径坐标，area[i] = 4πr_i^2 为球面积（m^2），
//       volume[i] = 4π/3 (r_{i+1}^3 - r_i^3) 为第 i 层球壳体积（m^3）。
Grid1D make_uniform_sphere_grid(int n, double radius) {
  if (n <= 0 || radius <= 0.0) {
    throw std::invalid_argument("sphere grid needs positive counts and radius");
  }

  Grid1D grid;
  grid.n = n;
  grid.face.resize(static_cast<std::size_t>(n + 1));
  grid.center.resize(static_cast<std::size_t>(n));
  grid.volume.resize(static_cast<std::size_t>(n));
  grid.area.resize(static_cast<std::size_t>(n + 1));
  grid.spacing.resize(static_cast<std::size_t>(n + 1));

  // 径向网格里 face[i] 是半径而不是直角坐标；第 i 个控制体是半径
  // face[i] 与 face[i+1] 之间的球壳，其体积必须用两球体积之差计算。
  const double dr = radius / static_cast<double>(n);
  for (int i = 0; i <= n; ++i) {
    grid.face[i] = i * dr;
    grid.area[i] = 4.0 * kPi * grid.face[i] * grid.face[i];
    grid.spacing[i] = dr;
  }
  for (int i = 0; i < n; ++i) {
    grid.center[i] = (i + 0.5) * dr;
    const double r0 = grid.face[i];
    const double r1 = grid.face[i + 1];
    grid.volume[i] = 4.0 * kPi / 3.0 * (r1 * r1 * r1 - r0 * r0 * r0);
  }
  return grid;
}

}  // namespace dfn
