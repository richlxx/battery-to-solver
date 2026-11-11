// 本文件实现负极与正极的开路电位（OCP）经验拟合及其对化学计量数的导数。
// 曲线形状（单调性、斜率符号）符合真实材料趋势，数值为简化拟合，仅用于
// 教学与演示，后续会由 Chen2020 参数表替换。

#include "dfn/thermo.hpp"

#include <algorithm>
#include <cmath>

namespace dfn {

// 负极（石墨类）半电池开路电位，随嵌锂程度 x = c_s / c_max 增大而降低。
// 参数 stoich：负极颗粒表面的无量纲化学计量数 c_s / c_max，内部钳制到 [0, 1]。
// 返回：value 为电位（V），slope 为 dU/dx（V），供 Newton 迭代构造
//       dU/dc_s = slope / c_max 时使用。x→0 时约 0.95 V，x→1 时约 0.05 V。
OcpValue negative_ocp(double stoich) {
  const double x = std::clamp(stoich, 0.0, 1.0);
  const double one_minus_x = 1.0 - x;
  OcpValue result;
  // 三项经验式：常数平台 + 低锂化端的指数上升 + 高锂化端的二次项。
  // 0.05、0.8、0.1 与指数衰减率 10 都是曲线形状参数，不是物理常数。
  result.value =
      0.05 + 0.8 * std::exp(-10.0 * x) + 0.1 * one_minus_x * one_minus_x;
  // 对 x 的解析导数：指数项求导得 -8*exp(-10x)，二次项求导得 -0.2*(1-x)。
  result.slope = -8.0 * std::exp(-10.0 * x) - 0.2 * one_minus_x;
  return result;
}

// 正极（层状氧化物类）半电池开路电位，随嵌锂程度 y = c_s / c_max 增大而降低。
// 参数 stoich：正极颗粒表面的无量纲化学计量数，钳制到 [0, 1]。
// 返回：value 为电位（V，y=0 时 4.1 V，y=1 时 3.9 V），slope 恒为 -0.2 V。
// 这里取线性形式，便于教学时手算核对 Jacobian。
OcpValue positive_ocp(double stoich) {
  const double y = std::clamp(stoich, 0.0, 1.0);
  OcpValue result;
  result.value = 3.9 + 0.2 * (1.0 - y);
  result.slope = -0.2;
  return result;
}

}  // namespace dfn
