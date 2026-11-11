#pragma once

namespace dfn {

// 电极开路电位（OCP）的经验拟合，供 SPM/SPMe/DFN 与拟合实验共用。
// 负、正极各一条曲线，自变量为化学计量数 stoich = c_s / c_max。

// 开路电位及其对化学计量数的导数。
struct OcpValue {
  double value = 0.0;  // V，开路电位
  double slope = 0.0;  // V / 单位化学计量数，即 dU/dstoich
};

// 半电池 OCP 的示意性拟合（斜率符号与取值范围是重点）。配置加载器支持
// Chen2020 参数表后，这些曲线会被真实数据替换。
OcpValue negative_ocp(double stoich);
OcpValue positive_ocp(double stoich);

}  // namespace dfn
