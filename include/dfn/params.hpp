#pragma once

#include <string>
#include <vector>

namespace dfn {

// 电池模型的物理常数与参数结构，全部使用 SI 单位。配置加载
// （sim/config.hpp）与各模型构造函数都以 Params 作为唯一输入。

// 法拉第常数，C/mol。
inline constexpr double kFaraday = 96485.33212;
// 理想气体常数，J/(mol K)。
inline constexpr double kGasConstant = 8.314462618;

// 单个多孔电极（负极或正极）的几何、传输与反应动力学参数。
// 厚度沿 x 方向，radius 是球形活性颗粒半径；"-" 表示无量纲量。
struct ElectrodeParams {
  double thickness = 0.0;    // m，电极厚度
  double radius = 0.0;       // m，活性颗粒半径
  double porosity = 0.0;     // -，孔隙率（电解液体积分数）
  double bruggeman = 1.5;    // -，Bruggeman 指数，修正有效传输系数
  double c_max = 0.0;        // mol/m^3，固相最大锂浓度
  double c0 = 0.0;           // mol/m^3，固相初始锂浓度
  double d_s = 0.0;          // m^2/s，固相扩散系数
  double k = 0.0;            // 反应速率常数（Butler-Volmer 前因子）
  double alpha_a = 0.5;      // 阳极传递系数
  double alpha_c = 0.5;      // 阴极传递系数
  double sigma = 0.0;        // S/m，固相电子电导率
};

// 隔膜参数：只参与电解液传输，没有活性颗粒。
struct SeparatorParams {
  double thickness = 0.0;    // m，隔膜厚度
  double porosity = 0.0;     // -，孔隙率
  double bruggeman = 1.5;    // -，Bruggeman 指数
};

// 电解液参数，整条 x 方向（负极、隔膜、正极）共用同一组数值。
struct ElectrolyteParams {
  double c0 = 0.0;           // mol/m^3，电解液初始锂离子浓度
  double d_e = 0.0;          // m^2/s，电解液扩散系数
  double t_plus = 0.0;       // -，锂离子迁移数
  double kappa = 0.0;        // S/m，电解液离子电导率
};

// 全电池参数：面积、温度，以及负极、正极、隔膜、电解液四部分信息。
struct Params {
  double area = 0.0;         // m^2，电极横截面积
  double temperature = 0.0;  // K，等温模型的温度
  ElectrodeParams neg;
  ElectrodeParams pos;
  SeparatorParams sep;
  ElectrolyteParams elyte;

  // 电池总厚度 L = L_neg + L_sep + L_pos，单位 m。
  double cell_length() const {
    return neg.thickness + sep.thickness + pos.thickness;
  }
};

// 用于示例与单元测试的示意参数；真实标定请替换为 Chen2020 参数表。
Params make_demo_params();

// 校验参数是否物理上自洽（厚度、浓度、扩散系数为正等），不合法时抛异常。
void validate(const Params& params);

}  // namespace dfn
