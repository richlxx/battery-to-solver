// 本文件实现演示参数集的构造（make_demo_params）与参数合法性校验（validate）。
// 参数统一使用 SI 单位：长度 m、浓度 mol/m^3、扩散系数 m^2/s、电导率 S/m、温度 K。

#include "dfn/params.hpp"

#include <stdexcept>

namespace dfn {

// 构造一套用于演示与测试的 P2D/DFN 参数：单层软包电芯、25 °C、1 M 电解液。
// 负极取石墨、正极取层状氧化物；数值为常见文献量级，用于跑通代码与画图，
// 并不对应某一商用型号电芯，后续由 Chen2020 参数表替换。
Params make_demo_params() {
  Params p;
  // 极片面积 0.0165 m^2，约合 150 mm × 110 mm 的单层极片；温度 298.15 K 即 25 °C。
  p.area = 0.0165;
  p.temperature = 298.15;

  // 负极：85 µm 厚多孔石墨，活性颗粒半径 5 µm，孔隙率 0.25；
  // 最大固相锂浓度 30000 mol/m^3，初始均匀荷锂 80%（对应满充附近的负极）。
  p.neg.thickness = 85.0e-6;
  p.neg.radius = 5.0e-6;
  p.neg.porosity = 0.25;
  p.neg.bruggeman = 1.5;
  p.neg.c_max = 30000.0;
  p.neg.c0 = 0.8 * p.neg.c_max;
  p.neg.d_s = 3.0e-14;
  p.neg.k = 1.0e-6;
  p.neg.sigma = 100.0;

  // 正极：75 µm 厚，孔隙率 0.33（压实更高）；最大固相浓度 60000 mol/m^3，
  // 初始荷锂 30%；固相电导率 10 S/m，比石墨低一个量级。
  p.pos.thickness = 75.0e-6;
  p.pos.radius = 5.0e-6;
  p.pos.porosity = 0.33;
  p.pos.bruggeman = 1.5;
  p.pos.c_max = 60000.0;
  p.pos.c0 = 0.3 * p.pos.c_max;
  p.pos.d_s = 4.0e-15;
  p.pos.k = 1.0e-6;
  p.pos.sigma = 10.0;

  // 隔膜：12 µm 厚，孔隙率 0.47。Bruggeman 指数 1.5 用于有效输运系数
  // D_eff = D*eps^1.5、kappa_eff = kappa*eps^1.5 的修正。
  p.sep.thickness = 12.0e-6;
  p.sep.porosity = 0.47;
  p.sep.bruggeman = 1.5;

  // 电解液：1 M（1000 mol/m^3）锂盐，Li+ 迁移数 0.26，
  // 本体扩散系数 5e-10 m^2/s，离子电导率 1 S/m。
  p.elyte.c0 = 1000.0;
  p.elyte.d_e = 5.0e-10;
  p.elyte.t_plus = 0.26;
  p.elyte.kappa = 1.0;

  return p;
}

// 校验参数是否满足物理与数值上的基本约束，不合法时抛出 std::invalid_argument。
// 只检查符号与区间（正厚度、孔隙率在 (0,1) 内等），不检查参数之间的耦合关系，
// 也不检查 Bruggeman 指数等不影响可解性的量。
void validate(const Params& params) {
  // require：条件为假时抛出异常；错误信息为英文，便于在日志中定位出错字段。
  const auto require = [](bool condition, const std::string& message) {
    if (!condition) {
      throw std::invalid_argument(message);
    }
  };

  // 单个电极的逐项校验，name 仅用于拼接错误信息（如 "negative electrode"）。
  // alpha_a、alpha_c 为正可保证 Butler-Volmer 的指数项不发散；
  // c0 落在 [0, c_max] 内，等价于初始化学计量数 c0/c_max 落在 [0, 1] 内。
  const auto check_electrode = [&require](const ElectrodeParams& e,
                                          const std::string& name) {
    require(e.thickness > 0.0, name + ": thickness must be positive");
    require(e.radius > 0.0, name + ": particle radius must be positive");
    require(e.porosity > 0.0 && e.porosity < 1.0,
            name + ": porosity must be in (0, 1)");
    require(e.c_max > 0.0, name + ": c_max must be positive");
    require(e.c0 >= 0.0 && e.c0 <= e.c_max,
            name + ": c0 must be in [0, c_max]");
    require(e.d_s > 0.0, name + ": d_s must be positive");
    require(e.k > 0.0, name + ": reaction rate must be positive");
    require(e.alpha_a > 0.0 && e.alpha_c > 0.0,
            name + ": transfer coefficients must be positive");
    require(e.sigma > 0.0, name + ": sigma must be positive");
  };

  // 电芯整体：面积为正，温度为正（绝对温度，K）。
  require(params.area > 0.0, "cell area must be positive");
  require(params.temperature > 0.0, "temperature must be positive");
  check_electrode(params.neg, "negative electrode");
  check_electrode(params.pos, "positive electrode");

  // 隔膜不含活性颗粒，只需检查厚度与孔隙率。
  require(params.sep.thickness > 0.0, "separator: thickness must be positive");
  require(params.sep.porosity > 0.0 && params.sep.porosity < 1.0,
          "separator: porosity must be in (0, 1)");

  // 电解液：c0、d_e、kappa 为正；t_plus 取 [0,1)（等于 1 意味着阴离子不迁移，
  // 会破坏本模型电荷守恒的推导）。
  require(params.elyte.c0 > 0.0, "electrolyte: c0 must be positive");
  require(params.elyte.d_e > 0.0, "electrolyte: d_e must be positive");
  require(params.elyte.t_plus >= 0.0 && params.elyte.t_plus < 1.0,
          "electrolyte: t_plus must be in [0, 1)");
  require(params.elyte.kappa > 0.0, "electrolyte: kappa must be positive");
}

}  // namespace dfn
