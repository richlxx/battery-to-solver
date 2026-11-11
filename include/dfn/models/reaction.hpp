#pragma once

#include "dfn/params.hpp"

namespace dfn {

// Butler-Volmer 反应动力学封装：计算界面反应通量及其偏导，供模型残差与
// Jacobian 共用。

// 单个界面的 Butler-Volmer 反应量及其偏导。符号约定：j 为正值表示
// 锂离开固相颗粒（脱锂，阳极方向）。
struct ReactionValue {
  double j = 0.0;        // mol/(m^2 s)，界面反应通量
  double dj_deta = 0.0;  // 对过电位 eta 的偏导
  double dj_dce = 0.0;   // 对电解液浓度 c_e 的偏导
  double dj_dcs = 0.0;   // 对表面固相浓度 c_s 的偏导
};

// 在给定局部状态（电解液浓度 ce、表面固相浓度 cs_surf、固相/电解液
// 电势 phis、phie）下计算反应通量与偏导；positive 选择正极的参数与符号。
ReactionValue reaction_at(const Params& params, const ElectrodeParams& electrode,
                          bool positive, double ce, double cs_surf, double phis,
                          double phie);

// 由给定通量 j 反解 Butler-Volmer 得到过电位。用于 SPM/SPMe 等
// 通量由外加电流直接确定的简化模型。
double overpotential_from_flux(const Params& params,
                               const ElectrodeParams& electrode, bool positive,
                               double ce, double cs_surf, double j);

}  // namespace dfn
