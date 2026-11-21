// 用 OpenFOAM 有限体积框架求解等温 DFN（P2D）模型。
//
// 与仓库自研求解器的物理方程一一对应，但离散与时间推进换成 OpenFOAM：
//   * x 方向（负极|隔膜|正极）建成 1D fvMesh，c_e、phi_e、phi_s 为 volScalarField；
//   * 颗粒径向扩散是每个 x 单元内的子网格，用三对角隐式推进（OpenFOAM
//     没有内建的“每单元一列径向网格”结构，这部分在求解器内部完成）；
//   * Butler-Volmer 在每个界面单元上做局部牛顿迭代，再以线性化源项进入
//     phi_s / phi_e 的 fvMatrix，整体用 Picard 迭代耦合。
//
// 单位全部为 SI：浓度 mol/m^3，电位 V，通量 mol/(m^2 s)，电流密度 A/m^2。
// 场的 dimensionSet 显式写出这些量纲，OpenFOAM 在矩阵装配时会检查一致性。
#include "fvCFD.H"
#include "IOdictionary.H"
#include "fixedValueFvPatchFields.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace Foam;

namespace {

// 物理常数。
constexpr scalar kFaraday = 96485.33212;   // C/mol
constexpr scalar kGasConstant = 8.314462618;  // J/(mol K)
constexpr scalar kPi = 3.14159265358979323846;

// OpenFOAM 的维度顺序：质量、长度、时间、温度、摩尔、电流、光强。
const dimensionSet dimConc(0, -3, 0, 0, 1, 0, 0);          // mol/m^3
const dimensionSet dimConcRate(0, -3, -1, 0, 1, 0, 0);     // mol/(m^3 s)
const dimensionSet dimPhi(1, 2, -3, 0, 0, -1, 0);          // V
const dimensionSet dimDiffusivity(0, 2, -1, 0, 0, 0, 0);   // m^2/s
const dimensionSet dimSigma(-1, -3, 3, 0, 0, 2, 0);        // S/m
const dimensionSet dimInvLength(0, -1, 0, 0, 0, 0, 0);     // 1/m
const dimensionSet dimVolCurrent(0, -3, 0, 0, 0, 1, 0);    // A/m^3
const dimensionSet dimVolCurrentPerVolt(-1, -5, 3, 0, 0, 2, 0);  // A/(m^3 V)

scalar clampValue(scalar value, scalar low, scalar high) {
  return std::max(low, std::min(value, high));
}

// 开路电位及其对化学计量数的导数，与 include/dfn/thermo.hpp 保持一致。
struct OcpValue {
  scalar value = 0.0;
  scalar slope = 0.0;
};

OcpValue negativeOcp(scalar stoich) {
  const scalar x = clampValue(stoich, 0.0, 1.0);
  const scalar oneMinusX = 1.0 - x;
  OcpValue result;
  result.value = 0.05 + 0.8 * std::exp(-10.0 * x) + 0.1 * oneMinusX * oneMinusX;
  result.slope = -8.0 * std::exp(-10.0 * x) - 0.2 * oneMinusX;
  return result;
}

OcpValue positiveOcp(scalar stoich) {
  const scalar y = clampValue(stoich, 0.0, 1.0);
  OcpValue result;
  result.value = 3.9 + 0.2 * (1.0 - y);
  result.slope = -0.2;
  return result;
}

// 单个电极的参数。nRadial 是颗粒径向壳层数，其余与 configs/demo.yaml 同名。
struct ElectrodeParams {
  scalar thickness = 0.0;   // m
  scalar radius = 0.0;      // m
  scalar porosity = 0.0;    // -
  scalar bruggeman = 1.5;   // -
  scalar cMax = 0.0;        // mol/m^3
  scalar c0 = 0.0;          // mol/m^3
  scalar dS = 0.0;          // m^2/s
  scalar k = 0.0;           // BV 前因子
  scalar alphaA = 0.5;
  scalar alphaC = 0.5;
  scalar sigma = 0.0;       // S/m
  label nRadial = 6;
  scalar sigmaEff = 0.0;    // S/m，Bruggeman 修正后
  scalar aS = 0.0;          // 1/m，比表面积
  scalar alphaSurf = 0.0;   // dr/(2 D_s)，表面浓度外推系数
};

struct ModelParams {
  scalar temperature = 298.15;
  ElectrodeParams neg;
  ElectrodeParams pos;
  scalar sepThickness = 0.0;
  scalar sepPorosity = 0.0;
  scalar sepBruggeman = 1.5;
  scalar ce0 = 0.0;
  scalar dE = 0.0;
  scalar tPlus = 0.0;
  scalar kappa = 0.0;

  scalar currentDensity = 0.0;  // A/m^2，正值放电
  label maxPicard = 30;
  scalar picardTol = 1e-4;
  scalar relaxPhi = 1.0;
  scalar relaxCe = 1.0;
  scalar voltageCutoff = 0.0;
};

ElectrodeParams readElectrode(const dictionary& dict) {
  ElectrodeParams electrode;
  electrode.thickness = readScalar(dict.lookup("thickness"));
  electrode.radius = readScalar(dict.lookup("radius"));
  electrode.porosity = readScalar(dict.lookup("porosity"));
  electrode.bruggeman = readScalar(dict.lookup("bruggeman"));
  electrode.cMax = readScalar(dict.lookup("c_max"));
  electrode.c0 = readScalar(dict.lookup("c0"));
  electrode.dS = readScalar(dict.lookup("d_s"));
  electrode.k = readScalar(dict.lookup("k"));
  electrode.alphaA = readScalar(dict.lookup("alpha_a"));
  electrode.alphaC = readScalar(dict.lookup("alpha_c"));
  electrode.sigma = readScalar(dict.lookup("sigma"));
  electrode.nRadial = dict.lookupOrDefault<label>("n_radial", 6);
  return electrode;
}

// 球形颗粒的径向网格：面积与体积按球坐标计算。
struct RadialGrid {
  label nR = 0;
  scalar radius = 0.0;
  scalar dr = 0.0;
  scalar dS = 0.0;
  scalarField area;
  scalarField volume;

  RadialGrid(label nR_, scalar radius_, scalar dS_)
      : nR(nR_), radius(radius_), dS(dS_) {
    dr = radius / static_cast<scalar>(nR);
    area.setSize(nR + 1);
    volume.setSize(nR);
    for (label i = 0; i <= nR; ++i) {
      const scalar r = static_cast<scalar>(i) * dr;
      area[i] = 4.0 * kPi * r * r;
    }
    for (label i = 0; i < nR; ++i) {
      const scalar r0 = static_cast<scalar>(i) * dr;
      const scalar r1 = static_cast<scalar>(i + 1) * dr;
      volume[i] = 4.0 * kPi / 3.0 * (r1 * r1 * r1 - r0 * r0 * r0);
    }
  }
};

// 三对角线性方程组的 Thomas 算法。lower[0] 与 upper[n-1] 不使用。
scalarField solveTridiagonal(const scalarField& lower, scalarField diag,
                             const scalarField& upper, scalarField rhs) {
  const label n = diag.size();
  for (label i = 1; i < n; ++i) {
    const scalar factor = lower[i] / diag[i - 1];
    diag[i] -= factor * upper[i - 1];
    rhs[i] -= factor * rhs[i - 1];
  }
  scalarField x(n);
  x[n - 1] = rhs[n - 1] / diag[n - 1];
  for (label i = n - 2; i >= 0; --i) {
    x[i] = (rhs[i] - upper[i] * x[i + 1]) / diag[i];
  }
  return x;
}

// Butler-Volmer 反应量及其偏导，符号约定与 src/models/reaction.cpp 相同：
// j > 0 表示锂离开固相进入电解液。
struct ReactionValue {
  scalar j = 0.0;
  scalar dj_dEta = 0.0;
  // 考虑 c_surf = c_last - alpha*j 后的约化导数：
  // dj/deta|_cLast = (dj/deta) / (1 + alpha*dj/dcs)。
  scalar dj_dEtaEff = 0.0;
  scalar dj_dCe = 0.0;
  scalar dj_dCs = 0.0;
  scalar eta = 0.0;
  scalar ocp = 0.0;
  scalar cSurf = 0.0;
};

scalar exchangeCurrent(const ModelParams& params,
                       const ElectrodeParams& electrode, scalar ce,
                       scalar cs) {
  const scalar cMax = electrode.cMax;
  const scalar csClamped = clampValue(cs, 1e-6 * cMax, cMax * (1.0 - 1e-6));
  const scalar ceClamped = std::max(ce, 1e-6 * params.ce0);
  return kFaraday * electrode.k * std::pow(ceClamped, electrode.alphaA) *
         std::pow(cMax - csClamped, electrode.alphaA) *
         std::pow(csClamped, electrode.alphaC);
}

ReactionValue reactionAt(const ModelParams& params,
                         const ElectrodeParams& electrode, bool positive,
                         scalar ce, scalar csSurf, scalar phiS, scalar phiE) {
  const scalar cMax = electrode.cMax;
  const scalar cs = clampValue(csSurf, 1e-6 * cMax, cMax * (1.0 - 1e-6));
  const scalar ceClamped = std::max(ce, 1e-6 * params.ce0);
  const scalar stoich = cs / cMax;
  const OcpValue ocp = positive ? positiveOcp(stoich) : negativeOcp(stoich);

  const scalar i0 = exchangeCurrent(params, electrode, ceClamped, cs);
  const scalar eta = phiS - phiE - ocp.value;
  const scalar a = electrode.alphaA * kFaraday / (kGasConstant * params.temperature);
  const scalar c = electrode.alphaC * kFaraday / (kGasConstant * params.temperature);
  const scalar ea = std::exp(clampValue(a * eta, -60.0, 60.0));
  const scalar ec = std::exp(clampValue(-c * eta, -60.0, 60.0));

  ReactionValue result;
  result.j = i0 / kFaraday * (ea - ec);
  result.dj_dEta = i0 / kFaraday * (a * ea + c * ec);
  result.dj_dCe = result.j * electrode.alphaA / ceClamped;
  result.dj_dCs =
      result.j * (-electrode.alphaA / (cMax - cs) + electrode.alphaC / cs) -
      result.dj_dEta * ocp.slope / cMax;
  const scalar reduced = 1.0 + electrode.alphaSurf * result.dj_dCs;
  result.dj_dEtaEff =
      result.dj_dEta / (std::abs(reduced) > 1e-12 ? reduced : 1e-12);
  result.eta = eta;
  result.ocp = ocp.value;
  result.cSurf = cs;
  return result;
}

// 给定颗粒最外层浓度、电解液浓度与两个电位，联立
//   j = BV(c_surf),  c_surf = c_last - alpha * j
// 求局部界面通量。外层用一个标量牛顿迭代，比把 c_surf 作为全局未知量更稳。
ReactionValue solveInterface(const ModelParams& params,
                             const ElectrodeParams& electrode, bool positive,
                             scalar ce, scalar cLast, scalar phiS, scalar phiE,
                             scalar jGuess) {
  const scalar cMin = 1e-6 * electrode.cMax;
  const scalar cMax = (1.0 - 1e-6) * electrode.cMax;
  const scalar jLow = (cLast - cMax) / electrode.alphaSurf;
  const scalar jHigh = (cLast - cMin) / electrode.alphaSurf;

  scalar j = clampValue(jGuess, jLow, jHigh);
  for (int iteration = 0; iteration < 50; ++iteration) {
    const scalar cSurf =
        clampValue(cLast - electrode.alphaSurf * j, cMin, cMax);
    const ReactionValue reaction =
        reactionAt(params, electrode, positive, ce, cSurf, phiS, phiE);
    const scalar residual = j - reaction.j;
    scalar derivative = 1.0 + electrode.alphaSurf * reaction.dj_dCs;
    if (std::abs(derivative) < 1e-14) {
      derivative = derivative >= 0.0 ? 1e-14 : -1e-14;
    }
    const scalar step = residual / derivative;
    j -= step;
    j = clampValue(j, jLow, jHigh);
    if (std::abs(step) < 1e-13 * (1.0 + std::abs(j))) {
      break;
    }
  }

  const scalar cSurf = clampValue(cLast - electrode.alphaSurf * j, cMin, cMax);
  ReactionValue reaction =
      reactionAt(params, electrode, positive, ce, cSurf, phiS, phiE);
  reaction.cSurf = cSurf;
  return reaction;
}

// 每个 x 单元独立求解球形颗粒扩散：
//   d c/dt = (1/r^2) d/dr (r^2 D dc/dr)，r=0 对称，r=R 处 -D dc/dr = j。
// 时间项用后向欧拉，径向用有限体积，最后用三对角求解。
// cOld 必须是时间步开始时的浓度：Picard 迭代中每次都用同一旧时间层重建，
// 若直接用上一轮迭代更新过的 c，等于在一步内重复推进 dt，浓度会被抽干。
void solveRadial(const RadialGrid& grid, const scalarField& cOld,
                 scalarField& cNew, const std::vector<scalar>& flux,
                 scalar dt) {
  const label nX = static_cast<label>(flux.size());
  const label nR = grid.nR;
  scalarField lower(nR, 0.0);
  scalarField upper(nR, 0.0);
  scalarField diag(nR, 0.0);
  scalarField rhs(nR, 0.0);

  for (label ix = 0; ix < nX; ++ix) {
    const label base = ix * nR;
    for (label ir = 0; ir < nR; ++ir) {
      const scalar area0 = grid.area[ir];
      const scalar area1 = grid.area[ir + 1];
      const scalar conductance = grid.dS / grid.dr;
      // r=R 的面用给定通量，不进入矩阵对角；r=0 的面积为 0，天然零通量。
      const scalar outerConductance = (ir + 1 < nR) ? area1 * conductance : 0.0;
      diag[ir] =
          grid.volume[ir] / dt + area0 * conductance + outerConductance;
      lower[ir] = (ir > 0) ? -area0 * conductance : 0.0;
      upper[ir] = (ir + 1 < nR) ? -area1 * conductance : 0.0;
      rhs[ir] = grid.volume[ir] * cOld[base + ir] / dt;
    }
    rhs[nR - 1] -= grid.area[nR] * flux[ix];

    const scalarField solved =
        solveTridiagonal(lower, diag, upper, rhs);
    for (label ir = 0; ir < nR; ++ir) {
      cNew[base + ir] = std::max(solved[ir], 0.0);
    }
  }
}

// 把内部单元值复制到边界，供 fvc::interpolate / laplacian 在集流体处使用。
void syncBoundary(volScalarField& field) {
  forAll(field.boundaryField(), patchi) {
    field.boundaryFieldRef()[patchi] ==
        field.boundaryField()[patchi].patchInternalField();
  }
}

// 体积平均化学计量数 <c_s>/c_max；c 按 (x, r) 展平。
scalar averageStoichiometry(const scalarField& c, const labelList& cells,
                            const RadialGrid& grid, const fvMesh& mesh,
                            scalar porosity, scalar cMax) {
  scalar particleVolume = 0.0;
  for (label ir = 0; ir < grid.nR; ++ir) {
    particleVolume += grid.volume[ir];
  }
  scalar numerator = 0.0;
  scalar denominator = 0.0;
  forAll(cells, ix) {
    scalar weighted = 0.0;
    for (label ir = 0; ir < grid.nR; ++ir) {
      weighted += c[ix * grid.nR + ir] * grid.volume[ir];
    }
    const scalar solidVolume = (1.0 - porosity) * mesh.V()[cells[ix]];
    numerator += (weighted / particleVolume) * solidVolume;
    denominator += solidVolume;
  }
  return denominator > 0.0 ? numerator / (denominator * cMax) : 0.0;
}

// 固相 + 液相总锂量（体积积分，单位 mol）。用于守恒漂移检查。
scalar totalLithium(const volScalarField& cE, const volScalarField& eps,
                    const scalarField& cNeg, const scalarField& cPos,
                    const labelList& negCells, const labelList& posCells,
                    const RadialGrid& negGrid, const RadialGrid& posGrid,
                    scalar negPorosity, scalar posPorosity,
                    const fvMesh& mesh) {
  scalar total = 0.0;
  const auto addSolid = [&](const scalarField& c, const labelList& cells,
                            const RadialGrid& grid, scalar porosity) {
    scalar particleVolume = 0.0;
    for (label ir = 0; ir < grid.nR; ++ir) {
      particleVolume += grid.volume[ir];
    }
    forAll(cells, ix) {
      scalar weighted = 0.0;
      for (label ir = 0; ir < grid.nR; ++ir) {
        weighted += c[ix * grid.nR + ir] * grid.volume[ir];
      }
      total += (weighted / particleVolume) * (1.0 - porosity) *
               mesh.V()[cells[ix]];
    }
  };
  addSolid(cNeg, negCells, negGrid, negPorosity);
  addSolid(cPos, posCells, posGrid, posPorosity);

  forAll(cE, cell) {
    total += eps[cell] * cE[cell] * mesh.V()[cell];
  }
  return total;
}

}  // namespace

int main(int argc, char* argv[]) {
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

  Info << "Reading constant/dfnProperties\n" << endl;
  IOdictionary dfnDict
  (
      IOobject
      (
          "dfnProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ_IF_MODIFIED,
          IOobject::NO_WRITE
      )
  );

  ModelParams params;
  params.temperature = readScalar(dfnDict.subDict("cell").lookup("temperature"));
  params.neg = readElectrode(dfnDict.subDict("negative"));
  params.pos = readElectrode(dfnDict.subDict("positive"));

  const dictionary& sepDict = dfnDict.subDict("separator");
  params.sepThickness = readScalar(sepDict.lookup("thickness"));
  params.sepPorosity = readScalar(sepDict.lookup("porosity"));
  params.sepBruggeman = readScalar(sepDict.lookup("bruggeman"));

  const dictionary& elyteDict = dfnDict.subDict("electrolyte");
  params.ce0 = readScalar(elyteDict.lookup("c0"));
  params.dE = readScalar(elyteDict.lookup("d_e"));
  params.tPlus = readScalar(elyteDict.lookup("t_plus"));
  params.kappa = readScalar(elyteDict.lookup("kappa"));

  const dictionary& controlDict = dfnDict.subDict("control");
  params.currentDensity = readScalar(controlDict.lookup("current_density"));
  params.maxPicard = controlDict.lookupOrDefault<label>("max_picard", 30);
  params.picardTol = controlDict.lookupOrDefault<scalar>("picard_tol", 1e-4);
  params.relaxPhi = controlDict.lookupOrDefault<scalar>("relax_phi", 1.0);
  params.relaxCe = controlDict.lookupOrDefault<scalar>("relax_ce", 1.0);
  params.voltageCutoff =
      controlDict.lookupOrDefault<scalar>("voltage_cutoff", 0.0);
  const bool debug = controlDict.lookupOrDefault<bool>("debug", false);

  // Bruggeman 修正与比表面积 a_s = 3(1-eps)/R。
  params.neg.sigmaEff = params.neg.sigma *
                        std::pow(1.0 - params.neg.porosity,
                                 params.neg.bruggeman);
  params.pos.sigmaEff = params.pos.sigma *
                        std::pow(1.0 - params.pos.porosity,
                                 params.pos.bruggeman);
  params.neg.aS =
      3.0 * (1.0 - params.neg.porosity) / params.neg.radius;
  params.pos.aS =
      3.0 * (1.0 - params.pos.porosity) / params.pos.radius;
  params.neg.alphaSurf = params.neg.radius /
                         static_cast<scalar>(params.neg.nRadial) /
                         (2.0 * params.neg.dS);
  params.pos.alphaSurf = params.pos.radius /
                         static_cast<scalar>(params.pos.nRadial) /
                         (2.0 * params.pos.dS);

  // 按单元中心所在区间划分三个区域。网格由 blockMesh 生成，顺序沿 x 递增。
  const scalar lNeg = params.neg.thickness;
  const scalar lSep = params.sepThickness;
  const scalar lPos = params.pos.thickness;
  labelList negCells;
  labelList sepCells;
  labelList posCells;
  labelList negLocal(mesh.nCells(), -1);
  labelList posLocal(mesh.nCells(), -1);
  forAll(mesh.C(), celli) {
    const scalar x = mesh.C()[celli].x();
    if (x < lNeg) {
      negLocal[celli] = negCells.size();
      negCells.append(celli);
    } else if (x < lNeg + lSep) {
      sepCells.append(celli);
    } else {
      posLocal[celli] = posCells.size();
      posCells.append(celli);
    }
  }
  if (negCells.empty() || sepCells.empty() || posCells.empty()) {
    FatalErrorInFunction
        << "Mesh does not contain all three regions; check blockMesh and "
        << "dfnProperties thicknesses" << exit(FatalError);
  }

  Info << "Regions: negative " << negCells.size() << " cells, separator "
       << sepCells.size() << " cells, positive " << posCells.size()
       << " cells" << endl;

  // 字段：读 0/ 中的边界条件，内部初值在 t=0 时按参数重设。
  volScalarField cE
  (
      IOobject("c_e", runTime.timeName(), mesh, IOobject::MUST_READ,
               IOobject::AUTO_WRITE),
      mesh
  );
  volScalarField phiE
  (
      IOobject("phi_e", runTime.timeName(), mesh, IOobject::MUST_READ,
               IOobject::AUTO_WRITE),
      mesh
  );
  volScalarField phiSNeg
  (
      IOobject("phi_s_neg", runTime.timeName(), mesh, IOobject::MUST_READ,
               IOobject::AUTO_WRITE),
      mesh
  );
  volScalarField phiSPos
  (
      IOobject("phi_s_pos", runTime.timeName(), mesh, IOobject::MUST_READ,
               IOobject::AUTO_WRITE),
      mesh
  );
  volScalarField cSurfNeg
  (
      IOobject("c_s_surface_neg", runTime.timeName(), mesh,
               IOobject::NO_READ, IOobject::AUTO_WRITE),
      mesh,
      dimensionedScalar("c_s_surface_neg", dimConc, 0.0)
  );
  volScalarField cSurfPos
  (
      IOobject("c_s_surface_pos", runTime.timeName(), mesh,
               IOobject::NO_READ, IOobject::AUTO_WRITE),
      mesh,
      dimensionedScalar("c_s_surface_pos", dimConc, 0.0)
  );

  // 材料属性场：孔隙率、比表面积、有效扩散系数与有效电导率。
  volScalarField eps
  (
      IOobject("eps", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("eps", dimless, 0.0)
  );
  volScalarField aS
  (
      IOobject("aS", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("aS", dimInvLength, 0.0)
  );
  volScalarField dEff
  (
      IOobject("dEff", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("dEff", dimDiffusivity, 0.0)
  );
  volScalarField kappaEff
  (
      IOobject("kappaEff", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("kappaEff", dimSigma, 0.0)
  );
  volScalarField sigmaNeg
  (
      IOobject("sigmaNeg", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("sigmaNeg", dimSigma, 0.0)
  );
  volScalarField sigmaPos
  (
      IOobject("sigmaPos", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("sigmaPos", dimSigma, 0.0)
  );
  // phi_s 只在各自电极内有方程；其余区域（隔膜与另一电极）用单位阵钉在 0。
  // 交界面上的面电导率也置 0，避免电子通路穿过隔膜。
  volScalarField maskNeg
  (
      IOobject("maskNeg", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("maskNeg", dimVolCurrentPerVolt, 1.0)
  );
  volScalarField maskPos
  (
      IOobject("maskPos", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("maskPos", dimVolCurrentPerVolt, 1.0)
  );

  forAll(negCells, i) {
    const label cell = negCells[i];
    eps[cell] = params.neg.porosity;
    aS[cell] = params.neg.aS;
    dEff[cell] = params.dE * std::pow(params.neg.porosity,
                                     params.neg.bruggeman);
    kappaEff[cell] = params.kappa *
                     std::pow(params.neg.porosity, params.neg.bruggeman);
    sigmaNeg[cell] = params.neg.sigmaEff;
    maskNeg[cell] = 0.0;
  }
  forAll(sepCells, i) {
    const label cell = sepCells[i];
    eps[cell] = params.sepPorosity;
    dEff[cell] = params.dE * std::pow(params.sepPorosity,
                                     params.sepBruggeman);
    kappaEff[cell] = params.kappa *
                     std::pow(params.sepPorosity, params.sepBruggeman);
  }
  forAll(posCells, i) {
    const label cell = posCells[i];
    eps[cell] = params.pos.porosity;
    aS[cell] = params.pos.aS;
    dEff[cell] = params.dE * std::pow(params.pos.porosity,
                                     params.pos.bruggeman);
    kappaEff[cell] = params.kappa *
                     std::pow(params.pos.porosity, params.pos.bruggeman);
    sigmaPos[cell] = params.pos.sigmaEff;
    maskPos[cell] = 0.0;
  }
  syncBoundary(eps);
  syncBoundary(aS);
  syncBoundary(dEff);
  syncBoundary(kappaEff);
  syncBoundary(sigmaNeg);
  syncBoundary(sigmaPos);
  syncBoundary(maskNeg);
  syncBoundary(maskPos);

  // phi_s 的有效面电导率：在电极/隔膜交界面上强制置 0，隔断电子通路。
  surfaceScalarField sigmaNegf
  (
      IOobject("sigmaNegf", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      fvc::interpolate(sigmaNeg)
  );
  surfaceScalarField sigmaPosf
  (
      IOobject("sigmaPosf", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      fvc::interpolate(sigmaPos)
  );
  for (label facei = 0; facei < mesh.nInternalFaces(); ++facei) {
    const label own = mesh.owner()[facei];
    const label nei = mesh.neighbour()[facei];
    if (negLocal[own] < 0 || negLocal[nei] < 0) {
      sigmaNegf.primitiveFieldRef()[facei] = 0.0;
    }
    if (posLocal[own] < 0 || posLocal[nei] < 0) {
      sigmaPosf.primitiveFieldRef()[facei] = 0.0;
    }
  }
  // 边界面的面电导率显式取相邻单元值，保证集流体上的固定通量真正生效。
  forAll(sigmaNegf.boundaryField(), patchi) {
    sigmaNegf.boundaryFieldRef()[patchi] ==
        sigmaNeg.boundaryField()[patchi].patchInternalField();
  }
  forAll(sigmaPosf.boundaryField(), patchi) {
    sigmaPosf.boundaryFieldRef()[patchi] ==
        sigmaPos.boundaryField()[patchi].patchInternalField();
  }

  // 正极集流体用固定电压边界；电压值在时间步内反解，使总电流等于目标值。
  const label posCollector = mesh.boundaryMesh().findPatchID("posCollector");
  const label negCollector = mesh.boundaryMesh().findPatchID("negCollector");
  if (posCollector < 0 || negCollector < 0) {
    FatalErrorInFunction
        << "Patches negCollector/posCollector are required"
        << exit(FatalError);
  }
  fixedValueFvPatchScalarField& posValue =
      refCast<fixedValueFvPatchScalarField>(
          phiSPos.boundaryFieldRef()[posCollector]);
  scalar posVoltage = 0.0;
  scalar collectorArea = 0.0;
  forAll(mesh.magSf().boundaryField()[posCollector], facei) {
    collectorArea += mesh.magSf().boundaryField()[posCollector][facei];
  }

  // 颗粒径向子网格与浓度数组（每个 x 单元一列）。
  const RadialGrid negGrid(params.neg.nRadial, params.neg.radius,
                           params.neg.dS);
  const RadialGrid posGrid(params.pos.nRadial, params.pos.radius,
                           params.pos.dS);
  scalarField cNeg(negCells.size() * negGrid.nR, params.neg.c0);
  scalarField cPos(posCells.size() * posGrid.nR, params.pos.c0);

  // t=0 的一致初值：c 均匀，phi_e=-U_neg，phi_s,pos=U_pos-U_neg。
  const scalar stoichNeg0 = params.neg.c0 / params.neg.cMax;
  const scalar stoichPos0 = params.pos.c0 / params.pos.cMax;
  const scalar uNeg0 = negativeOcp(stoichNeg0).value;
  const scalar uPos0 = positiveOcp(stoichPos0).value;
  if (runTime.timeIndex() == 0) {
    cE.primitiveFieldRef() = params.ce0;
    phiE.primitiveFieldRef() = -uNeg0;
    phiSNeg.primitiveFieldRef() = 0.0;
    phiSPos.primitiveFieldRef() = uPos0 - uNeg0;
    cE.correctBoundaryConditions();
    phiE.correctBoundaryConditions();
    phiSNeg.correctBoundaryConditions();
    phiSPos.correctBoundaryConditions();
    posVoltage = uPos0 - uNeg0;
    posValue == posVoltage;
    posValue.evaluate();
  }

  // 端电压：正极集流体处 φ_s 再沿梯度外推半个单元，与自研求解器定义一致。
  label lastPosCell = posCells[posCells.size() - 1];
  forAll(posCells, i) {
    if (mesh.C()[posCells[i]].x() > mesh.C()[lastPosCell].x()) {
      lastPosCell = posCells[i];
    }
  }
  // 端电压就是正极集流体的 Dirichlet 电位。
  const auto terminalVoltage = [&]() { return posVoltage; };

  const scalar lithium0 =
      totalLithium(cE, eps, cNeg, cPos, negCells, posCells, negGrid, posGrid,
                   params.neg.porosity, params.pos.porosity, mesh);

  // 输出时间序列，便于与 dfn-sim / labs 的结果直接比较。
  mkDir("postProcessing/dfnFoam");
  std::ofstream csv("postProcessing/dfnFoam/terminalVoltage.csv");
  csv << "t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,li_drift,"
         "picard_iter\n";
  csv << std::setprecision(10);

  const auto writeRow = [&](scalar t, scalar voltage, scalar stoichNeg,
                            scalar stoichPos, scalar drift, label picard) {
    csv << t << ',' << voltage << ',' << params.currentDensity << ','
        << stoichNeg << ',' << stoichPos << ',' << drift << ',' << picard
        << '\n';
  };

  const scalar stoichNeg = averageStoichiometry(
      cNeg, negCells, negGrid, mesh, params.neg.porosity, params.neg.cMax);
  const scalar stoichPos = averageStoichiometry(
      cPos, posCells, posGrid, mesh, params.pos.porosity, params.pos.cMax);
  Info << "Initial voltage " << terminalVoltage() << " V, stoich "
       << stoichNeg << "/" << stoichPos << endl;
  writeRow(0.0, terminalVoltage(), stoichNeg, stoichPos, 0.0, 0);

  // 每个时间步内用 Picard 迭代耦合：反应 -> 固相扩散 -> 电解液 -> 电势。
  scalarField srcCE(mesh.nCells(), 0.0);
  scalarField srcPhiNeg(mesh.nCells(), 0.0);
  scalarField srcPhiPos(mesh.nCells(), 0.0);
  scalarField srcPhiE(mesh.nCells(), 0.0);
  scalarField coefNeg(mesh.nCells(), 0.0);
  scalarField coefPos(mesh.nCells(), 0.0);
  scalarField coefE(mesh.nCells(), 0.0);
  volScalarField srcCEField
  (
      IOobject("srcCE", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("srcCE", dimConcRate, 0.0)
  );
  volScalarField srcPhiNegField
  (
      IOobject("srcPhiNeg", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("srcPhiNeg", dimVolCurrent, 0.0)
  );
  volScalarField srcPhiPosField
  (
      IOobject("srcPhiPos", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("srcPhiPos", dimVolCurrent, 0.0)
  );
  volScalarField srcPhiEField
  (
      IOobject("srcPhiE", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("srcPhiE", dimVolCurrent, 0.0)
  );
  volScalarField coefNegField
  (
      IOobject("coefNeg", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("coefNeg", dimVolCurrentPerVolt, 0.0)
  );
  volScalarField coefPosField
  (
      IOobject("coefPos", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("coefPos", dimVolCurrentPerVolt, 0.0)
  );
  volScalarField coefEField
  (
      IOobject("coefE", runTime.timeName(), mesh, IOobject::NO_READ,
               IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("coefE", dimVolCurrentPerVolt, 0.0)
  );

  std::vector<ReactionValue> reactionNeg(negCells.size());
  std::vector<ReactionValue> reactionPos(posCells.size());
  std::vector<scalar> jNeg(negCells.size(), 0.0);
  std::vector<scalar> jPos(posCells.size(), 0.0);

  // 假设正极固相电位均匀，计算给定集流体电压下的总电流。用于反解恒流电压。
  const auto positiveCurrentAtVoltage = [&](scalar voltage) {
    scalar current = 0.0;
    forAll(posCells, i) {
      const label cell = posCells[i];
      const scalar cLast = cPos[i * posGrid.nR + posGrid.nR - 1];
      const ReactionValue r = solveInterface(
          params, params.pos, true, cE[cell], cLast, voltage, phiE[cell],
          0.0);
      // 放电时正极嵌锂，j<0；取负号后总电流为正，与 currentDensity 同号。
      current -= params.pos.aS * kFaraday * r.j * mesh.V()[cell];
    }
    return current / collectorArea;
  };

  // 初值阶段先用均匀电位假设反解一次集流体电压，给后续 Newton 修正一个
  // 合理起点；进入时间步后直接复用上一迭代的电压并做实际电流修正。
  {
    scalar low = 0.0;
    scalar high = 6.0;
    for (int rootIter = 0; rootIter < 50; ++rootIter) {
      const scalar mid = 0.5 * (low + high);
      if (positiveCurrentAtVoltage(mid) > params.currentDensity) {
        low = mid;
      } else {
        high = mid;
      }
    }
    posVoltage = 0.5 * (low + high);
    posValue == posVoltage;
    posValue.evaluate();
  }

  bool stopped = false;
  scalarField cNegOld(cNeg);
  scalarField cPosOld(cPos);
  scalarField cEOld(cE.primitiveField());
  while (runTime.run() && !stopped) {
    ++runTime;
    const scalar dt = runTime.deltaTValue();
    cNegOld = cNeg;
    cPosOld = cPos;
    cEOld = cE.primitiveField();
    label picardUsed = 0;
    scalar change = 1.0;
    scalar picardResidual = 1.0;

    for (label iter = 0; iter < params.maxPicard; ++iter) {
      picardUsed = iter + 1;

      // 把当前集流体电压写入边界；电压本身在迭代末尾按实际正极电流修正。
      posValue == posVoltage;
      posValue.evaluate();

      std::vector<scalar> jNewNeg(negCells.size(), 0.0);
      std::vector<scalar> jNewPos(posCells.size(), 0.0);
      forAll(negCells, i) {
        const label cell = negCells[i];
        const scalar cLast = cNeg[i * negGrid.nR + negGrid.nR - 1];
        reactionNeg[i] = solveInterface(
            params, params.neg, false, cE[cell], cLast, phiSNeg[cell],
            phiE[cell], jNeg[i]);
        jNewNeg[i] = reactionNeg[i].j;
      }
      forAll(posCells, i) {
        const label cell = posCells[i];
        const scalar cLast = cPos[i * posGrid.nR + posGrid.nR - 1];
        reactionPos[i] = solveInterface(
            params, params.pos, true, cE[cell], cLast, phiSPos[cell],
            phiE[cell], jPos[i]);
        jNewPos[i] = reactionPos[i].j;
      }

      scalar deltaJ = 0.0;
      scalar scaleJ = 0.0;
      forAll(negCells, i) {
        deltaJ = std::max(deltaJ, std::abs(jNewNeg[i] - jNeg[i]));
        scaleJ = std::max(scaleJ, std::abs(jNewNeg[i]));
      }
      forAll(posCells, i) {
        deltaJ = std::max(deltaJ, std::abs(jNewPos[i] - jPos[i]));
        scaleJ = std::max(scaleJ, std::abs(jNewPos[i]));
      }
      change = deltaJ / (scaleJ + 1e-10);
      jNeg.swap(jNewNeg);
      jPos.swap(jNewPos);

      solveRadial(negGrid, cNegOld, cNeg, jNeg, dt);
      solveRadial(posGrid, cPosOld, cPos, jPos, dt);

      if (debug && runTime.timeIndex() <= 1 && iter == 0) {
        Info << "  radial check: cNeg first/last=" << cNeg[0] << "/"
             << cNeg[negGrid.nR - 1] << " cPos first/last="
             << cPos[0] << "/" << cPos[posGrid.nR - 1]
             << " negR=" << negGrid.nR << " negDr=" << negGrid.dr
             << " negV0=" << negGrid.volume[0]
             << " negA1=" << negGrid.area[1] << endl;
      }

      // 电解液质量守恒：eps dc/dt - div(D_eff grad c) = (1-t+) a_s j。
      for (label cell = 0; cell < mesh.nCells(); ++cell) {
        srcCE[cell] = 0.0;
      }
      forAll(negCells, i) {
        srcCE[negCells[i]] =
            (1.0 - params.tPlus) * params.neg.aS * jNeg[i];
      }
      forAll(posCells, i) {
        srcCE[posCells[i]] =
            (1.0 - params.tPlus) * params.pos.aS * jPos[i];
      }
      srcCEField.primitiveFieldRef() = srcCE;
      syncBoundary(srcCEField);
      // fvm::ddt 的旧时间层会被 fvMatrix::solve 内部刷新；Picard 迭代中
      // 每次装配前恢复步初值，使每次求解都是同一后向欧拉方程的不动点迭代。
      cE.oldTime().primitiveFieldRef() = cEOld;
      syncBoundary(cE.oldTime());
      fvScalarMatrix ceEqn
      (
          fvm::ddt(eps, cE) - fvm::laplacian(dEff, cE) == srcCEField
      );
      ceEqn.relax(params.relaxCe);
      ceEqn.solve();
      cE.correctBoundaryConditions();

      // 固相电荷守恒：div(sigma_eff grad phi_s) = a_s F j。
      // 用 j 对 eta 的偏导做隐式线性化，保证大 i0 下的稳定性。
      for (label cell = 0; cell < mesh.nCells(); ++cell) {
        srcPhiNeg[cell] = 0.0;
        srcPhiPos[cell] = 0.0;
        coefNeg[cell] = 0.0;
        coefPos[cell] = 0.0;
        srcPhiE[cell] = 0.0;
        coefE[cell] = 0.0;
      }
      forAll(negCells, i) {
        const label cell = negCells[i];
        const ReactionValue& r = reactionNeg[i];
        const scalar coupling = params.neg.aS * kFaraday;
        const scalar dEta = r.dj_dEtaEff;
        coefNeg[cell] = coupling * dEta;
        srcPhiNeg[cell] =
            coupling *
            (r.j - dEta * (phiE[cell] + r.ocp + r.eta));
        coefE[cell] = coupling * dEta;
        srcPhiE[cell] =
            -coupling *
            (r.j + dEta * (phiSNeg[cell] - r.ocp - r.eta));
      }
      forAll(posCells, i) {
        const label cell = posCells[i];
        const ReactionValue& r = reactionPos[i];
        const scalar coupling = params.pos.aS * kFaraday;
        const scalar dEta = r.dj_dEtaEff;
        coefPos[cell] = coupling * dEta;
        srcPhiPos[cell] =
            coupling *
            (r.j - dEta * (phiE[cell] + r.ocp + r.eta));
        coefE[cell] += coupling * dEta;
        srcPhiE[cell] +=
            -coupling *
            (r.j + dEta * (phiSPos[cell] - r.ocp - r.eta));
      }
      srcPhiNegField.primitiveFieldRef() = srcPhiNeg;
      srcPhiPosField.primitiveFieldRef() = srcPhiPos;
      srcPhiEField.primitiveFieldRef() = srcPhiE;
      coefNegField.primitiveFieldRef() = coefNeg;
      coefPosField.primitiveFieldRef() = coefPos;
      coefEField.primitiveFieldRef() = coefE;
      for (volScalarField* field : {&srcPhiNegField, &srcPhiPosField,
                                    &srcPhiEField, &coefNegField,
                                    &coefPosField, &coefEField}) {
        syncBoundary(*field);
      }

      // 负极：集流体 phi_s=0；隔膜内用单位阵钉住 0。
      phiSNeg.storePrevIter();
      fvScalarMatrix phiSNegEqn
      (
          fvm::laplacian(sigmaNegf, phiSNeg)
          - fvm::Sp(coefNegField, phiSNeg)
          - fvm::Sp(maskNeg, phiSNeg)
          == srcPhiNegField
      );
      phiSNegEqn.solve();
      phiSNeg.correctBoundaryConditions();
      phiSNeg.relax(params.relaxPhi);

      // 正极：集流体 fixedValue，电压由恒流反解给出；隔膜同样钉住 0。
      phiSPos.storePrevIter();
      fvScalarMatrix phiSPosEqn
      (
          fvm::laplacian(sigmaPosf, phiSPos)
          - fvm::Sp(coefPosField, phiSPos)
          - fvm::Sp(maskPos, phiSPos)
          == srcPhiPosField
      );
      if (debug && runTime.timeIndex() <= 2) {
        Info << "  phiSPos BC: sigmaf="
             << sigmaPosf.boundaryField()[posCollector][0]
             << " patchValue=" << posValue[0] << endl;
      }
      phiSPosEqn.solve();
      phiSPos.correctBoundaryConditions();
      phiSPos.relax(params.relaxPhi);

      // 液相电荷守恒：div(kappa_eff grad phi_e) + a_s F j = 0。
      phiE.storePrevIter();
      fvScalarMatrix phiEEqn
      (
          fvm::laplacian(kappaEff, phiE)
          - fvm::Sp(coefEField, phiE)
          == srcPhiEField
      );
      phiEEqn.solve();
      phiE.correctBoundaryConditions();
      phiE.relax(params.relaxPhi);

      // phi_e 的 Neumann 问题只确定梯度，整体平移量由正负极电流平衡决定。
      // 对平移量做内层 Newton 迭代，直到全电芯反应电流之和为零。
      scalar currentNeg = 0.0;
      scalar currentPos = 0.0;
      scalar dNegDPhiE = 0.0;
      scalar dPosDPhiS = 0.0;
      const scalar balanceScale =
          std::max(collectorArea * std::abs(params.currentDensity), 1e-30);
      for (int balanceIter = 0; balanceIter < 6; ++balanceIter) {
        currentNeg = 0.0;
        currentPos = 0.0;
        dNegDPhiE = 0.0;
        dPosDPhiS = 0.0;
        forAll(negCells, i) {
          const label cell = negCells[i];
          const scalar cLast = cNeg[i * negGrid.nR + negGrid.nR - 1];
          const ReactionValue r = solveInterface(
              params, params.neg, false, cE[cell], cLast, phiSNeg[cell],
              phiE[cell], jNeg[i]);
          currentNeg += params.neg.aS * kFaraday * r.j * mesh.V()[cell];
          dNegDPhiE -=
              params.neg.aS * kFaraday * r.dj_dEtaEff * mesh.V()[cell];
        }
        forAll(posCells, i) {
          const label cell = posCells[i];
          const scalar cLast = cPos[i * posGrid.nR + posGrid.nR - 1];
          const ReactionValue r = solveInterface(
              params, params.pos, true, cE[cell], cLast, phiSPos[cell],
              phiE[cell], jPos[i]);
          currentPos += params.pos.aS * kFaraday * r.j * mesh.V()[cell];
          dPosDPhiS +=
              params.pos.aS * kFaraday * r.dj_dEtaEff * mesh.V()[cell];
        }
        const scalar imbalance = currentNeg + currentPos;
        if (std::abs(imbalance) < 1e-10 * balanceScale) {
          break;
        }
        // d(imbalance)/d(phi_e) = dNegDPhiE - dPosDPhiS。
        const scalar derivative = dNegDPhiE - dPosDPhiS;
        if (std::abs(derivative) < 1e-30) {
          break;
        }
        phiE.primitiveFieldRef() -= imbalance / derivative;
        phiE.correctBoundaryConditions();
      }
      // 用平衡后的 phi_e 重新评估正极电流，供集流体电压修正使用。
      currentPos = 0.0;
      dPosDPhiS = 0.0;
      forAll(posCells, i) {
        const label cell = posCells[i];
        const scalar cLast = cPos[i * posGrid.nR + posGrid.nR - 1];
        const ReactionValue r = solveInterface(
            params, params.pos, true, cE[cell], cLast, phiSPos[cell],
            phiE[cell], jPos[i]);
        currentPos += params.pos.aS * kFaraday * r.j * mesh.V()[cell];
        dPosDPhiS +=
            params.pos.aS * kFaraday * r.dj_dEtaEff * mesh.V()[cell];
      }
      // 再用解出的正极总反应电流修正集流体电压，使实际电流（而非
      // 均匀电位假设下的电流）等于 targetPos。放电时 targetPos < 0。
      const scalar targetPos = -params.currentDensity * collectorArea;
      if (dPosDPhiS > 1e-30) {
        posVoltage += 0.7 * (targetPos - currentPos) / dPosDPhiS;
        posValue == posVoltage;
        posValue.evaluate();
      }

      // 电极之外的 φ_s 只是隔膜/另一电极里的占位未知量；显式清零，
      // 避免欠松弛留下的非物理残余进入诊断输出。
      forAll(mesh.C(), celli) {
        if (negLocal[celli] < 0) {
          phiSNeg[celli] = 0.0;
        }
        if (posLocal[celli] < 0) {
          phiSPos[celli] = 0.0;
        }
      }

      // 电位的更新量（相对场欠松弛前的值）作为 Picard 收敛判据；j 的相对
      // 变化在 j 接近零时会过于苛刻，不能反映场的实际收敛程度。
      scalar potentialChange = 0.0;
      forAll(mesh.C(), celli) {
        potentialChange = std::max(
            potentialChange, std::abs(phiSNeg[celli] -
                                      phiSNeg.prevIter()[celli]));
        potentialChange = std::max(
            potentialChange, std::abs(phiSPos[celli] -
                                      phiSPos.prevIter()[celli]));
        potentialChange = std::max(
            potentialChange,
            std::abs(phiE[celli] - phiE.prevIter()[celli]));
      }
      picardResidual = potentialChange;

      if (debug && runTime.timeIndex() <= 2) {
        scalar jMax = 0.0;
        forAll(negCells, i) {
          jMax = std::max(jMax, std::abs(jNeg[i]));
        }
        forAll(posCells, i) {
          jMax = std::max(jMax, std::abs(jPos[i]));
        }
        auto fieldRange = [](const volScalarField& field) {
          scalar lo = 1e300;
          scalar hi = -1e300;
          forAll(field, celli) {
            lo = std::min(lo, field[celli]);
            hi = std::max(hi, field[celli]);
          }
          return std::make_pair(lo, hi);
        };
        scalar cLastMin = 1e300;
        scalar cLastMax = -1e300;
        forAll(negCells, i) {
          const scalar value = cNeg[i * negGrid.nR + negGrid.nR - 1];
          cLastMin = std::min(cLastMin, value);
          cLastMax = std::max(cLastMax, value);
        }
        forAll(posCells, i) {
          const scalar value = cPos[i * posGrid.nR + posGrid.nR - 1];
          cLastMin = std::min(cLastMin, value);
          cLastMax = std::max(cLastMax, value);
        }
        const auto ceRange = fieldRange(cE);
        const auto pnRange = fieldRange(phiSNeg);
        const auto ppRange = fieldRange(phiSPos);
        const auto peRange = fieldRange(phiE);
        Info << "  iter " << iter << " jMax=" << jMax
             << " dJ=" << change << " ce=[" << ceRange.first << ","
             << ceRange.second << "] cLast=[" << cLastMin << ","
             << cLastMax << "] phiSneg=[" << pnRange.first << ","
             << pnRange.second << "] phiSpos=[" << ppRange.first << ","
             << ppRange.second << "] phiE=[" << peRange.first << ","
             << peRange.second << "]" << endl;
        Info << "  pos0: cell=" << posCells[0]
             << " phiS=" << phiSPos[posCells[0]]
             << " phiE=" << phiE[posCells[0]]
             << " cLast=" << cPos[posGrid.nR - 1]
             << " eta=" << reactionPos[0].eta
             << " ocp=" << reactionPos[0].ocp
             << " j=" << reactionPos[0].j
             << " coef=" << coefPos[posCells[0]] << endl;
        const label lastIdx = posCells.size() - 1;
        Info << "  posN: cell=" << posCells[lastIdx]
             << " phiS=" << phiSPos[posCells[lastIdx]]
             << " phiE=" << phiE[posCells[lastIdx]]
             << " cLast=" << cPos[lastIdx * posGrid.nR + posGrid.nR - 1]
             << " eta=" << reactionPos[lastIdx].eta
             << " ocp=" << reactionPos[lastIdx].ocp
             << " j=" << reactionPos[lastIdx].j
             << " coef=" << coefPos[posCells[lastIdx]] << endl;
        Info << "  posProfile:";
        forAll(posCells, i) {
          Info << " " << phiSPos[posCells[i]];
        }
        Info << " Vcell=" << posVoltage << endl;
        Info << "  masks: negCell0 phiSPos=" << phiSPos[negCells[0]]
             << " maskPos=" << maskPos[negCells[0]]
             << " ; posCell0 phiSNeg=" << phiSNeg[posCells[0]]
             << " maskNeg=" << maskNeg[posCells[0]] << endl;
        Info << "  phiEcell: neg0 phiE=" << phiE[negCells[0]]
             << " src=" << srcPhiE[negCells[0]]
             << " coef=" << coefE[negCells[0]]
             << " dEta=" << reactionNeg[0].dj_dEtaEff
             << " dCs=" << reactionNeg[0].dj_dCs
             << " ; pos0 phiE=" << phiE[posCells[0]]
             << " src=" << srcPhiE[posCells[0]]
             << " coef=" << coefE[posCells[0]]
             << " dEta=" << reactionPos[0].dj_dEtaEff
             << " dCs=" << reactionPos[0].dj_dCs << endl;
      }

      if (iter > 0 && picardResidual < params.picardTol) {
        break;
      }
    }

    // 用最终电位重算一次表面浓度，写诊断场并对齐输出。
    forAll(negCells, i) {
      const label cell = negCells[i];
      const scalar cLast = cNeg[i * negGrid.nR + negGrid.nR - 1];
      reactionNeg[i] = solveInterface(
          params, params.neg, false, cE[cell], cLast, phiSNeg[cell],
          phiE[cell], jNeg[i]);
      cSurfNeg[cell] = reactionNeg[i].cSurf;
      jNeg[i] = reactionNeg[i].j;
    }
    forAll(posCells, i) {
      const label cell = posCells[i];
      const scalar cLast = cPos[i * posGrid.nR + posGrid.nR - 1];
      reactionPos[i] = solveInterface(
          params, params.pos, true, cE[cell], cLast, phiSPos[cell],
          phiE[cell], jPos[i]);
      cSurfPos[cell] = reactionPos[i].cSurf;
      jPos[i] = reactionPos[i].j;
    }

    const scalar voltage = terminalVoltage();
    if (debug) {
      Info << "  terminal: lastCell=" << lastPosCell
           << " phiLast=" << phiSPos[lastPosCell]
           << " posVoltage=" << posVoltage
           << " posPatchValue="
           << phiSPos.boundaryField()[posCollector][0] << endl;
    }
    const scalar sNeg = averageStoichiometry(
        cNeg, negCells, negGrid, mesh, params.neg.porosity, params.neg.cMax);
    const scalar sPos = averageStoichiometry(
        cPos, posCells, posGrid, mesh, params.pos.porosity, params.pos.cMax);
    const scalar lithium = totalLithium(
        cE, eps, cNeg, cPos, negCells, posCells, negGrid, posGrid,
        params.neg.porosity, params.pos.porosity, mesh);
    const scalar drift = std::abs(lithium - lithium0) / lithium0;

    Info << "t=" << runTime.timeName() << " V=" << voltage
         << " stoich=" << sNeg << "/" << sPos << " picard=" << picardUsed
         << " dPhi=" << picardResidual << " dJ=" << change
         << " LiDrift=" << drift << endl;
    writeRow(runTime.value(), voltage, sNeg, sPos, drift, picardUsed);

    if (picardUsed >= params.maxPicard &&
        picardResidual >= params.picardTol) {
      WarningInFunction << "Picard did not converge at t="
                        << runTime.value() << " (dPhi=" << picardResidual
                        << ")" << endl;
    }

    runTime.write();

    if (params.voltageCutoff > 0.0 && params.currentDensity > 0.0 &&
        voltage <= params.voltageCutoff) {
      Info << "Voltage cutoff " << params.voltageCutoff << " V reached at t="
           << runTime.value() << " s" << endl;
      stopped = true;
    }
  }

  Info << "End" << nl << endl;
  return 0;
}
