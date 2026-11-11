#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string_view>

#include "dfn/thermo.hpp"

// 实验（第 0 章 入门 / 第 3 章 开路电压）：按 SOC 等距采样，把正负极化学计量比
// 映射到半电池开路电位并作差，得到整电池 OCV 曲线。
// 输出 CSV 列：soc（0..1）, sto_neg x, sto_pos y, ocv_v（V）。
namespace {

// 由 SOC 线性映射到负极锂化程度 x 与正极锂化程度 y，返回整电池 OCV：
// U_ocv = U_pos(y) - U_neg(x)，单位 V（正极电位减负极电位）。
double open_circuit_voltage(double soc) {
  const double x = 0.10 + 0.75 * soc;
  const double y = 0.90 - 0.80 * soc;
  return dfn::positive_ocp(y).value - dfn::negative_ocp(x).value;
}

// 打印命令行用法说明。
void print_usage() {
  std::cout << "Usage: ocv_curve [--n POINTS]\n"
               "Print soc, electrode stoichiometries and OCV as CSV.\n";
}

}  // namespace

int main(int argc, char** argv) {
  // 默认采样 101 个点（SOC 间隔 0.01）。
  int n = 101;

  // 解析命令行：--help 打印用法；--n 指定采样点数，必须 >= 2。
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      print_usage();
      return 0;
    }
    if (arg == "--n") {
      if (i + 1 >= argc) {
        std::cerr << "error: --n requires a value\n";
        return 2;
      }
      n = std::atoi(argv[++i]);
      if (n < 2) {
        std::cerr << "error: --n must be at least 2\n";
        return 2;
      }
      continue;
    }
    std::cerr << "error: unknown option " << arg << "\n";
    return 2;
  }

  std::cout << "soc,sto_neg,sto_pos,ocv_v\n"
            << std::fixed << std::setprecision(6);

  // SOC 在 [0, 1] 上均匀取 n 个点：soc_k = k / (n - 1)。
  for (int k = 0; k < n; ++k) {
    const double soc = static_cast<double>(k) / static_cast<double>(n - 1);
    // 随 SOC 升高负极嵌锂（x 增大），正极脱锂（y 减小），这里取线性近似。
    const double x = 0.10 + 0.75 * soc;
    const double y = 0.90 - 0.80 * soc;
    std::cout << soc << ',' << x << ',' << y << ','
              << open_circuit_voltage(soc) << '\n';
  }

  return 0;
}
