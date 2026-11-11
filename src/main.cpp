// 命令行主程序 dfn-sim：解析参数、覆盖 SimulationConfig 的默认值，
// 调用 run_simulation 并把时间序列结果写成 CSV（stdout 或文件），
// 最后在 stderr 上打印一行摘要。退出码：0 成功，1 仿真失败，2 参数或 I/O 错误。

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

#include "dfn/sim/config.hpp"
#include "dfn/sim/simulation.hpp"
#include "dfn/version.hpp"

namespace {

// 打印用法说明：逐条列出选项含义与默认值（这些默认值与 SimulationConfig 的
// 成员初值一致，改动其中一处时不要忘记同步另一处）。
void print_help() {
  std::cout
      << "dfn-sim " << dfn::kVersion << "\n\n"
      << "Usage: dfn-sim [options]\n\n"
      << "Options:\n"
      << "  -h, --help              Show this message\n"
      << "      --version           Show the version\n"
      << "      --config PATH       Load parameters (see configs/demo.yaml)\n"
      << "      --model NAME        spm, spme or dfn (default)\n"
      << "      --mode NAME         cc, cv, cccv (default) or pulse\n"
      << "                          cv holds --hold-voltage and needs an\n"
      << "                          initial state at that voltage\n"
      << "      --crate R           C-rate, 1C = 36 A/m^2 in this demo cell\n"
      << "      --current J         Current density in A/m^2 (overrides --crate)\n"
      << "      --charge            Charge instead of discharge (negative current)\n"
      << "      --v-min V           Discharge cutoff (default 3.5 V)\n"
      << "      --v-max V           Charge cutoff (default 4.0 V)\n"
      << "      --hold-voltage V    CV setpoint (default: the cutoff)\n"
      << "      --t-max S           End of the CC phase (default 20000 s)\n"
      << "      --hold S            CV hold time after the cutoff (default 600 s)\n"
      << "      --pulse-period S    Pulse period (default 600 s)\n"
      << "      --pulse-duty D      Fraction of the period on (default 0.5)\n"
      << "      --dt S              Output interval (default 60 s)\n"
      << "      --nx NEG SEP POS    x cells per region (default 8 4 8)\n"
      << "      --nr NEG POS        radial cells (default 6 6)\n"
      << "      --rtol R            Relative tolerance (default 1e-6)\n"
      << "      --out PATH          Write CSV to PATH instead of stdout\n"
      << "      --quiet             Do not print the summary\n\n"
      << "Example:\n"
      << "  dfn-sim --model dfn --mode cccv --crate 1 --dt 60 --out results/cccv.csv\n"
      << "  dfn-sim --model spme --mode cc --crate 1 --dt 60\n"
      << "  dfn-sim --model spm --mode cc --crate 1 --charge --v-max 4.0 --dt 10\n";
}

// 解析完整的 double：要求整段文本都被消费，否则返回 false（拒绝 "1.0x" 之类）。
bool parse_double(const char* text, double* value) {
  char* end = nullptr;
  const double parsed = std::strtod(text, &end);
  if (end == text || *end != '\0') {
    return false;
  }
  *value = parsed;
  return true;
}

// 解析 int：先按 double 解析再取整，若取整后与原值不等（如 "1.5"）则视为失败。
bool parse_int(const char* text, int* value) {
  double parsed = 0.0;
  if (!parse_double(text, &parsed)) {
    return false;
  }
  *value = static_cast<int>(parsed);
  return static_cast<double>(*value) == parsed;
}

// 把字符串映射为工况模式枚举：cc 恒流、cv 恒压、cccv 恒流恒压、pulse 脉冲。
bool parse_mode(std::string_view text, dfn::SimulationMode* mode) {
  if (text == "cc") {
    *mode = dfn::SimulationMode::cc;
  } else if (text == "cv") {
    *mode = dfn::SimulationMode::cv;
  } else if (text == "cccv") {
    *mode = dfn::SimulationMode::cccv;
  } else if (text == "pulse") {
    *mode = dfn::SimulationMode::pulse;
  } else {
    return false;
  }
  return true;
}

// 把字符串映射为模型种类枚举：spm / spme / dfn（分辨率由低到高）。
bool parse_model(std::string_view text, dfn::ModelKind* model) {
  if (text == "spm") {
    *model = dfn::ModelKind::spm;
  } else if (text == "spme") {
    *model = dfn::ModelKind::spme;
  } else if (text == "dfn") {
    *model = dfn::ModelKind::dfn;
  } else {
    return false;
  }
  return true;
}

// 模型枚举转字符串，用于结尾的摘要输出。
const char* model_name(dfn::ModelKind model) {
  switch (model) {
    case dfn::ModelKind::spm:
      return "spm";
    case dfn::ModelKind::spme:
      return "spme";
    case dfn::ModelKind::dfn:
      return "dfn";
  }
  return "?";
}

// 工况模式枚举转字符串，用于结尾的摘要输出。
const char* mode_name(dfn::SimulationMode mode) {
  switch (mode) {
    case dfn::SimulationMode::cc:
      return "cc";
    case dfn::SimulationMode::cv:
      return "cv";
    case dfn::SimulationMode::cccv:
      return "cccv";
    case dfn::SimulationMode::pulse:
      return "pulse";
  }
  return "?";
}

}  // namespace

// 程序入口：从左到右扫描 argv，每个选项就地覆盖 config 的对应字段，
// 因此后出现的选项优先；未给出的选项保持 SimulationConfig 的默认值。
// 参数缺失或非法时打印 error 并以退出码 2 结束。
int main(int argc, char** argv) {
  dfn::SimulationConfig config;
  std::string output_path;
  bool quiet = false;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    // 取下一个 argv 作为当前选项的值并把下标 i 前移；若已到末尾则报错退出。
    const auto value = [&](const char* name) -> const char* {
      if (i + 1 >= argc) {
        std::cerr << "error: " << name << " requires a value\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (arg == "-h" || arg == "--help") {
      print_help();
      return 0;
    }
    if (arg == "--version") {
      std::cout << dfn::kVersion << "\n";
      return 0;
    }
    // --config PATH：加载 YAML 参数文件（电芯、负极、正极、隔膜、电解液各段），
    // 建议放在最前面，其后的数值选项再逐项覆盖。
    if (arg == "--config") {
      try {
        config.params = dfn::load_params_yaml(value("--config"));
      } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 2;
      }
      continue;
    }
    // --model NAME：模型种类，默认 dfn（spm 最快、spme 居中、dfn 最完整）。
    if (arg == "--model") {
      if (!parse_model(value("--model"), &config.model)) {
        std::cerr << "error: --model must be spm, spme or dfn\n";
        return 2;
      }
      continue;
    }
    // --mode NAME：cc 恒流、cv 恒压、cccv 恒流后转恒压（默认）、pulse 脉冲；
    // cv 要求初始状态已处在 --hold-voltage 附近，否则会立刻判定截止。
    if (arg == "--mode") {
      if (!parse_mode(value("--mode"), &config.mode)) {
        std::cerr << "error: --mode must be cc, cv, cccv or pulse\n";
        return 2;
      }
      continue;
    }
    // --crate R：C 倍率，按 demo 电芯的 1C = 36 A/m^2 换算成电流密度。
    if (arg == "--crate") {
      if (!parse_double(value("--crate"), &config.c_rate)) {
        std::cerr << "error: --crate needs a number\n";
        return 2;
      }
      continue;
    }
    // --current J：直接指定电流密度（A/m^2），非零时优先于 --crate。
    if (arg == "--current") {
      if (!parse_double(value("--current"), &config.current_density)) {
        std::cerr << "error: --current needs a number\n";
        return 2;
      }
      continue;
    }
    // --charge：改为充电，施加电流取负号，截止条件切换到 --v-max。
    if (arg == "--charge") {
      config.charge = true;
      continue;
    }
    // --v-max V：充电截止电压，默认 4.0 V。
    if (arg == "--v-max") {
      if (!parse_double(value("--v-max"), &config.v_max)) {
        std::cerr << "error: --v-max needs a number\n";
        return 2;
      }
      continue;
    }
    // --v-min V：放电截止电压，默认 3.5 V。
    if (arg == "--v-min") {
      if (!parse_double(value("--v-min"), &config.v_min)) {
        std::cerr << "error: --v-min needs a number\n";
        return 2;
      }
      continue;
    }
    // --hold-voltage V：CV 阶段保持的电压，默认取当前工况对应的截止电压。
    if (arg == "--hold-voltage") {
      if (!parse_double(value("--hold-voltage"), &config.hold_voltage)) {
        std::cerr << "error: --hold-voltage needs a number\n";
        return 2;
      }
      continue;
    }
    // --t-max S：仿真时长上限，默认 20000 s；CC 阶段与 pulse 模式都用它封顶。
    if (arg == "--t-max") {
      if (!parse_double(value("--t-max"), &config.t_max)) {
        std::cerr << "error: --t-max needs a number\n";
        return 2;
      }
      continue;
    }
    // --hold S：达到截止条件后继续恒压的时间，默认 600 s。
    if (arg == "--hold") {
      if (!parse_double(value("--hold"), &config.t_hold)) {
        std::cerr << "error: --hold needs a number\n";
        return 2;
      }
      continue;
    }
    // --pulse-period S：pulse 模式的周期，默认 600 s。
    if (arg == "--pulse-period") {
      if (!parse_double(value("--pulse-period"), &config.t_pulse)) {
        std::cerr << "error: --pulse-period needs a number\n";
        return 2;
      }
      continue;
    }
    // --pulse-duty D：一个周期内施加电流的时间占比，默认 0.5。
    if (arg == "--pulse-duty") {
      if (!parse_double(value("--pulse-duty"), &config.pulse_duty)) {
        std::cerr << "error: --pulse-duty needs a number\n";
        return 2;
      }
      continue;
    }
    // --dt S：输出采样间隔，默认 60 s；不影响 IDA 内部的自适应步长。
    if (arg == "--dt") {
      if (!parse_double(value("--dt"), &config.output_interval)) {
        std::cerr << "error: --dt needs a number\n";
        return 2;
      }
      continue;
    }
    // --nx NEG SEP POS：负极 / 隔膜 / 正极沿 x 方向的网格数，默认 8 4 8。
    if (arg == "--nx") {
      if (!parse_int(value("--nx"), &config.n_x_neg) ||
          !parse_int(value("--nx"), &config.n_x_sep) ||
          !parse_int(value("--nx"), &config.n_x_pos)) {
        std::cerr << "error: --nx needs three integers\n";
        return 2;
      }
      continue;
    }
    // --nr NEG POS：负极 / 正极颗粒沿半径方向的网格数，默认 6 6。
    if (arg == "--nr") {
      if (!parse_int(value("--nr"), &config.n_r_neg) ||
          !parse_int(value("--nr"), &config.n_r_pos)) {
        std::cerr << "error: --nr needs two integers\n";
        return 2;
      }
      continue;
    }
    // --rtol R：传给 IDA 的相对误差容差，默认 1e-6；收紧会显著增加计算量。
    if (arg == "--rtol") {
      if (!parse_double(value("--rtol"), &config.rtol)) {
        std::cerr << "error: --rtol needs a number\n";
        return 2;
      }
      continue;
    }
    // --out PATH：把 CSV 写入文件；不给该选项时 CSV 写到 stdout。
    if (arg == "--out") {
      output_path = value("--out");
      continue;
    }
    // --quiet：关闭结尾的摘要输出（摘要本身写到 stderr）。
    if (arg == "--quiet") {
      quiet = true;
      continue;
    }
    // 未识别的选项直接报错，避免拼写错误被静默忽略。
    std::cerr << "error: unknown option " << arg << "\n";
    return 2;
  }

  // 运行仿真：内部按工况（cc/cv/cccv/pulse）分阶段调用 IDA 求解。
  const dfn::SimulationReport report = dfn::run_simulation(config);
  if (!report.ok) {
    std::cerr << "error: " << report.message << "\n";
    return 1;
  }

  // 写 CSV：默认 stdout（便于管道处理），给了 --out 则直接写文件；
  // 写文件失败（路径不存在、无权限等）按 I/O 错误返回退出码 2。
  try {
    if (output_path.empty()) {
      dfn::write_csv(report, std::cout);
    } else {
      dfn::write_csv(report, output_path);
    }
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 2;
  }

  // 摘要写到 stderr，保证 stdout 上的 CSV 不被混入文字：
  // 依次给出模型/工况/网格规模、未知量与非零元数、结束时电压电流、容量、
  // 锂守恒漂移、IDA 步数与 Jacobian 求值次数、实际耗时（毫秒）。
  if (!quiet && !report.samples.empty()) {
    const dfn::SimulationSample& last = report.samples.back();
    std::cerr << std::fixed << std::setprecision(4)
              << "model=" << model_name(config.model) << " mode="
              << mode_name(config.mode)
              << (config.charge ? "/charge" : "/discharge") << " grid="
              << config.n_x_neg << "/"
              << config.n_x_sep << "/" << config.n_x_pos << " nr="
              << config.n_r_neg << "/" << config.n_r_pos
              << " n=" << report.unknowns << " nnz=" << report.nnz
              << " samples=" << report.samples.size() << " t_end=" << report.t_end
              << " s V_end=" << last.voltage << " V I_end=" << last.current
              << " A/m^2 cap=" << report.capacity / 3600.0
              << " Ah/m^2 Li_drift=" << std::scientific << report.lithium_drift
              << std::defaultfloat << " steps=" << report.steps
              << " jac=" << report.jacobian_evals << " wall=" << report.wall_ms
              << " ms\n";
  }
  return 0;
}
