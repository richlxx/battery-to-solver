#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

#include "dfn/sim/config.hpp"
#include "dfn/sim/simulation.hpp"
#include "dfn/version.hpp"

namespace {

void print_help() {
  std::cout
      << "dfn-sim " << dfn::kVersion << "\n\n"
      << "Usage: dfn-sim [options]\n\n"
      << "Options:\n"
      << "  -h, --help              Show this message\n"
      << "      --version           Show the version\n"
      << "      --config PATH       Load parameters (see configs/demo.yaml)\n"
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
      << "  dfn-sim --mode cccv --crate 1 --dt 60 --out results/cccv.csv\n"
      << "  dfn-sim --mode cc --crate 1 --charge --v-max 4.0 --dt 10\n";
}

bool parse_double(const char* text, double* value) {
  char* end = nullptr;
  const double parsed = std::strtod(text, &end);
  if (end == text || *end != '\0') {
    return false;
  }
  *value = parsed;
  return true;
}

bool parse_int(const char* text, int* value) {
  double parsed = 0.0;
  if (!parse_double(text, &parsed)) {
    return false;
  }
  *value = static_cast<int>(parsed);
  return static_cast<double>(*value) == parsed;
}

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

int main(int argc, char** argv) {
  dfn::SimulationConfig config;
  std::string output_path;
  bool quiet = false;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
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
    if (arg == "--config") {
      try {
        config.params = dfn::load_params_yaml(value("--config"));
      } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 2;
      }
      continue;
    }
    if (arg == "--mode") {
      if (!parse_mode(value("--mode"), &config.mode)) {
        std::cerr << "error: --mode must be cc, cv, cccv or pulse\n";
        return 2;
      }
      continue;
    }
    if (arg == "--crate") {
      if (!parse_double(value("--crate"), &config.c_rate)) {
        std::cerr << "error: --crate needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--current") {
      if (!parse_double(value("--current"), &config.current_density)) {
        std::cerr << "error: --current needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--charge") {
      config.charge = true;
      continue;
    }
    if (arg == "--v-max") {
      if (!parse_double(value("--v-max"), &config.v_max)) {
        std::cerr << "error: --v-max needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--v-min") {
      if (!parse_double(value("--v-min"), &config.v_min)) {
        std::cerr << "error: --v-min needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--hold-voltage") {
      if (!parse_double(value("--hold-voltage"), &config.hold_voltage)) {
        std::cerr << "error: --hold-voltage needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--t-max") {
      if (!parse_double(value("--t-max"), &config.t_max)) {
        std::cerr << "error: --t-max needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--hold") {
      if (!parse_double(value("--hold"), &config.t_hold)) {
        std::cerr << "error: --hold needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--pulse-period") {
      if (!parse_double(value("--pulse-period"), &config.t_pulse)) {
        std::cerr << "error: --pulse-period needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--pulse-duty") {
      if (!parse_double(value("--pulse-duty"), &config.pulse_duty)) {
        std::cerr << "error: --pulse-duty needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--dt") {
      if (!parse_double(value("--dt"), &config.output_interval)) {
        std::cerr << "error: --dt needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--nx") {
      if (!parse_int(value("--nx"), &config.n_x_neg) ||
          !parse_int(value("--nx"), &config.n_x_sep) ||
          !parse_int(value("--nx"), &config.n_x_pos)) {
        std::cerr << "error: --nx needs three integers\n";
        return 2;
      }
      continue;
    }
    if (arg == "--nr") {
      if (!parse_int(value("--nr"), &config.n_r_neg) ||
          !parse_int(value("--nr"), &config.n_r_pos)) {
        std::cerr << "error: --nr needs two integers\n";
        return 2;
      }
      continue;
    }
    if (arg == "--rtol") {
      if (!parse_double(value("--rtol"), &config.rtol)) {
        std::cerr << "error: --rtol needs a number\n";
        return 2;
      }
      continue;
    }
    if (arg == "--out") {
      output_path = value("--out");
      continue;
    }
    if (arg == "--quiet") {
      quiet = true;
      continue;
    }
    std::cerr << "error: unknown option " << arg << "\n";
    return 2;
  }

  const dfn::SimulationReport report = dfn::run_simulation(config);
  if (!report.ok) {
    std::cerr << "error: " << report.message << "\n";
    return 1;
  }

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

  if (!quiet && !report.samples.empty()) {
    const dfn::SimulationSample& last = report.samples.back();
    std::cerr << std::fixed << std::setprecision(4) << "mode="
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
