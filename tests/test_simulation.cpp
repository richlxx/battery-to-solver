#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "dfn/sim/config.hpp"
#include "dfn/sim/simulation.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

bool nearly(double a, double b, double tolerance) {
  return std::abs(a - b) <= tolerance;
}

}  // namespace

int main() {
  const dfn::Params builtin = dfn::make_demo_params();
  const dfn::Params loaded =
      dfn::load_params_yaml(std::string(DFN_SOURCE_DIR) + "/configs/demo.yaml");
  if (!nearly(loaded.area, builtin.area, 0.0) ||
      !nearly(loaded.neg.thickness, builtin.neg.thickness, 0.0) ||
      !nearly(loaded.pos.c0, builtin.pos.c0, 0.0) ||
      !nearly(loaded.sep.porosity, builtin.sep.porosity, 0.0) ||
      !nearly(loaded.elyte.t_plus, builtin.elyte.t_plus, 0.0)) {
    return fail("YAML loader does not reproduce the demo parameters");
  }

#ifndef DFN_USE_KLU
  std::cout << "dfn sim ok (parameter loader only, no SUNDIALS)\n";
  return 0;
#else
  dfn::SimulationConfig config;
  config.mode = dfn::SimulationMode::cccv;
  config.n_x_neg = 4;
  config.n_x_sep = 2;
  config.n_x_pos = 4;
  config.n_r_neg = 4;
  config.n_r_pos = 4;
  config.t_hold = 400.0;
  config.output_interval = 300.0;

  const dfn::SimulationReport report = dfn::run_simulation(config);
  if (!report.ok) {
    std::cerr << report.message << "\n";
    return fail("CC-CV simulation failed");
  }
  if (!report.cutoff_reached) {
    return fail("CC phase did not reach the voltage cutoff");
  }
  if (report.samples.size() < 10) {
    return fail("too few output samples");
  }
  for (std::size_t k = 1; k < report.samples.size(); ++k) {
    if (report.samples[k].t <= report.samples[k - 1].t) {
      return fail("sample times are not increasing");
    }
  }
  if (!nearly(report.samples.front().t, 0.0, 1e-12) ||
      !nearly(report.samples.front().voltage, 3.98572, 2e-3)) {
    return fail("initial sample is wrong");
  }

  bool saw_cv = false;
  for (const dfn::SimulationSample& sample : report.samples) {
    if (sample.phase == "cv") {
      saw_cv = true;
      if (!nearly(sample.voltage, 3.5, 1e-6)) {
        return fail("CV phase does not hold the cutoff voltage");
      }
    }
  }
  if (!saw_cv) {
    return fail("CC-CV did not enter the CV phase");
  }
  const dfn::SimulationSample& last = report.samples.back();
  if (std::abs(last.current) >= 2.5 || std::abs(last.voltage - 3.5) > 1e-6) {
    return fail("CV current did not decay at the cutoff voltage");
  }
  if (report.capacity < 1.2e5 || report.capacity > 1.45e5) {
    std::cerr << "capacity " << report.capacity << "\n";
    return fail("discharged capacity is outside the plausible range");
  }
  if (report.lithium_drift > 1e-6) {
    return fail("lithium drifted during the CC-CV run");
  }

  std::ostringstream csv;
  dfn::write_csv(report, csv);
  std::size_t lines = 0;
  std::string first_line;
  std::istringstream reader(csv.str());
  std::string line;
  while (std::getline(reader, line)) {
    if (lines == 0) {
      first_line = line;
    }
    ++lines;
  }
  if (first_line !=
          "t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,phase" ||
      lines != report.samples.size() + 1) {
    return fail("CSV output does not match the sample list");
  }

  dfn::SimulationConfig pulse;
  pulse.mode = dfn::SimulationMode::pulse;
  pulse.n_x_neg = 4;
  pulse.n_x_sep = 2;
  pulse.n_x_pos = 4;
  pulse.n_r_neg = 3;
  pulse.n_r_pos = 3;
  pulse.t_max = 1200.0;
  pulse.t_pulse = 600.0;
  pulse.pulse_duty = 0.5;
  pulse.output_interval = 100.0;
  const dfn::SimulationReport pulse_report = dfn::run_simulation(pulse);
  if (!pulse_report.ok) {
    std::cerr << pulse_report.message << "\n";
    return fail("pulse simulation failed");
  }
  bool saw_on = false;
  bool saw_rest = false;
  for (const dfn::SimulationSample& sample : pulse_report.samples) {
    if (sample.phase == "pulse-on" && nearly(sample.current, 36.0, 1e-9)) {
      saw_on = true;
    }
    if (sample.phase == "pulse-rest" && nearly(sample.current, 0.0, 1e-12)) {
      saw_rest = true;
    }
  }
  if (!saw_on || !saw_rest) {
    return fail("pulse phases do not alternate current and rest");
  }
  if (!nearly(pulse_report.t_end, 1200.0, 1e-6)) {
    return fail("pulse run did not reach t_max");
  }

  std::cout << "dfn sim ok: t_event " << report.t_end << " s, capacity "
            << report.capacity / 3600.0 << " Ah/m^2, CV end " << last.current
            << " A/m^2, pulse t_end " << pulse_report.t_end << " s\n";
  return 0;
#endif
}
