#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "dfn/models/dfn.hpp"
#include "dfn/params.hpp"

namespace dfn {

enum class SimulationMode { cc, cv, cccv, pulse };

struct SimulationConfig {
  Params params = make_demo_params();
  SimulationMode mode = SimulationMode::cccv;

  // The demo cell stores about 36 Ah/m^2, so 36 A/m^2 is roughly 1C.
  // Current densities are A/m^2, positive = discharge. Charge flips the
  // applied current sign and turns v_max into the cutoff.
  double i_1c = 36.0;
  double c_rate = 1.0;
  double current_density = 0.0;  // if nonzero, overrides c_rate
  bool charge = false;

  double v_min = 3.5;
  double v_max = 4.0;
  double hold_voltage = 0.0;  // CV mode; zero means "use the cutoff voltage"

  double t_max = 20000.0;
  double t_hold = 600.0;
  double t_pulse = 600.0;
  double pulse_duty = 0.5;
  double output_interval = 60.0;
  double rtol = 1e-6;

  int n_x_neg = 8;
  int n_x_sep = 4;
  int n_x_pos = 8;
  int n_r_neg = 6;
  int n_r_pos = 6;
};

struct SimulationSample {
  double t = 0.0;
  double voltage = 0.0;
  double current = 0.0;  // A/m^2, positive = discharge
  double stoich_neg = 0.0;
  double stoich_pos = 0.0;
  std::string phase;
};

struct SimulationReport {
  bool ok = false;
  std::string message;
  std::vector<SimulationSample> samples;

  double t_end = 0.0;
  // A s/m^2, from the change of lithium stored in the negative electrode.
  // Independent of the output interval.
  double capacity = 0.0;
  bool cutoff_reached = false;
  double lithium_drift = 0.0;  // relative

  int steps = 0;
  int residual_evals = 0;
  int jacobian_evals = 0;
  int linear_setups = 0;
  int unknowns = 0;
  int nnz = 0;
  double wall_ms = 0.0;
};

SimulationReport run_simulation(const SimulationConfig& config);

// Write t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,phase to a stream.
void write_csv(const SimulationReport& report, std::ostream& out);
void write_csv(const SimulationReport& report, const std::string& path);

}  // namespace dfn
