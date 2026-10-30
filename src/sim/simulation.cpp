#include "dfn/sim/simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

#ifdef DFN_USE_KLU
#include "dfn/solvers/idaklu.hpp"
#endif

namespace dfn {

namespace {

using Clock = std::chrono::steady_clock;

#ifdef DFN_USE_KLU

struct PhaseTotals {
  int steps = 0;
  int residual_evals = 0;
  int jacobian_evals = 0;
  int linear_setups = 0;
};

// Grid points k * dt strictly between t0 and t1. ida_integrate appends t1
// itself, so the endpoint is not repeated here.
std::vector<double> sample_times(double t0, double t1, double dt) {
  std::vector<double> times;
  if (dt <= 0.0 || t1 <= t0) {
    return times;
  }
  const double first = std::floor(t0 / dt) + 1.0;
  for (double k = first;; k += 1.0) {
    const double t = k * dt;
    if (t >= t1 - 1e-9) {
      break;
    }
    times.push_back(t);
  }
  return times;
}

double controlled_current(const SimulationConfig& config) {
  const double magnitude =
      config.current_density != 0.0
          ? std::abs(config.current_density)
          : std::abs(config.c_rate) * config.i_1c;
  return config.charge ? -magnitude : magnitude;
}

void append_samples(const DfnModel& model, const SolveReport& report,
                    const std::string& phase, SimulationReport& out) {
  for (std::size_t k = 0; k < report.output_t.size(); ++k) {
    const double t = report.output_t[k];
    if (!out.samples.empty() && t <= out.samples.back().t + 1e-9) {
      continue;
    }
    const double* y = report.output_y[k].data();
    SimulationSample sample;
    sample.t = t;
    sample.voltage = model.terminal_voltage(y);
    sample.current = model.terminal_current(y);
    sample.stoich_neg = model.average_stoichiometry(y, false);
    sample.stoich_pos = model.average_stoichiometry(y, true);
    sample.phase = phase;
    out.samples.push_back(std::move(sample));
  }
}

bool run_phase(DfnModel& model, std::vector<double>& y, std::vector<double>& yp,
               double t0, double t1, const IdaOptions& options, double dt,
               const std::string& phase, SimulationReport& out,
               PhaseTotals& totals, bool* stopped_by_event, double* t_event) {
  IdaOptions phase_options = options;
  phase_options.rtol = options.rtol;
  phase_options.output_times = sample_times(t0, t1, dt);
  const SolveReport report =
      ida_integrate(model, y, yp, t0, t1, phase_options);
  *stopped_by_event = report.event_index >= 0;
  *t_event = report.t_event;

  totals.steps += report.steps;
  totals.residual_evals += report.residual_evals;
  totals.jacobian_evals += report.jacobian_evals;
  totals.linear_setups += report.linear_setups;

  if (!report.ok) {
    out.message = phase + ": " + report.message;
    return false;
  }
  append_samples(model, report, phase, out);
  y = report.y;
  yp = report.yp;
  return true;
}

// Lithium stored in one electrode per unit separator area, mol/m^2.
double electrode_lithium_per_area(const DfnModel& model, const double* y,
                                  bool positive) {
  const Params& params = model.params();
  const ElectrodeParams& electrode = positive ? params.pos : params.neg;
  return model.average_stoichiometry(y, positive) * electrode.c_max *
         (1.0 - electrode.porosity) * electrode.thickness;
}

#endif  // DFN_USE_KLU

}  // namespace

SimulationReport run_simulation(const SimulationConfig& config) {
  SimulationReport report;
  const auto start = Clock::now();

  try {
    validate(config.params);
  } catch (const std::exception& error) {
    report.message = error.what();
    return report;
  }
  if (config.n_x_neg <= 0 || config.n_x_sep <= 0 || config.n_x_pos <= 0 ||
      config.n_r_neg <= 0 || config.n_r_pos <= 0) {
    report.message = "grid sizes must be positive";
    return report;
  }
  if (config.t_max <= 0.0) {
    report.message = "t_max must be positive";
    return report;
  }
  if (config.mode == SimulationMode::pulse &&
      (config.t_pulse <= 0.0 || config.pulse_duty <= 0.0 ||
       config.pulse_duty >= 1.0)) {
    report.message = "pulse period must be positive and duty in (0, 1)";
    return report;
  }

#ifdef DFN_USE_KLU
  DfnModel model(config.params,
                 make_layout(config.n_x_neg, config.n_x_sep, config.n_x_pos,
                             config.n_r_neg, config.n_r_pos),
                 DfnOptions{});
  std::vector<double> y = model.initial_state();
  std::vector<double> yp(static_cast<std::size_t>(model.size()), 0.0);
  report.unknowns = model.size();
  report.nnz = model.jacobian_pattern().nnz();
  const double lithium0 = model.total_lithium(y.data());
  const double neg_lithium0 = electrode_lithium_per_area(model, y.data(), false);

  IdaOptions options;
  options.rtol = config.rtol;
  options.atol = model.suggested_atol();
  options.calculate_consistent_initial = true;

  PhaseTotals totals;
  const double i_cc = controlled_current(config);
  const double v_cut = config.charge ? config.v_max : config.v_min;
  const double i_cut = 0.05 * std::abs(i_cc);

  auto voltage_cutoff = [&model, v_cut, charge = config.charge](
                            double, const double* state, const double*) {
    const double voltage = model.terminal_voltage(state);
    return charge ? v_cut - voltage : voltage - v_cut;
  };

  double t = 0.0;
  bool stopped = false;
  double t_event = 0.0;

  switch (config.mode) {
    case SimulationMode::cc:
    case SimulationMode::cccv: {
      const double v0 = model.terminal_voltage(y.data());
      if (config.charge ? v0 >= v_cut : v0 <= v_cut) {
        report.message = config.charge
                             ? "initial voltage is already above the cutoff"
                             : "initial voltage is already below the cutoff";
        return report;
      }
      model.set_current_density(i_cc);
      options.events = {voltage_cutoff};
      if (!run_phase(model, y, yp, 0.0, config.t_max, options,
                     config.output_interval, "cc", report, totals, &stopped,
                     &t_event)) {
        return report;
      }
      t = stopped ? t_event : config.t_max;
      report.cutoff_reached = stopped;

      if (config.mode == SimulationMode::cccv && stopped &&
          config.t_hold > 0.0) {
        model.set_applied_voltage(v_cut);
        IdaOptions cv_options = options;
        cv_options.events.clear();
        if (i_cut > 0.0) {
          cv_options.events.push_back(
              [&model, i_cut](double, const double* state, const double*) {
                return std::abs(model.terminal_current(state)) - i_cut;
              });
        }
        bool cv_stopped = false;
        double t_cv_event = 0.0;
        if (!run_phase(model, y, yp, t, t + config.t_hold, cv_options,
                       config.output_interval, "cv", report, totals,
                       &cv_stopped, &t_cv_event)) {
          return report;
        }
        t = cv_stopped ? t_cv_event : t + config.t_hold;
      }
      break;
    }
    case SimulationMode::cv: {
      const double v_hold = config.hold_voltage > 0.0 ? config.hold_voltage
                                                      : v_cut;
      if (std::abs(model.terminal_voltage(y.data()) - v_hold) > 1e-4) {
        report.message =
            "cv mode must start at the hold voltage; run cccv instead of "
            "stepping the voltage on a fresh cell";
        return report;
      }
      model.set_applied_voltage(v_hold);
      if (i_cut > 0.0) {
        options.events.push_back(
            [&model, i_cut](double, const double* state, const double*) {
              return std::abs(model.terminal_current(state)) - i_cut;
            });
      }
      bool cv_stopped = false;
      double t_cv_event = 0.0;
      if (!run_phase(model, y, yp, 0.0, config.t_max, options,
                     config.output_interval, "cv", report, totals, &cv_stopped,
                     &t_cv_event)) {
        return report;
      }
      t = cv_stopped ? t_cv_event : config.t_max;
      break;
    }
    case SimulationMode::pulse: {
      const double v0 = model.terminal_voltage(y.data());
      if (config.charge ? v0 >= v_cut : v0 <= v_cut) {
        report.message = config.charge
                             ? "initial voltage is already above the cutoff"
                             : "initial voltage is already below the cutoff";
        return report;
      }
      options.events = {voltage_cutoff};
      const double on_time = config.t_pulse * config.pulse_duty;
      const double off_time = config.t_pulse - on_time;
      while (t < config.t_max - 1e-9) {
        const double t_on_end = std::min(t + on_time, config.t_max);
        model.set_current_density(i_cc);
        if (!run_phase(model, y, yp, t, t_on_end, options,
                       config.output_interval, "pulse-on", report, totals,
                       &stopped, &t_event)) {
          return report;
        }
        t = stopped ? t_event : t_on_end;
        if (stopped) {
          report.cutoff_reached = true;
          break;
        }
        const double t_off_end = std::min(t + off_time, config.t_max);
        model.set_current_density(0.0);
        if (!run_phase(model, y, yp, t, t_off_end, options,
                       config.output_interval, "pulse-rest", report, totals,
                       &stopped, &t_event)) {
          return report;
        }
        t = stopped ? t_event : t_off_end;
        if (stopped) {
          report.cutoff_reached = true;
          break;
        }
      }
      break;
    }
  }

  const double lithium1 = model.total_lithium(y.data());
  report.lithium_drift = std::abs(lithium1 - lithium0) / lithium0;
  const double moved =
      neg_lithium0 - electrode_lithium_per_area(model, y.data(), false);
  report.capacity = kFaraday * (config.charge ? -moved : moved);
  report.t_end = report.samples.empty() ? t : report.samples.back().t;
  report.steps = totals.steps;
  report.residual_evals = totals.residual_evals;
  report.jacobian_evals = totals.jacobian_evals;
  report.linear_setups = totals.linear_setups;
  report.ok = true;
#else
  (void)config;
  report.message =
      "dfn-sim was built without SUNDIALS; configure with -DDFN_USE_KLU=ON";
#endif

  report.wall_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  return report;
}

void write_csv(const SimulationReport& report, std::ostream& out) {
  out << "t_s,voltage_v,current_a_m2,stoich_neg,stoich_pos,phase\n";
  out << std::setprecision(10);
  for (const SimulationSample& sample : report.samples) {
    out << sample.t << ',' << sample.voltage << ',' << sample.current << ','
        << sample.stoich_neg << ',' << sample.stoich_pos << ','
        << sample.phase << '\n';
  }
}

void write_csv(const SimulationReport& report, const std::string& path) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("cannot open " + path + " for writing");
  }
  write_csv(report, file);
}

}  // namespace dfn
