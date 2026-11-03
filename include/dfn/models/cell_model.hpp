#pragma once

#include <stdexcept>
#include <vector>

#include "dfn/params.hpp"
#include "dfn/solvers/dae.hpp"

namespace dfn {

// Common interface of the SPM -> SPMe -> DFN ladder. The simulation driver
// only depends on this interface; each model owns its own state layout.
class CellModel : public DaeModel {
 public:
  virtual std::vector<double> initial_state() const = 0;
  virtual std::vector<double> suggested_atol() const = 0;

  // Terminal voltage and current density, A/m^2 (positive = discharge).
  virtual double terminal_voltage(const double* y) const = 0;
  virtual double terminal_current(const double* y) const = 0;

  virtual double average_stoichiometry(const double* y, bool positive) const = 0;
  virtual double total_lithium(const double* y) const = 0;

  virtual void set_current_density(double value) = 0;
  virtual void set_applied_voltage(double) {
    throw std::logic_error("this model does not support voltage control");
  }
  virtual bool supports_voltage_control() const { return false; }

  virtual const Params& params() const = 0;
  virtual const char* name() const = 0;
};

}  // namespace dfn
