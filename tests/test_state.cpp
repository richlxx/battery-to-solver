#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "dfn/params.hpp"
#include "dfn/state.hpp"

namespace {

int fail(const char* message) {
  std::cerr << "FAIL: " << message << "\n";
  return 1;
}

}  // namespace

int main() {
  const dfn::Params params = dfn::make_demo_params();
  dfn::validate(params);

  const double expected_length = 85.0e-6 + 12.0e-6 + 75.0e-6;
  if (std::abs(params.cell_length() - expected_length) > 1e-18) {
    return fail("cell length");
  }

  const dfn::StateLayout layout = dfn::make_layout(3, 2, 4, 5, 6);
  const int expected_size = 3 * 5 + 4 * 6 + 9 + 3 + 4 + 9;
  if (layout.size() != expected_size) {
    return fail("layout size");
  }

  std::vector<int> hits(static_cast<std::size_t>(layout.size()), 0);
  for (int i_x = 0; i_x < layout.n_x_neg; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_neg; ++i_r) {
      ++hits[layout.cs_neg(i_x, i_r)];
    }
  }
  for (int i_x = 0; i_x < layout.n_x_pos; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_pos; ++i_r) {
      ++hits[layout.cs_pos(i_x, i_r)];
    }
  }
  for (int i_x = 0; i_x < layout.n_x_total(); ++i_x) {
    ++hits[layout.ce(i_x)];
    ++hits[layout.phie(i_x)];
  }
  for (int i_x = 0; i_x < layout.n_x_neg; ++i_x) {
    ++hits[layout.phis_neg(i_x)];
  }
  for (int i_x = 0; i_x < layout.n_x_pos; ++i_x) {
    ++hits[layout.phis_pos(i_x)];
  }
  for (int count : hits) {
    if (count != 1) {
      return fail("index coverage");
    }
  }

  const std::vector<double> y = dfn::make_initial_state(params, layout);
  if (y.size() != static_cast<std::size_t>(layout.size())) {
    return fail("initial state size");
  }
  if (std::abs(y[layout.cs_neg(2, 4)] - params.neg.c0) > 1e-12) {
    return fail("negative concentration");
  }
  if (std::abs(y[layout.cs_pos(3, 5)] - params.pos.c0) > 1e-12) {
    return fail("positive concentration");
  }
  if (std::abs(y[layout.ce(4)] - params.elyte.c0) > 1e-12) {
    return fail("electrolyte concentration");
  }
  if (y[layout.phis_neg(0)] != 0.0 || y[layout.phie(0)] != 0.0) {
    return fail("potentials should start at zero");
  }

  dfn::Params bad = params;
  bad.neg.thickness = -1.0;
  bool threw = false;
  try {
    dfn::validate(bad);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  if (!threw) {
    return fail("validation should reject negative thickness");
  }

  std::cout << "state layout ok (" << layout.size() << " unknowns)\n";
  return 0;
}
