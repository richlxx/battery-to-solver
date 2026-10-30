#include "dfn/state.hpp"

#include <stdexcept>

namespace dfn {

StateLayout make_layout(int n_x_neg, int n_x_sep, int n_x_pos,
                        int n_r_neg, int n_r_pos) {
  if (n_x_neg <= 0 || n_x_sep <= 0 || n_x_pos <= 0 ||
      n_r_neg <= 0 || n_r_pos <= 0) {
    throw std::invalid_argument("grid counts must be positive");
  }
  StateLayout layout;
  layout.n_x_neg = n_x_neg;
  layout.n_x_sep = n_x_sep;
  layout.n_x_pos = n_x_pos;
  layout.n_r_neg = n_r_neg;
  layout.n_r_pos = n_r_pos;
  return layout;
}

std::vector<double> make_initial_state(const Params& params,
                                       const StateLayout& layout) {
  std::vector<double> y(static_cast<std::size_t>(layout.size()), 0.0);

  for (int i_x = 0; i_x < layout.n_x_neg; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_neg; ++i_r) {
      y[layout.cs_neg(i_x, i_r)] = params.neg.c0;
    }
  }
  for (int i_x = 0; i_x < layout.n_x_pos; ++i_x) {
    for (int i_r = 0; i_r < layout.n_r_pos; ++i_r) {
      y[layout.cs_pos(i_x, i_r)] = params.pos.c0;
    }
  }
  for (int i_x = 0; i_x < layout.n_x_total(); ++i_x) {
    y[layout.ce(i_x)] = params.elyte.c0;
  }

  return y;
}

}  // namespace dfn
