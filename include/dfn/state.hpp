#pragma once

#include <vector>

#include "dfn/params.hpp"

namespace dfn {

struct StateLayout {
  int n_x_neg = 0;
  int n_x_sep = 0;
  int n_x_pos = 0;
  int n_r_neg = 0;
  int n_r_pos = 0;

  int n_x_total() const { return n_x_neg + n_x_sep + n_x_pos; }

  int cs_neg_offset() const { return 0; }
  int cs_pos_offset() const { return n_x_neg * n_r_neg; }
  int ce_offset() const { return cs_pos_offset() + n_x_pos * n_r_pos; }
  int phis_neg_offset() const { return ce_offset() + n_x_total(); }
  int phis_pos_offset() const { return phis_neg_offset() + n_x_neg; }
  int phie_offset() const { return phis_pos_offset() + n_x_pos; }
  int size() const { return phie_offset() + n_x_total(); }

  int cs_neg(int i_x, int i_r) const {
    return cs_neg_offset() + i_x * n_r_neg + i_r;
  }
  int cs_pos(int i_x, int i_r) const {
    return cs_pos_offset() + i_x * n_r_pos + i_r;
  }
  int ce(int i_x) const { return ce_offset() + i_x; }
  int phis_neg(int i_x) const { return phis_neg_offset() + i_x; }
  int phis_pos(int i_x) const { return phis_pos_offset() + i_x; }
  int phie(int i_x) const { return phie_offset() + i_x; }
};

StateLayout make_layout(int n_x_neg, int n_x_sep, int n_x_pos,
                        int n_r_neg, int n_r_pos);

std::vector<double> make_initial_state(const Params& params,
                                       const StateLayout& layout);

}  // namespace dfn
