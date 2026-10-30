#pragma once

#include <vector>

namespace dfn {

struct Grid1D {
  int n = 0;
  std::vector<double> face;     // n + 1
  std::vector<double> center;   // n
  std::vector<double> volume;   // n
  std::vector<double> area;     // n + 1
  std::vector<double> spacing;  // n + 1, center-to-center distance at faces
};

Grid1D make_uniform_line_grid(int n, double length, double cross_section);

// Concatenates uniform cell blocks with exact interfaces. counts[k] cells of
// length lengths[k] are placed back to back; faces at the block boundaries.
Grid1D make_segmented_line_grid(const std::vector<int>& counts,
                                const std::vector<double>& lengths,
                                double cross_section);

Grid1D make_uniform_sphere_grid(int n, double radius);

}  // namespace dfn
