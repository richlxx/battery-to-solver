#include "dfn/grid.hpp"

#include <cmath>
#include <stdexcept>

namespace dfn {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

Grid1D make_uniform_line_grid(int n, double length, double cross_section) {
  if (n <= 0 || length <= 0.0 || cross_section <= 0.0) {
    throw std::invalid_argument("line grid needs positive counts and sizes");
  }

  Grid1D grid;
  grid.n = n;
  grid.face.resize(static_cast<std::size_t>(n + 1));
  grid.center.resize(static_cast<std::size_t>(n));
  grid.volume.resize(static_cast<std::size_t>(n));
  grid.area.assign(static_cast<std::size_t>(n + 1), cross_section);
  grid.spacing.resize(static_cast<std::size_t>(n + 1));

  const double dx = length / static_cast<double>(n);
  for (int i = 0; i <= n; ++i) {
    grid.face[i] = i * dx;
    grid.spacing[i] = dx;
  }
  for (int i = 0; i < n; ++i) {
    grid.center[i] = (i + 0.5) * dx;
    grid.volume[i] = cross_section * dx;
  }
  return grid;
}

Grid1D make_segmented_line_grid(const std::vector<int>& counts,
                                const std::vector<double>& lengths,
                                double cross_section) {
  if (counts.empty() || counts.size() != lengths.size() ||
      cross_section <= 0.0) {
    throw std::invalid_argument("segmented grid needs matching, non-empty blocks");
  }
  int total = 0;
  for (std::size_t k = 0; k < counts.size(); ++k) {
    if (counts[k] <= 0 || lengths[k] <= 0.0) {
      throw std::invalid_argument("segmented grid needs positive blocks");
    }
    total += counts[k];
  }

  Grid1D grid;
  grid.n = total;
  grid.face.assign(static_cast<std::size_t>(total) + 1, 0.0);
  grid.center.assign(static_cast<std::size_t>(total), 0.0);
  grid.volume.assign(static_cast<std::size_t>(total), 0.0);
  grid.area.assign(static_cast<std::size_t>(total) + 1, cross_section);
  grid.spacing.assign(static_cast<std::size_t>(total) + 1, 0.0);

  std::vector<double> boundary(counts.size() + 1, 0.0);
  int offset = 0;
  for (std::size_t r = 0; r < counts.size(); ++r) {
    const double dx = lengths[r] / static_cast<double>(counts[r]);
    boundary[r + 1] = boundary[r] + lengths[r];
    for (int k = 0; k < counts[r]; ++k) {
      const int i = offset + k;
      grid.center[static_cast<std::size_t>(i)] =
          boundary[r] + (static_cast<double>(k) + 0.5) * dx;
      grid.volume[static_cast<std::size_t>(i)] = cross_section * dx;
    }
    offset += counts[r];
  }

  grid.face[0] = 0.0;
  grid.face[static_cast<std::size_t>(total)] = boundary.back();
  offset = counts[0];
  for (std::size_t r = 1; r < counts.size(); ++r) {
    grid.face[static_cast<std::size_t>(offset)] = boundary[r];
    offset += counts[r];
  }
  for (int i = 1; i < total; ++i) {
    if (grid.face[static_cast<std::size_t>(i)] == 0.0) {
      grid.face[static_cast<std::size_t>(i)] =
          0.5 * (grid.center[static_cast<std::size_t>(i - 1)] +
                 grid.center[static_cast<std::size_t>(i)]);
    }
    grid.spacing[static_cast<std::size_t>(i)] =
        grid.center[static_cast<std::size_t>(i)] -
        grid.center[static_cast<std::size_t>(i - 1)];
  }
  grid.spacing[0] = 0.5 * (grid.center[0] - grid.face[0]);
  grid.spacing[static_cast<std::size_t>(total)] =
      grid.face[static_cast<std::size_t>(total)] -
      grid.center[static_cast<std::size_t>(total - 1)];
  return grid;
}

Grid1D make_uniform_sphere_grid(int n, double radius) {
  if (n <= 0 || radius <= 0.0) {
    throw std::invalid_argument("sphere grid needs positive counts and radius");
  }

  Grid1D grid;
  grid.n = n;
  grid.face.resize(static_cast<std::size_t>(n + 1));
  grid.center.resize(static_cast<std::size_t>(n));
  grid.volume.resize(static_cast<std::size_t>(n));
  grid.area.resize(static_cast<std::size_t>(n + 1));
  grid.spacing.resize(static_cast<std::size_t>(n + 1));

  const double dr = radius / static_cast<double>(n);
  for (int i = 0; i <= n; ++i) {
    grid.face[i] = i * dr;
    grid.area[i] = 4.0 * kPi * grid.face[i] * grid.face[i];
    grid.spacing[i] = dr;
  }
  for (int i = 0; i < n; ++i) {
    grid.center[i] = (i + 0.5) * dr;
    const double r0 = grid.face[i];
    const double r1 = grid.face[i + 1];
    grid.volume[i] = 4.0 * kPi / 3.0 * (r1 * r1 * r1 - r0 * r0 * r0);
  }
  return grid;
}

}  // namespace dfn
