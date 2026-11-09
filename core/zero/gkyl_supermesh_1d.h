#pragma once

#include <math.h>

#include <gkyl_util.h>

// One interval of the common refinement of two 1D grids: [lo, up] lies in the cell k_l of the
// left grid and in the cell k_r of the right grid (-1: outside that grid).
struct gkyl_supermesh_1d_interval {
  double lo, up; // ends of the interval in the common coordinate x
  int k_l, k_r; // cell of each grid containing the interval (-1: none)
};

/**
 * Largest number of intervals of gkyl_supermesh_1d_merge.
 *
 * @param num_edges_l Number of edges of the left grid.
 * @param num_edges_r Number of edges of the right grid.
 * @param num_kinks Number of kinks.
 * @return Capacity the output array needs.
 */
GKYL_CU_DH static inline int
gkyl_supermesh_1d_max_intervals(int num_edges_l, int num_edges_r, int num_kinks)
{
  return num_edges_l + num_edges_r + num_kinks - 1;
}

/**
 * The supermesh of two 1D grids and a set of kinks: the common refinement of the two grids, each
 * mapped to a common coordinate x through its own affine map x_s = r_s x + c_s (r_s > 0), further
 * split at the kinks. Every interval lies in at most one cell of each grid. Intervals not longer
 * than tol are dropped (coincident edges, kinks on edges), and so are the gaps outside both
 * grids, which exist when the two mapped grids are disjoint. Kinks outside both grids split
 * nothing.
 *
 * @param num_edges_l Number of edges of the left grid.
 * @param edges_l Sorted edges of the left grid, in its own coordinate.
 * @param r_l Dilation of the left grid's map.
 * @param c_l Shift of the left grid's map.
 * @param num_edges_r Number of edges of the right grid.
 * @param edges_r Sorted edges of the right grid, in its own coordinate.
 * @param r_r Dilation of the right grid's map.
 * @param c_r Shift of the right grid's map.
 * @param num_kinks Number of kinks.
 * @param kinks Sorted kinks, in the common coordinate.
 * @param tol Length at or below which an interval is dropped.
 * @param intervals Output intervals, sorted (capacity: gkyl_supermesh_1d_max_intervals).
 * @return Number of intervals written.
 */
GKYL_CU_DH static inline int
gkyl_supermesh_1d_merge(
  int num_edges_l, const double *edges_l, double r_l, double c_l, int num_edges_r,
  const double *edges_r, double r_r, double c_r, int num_kinks, const double *kinks, double tol,
  struct gkyl_supermesh_1d_interval *intervals
)
{
  const double past_end = 1e300; // finite sentinel (no infinities under -ffast-math)
  int next_l = 1, next_r = 1, next_k = 0, count = 0;
  double x_l = (edges_l[0] - c_l) / r_l, x_r = (edges_r[0] - c_r) / r_r;
  double pos = fmin(x_l, x_r);
  while (next_k < num_kinks && kinks[next_k] <= pos) {
    ++next_k;
  }
  double x_k = next_k < num_kinks ? kinks[next_k] : past_end;

  // The grid whose first edge is pos is entered at the start; a grid starting later is entered
  // when the sweep reaches its first edge, which the event loop below does as any other edge.
  int k_l = x_l <= pos ? 0 : -1, k_r = x_r <= pos ? 0 : -1;
  if (k_l == 0) {
    x_l = next_l < num_edges_l ? (edges_l[next_l] - c_l) / r_l : past_end;
  } else {
    next_l = 0;
  }
  if (k_r == 0) {
    x_r = next_r < num_edges_r ? (edges_r[next_r] - c_r) / r_r : past_end;
  } else {
    next_r = 0;
  }

  while (next_l < num_edges_l || next_r < num_edges_r) {
    double next = fmin(x_k, fmin(x_l, x_r));
    if (next - pos > tol && (k_l >= 0 || k_r >= 0)) {
      intervals[count].lo = pos;
      intervals[count].up = next;
      intervals[count].k_l = k_l;
      intervals[count].k_r = k_r;
      ++count;
    }
    pos = next;
    if (x_l <= x_r && x_l <= x_k) {
      k_l = next_l < num_edges_l - 1 ? next_l : -1;
      ++next_l;
      x_l = next_l < num_edges_l ? (edges_l[next_l] - c_l) / r_l : past_end;
    } else if (x_r <= x_k) {
      k_r = next_r < num_edges_r - 1 ? next_r : -1;
      ++next_r;
      x_r = next_r < num_edges_r ? (edges_r[next_r] - c_r) / r_r : past_end;
    } else {
      ++next_k;
      x_k = next_k < num_kinks ? kinks[next_k] : past_end;
    }
  }
  return count;
}
