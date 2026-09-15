#pragma once

#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Construct a node from its global cell index. Avoid accumulating quadrature
// increments from a rank-local origin: the same node must have the same value
// on either side of an MPI face. num_nodes is the number of Gauss nodes per
// cell, or the number of equal intervals per cell for corner/face nodes.
static inline double
gk_geometry_node_coord(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_range *nrange, int dir, int node, int num_nodes, bool quadrature)
{
  int global_node = node-nrange->lower[dir]
    + num_nodes*(local->lower[dir]-global->lower[dir]);
  double cell;
  if (quadrature)
    cell = global_node/num_nodes + 0.5*(1.0 + gkyl_gauss_ordinates[num_nodes][global_node%num_nodes]);
  else {
    if (global_node == 0) return grid->lower[dir];
    if (global_node == num_nodes*grid->cells[dir]) return grid->upper[dir];
    cell = (double) global_node/num_nodes;
  }
  return grid->lower[dir] + cell*grid->dx[dir];
}
