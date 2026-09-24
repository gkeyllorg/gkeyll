#pragma once

#include <gkyl_rect_grid.h>

/* Find upper and lower boundaries of given cell
 * @params grid: grid struct which contains the cells
 * @params cell_in: given cell
 * @params dim_trans: The dimensions to check
 * @params known_index: Any already known indices
 * @params lower_boundaries: lower sides of given cell (output)
 * @params upper_boundaries: upper sides of given cell (output)
 */
GKYL_CU_DH
void in_dir(const struct gkyl_rect_grid *grid, int *cell_in, const int *dim_trans,
  const int *known_index, double lower_boundaries[], double upper_boundaries[])
{
  int ndim, check_index;
  const double *dx, *lower;
  ndim = grid -> ndim;
  dx = grid -> dx;
  lower = grid -> lower;
  for (int d=0; d<ndim; d++) {
    check_index = known_index[d] < 0 ? cell_in[dim_trans[d]] : known_index[d];
    lower_boundaries[d] = lower[d]+(check_index-1)*dx[d];
    upper_boundaries[d] = lower[d]+(check_index)*dx[d];
  }
}

/* Checks if given point is in given cell
 * @params grid: grid struct which contains the cells
 * @params point: coordinates of given point
 * @params cell_in: given cell
 * @params dim_trans: The dimensions to check
 * @params known_index: Any already known indices
 * @returns bool: true if point is in given cell
 */
GKYL_CU_DH
bool is_in_cell(const struct gkyl_rect_grid *grid, const double *point,
  int *cell_in, const int *dim_trans, const int *known_index)
{
  int ndim;
  ndim = grid -> ndim;
  double lower_boundaries[ndim], upper_boundaries[ndim];
  for (int d=0; d<ndim; d++) {
    lower_boundaries[d]=0;
    upper_boundaries[d]=0;
  }
  in_dir(grid, cell_in, dim_trans, known_index, lower_boundaries, upper_boundaries);
  bool in_cell = true;
  double eps = 1.0e-6;
  for (int d=0; d<ndim; d++) {
    if (lower_boundaries[d]-eps*grid->dx[d]>point[d] || upper_boundaries[d]+eps*grid->dx[d]<point[d]) {
      in_cell = false;
      break;
    }
  }
  return in_cell;
}

