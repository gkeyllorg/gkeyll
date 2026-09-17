#pragma once

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include <gkyl_util.h>

/**
 * Rectangular grid object.
 */
struct gkyl_rect_grid {
  int ndim; // number of dimensions
  double lower[GKYL_MAX_DIM]; // lower-left corner
  double upper[GKYL_MAX_DIM]; // upper-right corner
  int cells[GKYL_MAX_DIM]; // number of cells    
  double dx[GKYL_MAX_DIM]; // cell spacing
  double cellVolume; // cell volume
};

/**
 * Initialize a grid object.
 *
 * @param grid Grid object to initialize.
 * @param ndim Dimension of grid
 * @param lower Coordinates of lower-left corner of grid
 * @param upper Coordinates of upper-right corner of grid
 * @param cells Number of cells in each direction
 */
void gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells);

/**
 * Create and initialize a new grid object. Release with gkyl_rect_grid_release.
 *
 * @param grid Grid object to initialize.
 * @param ndim Dimension of grid
 * @param lower Coordinates of lower-left corner of grid
 * @param upper Coordinates of upper-right corner of grid
 * @param cells Number of cells in each direction
 */
struct gkyl_rect_grid* gkyl_rect_grid_new(int ndim,
  const double *lower, const double *upper, const int *cells);

/**
 * Find cell indices of point
 *
 * @param grid Grid object.
 * @param point The point to find the cell indices at.
 * @param pick_lower If point is on a cell boundary, pick the lower cell if true (per direction).
 * @param known_index Any known indices of where the point is (<0 if not known).
 * @param cell_index Pointer to cell indices.
 * Asserts: point lies within cell(s) specified by known_index (if specified). 
 */
GKYL_CU_DH
void gkyl_rect_grid_find_cell(const struct gkyl_rect_grid *grid, const double *point,
  const bool *pick_lower, const int *known_index, int *cell_index);

/**
 * Get cell-center coordinates. Note that idx is a 1-based cell index,
 * i.e. the lower-left corner is (1,1,...).
 *
 * @param grid Grid object
 * @param idx Index of cell (lower-left corner has all index (1,1,...) )
 * @param xc On output, cell-center coordinates of cell 'idx'
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_cell_center(const struct gkyl_rect_grid *grid,
  const int *idx, double *xc)
{
  for (int i=0; i<grid->ndim; ++i)
    xc[i] = grid->lower[i]+(idx[i]-0.5)*grid->dx[i];
}

/**
 * Get coordinate of lower-left node. Note that idx is a 1-based cell
 * index, i.e. the lower-left corner is (1,1,...).
 *
 * @param grid Grid object
 * @param idx Index of cell (lower-left corner has all index (1,1,...) )
 * @param xn On output, coordinates of lower-left node
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_ll_node(const struct gkyl_rect_grid *grid,
  const int *idx, double *xc)
{
  for (int i=0; i<grid->ndim; ++i)
    xc[i] = grid->lower[i]+(idx[i]-1)*grid->dx[i];
}

/**
 * Get index extents in direction @a dir. The extents are inclusive.
 *
 * @param grid Grid object
 * @param dir Direction in which to get extents
 * @param ext On output, inclusive extents in direction @a dir.
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_extents(const struct gkyl_rect_grid *grid, int dir, int ext[2])
{
  ext[0] = 1; ext[1] = grid->cells[dir];
}

/**
 * Get index of point with coordinate @a xn
 *
 * @param grid Grid object
 * @param xn Coordinate of point in grid
 * @param idx On output, index of point in grid
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_coord_idx(const struct gkyl_rect_grid *grid,
  const double *xn, int *idx)
{
  for (int d=0; d<grid->ndim; ++d) {
    int ext[2]; gkyl_rect_grid_extents(grid, d, ext);
    double xlower = grid->lower[d], dx = grid->dx[d];
    idx[d] = ext[0] + (int) floor((xn[d]-xlower)/dx);
  }
}

/**
 * Compare grids
 *
 * @param grid1 Grid object to compare
 * @param grid2 Grid object to compare
 * @return true if the grids are the same, false otherwise
 */
bool gkyl_rect_grid_cmp(const struct gkyl_rect_grid *grid1, struct gkyl_rect_grid *grid2);

/**
 * Write grid data to given file pointer. If file pointer is not stdout or
 * stderr, data is written to previously opened file in binary format.
 *
 * @param grid Grid object to write
 * @param nm Name of grid (for fp=stdout or stderr);
 * @param fp File handle to write to.
 */
void gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, const char *nm, FILE *fp);

/**
 * Read data from file and initialize grid. File must be opened by
 * caller of this function.
 *
 * @param grid Grid object to initialize
 * @param fp File handle to read fromx.
 * @return True if read succeeded, false otherwise
 */
bool gkyl_rect_grid_read(struct gkyl_rect_grid *grid, FILE *fp);

/**
 * Release grid object that was created with gkyl_rect_grid_new.
 *
 * @param grid Grid object.
 */
void gkyl_rect_grid_release(struct gkyl_rect_grid *grid);
