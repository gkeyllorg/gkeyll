#pragma once

#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_ref_count.h>

// Decomposition object
struct gkyl_rect_decomp {
  int ndim; // dimension of decomposition
  int ndecomp; // number of sub-domains
  struct gkyl_range parent_range; // range that was decomposed
  struct gkyl_range *ranges; // decomposed ranges

  struct gkyl_ref_count ref_count;
};

// List of neighbors 
struct gkyl_rect_decomp_neigh {
  int num_neigh; // number of neighbors
  const int *neigh; // list of neighbors

  // following information is not typically needed
  const int *dir; // direction in which neigh[i] range is located
  const int *edge; // edge on which neigh[i] range is located
};

/**
 * Create a new decomposition of @a range, given @a cuts in each
 * direction. The total number of decomposed ranges are product of all
 * cuts. The decomposed ranges index independent of @a range,
 * i.e. decomposed ranges are NOT sub-ranges of @a range.
 *
 * @param ndim Number of dimensions
 * @param cuts Cuts in each direction.
 * @param range Range to decompose
 * @return Decomposition of @a range
 */
struct gkyl_rect_decomp* gkyl_rect_decomp_new_from_cuts(int ndim, const int cuts[],
  const struct gkyl_range *range);

/**
 * Create a new decomposition given @a cuts and cells in each
 * direction. The total number of decomposed ranges are product of all
 * cuts.
 *
 * @param ndim Number of dimensions
 * @param cuts Cuts in each direction.
 * @param cells Number of cells in each direction
 * @return Decomposition of range based on cuts
 */
struct gkyl_rect_decomp *gkyl_rect_decomp_new_from_cuts_and_cells(int ndim,
  const int cuts[], const int cells[]);

/**
 * Create a new decomposition from a given decomposition. The new
 * decomposition extends each region by a tensor product with @a
 * arange.
 *
 * @param arange Range to extend by
 * @return New extended decomposition
 */
struct gkyl_rect_decomp *gkyl_rect_decomp_extended_new(const struct gkyl_range *arange,
  const struct gkyl_rect_decomp *decomp);

/**
 * Acquire a pointer to the decomposition.
 *
 * @param decomp Decom to acquire pointer to
 * @return New decomposition
 */
struct gkyl_rect_decomp* gkyl_rect_decomp_acquire(const struct gkyl_rect_decomp *decomp);

/**
 * Check if decomposition is  a valid covering of the range.
 *
 * NOTE: This function internally allocates memory over the complete
 * parent range. This can be a problem if the parent range is huge.
 *
 * @param decomp Demposition to check
 * @return true if this is a valid covering
 */
bool gkyl_rect_decomp_check_covering(const struct gkyl_rect_decomp *decomp);

/**
 * Compute the neighbor of range @a nidx. The returned object must be
 * freed using the gkyl_rect_decomp_neigh_release call.
 *
 * @param decomp Decomposition object
 * @param inc_corners If true, corner neighbors are also included
 * @param nidx Index of range for which neighbor data is needed
 * @return Neighbor list for range nidx
 */
struct gkyl_rect_decomp_neigh* gkyl_rect_decomp_calc_neigh(
  const struct gkyl_rect_decomp *decomp, bool inc_corners, int nidx);

/**
 * Compute the periodic neighbor of range @a nidx in the specified
 * direction. The returned object must be freed using the
 * gkyl_rect_decomp_neigh_release call.
 *
 * @param decomp Decomposition object
 * @param dir Direction to compute periodic neighbors
 * @param inc_corners If true, corner neighbors are also included
 * @param nidx Index of range for which neighbor data is needed
 * @return Periodic neighbor list for range nidx
 */
struct gkyl_rect_decomp_neigh* gkyl_rect_decomp_calc_periodic_neigh(
  const struct gkyl_rect_decomp *decomp, int dir, bool inc_corners, int nidx);

/**
 * Free neighbor memory
 *
 * @param ng Neighbor data to free
 */
void gkyl_rect_decomp_neigh_release(struct gkyl_rect_decomp_neigh *ng);

/**
 * Compute cumulative offet of @a nidx range in the decomp. The
 * cumulative offset is the global linear index of the first cell in
 * the local range.
 *
 * @param decomp Decomposition object
 * @param nidx Index of range for which offset is needed
 * @return Offest of first cell in range[nidx]
 */
long gkyl_rect_decomp_calc_offset(const struct gkyl_rect_decomp *decomp, int nidx);

/**
 * Free decomposition.
 *
 * @param decomp Decomposition to free
 */
void gkyl_rect_decomp_release(struct gkyl_rect_decomp *decomp);

// The functions below are utility functions to construct properly
// nested ranges that extend over the grid or over local ranges, given
// ghost cells.

/**
 * Create range over global region given cells in each direction.
 *
 * @param ndim Grid dimension
 * @param cells Number of cells in each direction
 * @param range On output, global range
 */
void gkyl_create_global_range(int ndim, const int *cells, struct gkyl_range *range);

/**
 * Create range and extended ranges from grid and ghost-cell data. The
 * range is a sub-range of the extended range.
 *
 * @param grid Grid to compute ranges for
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning grid+ghost-cells
 * @param range On output, range spanning grid. Sub-range of ext_range.
 */
void gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range,
  struct gkyl_range *range);

/**
 * Create range and extended ranges from given range and ghost-cell
 * data. The range is a sub-range of the extended range.
 *
 * @param inrange Input range to use
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning inrange+ghost-cells
 * @param range On output, range same as inrange, but sub-range of ext_range.
 */
void gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);

/**
 * Create range and extended vertex ranges from given range and ghost-cell
 * data. The range is a sub-range of the extended range.
 *
 * @param inrange Input range to use
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning inrange+ghost-cells
 * @param range On output, range same as inrange, but sub-range of ext_range.
 */
void gkyl_create_vertex_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);

/**
 * Return the cuts used to create the the decomposition object.
 * 
 * @param decomp Decomposition object.
 * @param cuts Output cuts in each direction.
 */
void gkyl_rect_decomp_get_cuts(struct gkyl_rect_decomp* decomp, int* cuts);
