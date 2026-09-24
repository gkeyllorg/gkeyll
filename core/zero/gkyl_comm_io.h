#include <gkyl_comm.h>

/**
 * Write out grid and array data to file in .gkyl format so postgkyl
 * can understand it.
 *
 * @param comm Communicator
 * @param grid Grid object to write
 * @param range Range describing portion of the array to output.
 * @param meta Meta-data to write. Set to NULL or 0 if no metadata
 * @param arr Array object to write
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_comm_array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname);

/**
 * Read array data from .gkyl format. The input grid must be
 * pre-computed and must match the grid in the array. An error is
 * returned if this is not the case.
 *
 * @param comm Communicator
 * @param grid Grid object for read
 * @param range Range describing portion of the array to read.
 * @param arr Array object to read
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_comm_array_read(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname);
