#pragma once

#include <stdio.h>

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Read status flags
enum gkyl_array_rio_status {
  GKYL_ARRAY_RIO_SUCCESS = 0,
  GKYL_ARRAY_RIO_BAD_VERSION,
  GKYL_ARRAY_RIO_FOPEN_FAILED,
  GKYL_ARRAY_RIO_FREAD_FAILED,
  GKYL_ARRAY_RIO_DATA_MISMATCH,
  GKYL_ARRAY_RIO_META_FAILED
};

/**
 * Return character string corresponding to the status enum flag.
 *
 * @param status Status flag
 * @return string corresponding to flag
 */
const char* gkyl_array_rio_status_msg(enum gkyl_array_rio_status status);

// Array header data to write: this is for low-level control and is
// typically not something most users would ever encounter
struct gkyl_array_header_info {
  uint64_t file_type; // file type
  enum gkyl_elem_type etype; // element type
  uint64_t esznc; // elem sz * number of components
  uint64_t tot_cells; // total number of cells in grid
  uint64_t meta_size; // size in bytes of meta-data embedded in header
  char *meta; // meta-data as byte array
  uint64_t nrange; // number of ranges
};

/**
 * Read grid and array data header data from file. Note that only the
 * header is read and NOT the array data itself. If the header has
 * meta-data (meta_size > 0) then the meta char array must be freed
 * using gkyl_free.
 *
 * @param grid Grid object to read
 * @param hrd On output, Header data.
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag
 */
enum gkyl_array_rio_status gkyl_grid_sub_array_header_read(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const char *fname);

/**
 * Free header info if needed (only of meta_size > 0) does this call
 * actually free anything.
 *
 * @param info Header info to free
 */
void gkyl_array_header_info_release(struct gkyl_array_header_info *info);

/**
 * Write out grid and array data to file in .gkyl format so postgkyl
 * can understand it.
 *
 * @param grid Grid object to write
 * @param range Range describing portion of the array to output.
 * @param meta Meta-data to write. Set to NULL or 0 if no metadata
 * @param arr Array object to write
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag
 */
enum gkyl_array_rio_status gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname);

/**
 * Read grid and array data from file. The input array must be
 * pre-allocated and must be big enough to hold the read data.
 * 
 * @param grid Grid object to read
 * @param range Range describing portion of the array.
 * @param arr Array object to read
 * @param fname Name of input file
 * @return Status flag
 */
enum gkyl_array_rio_status gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  struct gkyl_array *arr, const char* fname);

/**
 * Read grid and array data from file, creating a new array.
 * 
 * @param grid On outout, grid on which array is defined.
 * @param fname Name of input file
 * @return Newly created array object. NULL if failed
 */
struct gkyl_array *gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid,
  const char* fname);
