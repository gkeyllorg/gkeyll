#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_eval_offset_fd gkyl_eval_offset_fd;

// Struct to describe offset with respect to the cell-center: 0.0 is
// cell center, -0.5 the left edge and 0.5 the right edge.
struct gkyl_offset_descr {
  double od_off[GKYL_MAX_DIM];
};

// input packaged as a struct
struct gkyl_eval_offset_fd_inp {
  const struct gkyl_rect_grid *grid; // grid on which to project

  int num_ret_vals; // number of return values in eval function
  struct gkyl_offset_descr *offsets; // size num_ret_vals
  
  evalf_t eval; // function to project
  void *ctx; // function context
};

/**
 * Create new updater to evaluate function on a finite-difference
 * grid. Each returned value is evaluated on an internal node that is
 * potentially offset from the cell-center. This updater is useful for
 * initializing finite-difference grids.
 *
 * @param inp Input parameters
 * @return New updater pointer.
 */
gkyl_eval_offset_fd* gkyl_eval_offset_fd_new(const struct gkyl_eval_offset_fd_inp *inp);

/**
 * Compute function of nodes. The update_rng MUST be a sub-range of
 * the range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Updater to run
 * @param tm Time at which projection must be computed
 * @param update_rng Range on which to run projection.
 * @param out Output array
 */
void gkyl_eval_offset_fd_advance(const gkyl_eval_offset_fd *up,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_eval_offset_fd_release(gkyl_eval_offset_fd *up);
