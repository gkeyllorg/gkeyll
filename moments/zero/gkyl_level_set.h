#include <float.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_null_comm.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_prop.h>

/**
 * Reinitialize the level set function in the Riemann ghost fluid method for the Euler equations.
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
euler_rgfm_reinit_level_set(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic Maxwell equations.
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_maxwell_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic Maxwell equations in the tetrad basis.
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_maxwell_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic Euler equations (general equation of state).
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_euler_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic Euler equations in the tetrad basis (general equation of state).
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_euler_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic Euler equations (ultra-relativistic equation of state).
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_ultra_rel_euler_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic Euler equations in the tetrad basis (ultra-relativistic equation of state).
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_ultra_rel_euler_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic two-fluid equations.
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
 void
 gr_twofluid_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
   struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic two-fluid equations in the tetrad basis.
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
 void
 gr_twofluid_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
   struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic magnetohydrodynamics equations.
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_mhd_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);

/**
 * Impose gauge conditions for the general relativistic magnetohydrodynamics equations in the tetrad basis.
 *
 * @param wv Wave propagation object.
 * @param update_range Range of cells to be updated.
 * @param idxl Index of cell(s) to update.
 * @param loidx_c Lower index of cells to update.
 * @param upidx_c Upper index of cells to update.
 * @param qout Output array of fluid variables.
 * @param dir Direction in which to perform the update.
 */
void
gr_mhd_tetrad_impose_gauge(gkyl_wave_prop *wv, const struct gkyl_range *update_range, int idxl[GKYL_MAX_DIM], int loidx_c, int upidx_c,
  struct gkyl_array *qout, int dir);
