#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_evalf_def.h>
#include <gkyl_eqn_type.h>
#include <assert.h>

// Object type.
typedef struct gkyl_bc_twistshift gkyl_bc_twistshift;

struct gkyl_bc_twistshift_inp {
  int bc_dir; // Direction in which to apply this BC.
  int shift_dir; // Direction of the shift.
  int shear_dir; // Direction in which the shift varies (shear).
  enum gkyl_edge_loc edge; // Edge to apply this BC at (lower/upper).
  int cdim; // Configuration space dimensions.
  const struct gkyl_range *bcdir_ext_update_r; // Local range where to apply BC, extended in bc_dir.
  const int *num_ghost; // Number of ghost cells in each direction.
  const struct gkyl_basis *basis; // Basis of the field shifted.
  const struct gkyl_rect_grid *grid; // Grid the field shifted is defined on.
  evalf_t shift_func; // Function defining the shift.
  void *shift_func_ctx; // Context for shift_func.
  struct gkyl_array *shift_dg; // Discretized shift.
  bool use_gpu; // Whether to apply the BC using the GPU.
  // Optional inputs:
  int shift_poly_order; // Basis order for the DG representation of the shift.
  enum gkyl_closed_flux_bc_type type; // Filtered (default) or unfiltered twist-shift.
  int filter_half_width; // Filter stencil half-width M in cells of grid (default 1).
  double filter_cutoff_wavelength; // Filter cutoff wavelength (default 2*dx along shear_dir).
  int upsample_factor; // Supersampling factor along shear_dir (default 4).
};

/**
 * Create a new updater to apply the twist-shift BC, optionally with
 * supersampling and low-pass filtering to de-alias the shifted field. It
 * combines a gkyl_twistshift_dg updater with, when requested, a
 * gkyl_dg_lowpass_filter and gkyl_dg_interpolate operators.
 *
 * @param inp bc_twistshift_inp struct containing the inputs to the updater.
 * @return New updater pointer.
 */
struct gkyl_bc_twistshift* gkyl_bc_twistshift_inew(const struct gkyl_bc_twistshift_inp *inp);

/**
 * Create a new updater to apply twist-shift BCs, passing each argument separately.
 *
 * @param bc_dir Direction in which to apply this BC.
 * @param shift_dir Direction of the shift.
 * @param shear_dir Direction in which the shift varies (shear).
 * @param edge Edge to apply this BC at (lower/upper).
 * @param cdim Configuration space dimensions.
 * @param bcdir_ext_update_r Local range where to apply BC, extended in bc_dir.
 * @param num_ghost Number of ghost cells in each direction.
 * @param basis Basis of the field shifted.
 * @param grid Grid the field shifted is defined on.
 * @param shift_func Function defining the shift.
 * @param shift_func_ctx Context for shift_func.
 * @param shift_dg Discretized shift.
 * @param shift_poly_order Basis order for the DG representation of the shift (optional).
 * @param type Filtered (default) or unfiltered twist-shift.
 * @param filter_half_width Filter stencil half-width M in cells of grid (0 = default 1).
 * @param filter_cutoff_wavelength Filter cutoff wavelength (0 = default 2*dx along shear_dir).
 * @param upsample_factor Supersampling factor along shear_dir (0 = default 4).
 * @param use_gpu Whether to apply the BC using the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_twistshift* gkyl_bc_twistshift_new(int bc_dir, int shift_dir, int shear_dir,
  enum gkyl_edge_loc edge, int cdim, const struct gkyl_range *bcdir_ext_update_r, const int *num_ghost,
  const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, evalf_t shift_func, void *shift_func_ctx,
  struct gkyl_array *shift_dg, int shift_poly_order, enum gkyl_closed_flux_bc_type type,
  int filter_half_width, double filter_cutoff_wavelength, int upsample_factor, bool use_gpu);

/**
 * Apply the twist-shift periodic BC. Expects periodicity along bc_dir to have
 * been applied to the donor field beforehand. Can be used in-place.
 *
 * @param up Twist-shift BC updater object.
 * @param fdo Donor field.
 * @param ftar Target field.
 */
void gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Free memory associated with the bc_twistshift updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up);
