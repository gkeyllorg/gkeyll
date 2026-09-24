#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_eqn_type.h>

struct gkyl_dg_cx_inp {
  struct gkyl_basis* cbasis; // Configuration-space basis-functions
  const struct gkyl_range *conf_rng; // Configuration-space range
  const struct gkyl_range *conf_rng_ext; // Configuration-space extended range
  double vt_sq_ion_min; // Min vtSq that can be represented on ion grid
  double vt_sq_neut_min; // Min vtSq that can be represented on neut grid
  enum gkyl_ion_type type_ion; // Enum for type of ion for CX (H,D,HE,NE)
};

// Object type
typedef struct gkyl_dg_cx gkyl_dg_cx;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param gkyl_dg_cx_inp
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_cx* gkyl_dg_cx_new(struct gkyl_dg_cx_inp *inp, bool use_gpu); 


/**
 * Compute CX reaction rate coefficient for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param cx charge exchange object.
 * @param maxwellian_moms_ion  Ion Maxwellian moments (n, upar, T/m).
 * @param maxwellian_moms_neut Neutral Maxwellian moments (n, ux, uy, uz, T/m).
 * @param upar_b_i Ion drift velocity vector (upar b_x, upar b_y, upar b_z).
 * @param coef_cx Output reaction rate coefficient
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T]) 
 */
void gkyl_dg_cx_coll(const struct gkyl_dg_cx *up, 
  struct gkyl_array *maxwellian_moms_ion, struct gkyl_array *maxwellian_moms_neut,
  struct gkyl_array *upar_b_i, struct gkyl_array *coef_cx, struct gkyl_array *cflrate);

/**
 * Delete updater.
 *
 * @param cx Updater to delete.
 */
void gkyl_dg_cx_release(gkyl_dg_cx *cx);
