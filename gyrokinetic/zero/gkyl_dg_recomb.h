#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_dg_recomb_inp {
  const struct gkyl_rect_grid* grid; // Grid object needed for fmax
  struct gkyl_basis* cbasis; // Configuration-space basis-functions
  struct gkyl_basis* pbasis; // Phase-space basis-functions
  const struct gkyl_range *conf_rng; // Configuration-space range
  const struct gkyl_range *conf_rng_ext; // Configuration-space extended range
  const struct gkyl_range *phase_rng; // Phase-space range
  double mass_self; // Mass of the species
  enum gkyl_ion_type type_ion; // Enum for type of ion for ionization (H thru 0)
  int charge_state; // Ion charge state
  enum gkyl_react_self_type type_self; // Species type (ion, electron, or receiver)
};


// Object type
typedef struct gkyl_dg_recomb gkyl_dg_recomb;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param gkyl_dg_recomb_inp
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_recomb* gkyl_dg_recomb_new(struct gkyl_dg_recomb_inp *inp, bool use_gpu); 

/**
 * Compute recombination collision term for use in neutral reactions. 
 *
 * @param recomb Recombation object.
 * @param prim_vars_elc (n, upar, T/m) for the electrons
 * @param coef_recomb Output reaction rate coefficient
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T]) 
 */
void gkyl_dg_recomb_coll(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *prim_vars_elc, 
  struct gkyl_array *coef_recomb, struct gkyl_array *cflrate);

void gkyl_dg_recomb_coll_cu(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *prim_vars_elc, 
  struct gkyl_array *coef_recomb, struct gkyl_array *cflrate);

/**
 * Delete updater.
 *
 * @param recomb Updater to delete.
 */
void gkyl_dg_recomb_release(gkyl_dg_recomb *recomb);
