#pragma once

#include <gkyl_eqn_type.h>
#include <read_adas.h>

// Private header, not for direct use in user code

// Primary struct in this updater.
struct gkyl_dg_iz {
  int cdim; // number of configuration space dimensions
  int poly_order; // polynomial order of DG basis

  const struct gkyl_range *conf_rng; // Configuration-space range
  bool use_gpu;
  
  double elem_charge; // elementary charge value
  double mass_elc; // mass of the electron

  int resM0;
  int resTe;
  double dlogM0;
  double dlogTe;
  double E;

  double minLogM0, minLogTe, maxLogM0, maxLogTe;

  enum gkyl_react_self_type type_self;

  struct gkyl_array *ioniz_data;
  struct gkyl_range adas_rng;
  struct gkyl_basis adas_basis;
  struct gkyl_basis *basis_on_dev;

  struct gkyl_dg_iz *on_dev; // pointer to itself or device data
};

#ifdef GKYL_HAVE_CUDA
/**
 * Create new ionization updater type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_iz* gkyl_dg_iz_cu_dev_new(struct gkyl_dg_iz_inp *inp);

/**
 * Compute ionization collision term for use in neutral reactions. 
 * 
 *
 * @param iz Ionization object.
 * @param maxwellian_moms_elc Electron Maxwellian moments (n, upar, T/m).
 * @param vtSq_iz1 First thermal Speed for ionization fmax (primary electrons).
 * @param vtSq_iz2 Second thermal Speed for ionization fmax (secondary electrons).
 * @param coef_iz Output reaction rate coefficient.
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T]).
 */
void gkyl_dg_iz_coll_cu(const struct gkyl_dg_iz *up, 
  const struct gkyl_array *maxwellian_moms_elc, 
  struct gkyl_array *vtSq_iz1, struct gkyl_array *vtSq_iz2,
  struct gkyl_array *coef_iz, struct gkyl_array *cflrate);
#endif
