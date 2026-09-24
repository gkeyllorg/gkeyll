#pragma once

#include <gkyl_eqn_type.h>
#include <read_adas.h>

// Private header, not for direct use in user code

// Primary struct in this updater.
struct gkyl_dg_recomb {
  const struct gkyl_rect_grid *grid; // conf grid object
  int cdim; // number of configuration space dimensions
  int poly_order; // polynomial order of DG basis

  const struct gkyl_range *conf_rng; // Configuration-space range
  const struct gkyl_range *conf_rng_ext; // Configuration-space extended range
  const struct gkyl_range *phase_rng; // Phase-space range
  bool use_gpu;
  
  double elem_charge; // elementary charge value
  double mass_elc; // mass of the electron

  int vdim;
  int resM0;
  int resTe;
  double dlogM0;
  double dlogTe;
  double minLogM0, minLogTe, maxLogM0, maxLogTe;

  enum gkyl_react_self_type type_self;

  struct gkyl_basis *cbasis;
  struct gkyl_basis *pbasis;

  struct gkyl_range adas_rng;
  struct gkyl_basis adas_basis;
  struct gkyl_basis *basis_on_dev;
  
  struct gkyl_array *recomb_data;
  struct gkyl_array *vtSq_elc;

  struct gkyl_dg_recomb *on_dev; // pointer to itself or device data
};


