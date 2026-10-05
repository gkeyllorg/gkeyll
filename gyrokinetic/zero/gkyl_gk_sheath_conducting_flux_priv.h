#pragma once

#include <gkyl_bc_sheath_gyrokinetic_priv.h>
#include <gkyl_gk_collisionless_flux_priv.h>
#include <gkyl_gk_sheath_conducting_flux.h>

struct gkyl_gk_sheath_conducting_flux {
  int dir;
  int cdim;
  enum gkyl_edge_loc edge;
  double q2Dm;
  bool use_gpu;
  struct gkyl_rect_grid phase_grid;
  const struct gkyl_basis *phase_basis;
  const struct gkyl_basis *surf_basis;
  struct gkyl_range conf_range;
  struct gkyl_range conf_ext_range;
  struct gkyl_range phase_ext_range;
  struct gkyl_range skin_range;
  const struct gkyl_gk_collisionless_flux *flux_op;
  sheath_reflectedf_t reflectedf;
  struct gkyl_gk_sheath_conducting_flux *on_dev;
};

#ifdef GKYL_HAVE_CUDA
void gkyl_gk_sheath_conducting_flux_advance_cu(
  const struct gkyl_gk_sheath_conducting_flux *up, const struct gkyl_array *phi,
  const struct gkyl_array *gyro_phi, const struct gkyl_array *phi_wall,
  const struct gkyl_array *fin, struct gkyl_array *flux_surf
);

void gkyl_gk_sheath_conducting_flux_set_cu_ptrs(
  struct gkyl_gk_sheath_conducting_flux *up, const struct gkyl_basis *phase_basis
);
#endif
