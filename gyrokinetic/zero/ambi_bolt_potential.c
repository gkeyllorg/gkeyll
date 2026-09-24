#include <gkyl_ambi_bolt_potential.h>
#include <gkyl_ambi_bolt_potential_priv.h>
#include <gkyl_alloc.h>

gkyl_ambi_bolt_potential*
gkyl_ambi_bolt_potential_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  double mass_e, double charge_e, double temp_e, bool use_gpu)
{
  struct gkyl_ambi_bolt_potential *up = gkyl_malloc(sizeof(struct gkyl_ambi_bolt_potential));

  up->cdim = basis->ndim;
  up->num_basis = basis->num_basis;
  up->use_gpu = use_gpu;

  up->dz = grid->dx[up->cdim-1];
  up->mass_e = mass_e;
  up->charge_e = charge_e;
  up->temp_e = temp_e;

  up->kernels = gkyl_malloc(sizeof(struct gkyl_ambi_bolt_potential_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_ambi_bolt_potential_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
#endif

  // Select sheath_calc and phi_calc kernels.
  ambi_bolt_potential_choose_kernels(basis, up->kernels);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    ambi_bolt_potential_choose_kernels_cu(basis, up->kernels_cu);
#endif

  return up;
}

void gkyl_ambi_bolt_potential_sheath_calc(struct gkyl_ambi_bolt_potential *up, enum gkyl_edge_loc edge,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_array *cmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *gammai, const struct gkyl_array *m0i, const struct gkyl_array *Jm0i,
  struct gkyl_array *sheath_vals)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_ambi_bolt_potential_sheath_calc_cu(up, edge, skin_r, ghost_r,
      cmag, jacobtot_inv, gammai, m0i, Jm0i, sheath_vals);
#endif

  unsigned int keridx = (edge == GKYL_LOWER_EDGE) ? 0 : 1;

  int idx_s[GKYL_MAX_CDIM]; // Skin index.

  // MF 2022/02/27: Although we wish to obtain skin-cell quantities (sheath
  // potential and density) we loop over the ghost cell because it is easier
  // given that ghost grids exist in the app, but we have not created skin
  // grids. And gammai and sheath_vals are defined on the ghost grid.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(up->cdim, iter.idx, idx_s);
    // Assume only 1 ghost cell on either side along the field line.
    idx_s[up->cdim-1] = edge == GKYL_LOWER_EDGE ? iter.idx[up->cdim-1]+1 : iter.idx[up->cdim-1]-1;

    long ghost_loc = gkyl_range_idx(ghost_r, iter.idx);
    long skin_loc = gkyl_range_idx(skin_r, idx_s);

    const double *cmag_p = (const double*) gkyl_array_cfetch(cmag, skin_loc);
    const double *jacobtotinv_p = (const double*) gkyl_array_cfetch(jacobtot_inv, skin_loc);
    const double *m0i_p = (const double*) gkyl_array_cfetch(m0i, skin_loc);
    const double *Jm0i_p = (const double*) gkyl_array_cfetch(Jm0i, skin_loc);
    const double *gammai_p = (const double*) gkyl_array_cfetch(gammai, ghost_loc);
    double *out_p = (double*) gkyl_array_cfetch(sheath_vals, ghost_loc);

    up->kernels->sheath_calc[keridx](up->dz, up->charge_e, up->mass_e, up->temp_e,
      cmag_p, jacobtotinv_p, gammai_p, m0i_p, Jm0i_p, out_p);
  }
}

void gkyl_ambi_bolt_potential_phi_calc(struct gkyl_ambi_bolt_potential *up,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  const struct gkyl_array *m0i, const struct gkyl_array *sheath_vals, struct gkyl_array *phi)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_ambi_bolt_potential_phi_calc_cu(up, local, local_ext, m0i, sheath_vals, phi);
#endif

  int idx_g[GKYL_MAX_CDIM]; // Index in ghost grid sheath_vals is defined on.

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, local);
  while (gkyl_range_iter_next(&iter)) {
    // We assume each MPI rank calls this over the local range and
    // that the sheath values are defined on the lower local ghost range.
    gkyl_copy_int_arr(up->cdim, iter.idx, idx_g);
    idx_g[up->cdim-1] = local_ext->lower[up->cdim-1];

    long loc = gkyl_range_idx(local, iter.idx);
    long ghost_loc = gkyl_range_idx(local_ext, idx_g);

    const double *m0i_p = (const double*) gkyl_array_cfetch(m0i, loc);
    const double *sheathvals_p = (const double*) gkyl_array_cfetch(sheath_vals, ghost_loc);
    double *phi_p = (double*) gkyl_array_cfetch(phi, loc);

    up->kernels->phi_calc(up->charge_e, up->temp_e, m0i_p, sheathvals_p, phi_p);
  }
}

void gkyl_ambi_bolt_potential_release(gkyl_ambi_bolt_potential *up)
{
  gkyl_free(up->kernels);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels_cu);
#endif
  gkyl_free(up);
}
