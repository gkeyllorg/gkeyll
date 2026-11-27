#include <gkyl_alloc.h>
#include <gkyl_boltz_potential_priv.h>

void
boltz_potential_init(
  struct boltz_potential *up, const struct gkyl_basis *basis, double mass_e, double charge_e,
  double temp_e, bool use_gpu
)
{
  up->cdim = basis->ndim;
  up->mass_e = mass_e;
  up->charge_e = charge_e;
  up->temp_e = temp_e;
  up->use_gpu = use_gpu;
  up->on_dev = up;
}

void
boltz_potential_boundary(
  const struct boltz_potential *up, enum gkyl_edge_loc edge, const struct gkyl_range *skin,
  const struct gkyl_array *flux, const struct gkyl_array *area, const struct gkyl_array *m0_elc,
  const struct gkyl_array *wall, struct gkyl_array *end_vals
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return boltz_potential_boundary_cu(up, edge, skin, flux, area, m0_elc, wall, end_vals);
  }
#endif

  long linc = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, skin);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(skin, iter.idx);

    up->boundary[edge == GKYL_UPPER_EDGE](
      gkyl_array_cfetch(flux, linc), gkyl_array_cfetch(area, linc),
      gkyl_array_cfetch(m0_elc, linidx), wall ? gkyl_array_cfetch(wall, linidx) : NULL,
      gkyl_array_fetch(end_vals, linc)
    );
    ++linc;
  }
}

void
boltz_potential_reference(
  const struct boltz_potential *up, const struct gkyl_range *ghost, const struct gkyl_array *lower,
  const struct gkyl_array *upper, struct gkyl_array *nphi_ref
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return boltz_potential_reference_cu(up, ghost, lower, upper, nphi_ref);
  }
#endif

  long linc = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost);
  while (gkyl_range_iter_next(&iter)) {
    up->reference(
      up->charge_e, up->temp_e, up->mass_e, gkyl_array_cfetch(lower, linc),
      gkyl_array_cfetch(upper, linc), gkyl_array_fetch(nphi_ref, gkyl_range_idx(ghost, iter.idx))
    );
    ++linc;
  }
}

void
boltz_potential_phi_calc(
  const struct boltz_potential *up, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_array *m0_elc,
  const struct gkyl_array *nphi_ref, struct gkyl_array *phi
)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return boltz_potential_phi_calc_cu(up, local, local_ext, m0_elc, nphi_ref, phi);
  }
#endif

  int idx_g[GKYL_MAX_CDIM]; // Index in ghost grid nphi_ref is defined on.

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, local);
  while (gkyl_range_iter_next(&iter)) {
    // We assume each MPI rank calls this over the local range and
    // that the reference values are defined on the lower local ghost range.
    gkyl_copy_int_arr(up->cdim, iter.idx, idx_g);
    idx_g[up->cdim - 1] = local_ext->lower[up->cdim - 1];

    long loc = gkyl_range_idx(local, iter.idx);
    long ghost_loc = gkyl_range_idx(local_ext, idx_g);

    const double *m0_elc_p = (const double *)gkyl_array_cfetch(m0_elc, loc);
    const double *nphi_ref_p = (const double *)gkyl_array_cfetch(nphi_ref, ghost_loc);
    double *phi_p = (double *)gkyl_array_cfetch(phi, loc);

    up->phi(up->charge_e, up->temp_e, m0_elc_p, nphi_ref_p, phi_p);
  }
}

void
boltz_potential_release(struct boltz_potential *up)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->on_dev);
  }
#endif
}
