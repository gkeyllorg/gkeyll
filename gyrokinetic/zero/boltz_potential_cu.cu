extern "C" {
#include <gkyl_boltz_potential_priv.h>
}

// Evaluate the selected face and pack compact surface data.
__global__ static void
boltz_potential_boundary_cu_ker(
  const struct boltz_potential *up, int edge, struct gkyl_range skin, const struct gkyl_array *flux,
  const struct gkyl_array *area, const struct gkyl_array *m0_elc, const struct gkyl_array *wall,
  struct gkyl_array *end_vals
)
{
  int idx[GKYL_MAX_CDIM];
  for (long i = threadIdx.x + blockIdx.x * blockDim.x; i < skin.volume;
       i += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&skin, i, idx);
    long loc = gkyl_range_idx(&skin, idx);
    up->boundary[edge](
      (const double *)gkyl_array_cfetch(flux, i), (const double *)gkyl_array_cfetch(area, i),
      (const double *)gkyl_array_cfetch(m0_elc, loc),
      wall ? (const double *)gkyl_array_cfetch(wall, loc) : NULL,
      (double *)gkyl_array_fetch(end_vals, i)
    );
  }
}

void
boltz_potential_boundary_cu(
  const struct boltz_potential *up, enum gkyl_edge_loc edge, const struct gkyl_range *skin,
  const struct gkyl_array *flux, const struct gkyl_array *area, const struct gkyl_array *m0_elc,
  const struct gkyl_array *wall, struct gkyl_array *end_vals
)
{
  boltz_potential_boundary_cu_ker<<<skin->nblocks, skin->nthreads>>>(
    up->on_dev, edge == GKYL_UPPER_EDGE, *skin, flux->on_dev, area->on_dev, m0_elc->on_dev,
    wall ? wall->on_dev : NULL, end_vals->on_dev
  );
}

// Calculate model-specific reference data in the local lower ghost slab.
__global__ static void
boltz_potential_reference_cu_ker(
  const struct boltz_potential *up, struct gkyl_range ghost, const struct gkyl_array *lower,
  const struct gkyl_array *upper, struct gkyl_array *nphi_ref
)
{
  int idx[GKYL_MAX_CDIM];
  for (long i = threadIdx.x + blockIdx.x * blockDim.x; i < ghost.volume;
       i += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&ghost, i, idx);
    up->reference(
      up->charge_e, up->temp_e, up->mass_e, (const double *)gkyl_array_cfetch(lower, i),
      (const double *)gkyl_array_cfetch(upper, i),
      (double *)gkyl_array_fetch(nphi_ref, gkyl_range_idx(&ghost, idx))
    );
  }
}

void
boltz_potential_reference_cu(
  const struct boltz_potential *up, const struct gkyl_range *ghost, const struct gkyl_array *lower,
  const struct gkyl_array *upper, struct gkyl_array *nphi_ref
)
{
  boltz_potential_reference_cu_ker<<<ghost->nblocks, ghost->nthreads>>>(
    up->on_dev, *ghost, lower->on_dev, upper->on_dev, nphi_ref->on_dev
  );
}

// Apply the common volume response using reference data at the same transverse index.
__global__ static void
boltz_potential_phi_cu_ker(
  const struct boltz_potential *up, struct gkyl_range local, struct gkyl_range local_ext,
  const struct gkyl_array *m0_elc, const struct gkyl_array *nphi_ref, struct gkyl_array *phi
)
{
  int idx[GKYL_MAX_CDIM];
  for (long i = threadIdx.x + blockIdx.x * blockDim.x; i < local.volume;
       i += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&local, i, idx);
    long loc = gkyl_range_idx(&local, idx);
    idx[up->cdim - 1] = local_ext.lower[up->cdim - 1];
    long ghost_loc = gkyl_range_idx(&local_ext, idx);
    up->phi(
      up->charge_e, up->temp_e, (const double *)gkyl_array_cfetch(m0_elc, loc),
      (const double *)gkyl_array_cfetch(nphi_ref, ghost_loc), (double *)gkyl_array_fetch(phi, loc)
    );
  }
}

void
boltz_potential_phi_calc_cu(
  const struct boltz_potential *up, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_array *m0_elc,
  const struct gkyl_array *nphi_ref, struct gkyl_array *phi
)
{
  boltz_potential_phi_cu_ker<<<local->nblocks, local->nthreads>>>(
    up->on_dev, *local, *local_ext, m0_elc->on_dev, nphi_ref->on_dev, phi->on_dev
  );
}
