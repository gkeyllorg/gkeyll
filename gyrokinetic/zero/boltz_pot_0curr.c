#include <gkyl_alloc.h>
#include <gkyl_boltz_pot_0curr_priv.h>

gkyl_boltz_pot_0curr *
gkyl_boltz_pot_0curr_new(
  const struct gkyl_basis *basis, double mass_e, double charge_e, double temp_e, bool use_gpu
)
{
  assert(basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
  assert(basis->ndim >= 1 && basis->ndim <= 3);
  assert(basis->poly_order >= 1 && basis->poly_order <= 2);
  gkyl_boltz_pot_0curr *up = gkyl_malloc(sizeof(*up));
  boltz_potential_init(&up->common, basis, mass_e, charge_e, temp_e, use_gpu);
  boltz_potential_set_kernels(
    &up->common, basis->poly_order,
    boltz_pot_0curr_reference_kernels[basis->ndim - 1][basis->poly_order - 1]
  );
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    boltz_pot_0curr_cu_dev_new(up, basis->poly_order);
  }
#else
  assert(!use_gpu);
#endif
  return up;
}

void
gkyl_boltz_pot_0curr_boundary(
  const gkyl_boltz_pot_0curr *up, enum gkyl_edge_loc edge, const struct gkyl_range *skin,
  const struct gkyl_array *flux, const struct gkyl_array *area, const struct gkyl_array *m0_elc,
  const struct gkyl_array *wall, struct gkyl_array *end_vals
)
{
  boltz_potential_boundary(&up->common, edge, skin, flux, area, m0_elc, wall, end_vals);
}

void
gkyl_boltz_pot_0curr_reference(
  const gkyl_boltz_pot_0curr *up, const struct gkyl_range *ghost, const struct gkyl_array *lower,
  const struct gkyl_array *upper, struct gkyl_array *nphi_ref
)
{
  boltz_potential_reference(&up->common, ghost, lower, upper, nphi_ref);
}

void
gkyl_boltz_pot_0curr_phi_calc(
  const gkyl_boltz_pot_0curr *up, const struct gkyl_range *local,
  const struct gkyl_range *local_ext, const struct gkyl_array *m0_elc,
  const struct gkyl_array *nphi_ref, struct gkyl_array *phi
)
{
  boltz_potential_phi_calc(&up->common, local, local_ext, m0_elc, nphi_ref, phi);
}

void
gkyl_boltz_pot_0curr_release(gkyl_boltz_pot_0curr *up)
{
  boltz_potential_release(&up->common);
  gkyl_free(up);
}
