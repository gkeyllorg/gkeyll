#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_fluid_em_coupling.h>
#include <gkyl_dg_calc_fluid_em_coupling_priv.h>
#include <gkyl_util.h>

struct gkyl_dg_calc_fluid_em_coupling *
gkyl_dg_calc_fluid_em_coupling_inew(const struct gkyl_dg_calc_fluid_em_coupling_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_calc_fluid_em_coupling_cu_dev_inew(inp);
  }
#endif
  gkyl_dg_calc_fluid_em_coupling *up = gkyl_malloc(sizeof(gkyl_dg_calc_fluid_em_coupling));

  int nc = inp->conf_basis->num_basis;
  int cdim = inp->conf_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;
  enum gkyl_basis_type b_type = inp->conf_basis->b_type;
  up->mem_range = *inp->mem_range;

  up->fluid_em_coupling_set = choose_fluid_em_coupling_set_kern(b_type, cdim, poly_order);
  up->fluid_em_coupling_copy = choose_fluid_em_coupling_copy_kern(b_type, cdim, poly_order);

  // Linear system size is nc*(3*num_species + 3)
  up->num_species = inp->num_species;
  up->As = gkyl_nmat_new(
    inp->mem_range->volume, nc * (3 * up->num_species + 3), nc * (3 * up->num_species + 3)
  );
  up->xs = gkyl_nmat_new(inp->mem_range->volume, nc * (3 * up->num_species + 3), 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  // Boolean for whether or not self-consistent EM fields are static
  up->is_field_static = inp->is_field_static;

  // Needed constants for the source solve
  up->epsilon0 = inp->epsilon0;
  for (int n = 0; n < inp->num_species; ++n) {
    up->qbym[n] = inp->qbym[n];
    up->has_energy[n] = inp->has_energy ? inp->has_energy[n] : true;
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_dg_calc_fluid_em_coupling_advance(
  struct gkyl_dg_calc_fluid_em_coupling *up, double dt,
  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *ext_em,
  const struct gkyl_array *app_current, struct gkyl_array *fluid[GKYL_MAX_SPECIES],
  struct gkyl_array *em
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(em)) {
    return gkyl_dg_calc_fluid_em_coupling_advance_cu(
      up, dt, app_accel, ext_em, app_current, fluid, em
    );
  }
#endif
  int num_species = up->num_species;
  double *fluids[GKYL_MAX_SPECIES];
  const double *fluids_c[GKYL_MAX_SPECIES];
  const double *app_accels[GKYL_MAX_SPECIES];

  // Set the block-diagonal (nodal) systems of every cell, solve them in a batch, copy back.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);
    for (int n = 0; n < num_species; ++n) {
      app_accels[n] = gkyl_array_cfetch(app_accel[n], loc);
      fluids_c[n] = gkyl_array_cfetch(fluid[n], loc);
    }
    up->fluid_em_coupling_set(
      count, up->num_species, up->qbym, up->epsilon0, up->is_field_static, dt, up->As, up->xs,
      app_accels, gkyl_array_cfetch(ext_em, loc), gkyl_array_cfetch(app_current, loc), fluids_c,
      gkyl_array_cfetch(em, loc)
    );
    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);
    for (int n = 0; n < num_species; ++n) {
      fluids[n] = gkyl_array_fetch(fluid[n], loc);
    }
    up->fluid_em_coupling_copy(
      count, up->num_species, up->qbym, up->has_energy, up->epsilon0, up->xs, fluids,
      gkyl_array_fetch(em, loc)
    );
    count += 1;
  }
}

void
gkyl_dg_calc_fluid_em_coupling_release(gkyl_dg_calc_fluid_em_coupling *up)
{
  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);

  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }

  gkyl_free(up);
}
