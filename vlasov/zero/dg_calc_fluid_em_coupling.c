#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_fluid_em_coupling.h>
#include <gkyl_dg_calc_fluid_em_coupling_priv.h>
#include <gkyl_util.h>

struct gkyl_dg_calc_fluid_em_coupling*
gkyl_dg_calc_fluid_em_coupling_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range *mem_range, 
  int num_fluids, double qbym[GKYL_MAX_SPECIES], double epsilon0,
  bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_fluid_em_coupling_cu_dev_new(cbasis, mem_range, 
      num_fluids, qbym, epsilon0);
  } 
#endif     
  gkyl_dg_calc_fluid_em_coupling *up = gkyl_malloc(sizeof(gkyl_dg_calc_fluid_em_coupling));

  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->mem_range = *mem_range;

  up->fluid_em_coupling_set = choose_fluid_em_coupling_set_kern(b_type, cdim, poly_order);
  up->fluid_em_coupling_copy = choose_fluid_em_coupling_copy_kern(b_type, cdim, poly_order);
  up->fluid_em_coupling_energy = choose_fluid_em_coupling_energy_kern(b_type, cdim, poly_order);

  // Linear system size is nc*(3*num_fluids + 3)
  up->num_fluids = num_fluids;
  up->As = gkyl_nmat_new(mem_range->volume, nc*(3*up->num_fluids + 3), nc*(3*up->num_fluids + 3));
  up->xs = gkyl_nmat_new(mem_range->volume, nc*(3*up->num_fluids + 3), 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  // Needed constants for the source solve
  up->epsilon0 = epsilon0;
  for (int n = 0; n < num_fluids; ++n) {
    up->qbym[n] = qbym[n];
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void 
gkyl_dg_calc_fluid_em_coupling_advance(struct gkyl_dg_calc_fluid_em_coupling *up, double dt, 
  const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], struct gkyl_array* em)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(em)) {
    return gkyl_dg_calc_fluid_em_coupling_advance_cu(up, 
      dt, app_accel, ext_em, app_current, fluid, em);
  }
#endif
  int num_fluids = up->num_fluids;
  double *fluids[GKYL_MAX_SPECIES];
  const double *app_accels[GKYL_MAX_SPECIES];

  // First loop over mem_range for solving linear systems to compute primitive moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    for (int n=0; n<num_fluids; ++n) {
      fluids[n] = gkyl_array_fetch(fluid[n], loc);
      app_accels[n] = gkyl_array_cfetch(app_accel[n], loc);
    }
    const double *ext_em_d = gkyl_array_cfetch(ext_em, loc);
    const double *app_current_d = gkyl_array_cfetch(app_current, loc);
    double *em_d = gkyl_array_fetch(em, loc);

    up->fluid_em_coupling_set(count, up->num_fluids, up->qbym, up->epsilon0, dt, 
      up->As, up->xs, 
      app_accels, ext_em_d, app_current_d,
      fluids, em_d);

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    for (int n=0; n<num_fluids; ++n) {
      fluids[n] = gkyl_array_fetch(fluid[n], loc);
    }
    double *em_d = gkyl_array_fetch(em, loc);

    up->fluid_em_coupling_copy(count, up->num_fluids, up->qbym, up->epsilon0, 
      up->xs, fluids, em_d);

    count += 1;
  }
}

void 
gkyl_dg_calc_fluid_em_coupling_energy(struct gkyl_dg_calc_fluid_em_coupling *up, 
  const struct gkyl_array* ke_old, const struct gkyl_array* ke_new, struct gkyl_array* fluid)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid)) {
    return gkyl_dg_calc_fluid_em_coupling_energy_cu(up, 
      ke_old, ke_new, fluid);
  }
#endif  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *ke_old_d = gkyl_array_cfetch(ke_old, loc);
    const double *ke_new_d = gkyl_array_cfetch(ke_new, loc);

    double *fluid_d = gkyl_array_fetch(fluid, loc);

    up->fluid_em_coupling_energy(ke_old_d, ke_new_d, fluid_d);
  }
}

void 
gkyl_dg_calc_fluid_em_coupling_release(gkyl_dg_calc_fluid_em_coupling *up)
{
  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);
  
  if (GKYL_IS_CU_ALLOC(up->flags)) 
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
