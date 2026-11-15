/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_em_coupling.h>
#include <gkyl_dg_calc_pkpm_em_coupling_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_pkpm_em_coupling_set_cu_kernel(
  gkyl_dg_calc_pkpm_em_coupling *up, struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_range conf_range, double dt, const struct gkyl_array *app_accel0,
  const struct gkyl_array *app_accel1, const struct gkyl_array *ext_em, const struct gkyl_array *app_current,
  const struct gkyl_array *fluid0, const struct gkyl_array *fluid1, const struct gkyl_array *em
)
{
  const double *app_accels[GKYL_MAX_SPECIES];
  const double *fluids[GKYL_MAX_SPECIES];
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    app_accels[0] = (const double *)gkyl_array_cfetch(app_accel0, loc);
    fluids[0] = (const double *)gkyl_array_cfetch(fluid0, loc);
    if (up->num_species == 2) {
      app_accels[1] = (const double *)gkyl_array_cfetch(app_accel1, loc);
      fluids[1] = (const double *)gkyl_array_cfetch(fluid1, loc);
    }
    up->pkpm_em_coupling_set(
      linc1, up->num_species, up->qbym, up->epsilon0, up->pkpm_field_static, dt, As, xs, app_accels,
      (const double *)gkyl_array_cfetch(ext_em, loc), (const double *)gkyl_array_cfetch(app_current, loc),
      fluids, (const double *)gkyl_array_cfetch(em, loc)
    );
  }
}

__global__ static void
gkyl_dg_calc_pkpm_em_coupling_copy_cu_kernel(
  gkyl_dg_calc_pkpm_em_coupling *up, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array *fluid0, struct gkyl_array *fluid1, struct gkyl_array *em
)
{
  double *fluids[GKYL_MAX_SPECIES];
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    fluids[0] = (double *)gkyl_array_fetch(fluid0, loc);
    if (up->num_species == 2) {
      fluids[1] = (double *)gkyl_array_fetch(fluid1, loc);
    }
    up->pkpm_em_coupling_copy(
      linc1, up->num_species, up->qbym, up->epsilon0, xs, fluids, (double *)gkyl_array_fetch(em, loc)
    );
  }
}

// Host-side wrapper for the nodal implicit source update (one or two species)
void
gkyl_dg_calc_pkpm_em_coupling_advance_cu(
  struct gkyl_dg_calc_pkpm_em_coupling *up, double dt,
  const struct gkyl_array *app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *ext_em,
  const struct gkyl_array *app_current, struct gkyl_array *fluid[GKYL_MAX_SPECIES], struct gkyl_array *em
)
{
  struct gkyl_range conf_range = up->mem_range;
  int ns = up->num_species;
  assert(ns == 1 || ns == 2);
  int s1 = ns == 2 ? 1 : 0;
  gkyl_dg_calc_pkpm_em_coupling_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(
    up->on_dev, up->As->on_dev, up->xs->on_dev, conf_range, dt, app_accel[0]->on_dev, app_accel[s1]->on_dev,
    ext_em->on_dev, app_current->on_dev, fluid[0]->on_dev, fluid[s1]->on_dev,
    em->on_dev
  );
  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);
  gkyl_dg_calc_pkpm_em_coupling_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(
    up->on_dev, up->xs->on_dev, conf_range, fluid[0]->on_dev, fluid[s1]->on_dev, em->on_dev
  );
}

__global__ static void
dg_calc_pkpm_em_coupling_set_cu_dev_ptrs(
  struct gkyl_dg_calc_pkpm_em_coupling *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  up->pkpm_em_coupling_set = choose_pkpm_em_coupling_set_kern(b_type, cdim, poly_order);
  up->pkpm_em_coupling_copy = choose_pkpm_em_coupling_copy_kern(b_type, cdim, poly_order);
}

gkyl_dg_calc_pkpm_em_coupling *
gkyl_dg_calc_pkpm_em_coupling_cu_dev_inew(const struct gkyl_dg_calc_pkpm_em_coupling_inp *inp)
{
  struct gkyl_dg_calc_pkpm_em_coupling *up =
    (struct gkyl_dg_calc_pkpm_em_coupling *)gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_em_coupling));

  int nc = inp->conf_basis->num_basis;
  int cdim = inp->conf_basis->ndim;
  int poly_order = inp->conf_basis->poly_order;
  enum gkyl_basis_type b_type = inp->conf_basis->b_type;
  up->mem_range = *inp->mem_range;

  // Linear system size is nc*(3*num_species + 3)
  up->num_species = inp->num_species;
  up->As = gkyl_nmat_cu_dev_new(
    inp->mem_range->volume, nc * (3 * up->num_species + 3), nc * (3 * up->num_species + 3)
  );
  up->xs = gkyl_nmat_cu_dev_new(inp->mem_range->volume, nc * (3 * up->num_species + 3), 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  // Boolean for whether or not self-consistent EM fields are static
  up->pkpm_field_static = inp->pkpm_field_static;

  // Needed constants for the source solve
  up->epsilon0 = inp->epsilon0;
  for (int n = 0; n < inp->num_species; ++n) {
    up->qbym[n] = inp->qbym[n];
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_em_coupling *up_cu =
    (struct gkyl_dg_calc_pkpm_em_coupling *)gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_em_coupling));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_em_coupling), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_em_coupling_set_cu_dev_ptrs<<<1, 1>>>(up_cu, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  return up;
}
