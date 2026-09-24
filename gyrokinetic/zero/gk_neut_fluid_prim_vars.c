#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_gk_neut_fluid_prim_vars.h>
#include <gkyl_gk_neut_fluid_prim_vars_priv.h>
#include <gkyl_util.h>

void gkyl_gk_neut_fluid_prim_vars_udrift_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = up->udrift_ncomp;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_udrift_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_set_prob_ker(count, up->As, up->xs, moms_d);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[up->udrift_ncomp*up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->udrift_get_sol_ker(count, up->xs, prim_vars_d);

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }

    count += nprob;
  }
}

void gkyl_gk_neut_fluid_prim_vars_pressure_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = up->udrift_ncomp;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_pressure_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_set_prob_ker(count, up->As, up->xs, moms_d);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    double udrift_d[up->udrift_ncomp*up->num_basis];
    for (int i=0; i<up->num_basis; i++)
      udrift_d[i] = 0.0;

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->udrift_get_sol_ker(count, up->xs, udrift_d);
    up->pressure_ker(up->gas_gamma, moms_d, udrift_d, prim_vars_d);

    for (int i=0; i<up->num_basis; i++)
      prim_vars_d[i] *= up->thermalE_fac;

    if (up->is_integrated) {
      out_d[out_coff] = up->integrated_fac * prim_vars_d[0];
    }

    count += nprob;
  }
}

void gkyl_gk_neut_fluid_prim_vars_temp_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = 1;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_temp_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->temp_set_prob_ker(count, up->As, up->xs, moms_d, up->gas_gamma, up->mass);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->temp_get_sol_ker(count, up->xs, prim_vars_d);

    if (up->is_integrated) {
      out_d[out_coff] = up->integrated_fac * prim_vars_d[0];
    }

    count += nprob;
  }
}

void gkyl_gk_neut_fluid_prim_vars_udrift_pressure_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = up->udrift_ncomp;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_udrift_pressure_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_set_prob_ker(count, up->As, up->xs, moms_d);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[(up->udrift_ncomp+1)*up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    double* udrift_d = prim_vars_d;
    double* pressure_d = &prim_vars_d[up->udrift_ncomp*up->num_basis];

    up->udrift_get_sol_ker(count, up->xs, udrift_d);
    up->pressure_ker(up->gas_gamma, moms_d, udrift_d, pressure_d);

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+1; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }
    count += nprob;
  }
}

void gkyl_gk_neut_fluid_prim_vars_udrift_temp_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = up->udrift_ncomp+1;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_udrift_temp_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_temp_set_prob_ker(count, up->As, up->xs, moms_d, up->gas_gamma, up->mass);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[(up->udrift_ncomp+1)*up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->udrift_temp_get_sol_ker(count, up->xs, prim_vars_d);

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+1; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }
    count += nprob;
  }
}

void gkyl_gk_neut_fluid_prim_vars_lte_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = up->udrift_ncomp+1;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_lte_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->udrift_temp_set_prob_ker(count, up->As, up->xs, moms_d, up->gas_gamma, up->mass);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[(up->udrift_ncomp+2)*up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];


    up->udrift_temp_get_sol_ker(count, up->xs, &prim_vars_d[up->num_basis]);

    // Scale rho and temp by 1/m.
    for (int i=0; i<up->num_basis; i++) {
      prim_vars_d[out_coff+i] = moms_d[i]/up->mass;
      prim_vars_d[out_coff+(up->udrift_ncomp+1)*up->num_basis+i] *= 1.0/up->mass;
    }

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+2; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }

    count += nprob;
  }
}

void gkyl_gk_neut_fluid_prim_vars_flow_energy_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = 1;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_flow_energy_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->flowE_set_prob_ker(count, up->As, up->xs, moms_d);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    up->flowE_get_sol_ker(count, up->xs, prim_vars_d);

    if (up->is_integrated) {
      out_d[out_coff] = up->integrated_fac * prim_vars_d[0];
    }

    count += nprob;
  }
}

void gkyl_gk_neut_fluid_prim_vars_mass_momentum_flow_thermal_energy_advance(struct gkyl_gk_neut_fluid_prim_vars *up,
  const struct gkyl_array* moms, struct gkyl_array *out, int out_coff)
{
  int nprob = 1;
  assert(up->As->num == nprob*up->mem_range.volume);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_gk_neut_fluid_prim_vars_mass_momentum_flow_thermal_energy_advance_cu(up, moms, out, out_coff);
  }
#endif

  // Loop over mem_range for solving linear systems to compute primitive moments.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);

    up->flowE_set_prob_ker(count, up->As, up->xs, moms_d);

    count += nprob;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  double prim_vars_buff[(up->udrift_ncomp+3)*up->num_basis]; // Buffer needed for integrated moms.
  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *moms_d = gkyl_array_cfetch(moms, linidx);
    double* out_d = gkyl_array_fetch(out, linidx);

    double* prim_vars_d = up->is_integrated? prim_vars_buff : &out_d[out_coff];

    int fourth_comp_off = 4*up->num_basis;
    for (int i=0; i<fourth_comp_off; i++)
      prim_vars_d[i] = moms_d[i];

    up->flowE_get_sol_ker(count, up->xs, &prim_vars_d[fourth_comp_off]);

    for (int i=0; i<up->num_basis; i++)
      prim_vars_d[5*up->num_basis+i] = moms_d[fourth_comp_off+i] - prim_vars_d[fourth_comp_off+i];

    if (up->is_integrated) {
      for (int i=0; i<up->udrift_ncomp+3; i++)
        out_d[out_coff+i] = up->integrated_fac * prim_vars_d[i*up->num_basis];
    }

    count += nprob;
  }
}

gkyl_gk_neut_fluid_prim_vars*
gkyl_gk_neut_fluid_prim_vars_new(double gas_gamma, double mass, const struct gkyl_basis* cbasis,
  struct gkyl_rect_grid *grid, const struct gkyl_range *mem_range,
  enum gkyl_gk_neut_fluid_prim_vars_type prim_vars_type, bool is_integrated, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_gk_neut_fluid_prim_vars_cu_dev_new(gas_gamma, mass, cbasis, grid, mem_range, prim_vars_type, is_integrated);
  }
#endif
  gkyl_gk_neut_fluid_prim_vars *up = gkyl_malloc(sizeof(gkyl_gk_neut_fluid_prim_vars));

  up->gas_gamma = gas_gamma;
  up->mass = mass;

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->num_basis = cbasis->num_basis;
  up->udrift_ncomp = 3;
  up->mem_range = *mem_range;
  up->is_integrated = is_integrated;
  
  up->integrated_fac = 0.0;
  if (up->is_integrated) {
    up->integrated_fac = 1.0;
    for (int d=0; d<up->cdim; d++)
      up->integrated_fac *= (grid->dx[d]/2.0)*sqrt(2.0);
  }

  int nprob;
  up->thermalE_fac = 0.0;
  if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT) {
    nprob = up->udrift_ncomp;
  }
  else if ( (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_PRESSURE) ||
            (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_THERMAL_ENERGY) ) {
    nprob = up->udrift_ncomp;
    up->thermalE_fac = prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_PRESSURE? 1.0 : 1.0/(up->gas_gamma-1.0);
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_TEMP) {
    nprob = 1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_PRESSURE) {
    nprob = up->udrift_ncomp;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_UDRIFT_TEMP) {
    nprob = up->udrift_ncomp+1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_LTE) {
    nprob = up->udrift_ncomp+1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_FLOW_ENERGY) {
    nprob = 1;
  }
  else if (prim_vars_type == GKYL_GK_NEUT_FLUID_PRIM_VARS_MASS_MOMENTUM_FLOW_THERMAL_ENERGY) {
    nprob = 1;
  }

  // There are udrift_ncomp*range->volume linear systems to be solved
  // for 3 components of u: ux, uy, uz.
  up->As = gkyl_nmat_new(nprob*mem_range->volume, up->num_basis, up->num_basis);
  up->xs = gkyl_nmat_new(nprob*mem_range->volume, up->num_basis, 1);
  if (up->poly_order > 1) {
    up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);
  }

  // Assign all kernels so one can use the updater for multiple calc types.
  up->udrift_set_prob_ker = choose_udrift_set_prob_ker(b_type, cdim, poly_order);
  up->udrift_get_sol_ker = choose_udrift_get_sol_ker(b_type, cdim, poly_order);
  up->pressure_ker = choose_pressure_ker(b_type, cdim, poly_order);
  up->temp_set_prob_ker = choose_temp_set_prob_ker(b_type, cdim, poly_order);
  up->temp_get_sol_ker = choose_temp_get_sol_ker(b_type, cdim, poly_order);
  up->udrift_temp_set_prob_ker = choose_udrift_temp_set_prob_ker(b_type, cdim, poly_order);
  up->udrift_temp_get_sol_ker = choose_udrift_temp_get_sol_ker(b_type, cdim, poly_order);
  up->flowE_set_prob_ker = choose_flowE_set_prob_ker(b_type, cdim, poly_order);
  up->flowE_get_sol_ker = choose_flowE_get_sol_ker(b_type, cdim, poly_order);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void gkyl_gk_neut_fluid_prim_vars_release(gkyl_gk_neut_fluid_prim_vars *up)
{
  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  if (up->poly_order > 1) {
    gkyl_nmat_linsolve_lu_release(up->mem);
  }

  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);

  gkyl_free(up);
}
