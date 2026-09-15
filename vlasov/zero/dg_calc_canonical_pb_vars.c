#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_canonical_pb_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_canonical_pb_vars*
gkyl_dg_calc_canonical_pb_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu)
{
  // The Serendipity canonical-PB moment tables have no 3x3v row (only the
  // tensor p=1 hybrid is generated in 3x3v); fail loudly instead of indexing
  // past the table, on both host and device.
  assert(!(phase_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY
    && conf_basis->ndim == 3 && phase_basis->ndim == 6));

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_canonical_pb_vars_cu_dev_new(phase_grid, 
      conf_basis, phase_basis);
  } 
#endif     
  gkyl_dg_calc_canonical_pb_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_canonical_pb_vars));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->canonical_pb_pressure = choose_canonical_pb_pressure_kern(phase_basis->b_type, cv_index[cdim].vdim[vdim], cdim, poly_order);
  up->canonical_pb_covariant_u_i = choose_canonical_pb_m1i_contra_to_cov_kern(phase_basis->b_type, cv_index[cdim].vdim[vdim], cdim, poly_order);
  assert(up->canonical_pb_pressure);
  assert(up->canonical_pb_covariant_u_i);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}


void gkyl_canonical_pb_contra_to_covariant_m1i(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij, const struct gkyl_array *V_drift, const struct gkyl_array *M1i,
  struct gkyl_array *V_drift_cov, struct gkyl_array *M1i_cov)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(V_drift_cov)) {
    return gkyl_canonical_pb_contra_to_covariant_m1i_cu(up, conf_range, h_ij, V_drift, M1i, V_drift_cov, M1i_cov);
  }
#endif
  int cdim = up->cdim;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *h_ij_d = gkyl_array_cfetch(h_ij, loc);
    const double *v_i_d = gkyl_array_cfetch(V_drift, loc);
    const double *nv_i_d = gkyl_array_cfetch(M1i, loc);
   
    double *v_i_cov_d = gkyl_array_fetch(V_drift_cov, loc);
    double *nv_i_cov_d = gkyl_array_fetch(M1i_cov, loc);

    up->canonical_pb_covariant_u_i(h_ij_d, v_i_d, nv_i_d, v_i_cov_d, nv_i_cov_d);
  }
}

void gkyl_canonical_pb_pressure(struct gkyl_dg_calc_canonical_pb_vars *up, const struct gkyl_range *conf_range,
 const struct gkyl_array *h_ij_inv, 
 const struct gkyl_array *MEnergy, const struct gkyl_array *V_drift, const struct gkyl_array *M1i,
 struct gkyl_array *pressure)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pressure)) {
    return gkyl_canonical_pb_pressure_cu(up, conf_range, h_ij_inv, MEnergy, V_drift, M1i, pressure);
  }
#endif
  int cdim = up->cdim;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *h_ij_inv_d = gkyl_array_cfetch(h_ij_inv, loc);
    const double *MEnergy_d = gkyl_array_cfetch(MEnergy, loc);
    const double *v_j_d = gkyl_array_cfetch(V_drift, loc);
    const double *nv_i_d = gkyl_array_cfetch(M1i, loc);

    double* d_Jv_P_d = gkyl_array_fetch(pressure, loc);

    up->canonical_pb_pressure(h_ij_inv_d, MEnergy_d, v_j_d, nv_i_d, d_Jv_P_d);
  }
}

void gkyl_dg_calc_canonical_pb_vars_release(gkyl_dg_calc_canonical_pb_vars *up)
{
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
