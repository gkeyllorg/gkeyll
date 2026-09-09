#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_canonical_pb_fluid_vars.h>
#include <gkyl_dg_calc_canonical_pb_fluid_vars_priv.h>
#include <gkyl_wv_canonical_pb_fluid.h>
#include <gkyl_util.h>

gkyl_dg_calc_canonical_pb_fluid_vars*
gkyl_dg_calc_canonical_pb_fluid_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_wv_eqn *wv_eqn, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_canonical_pb_fluid_vars_cu_dev_new(conf_grid, 
      conf_basis, conf_range, conf_ext_range, wv_eqn);
  } 
#endif     
  gkyl_dg_calc_canonical_pb_fluid_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_canonical_pb_fluid_vars));

  up->conf_grid = *conf_grid;
  up->conf_basis = *conf_basis; 
  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  up->cdim = cdim;
  up->alpha = 0.0;
  up->is_modified = 0; 

  up->adiabatic_coupling_phi_n = 0;
  if (wv_eqn->type == GKYL_EQN_CAN_PB_HASEGAWA_MIMA) {
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_hasegawa_mima_source_kern(conf_basis->b_type, cdim, poly_order);
  }
  else if (wv_eqn->type == GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI) {
    up->alpha = gkyl_wv_can_pb_hasegawa_wakatani_alpha(wv_eqn); 
    up->is_modified = gkyl_wv_can_pb_hasegawa_wakatani_is_modified(wv_eqn); 

    // Temporary array for holding the density and combined potential and density for computing the adiabatic coupling.
    // These are stored separately from the input phi and fluid arrays in case we are solving the
    // modified Hasegawa-Wakatani system and need to subtract the zonal components. 
    up->n = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_ext_range->volume);
    // Component 0 is n, Component 1 is phi. 
    up->adiabatic_coupling_phi_n = gkyl_array_new(GKYL_DOUBLE, 2*conf_basis->num_basis, conf_ext_range->volume);

    // Set up the array averaging to correctly subtact the zonal component of fluctuations
    // Currently assumes that updater has *no decomposition* in y and thus owns the whole y range.
    if (up->is_modified) {
      // Make the one-dimensional x basis, and ranges for constructing the average. 
      struct gkyl_basis basis_x;
      gkyl_cart_modal_serendip(&basis_x, 1, poly_order);
      gkyl_range_init(&up->x_local, 1, &conf_range->lower[0], &conf_range->upper[0]);
      gkyl_range_init(&up->x_local_ext, 1, &conf_ext_range->lower[0], &conf_ext_range->upper[0]);
      // Integration over y only, (x,y) to (x).
      int int_dim_y[] = {0,1,0};
      struct gkyl_array_average_inp inp_int_y = {
        .grid = &up->conf_grid,
        .basis = up->conf_basis,
        .basis_avg = basis_x,
        .local = conf_range,
        .local_avg = &up->x_local,
        .local_avg_ext = &up->x_local_ext,
        .weight = NULL,
        .avg_dim = int_dim_y,
        .use_gpu = use_gpu
      };
      up->int_y = gkyl_array_average_inew(&inp_int_y);
      up->phi_zonal = gkyl_array_new(GKYL_DOUBLE, basis_x.num_basis, up->x_local_ext.volume);
      up->n_zonal = gkyl_array_new(GKYL_DOUBLE, basis_x.num_basis, up->x_local_ext.volume);
      up->subtract_zonal = choose_canonical_pb_fluid_subtract_zonal_kern(conf_basis->b_type, cdim, poly_order);
    }
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_hasegawa_wakatani_source_kern(conf_basis->b_type, cdim, poly_order);
  }
  else {
    // Default source kernel; immediately returns and does not do anything. 
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_default_source_kern(conf_basis->b_type, cdim, poly_order);
  }
  for (int d=0; d<cdim; ++d) {
    up->alpha_surf[d] = choose_canonical_pb_fluid_alpha_surf_kern(conf_basis->b_type, d, cdim, poly_order);
    up->alpha_edge_surf[d] = choose_canonical_pb_fluid_alpha_edge_surf_kern(conf_basis->b_type, d, cdim, poly_order);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_array* phi,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(alpha_surf)) {
    return gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf_cu(up, conf_range, conf_ext_range, phi, 
      alpha_surf, sgn_alpha_surf, const_sgn_alpha);
  }
#endif
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx);
    long loc = gkyl_range_idx(conf_range, idx);
    gkyl_rect_grid_cell_center(&up->conf_grid, idx, xc);

    double* alpha_surf_d = gkyl_array_fetch(alpha_surf, loc);
    double* sgn_alpha_surf_d = gkyl_array_fetch(sgn_alpha_surf, loc);
    int* const_sgn_alpha_d = gkyl_array_fetch(const_sgn_alpha, loc);

    // Fill in the configuration space alpha_surf
    for (int dir = 0; dir<cdim; ++dir) {
      const_sgn_alpha_d[dir] = up->alpha_surf[dir](xc, up->conf_grid.dx, 
        (const double*) gkyl_array_cfetch(phi, loc),
        alpha_surf_d, sgn_alpha_surf_d);

      // If the configuration space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the potential in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_ext = gkyl_range_idx(conf_ext_range, idx_edge);

        double* alpha_surf_ext_d = gkyl_array_fetch(alpha_surf, loc_ext);
        double* sgn_alpha_surf_ext_d = gkyl_array_fetch(sgn_alpha_surf, loc_ext);
        int* const_sgn_alpha_ext_d = gkyl_array_fetch(const_sgn_alpha, loc_ext);
        const_sgn_alpha_ext_d[dir] = up->alpha_edge_surf[dir](xc, up->conf_grid.dx, 
          (const double*) gkyl_array_cfetch(phi, loc),
          alpha_surf_ext_d, sgn_alpha_surf_ext_d);
      }  
    }
  }
}

void gkyl_canonical_pb_fluid_vars_source(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array *phi, const struct gkyl_array *n0, 
  const struct gkyl_array *fluid, struct gkyl_array *rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_canonical_pb_fluid_vars_source_cu(up, conf_range, phi, n0, fluid, rhs);
  }
#endif
  // If alpha is specified, we are solving Hasegawa-Wakatani and need to check 
  // whether we are solving the modified version of Hasegawa-Wakatani which 
  // requires computing the zonal components of phi, n: f_zonal = 1/Ly int f dy
  // and subtracting the zonal components off the adiabatic coupling, f_tilde = f - f_zonal. 
  // Otherwise, we just copy n and phi into a temporary array for use in the updater.
  if (up->alpha > 0.0) {
    gkyl_array_set_offset(up->n, 1.0, fluid, up->conf_basis.num_basis);
    gkyl_array_set_offset(up->adiabatic_coupling_phi_n, 1.0, phi, 0);
    gkyl_array_set_offset(up->adiabatic_coupling_phi_n, 1.0, up->n, up->conf_basis.num_basis);
    if (up->is_modified) {
      // Compute the zonal components of phi and n. 
      gkyl_array_average_advance(up->int_y, phi, up->phi_zonal);
      gkyl_array_average_advance(up->int_y, up->n, up->n_zonal);
      // Iterate over the 2D grid and subtract off the 1D zonal component
      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, conf_range);
      while (gkyl_range_iter_next(&iter)) {
        long loc = gkyl_range_idx(conf_range, iter.idx);  
        long loc_1d = gkyl_range_idx(&up->x_local, iter.idx);  

        const double *phi_zonal_d = gkyl_array_cfetch(up->phi_zonal, loc_1d);
        const double *n_zonal_d = gkyl_array_cfetch(up->n_zonal, loc_1d);

        double *adiabatic_coupling_phi_n_d = gkyl_array_fetch(up->adiabatic_coupling_phi_n, loc); 
        up->subtract_zonal(phi_zonal_d, n_zonal_d, adiabatic_coupling_phi_n_d);   
      }
    }
  }

  // Loop over the grid and compute the source update. 
  // For equation systems such as incompressible Euler, this function returns immediately (no sources).
  // For Hasegawa-Mima, computes {phi, n0} where {., .} is the canonical Poisson bracket.
  // For Hasegawa-Wakatani, computes alpha*(phi - n) + {phi, n0}, where the first
  // term potentially has the zonal components subtracted off if we are solving modified Hasegawa-Wakatani. 
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *phi_d = gkyl_array_cfetch(phi, loc);
    const double *n0_d = gkyl_array_cfetch(n0, loc);

    double* rhs_d = gkyl_array_fetch(rhs, loc);

    up->canonical_pb_fluid_source(up->conf_grid.dx, up->alpha, phi_d, n0_d, 
      up->adiabatic_coupling_phi_n ? (const double*) gkyl_array_cfetch(up->adiabatic_coupling_phi_n, loc) : 0, 
      rhs_d);
  }
}

void gkyl_dg_calc_canonical_pb_fluid_vars_release(gkyl_dg_calc_canonical_pb_fluid_vars *up)
{
  // If alpha was specified, we were solving Hasegawa-Wakatani
  // and need to free specific Hasegawa-Wakatani allocated memory. 
  if (up->alpha > 0.0) {
    gkyl_array_release(up->n); 
    gkyl_array_release(up->adiabatic_coupling_phi_n); 
    if (up->is_modified) {
      gkyl_array_release(up->phi_zonal); 
      gkyl_array_release(up->n_zonal); 
      gkyl_array_average_release(up->int_y); 
    }
  }

  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}
