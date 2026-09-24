#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_vars_priv.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_util.h>

gkyl_dg_calc_pkpm_vars*
gkyl_dg_calc_pkpm_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_range *mem_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, double limiter_fac, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_pkpm_vars_cu_dev_new(conf_grid, cbasis, mem_range, wv_eqn, geom, limiter_fac);
  } 
#endif     
  gkyl_dg_calc_pkpm_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_vars));

  up->conf_grid = *conf_grid;
  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->Ncomp = 9;
  up->mem_range = *mem_range;

  up->wv_eqn = gkyl_wv_eqn_acquire(wv_eqn);
  up->geom = gkyl_wave_geom_acquire(geom);
  // Limiter factor for relationship between slopes and cell average differences
  // By default, this factor is 1/sqrt(3) because cell_avg(f) = f0/sqrt(2^cdim)
  // and a cell slope estimate from two adjacent cells is (for the x variation): 
  // integral(psi_1 [cell_avg(f_{i+1}) - cell_avg(f_{i})]*x) = sqrt(2^cdim)/sqrt(3)*[cell_avg(f_{i+1}) - cell_avg(f_{i})]
  // where psi_1 is the x cell slope basis in our orthonormal expansion psi_1 = sqrt(3)/sqrt(2^cdim)*x
  // This factor can be made smaller (larger) to increase (decrease) the diffusion from the slope limiter
  if (limiter_fac == 0.0) {
    up->limiter_fac = 0.5773502691896258;
  }
  else {
    up->limiter_fac = limiter_fac;
  }

  // Tolerance in mass density and average normal velocity at the interface
  // for switching to Lax fluxes in computing penalization of the momentum solve
  up->tol = 1.0e-12;
  up->force_lax = false;

  up->pkpm_set = choose_pkpm_set_kern(b_type, cdim, poly_order);
  up->pkpm_copy = choose_pkpm_copy_kern(b_type, cdim, poly_order);
  up->pkpm_u_set = choose_pkpm_u_set_kern(b_type, cdim, poly_order);
  up->pkpm_u_copy = choose_pkpm_u_copy_kern(b_type, cdim, poly_order);  
  up->pkpm_pressure = choose_pkpm_pressure_kern(b_type, cdim, poly_order);
  up->pkpm_p_force = choose_pkpm_p_force_kern(b_type, cdim, poly_order);
  up->pkpm_source = choose_pkpm_source_kern(b_type, cdim, poly_order);
  up->pkpm_int = choose_pkpm_int_kern(b_type, cdim, poly_order);
  up->pkpm_io = choose_pkpm_io_kern(b_type, cdim, poly_order);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) {
    up->pkpm_accel[d] = choose_pkpm_accel_kern(d, b_type, cdim, poly_order);
    up->pkpm_penalization[d] = choose_pkpm_penalization_kern(d, b_type, cdim, poly_order);
    up->pkpm_limiter[d] = choose_pkpm_limiter_kern(d, b_type, cdim, poly_order);
  }

  // There are Ncomp*range->volume linear systems to be solved 
  // 6 components: ux, uy, uz, div(p_par b)/rho, p_perp/rho, rho/p_perp
  up->As = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  // Linear system for just solving for ux, uy, uz
  up->As_u = gkyl_nmat_new(3*mem_range->volume, nc, nc);
  up->xs_u = gkyl_nmat_new(3*mem_range->volume, nc, 1);
  up->mem_u = gkyl_nmat_linsolve_lu_new(up->As_u->num, up->As_u->nr);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_pkpm_vars_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* p_ij, const struct gkyl_array* pkpm_div_ppar, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* prim, struct gkyl_array* prim_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(prim)) {
    return gkyl_dg_calc_pkpm_vars_advance_cu(up, 
      vlasov_pkpm_moms, euler_pkpm, 
      p_ij, pkpm_div_ppar, 
      cell_avg_prim, prim, prim_surf);
  }
#endif

  // First loop over mem_range for solving linear systems to compute primitive moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);
    const double *pkpm_div_ppar_d = gkyl_array_cfetch(pkpm_div_ppar, loc);

    int* cell_avg_prim_d = gkyl_array_fetch(cell_avg_prim, loc);

    cell_avg_prim_d[0] = up->pkpm_set(count, up->As, up->xs, 
      vlasov_pkpm_moms_d, euler_pkpm_d, p_ij_d, pkpm_div_ppar_d);

    count += up->Ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    double* prim_d = gkyl_array_fetch(prim, loc);
    double* prim_surf_d = gkyl_array_fetch(prim_surf, loc);

    up->pkpm_copy(count, up->xs, prim_d, prim_surf_d);

    count += up->Ncomp;
  }
}

void gkyl_dg_calc_pkpm_vars_u(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* pkpm_u)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_u)) {
    return gkyl_dg_calc_pkpm_vars_u_cu(up, 
      vlasov_pkpm_moms, euler_pkpm, 
      cell_avg_prim, pkpm_u);
  }
#endif

  // First loop over mem_range for solving linear systems to compute primitive moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);

    int* cell_avg_prim_d = gkyl_array_fetch(cell_avg_prim, loc);

    cell_avg_prim_d[0] = up->pkpm_u_set(count, up->As_u, up->xs_u, 
      vlasov_pkpm_moms_d, euler_pkpm_d);

    count += 3;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem_u, up->As_u, up->xs_u);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    double* pkpm_u_d = gkyl_array_fetch(pkpm_u, loc);

    up->pkpm_u_copy(count, up->xs_u, pkpm_u_d);

    count += 3;
  }
}

void gkyl_dg_calc_pkpm_vars_pressure(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(p_ij)) {
    return gkyl_dg_calc_pkpm_vars_pressure_cu(up, conf_range, 
      bvar, vlasov_pkpm_moms, p_ij);
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *bvar_d = gkyl_array_cfetch(bvar, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);

    double* p_ij_d = gkyl_array_fetch(p_ij, loc);

    up->pkpm_pressure(bvar_d, vlasov_pkpm_moms_d, p_ij_d);
  }
}

void gkyl_dg_calc_pkpm_vars_accel(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* prim_surf, const struct gkyl_array* prim, 
  const struct gkyl_array* bvar, const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_accel)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_accel)) {
    return gkyl_dg_calc_pkpm_vars_accel_cu(up, conf_range, 
      prim_surf, prim, bvar, div_b, nu, pkpm_accel);
  }
#endif

  // Loop over configuration space range to compute gradients and pkpm acceleration variables
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(conf_range, idxc);

    const double *prim_surf_c = gkyl_array_cfetch(prim_surf, linc);
  
    const double *prim_d = gkyl_array_cfetch(prim, linc);
    const double *bvar_d = gkyl_array_cfetch(bvar, linc);
    const double *div_b_d = gkyl_array_cfetch(div_b, linc);
    const double *nu_d = gkyl_array_cfetch(nu, linc);

    double *pkpm_accel_d = gkyl_array_fetch(pkpm_accel, linc);

    // Compute T_perp/m div(b) and p_force
    up->pkpm_p_force(prim_d, div_b_d, pkpm_accel_d);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      const double *prim_surf_l = gkyl_array_cfetch(prim_surf, linl);
      const double *prim_surf_r = gkyl_array_cfetch(prim_surf, linr);

      up->pkpm_accel[dir](up->conf_grid.dx, 
        prim_surf_l, prim_surf_c, prim_surf_r, 
        prim_d, bvar_d, nu_d,
        pkpm_accel_d);
    }
  }
}

void gkyl_dg_calc_pkpm_vars_penalization(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_penalization)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_penalization)) {
    return gkyl_dg_calc_pkpm_vars_penalization_cu(up, 
      conf_range, conf_range_ext, 
      vlasov_pkpm_moms, p_ij, prim, euler_pkpm, 
      pkpm_lax, pkpm_penalization);
  }
#endif

  // Loop over configuration space range to compute gradients and pkpm acceleration variables
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(conf_range, idxc);
    const struct gkyl_wave_cell_geom *geom = gkyl_wave_geom_get(up->geom, idxc);
  
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, linc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, linc);
    const double *prim_d = gkyl_array_cfetch(prim, linc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, linc);

    double *pkpm_lax_d = gkyl_array_fetch(pkpm_lax, linc);
    double *pkpm_penalization_d = gkyl_array_fetch(pkpm_penalization, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);

      // Each cell owns their *lower* edge surface evaluation so we need both 
      // the cell we are in (linc) and our lower neighbor (linl) to compute 
      // the penalization surface expansions at the lower edge surface. 
      idxl[dir] = idxl[dir]-1; 
      long linl = gkyl_range_idx(conf_range, idxl); 

      const double *vlasov_pkpm_moms_l = gkyl_array_cfetch(vlasov_pkpm_moms, linl);
      const double *p_ij_l = gkyl_array_cfetch(p_ij, linl);
      const double *prim_l = gkyl_array_cfetch(prim, linl);
      const double *euler_pkpm_l = gkyl_array_cfetch(euler_pkpm, linl);

      up->pkpm_penalization[dir](up->tol, up->force_lax, up->wv_eqn, geom, 
        vlasov_pkpm_moms_l, vlasov_pkpm_moms_d, p_ij_l, p_ij_d, 
        prim_l, prim_d, euler_pkpm_l, euler_pkpm_d, 
        pkpm_lax_d, pkpm_penalization_d);

      // If the configuration-space index is at the local configuration space upper value, 
      // we are at the configuration space upper edge and we also need to evaluate the 
      // penalization terms at the upper edge interface. We index into the ghost cells (linr)
      // and following the convention of each cell owning their lower surface expansion,
      // the upper edge surface expansions are store in the ghost cells of the array. 
      if (idxc[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, idxc, idxr);
        idxr[dir] = idxr[dir]+1; 
        long linr = gkyl_range_idx(conf_range_ext, idxr);

        const struct gkyl_wave_cell_geom *geom_r = gkyl_wave_geom_get(up->geom, idxr);

        const double *vlasov_pkpm_moms_r = gkyl_array_cfetch(vlasov_pkpm_moms, linr);
        const double *p_ij_r = gkyl_array_cfetch(p_ij, linr);
        const double *prim_r = gkyl_array_cfetch(prim, linr);
        const double *euler_pkpm_r = gkyl_array_cfetch(euler_pkpm, linr);

        double *pkpm_lax_r = gkyl_array_fetch(pkpm_lax, linr);
        double *pkpm_penalization_r = gkyl_array_fetch(pkpm_penalization, linr);

        up->pkpm_penalization[dir](up->tol, up->force_lax, up->wv_eqn, geom_r, 
          vlasov_pkpm_moms_d, vlasov_pkpm_moms_r, p_ij_d, p_ij_r, 
          prim_d, prim_r, euler_pkpm_d, euler_pkpm_r, 
          pkpm_lax_r, pkpm_penalization_r);
      }
    }
  }
}

void gkyl_dg_calc_pkpm_integrated_vars(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* prim, 
  struct gkyl_array* pkpm_int_vars)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_int_vars)) {
    return gkyl_dg_calc_pkpm_integrated_vars_cu(up, conf_range, 
      vlasov_pkpm_moms, euler_pkpm, prim, pkpm_int_vars);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *prim_d = gkyl_array_cfetch(prim, loc);
    
    double *pkpm_int_vars_d = gkyl_array_fetch(pkpm_int_vars, loc);
    up->pkpm_int(vlasov_pkpm_moms_d, euler_pkpm_d, prim_d, pkpm_int_vars_d);
  }
}

void gkyl_dg_calc_pkpm_vars_source(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* qmem, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_dg_calc_pkpm_vars_source_cu(up, conf_range, 
      qmem, vlasov_pkpm_moms, euler_pkpm, rhs);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *qmem_d = gkyl_array_cfetch(qmem, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);

    double *rhs_d = gkyl_array_fetch(rhs, loc);
    up->pkpm_source(qmem_d, vlasov_pkpm_moms_d, euler_pkpm_d, rhs_d);
  }
}

void gkyl_dg_calc_pkpm_vars_io(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_vars_io)) {
    return gkyl_dg_calc_pkpm_vars_io_cu(up, conf_range, 
      vlasov_pkpm_moms, euler_pkpm, p_ij, prim, pkpm_accel, 
      fluid_io, pkpm_vars_io);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);
    const double *prim_d = gkyl_array_cfetch(prim, loc);
    const double *pkpm_accel_d = gkyl_array_cfetch(pkpm_accel, loc);
    
    double *fluid_io_d = gkyl_array_fetch(fluid_io, loc);
    double *pkpm_vars_io_d = gkyl_array_fetch(pkpm_vars_io, loc);
    up->pkpm_io(vlasov_pkpm_moms_d, euler_pkpm_d, p_ij_d, prim_d, pkpm_accel_d, 
      fluid_io_d, pkpm_vars_io_d);
  }
}

void gkyl_dg_calc_pkpm_vars_limiter(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* prim, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  struct gkyl_array* fluid)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid)) {
    return gkyl_dg_calc_pkpm_vars_limiter_cu(up, conf_range, 
      prim, vlasov_pkpm_moms, p_ij, fluid);
  }
#endif
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(conf_range, idxc);
    const struct gkyl_wave_cell_geom *geom = gkyl_wave_geom_get(up->geom, idxc);

    const double *prim_c = gkyl_array_cfetch(prim, linc);
    const double *vlasov_pkpm_moms_c = gkyl_array_cfetch(vlasov_pkpm_moms, linc);
    const double *p_ij_c = gkyl_array_cfetch(p_ij, linc);

    double *fluid_c = gkyl_array_fetch(fluid, linc);
    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      const double *vlasov_pkpm_moms_l = gkyl_array_cfetch(vlasov_pkpm_moms, linl);
      const double *vlasov_pkpm_moms_r = gkyl_array_cfetch(vlasov_pkpm_moms, linr);
      const double *p_ij_l = gkyl_array_cfetch(p_ij, linl);
      const double *p_ij_r = gkyl_array_cfetch(p_ij, linr);

      double *fluid_l = gkyl_array_fetch(fluid, linl);
      double *fluid_r = gkyl_array_fetch(fluid, linr);

      up->pkpm_limiter[dir](up->limiter_fac, up->wv_eqn, geom, prim_c, 
        vlasov_pkpm_moms_l, vlasov_pkpm_moms_c, vlasov_pkpm_moms_r, 
        p_ij_l, p_ij_c, p_ij_r, 
        fluid_l, fluid_c, fluid_r); 
    }
  }
}

void gkyl_dg_calc_pkpm_vars_release(gkyl_dg_calc_pkpm_vars *up)
{
  gkyl_wv_eqn_release(up->wv_eqn);
  gkyl_wave_geom_release(up->geom);

  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);

  gkyl_nmat_release(up->As_u);
  gkyl_nmat_release(up->xs_u);
  gkyl_nmat_linsolve_lu_release(up->mem_u);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
