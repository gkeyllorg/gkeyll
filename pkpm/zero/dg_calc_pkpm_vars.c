#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_pkpm_vars *
gkyl_dg_calc_pkpm_vars_new(
  const struct gkyl_rect_grid *conf_grid, const struct gkyl_basis *cbasis,
  const struct gkyl_vlasov_position_map *pos_map, bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_dg_calc_pkpm_vars_cu_dev_new(conf_grid, cbasis, pos_map);
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
  assert(pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(pos_map);
  up->jacob_pos = pos_map->jacob_pos;

  for (int d = 0; d < cdim; ++d) {
    up->pkpm_in_cell_grad[d] = choose_pkpm_in_cell_grad_kern(b_type, cdim, poly_order, d);
  }
  up->pkpm_closure = choose_pkpm_closure_kern(b_type, cdim, poly_order);
  up->pkpm_prim = choose_pkpm_prim_kern(b_type, cdim, poly_order);
  up->pkpm_int = choose_pkpm_int_kern(b_type, cdim, poly_order);
  up->pkpm_grad_ln_n = choose_pkpm_grad_ln_n_kern(b_type, cdim, poly_order);
  up->pkpm_fluid_positivity = choose_pkpm_fluid_positivity_kern(b_type, cdim, poly_order);
  up->pkpm_normalize = choose_pkpm_normalize_kern(b_type, cdim, poly_order);
  up->pkpm_transverse = choose_pkpm_transverse_kern(b_type, cdim, poly_order);
  up->pkpm_source = choose_pkpm_source_kern(b_type, cdim, poly_order);
  up->pkpm_diag = choose_pkpm_diag_kern(b_type, cdim, poly_order);
  up->num_null_rho = 0;
  up->num_null_rho_dev = 0;
  up->closure_cond_min = 1.0;
  up->closure_cond_min_dev = 0;
  up->num_pos_dev = 0;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

int
gkyl_dg_calc_pkpm_vars_closure(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *wnode, const struct gkyl_array *theta,
  const struct gkyl_array *fin,
  struct gkyl_array *c_0, struct gkyl_array *alpha_0, struct gkyl_array *c_h, struct gkyl_array *a_perp,
  struct gkyl_array *vel_flux_surf, struct gkyl_array *cflrate
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(alpha_0)) {
    return gkyl_dg_calc_pkpm_vars_closure_cu(
      up, conf_range, phase_range, phase_grid, norm_moms, rhs_moms, bvar, wnode, theta, fin, c_0, alpha_0, c_h,
      a_perp, vel_flux_surf, cflrate
    );
  }
#endif
  assert(up->pkpm_closure);
  int cdim = conf_range->ndim;
  double w_lo = phase_grid->lower[cdim], w_up = phase_grid->upper[cdim], dv = phase_grid->dx[cdim];
  int idx[GKYL_MAX_DIM];
  double cfl[2], cond[1];
  int num_null = 0;
  double cond_min = 1.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    gkyl_copy_int_arr(cdim, iter.idx, idx);
    idx[cdim] = phase_range->lower[cdim];
    long plo = gkyl_range_idx(phase_range, idx);
    idx[cdim] = phase_range->upper[cdim];
    long pup = gkyl_range_idx(phase_range, idx);
    idx[cdim] = phase_range->upper[cdim] + 1; // velocity ghost cell: its slot is the upper edge face
    long pgh = gkyl_range_idx(phase_range, idx);
    num_null += up->pkpm_closure(
      w_lo, w_up, dv, gkyl_array_cfetch(norm_moms, loc), gkyl_array_cfetch(rhs_moms, loc),
      gkyl_array_cfetch(bvar, loc), gkyl_array_cfetch(wnode, loc), gkyl_array_cfetch(theta, loc),
      gkyl_array_cfetch(alpha_0, loc),
      gkyl_array_cfetch(c_h, loc), gkyl_array_cfetch(fin, plo), gkyl_array_cfetch(fin, pup),
      gkyl_array_fetch(c_0, loc), gkyl_array_fetch(alpha_0, loc), gkyl_array_fetch(c_h, loc),
      gkyl_array_fetch(a_perp, loc), gkyl_array_fetch(vel_flux_surf, plo), gkyl_array_fetch(vel_flux_surf, pgh),
      cfl, cond
    );
    cond_min = fmin(cond_min, cond[0]);
    // the column rate of the velocity-space terms, the same in every cell of the column
    for (int v = phase_range->lower[cdim]; v <= phase_range->upper[cdim]; ++v) {
      idx[cdim] = v;
      double *c = gkyl_array_fetch(cflrate, gkyl_range_idx(phase_range, idx));
      c[0] += cfl[0];
    }
  }
  up->num_null_rho = num_null;
  up->closure_cond_min = cond_min;
  return num_null;
}

int
gkyl_dg_calc_pkpm_vars_prim(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *u, struct gkyl_array *xi
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(xi)) {
    return gkyl_dg_calc_pkpm_vars_prim_cu(up, conf_range, fluid, u, xi);
  }
#endif
  assert(up->pkpm_prim);
  int num_null = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    num_null += up->pkpm_prim(gkyl_array_cfetch(fluid, loc), gkyl_array_fetch(u, loc), gkyl_array_fetch(xi, loc));
  }
  return num_null;
}

void
gkyl_dg_calc_pkpm_vars_int(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(int_pkpm_vars)) {
    return gkyl_dg_calc_pkpm_vars_int_cu(up, conf_range, fluid, vlasov_pkpm_moms, u, int_pkpm_vars);
  }
#endif
  assert(up->pkpm_int);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    up->pkpm_int(
      gkyl_array_cfetch(fluid, loc), gkyl_array_cfetch(vlasov_pkpm_moms, loc), gkyl_array_cfetch(u, loc),
      gkyl_array_fetch(int_pkpm_vars, loc)
    );
  }
}

// The three constraints of the normalized distributions [M0(F), M1(F), (M2(F) + 2 M0(G))/3]
// from the raw moments [M0, M2, M0_G, M1, M3, M1_G]; they are [1, 0, 1] in every mode of
// every cell when the closures hold.
void
gkyl_dg_calc_pkpm_vars_constraints(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, int num_basis,
  const struct gkyl_array *norm_moms, struct gkyl_array *constraints
)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *m = gkyl_array_cfetch(norm_moms, loc);
    double *c = gkyl_array_fetch(constraints, loc);
    for (int k = 0; k < num_basis; ++k) {
      c[k] = m[k];
      c[num_basis + k] = m[3 * num_basis + k];
      c[2 * num_basis + k] = (m[num_basis + k] + 2.0 * m[2 * num_basis + k]) / 3.0;
    }
  }
}

double
gkyl_dg_calc_pkpm_vars_closure_cond_min(const struct gkyl_dg_calc_pkpm_vars *up)
{
  return up->closure_cond_min;
}

int
gkyl_dg_calc_pkpm_vars_grad_ln_n(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rho, const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(grad_ln_n)) {
    return gkyl_dg_calc_pkpm_vars_grad_ln_n_cu(up, conf_range, rho, grad_rho, grad_ln_n);
  }
#endif
  assert(up->pkpm_grad_ln_n);
  int num_null = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    num_null += up->pkpm_grad_ln_n(
      gkyl_array_cfetch(rho, loc), gkyl_array_cfetch(grad_rho, loc), gkyl_array_fetch(grad_ln_n, loc)
    );
  }
  return num_null;
}

void
gkyl_dg_calc_pkpm_vars_positivity(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, double eps,
  struct gkyl_array *fluid, struct gkyl_array *count, int *num_rho, int *num_e, int *num_xi,
  int *num_avg
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid)) {
    return gkyl_dg_calc_pkpm_vars_positivity_cu(
      up, conf_range, eps, fluid, count, num_rho, num_e, num_xi, num_avg
    );
  }
#endif
  assert(up->pkpm_fluid_positivity);
  int nrho = 0, ne = 0, nxi = 0, navg = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    int flag = up->pkpm_fluid_positivity(eps, gkyl_array_fetch(fluid, loc));
    double *cnt = gkyl_array_fetch(count, loc);
    if (flag < 0) {
      cnt[3] += 1.0;
      navg += 1;
    }
    if (flag > 0 && (flag & 1)) {
      cnt[0] += 1.0;
      nrho += 1;
    }
    if (flag > 0 && (flag & 2)) {
      cnt[1] += 1.0;
      ne += 1;
    }
    if (flag > 0 && (flag & 4)) {
      cnt[2] += 1.0;
      nxi += 1;
    }
  }
  *num_rho = nrho;
  *num_e = ne;
  *num_xi = nxi;
  *num_avg = navg;
}

void
gkyl_dg_calc_pkpm_vars_transverse(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, struct gkyl_array *f
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(f)) {
    gkyl_dg_calc_pkpm_vars_transverse_cu(up, conf_range, phase_range, bvar, f);
    return;
  }
#endif
  assert(up->pkpm_transverse);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);
    long pidx = gkyl_range_idx(phase_range, iter.idx);
    up->pkpm_transverse(gkyl_array_cfetch(bvar, cidx), gkyl_array_fetch(f, pidx));
  }
}

int
gkyl_dg_calc_pkpm_vars_normalize_dist(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *n, const struct gkyl_array *f,
  struct gkyl_array *fbar
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fbar)) {
    return gkyl_dg_calc_pkpm_vars_normalize_dist_cu(up, conf_range, phase_range, n, f, fbar);
  }
#endif
  assert(up->pkpm_normalize);
  int num_null = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);
    long pidx = gkyl_range_idx(phase_range, iter.idx);
    num_null += up->pkpm_normalize(
      gkyl_array_cfetch(n, cidx), gkyl_array_cfetch(f, pidx), gkyl_array_fetch(fbar, pidx)
    );
  }
  return num_null;
}

void
gkyl_dg_calc_pkpm_vars_diag(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *alpha_0,
  const struct gkyl_array *xi, const struct gkyl_array *c_h, struct gkyl_array *fluid_io,
  struct gkyl_array *pkpm_vars_io
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid_io)) {
    return gkyl_dg_calc_pkpm_vars_diag_cu(
      up, conf_range, vlasov_pkpm_moms, euler_pkpm, bvar, grad_b, u, grad_u, theta, alpha_0, xi, c_h,
      fluid_io, pkpm_vars_io
    );
  }
#endif
  assert(up->pkpm_diag);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    up->pkpm_diag(
      up->conf_grid.dx, gkyl_array_cfetch(vlasov_pkpm_moms, loc), gkyl_array_cfetch(euler_pkpm, loc),
      gkyl_array_cfetch(bvar, loc), gkyl_array_cfetch(grad_b, loc), gkyl_array_cfetch(u, loc),
      gkyl_array_cfetch(grad_u, loc), gkyl_array_cfetch(theta, loc), gkyl_array_cfetch(alpha_0, loc),
      gkyl_array_cfetch(xi, loc), gkyl_array_cfetch(c_h, loc), gkyl_array_fetch(fluid_io, loc),
      gkyl_array_fetch(pkpm_vars_io, loc)
    );
  }
}

void
gkyl_dg_calc_pkpm_vars_in_cell_grad(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *blend, const struct gkyl_array *q, struct gkyl_array *grad
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(grad)) {
    return gkyl_dg_calc_pkpm_vars_in_cell_grad_cu(up, conf_range, blend, q, grad);
  }
#endif
  int cdim = up->cdim;
  int idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *q_c = gkyl_array_cfetch(q, loc);
    const double *blend_c = gkyl_array_cfetch(blend, loc);
    double *gd = gkyl_array_fetch(grad, loc);
    for (int d = 0; d < cdim; ++d) {
      assert(up->pkpm_in_cell_grad[d]);
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      gkyl_copy_int_arr(cdim, iter.idx, idx_r);
      idx_l[d] = idx_l[d] - 1;
      idx_r[d] = idx_r[d] + 1;
      long loc_r = gkyl_range_idx(conf_range, idx_r);
      const double *blend_r = gkyl_array_cfetch(blend, loc_r);
      up->pkpm_in_cell_grad[d](
        up->conf_grid.dx, gkyl_array_cfetch(up->jacob_pos, loc), blend_c[d], blend_r[d],
        gkyl_array_cfetch(q, gkyl_range_idx(conf_range, idx_l)), q_c, gkyl_array_cfetch(q, loc_r),
        gd
      );
    }
  }
}

void
gkyl_dg_calc_pkpm_vars_source(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *qmem, const struct gkyl_array *fluid, struct gkyl_array *rhs
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_dg_calc_pkpm_vars_source_cu(up, conf_range, qmem, fluid, rhs);
  }
#endif
  assert(up->pkpm_source);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    up->pkpm_source(gkyl_array_cfetch(qmem, loc), gkyl_array_cfetch(fluid, loc), gkyl_array_fetch(rhs, loc));
  }
}

void
gkyl_dg_calc_pkpm_vars_release(gkyl_dg_calc_pkpm_vars *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);

  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->num_null_rho_dev);
    gkyl_cu_free(up->closure_cond_min_dev);
    gkyl_cu_free(up->num_pos_dev);
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}
