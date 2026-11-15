#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_bvar_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_pkpm_bvar *
gkyl_dg_calc_pkpm_bvar_inew(const struct gkyl_dg_calc_pkpm_bvar_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_calc_pkpm_bvar_cu_dev_inew(inp);
  }
#endif
  gkyl_dg_calc_pkpm_bvar *up = gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_bvar));

  up->cdim = inp->conf_basis->ndim;
  for (int d = 0; d < inp->conf_basis->ndim; ++d) {
    up->dx[d] = inp->conf_grid->dx[d];
  }
  assert(inp->pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  up->pkpm_bvar = choose_pkpm_bvar_kern(inp->conf_basis->b_type, inp->conf_basis->ndim, inp->conf_basis->poly_order);
  up->pkpm_dtB_over_magB = choose_pkpm_dtB_over_magB_kern(inp->conf_basis->b_type, inp->conf_basis->ndim, inp->conf_basis->poly_order);
  for (int d = 0; d < inp->conf_basis->ndim; ++d) {
    up->pkpm_grad[d] = choose_pkpm_bvar_grad_kern(inp->conf_basis->b_type, inp->conf_basis->ndim, inp->conf_basis->poly_order, d);
  }
  // every kernel exists for the basis (p=1 serendipity)
  assert(up->pkpm_bvar && up->pkpm_dtB_over_magB);
  for (int d = 0; d < inp->conf_basis->ndim; ++d) {
    assert(up->pkpm_grad[d]);
  }
  up->num_null_b = 0;
  up->num_null_b_dev = 0;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_dg_calc_pkpm_bvar_advance(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *bvar
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(bvar)) {
    return gkyl_dg_calc_pkpm_bvar_advance_cu(up, conf_range, em, bvar);
  }
#endif
  int num_null_b = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    num_null_b += up->pkpm_bvar(gkyl_array_cfetch(em, loc), gkyl_array_fetch(bvar, loc));
  }
  up->num_null_b = num_null_b;
}

int
gkyl_dg_calc_pkpm_bvar_dtB_over_magB(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, const struct gkyl_array *dt_em, struct gkyl_array *dtB_over_magB
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(dtB_over_magB)) {
    return gkyl_dg_calc_pkpm_bvar_dtB_over_magB_cu(up, conf_range, em, dt_em, dtB_over_magB);
  }
#endif
  int num_null = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    num_null += up->pkpm_dtB_over_magB(
      gkyl_array_cfetch(em, loc), gkyl_array_cfetch(dt_em, loc), gkyl_array_fetch(dtB_over_magB, loc)
    );
  }
  return num_null;
}

int
gkyl_dg_calc_pkpm_bvar_num_null_b(const struct gkyl_dg_calc_pkpm_bvar *up)
{
  return up->num_null_b;
}

void
gkyl_dg_calc_pkpm_bvar_grad(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, struct gkyl_array *grad_b
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(grad_b)) {
    return gkyl_dg_calc_pkpm_bvar_grad_cu(up, conf_range, bvar, grad_b);
  }
#endif
  int cdim = up->cdim;
  int idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    const double *b_c = gkyl_array_cfetch(bvar, loc);
    double *gb = gkyl_array_fetch(grad_b, loc);
    for (int d = 0; d < cdim; ++d) {
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      gkyl_copy_int_arr(cdim, iter.idx, idx_r);
      idx_l[d] = idx_l[d] - 1;
      idx_r[d] = idx_r[d] + 1;
      up->pkpm_grad[d](
        up->dx, gkyl_array_cfetch(up->jacob_pos, loc), gkyl_array_cfetch(bvar, gkyl_range_idx(conf_range, idx_l)), b_c,
        gkyl_array_cfetch(bvar, gkyl_range_idx(conf_range, idx_r)), gb
      );
    }
  }
}

void
gkyl_dg_calc_pkpm_bvar_release(gkyl_dg_calc_pkpm_bvar *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->num_null_b_dev);
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

gkyl_dg_calc_pkpm_bvar *
gkyl_dg_calc_pkpm_bvar_cu_dev_inew(const struct gkyl_dg_calc_pkpm_bvar_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_calc_pkpm_bvar_advance_cu(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *bvar
)
{
  assert(false);
}

#endif

void
gkyl_dg_calc_pkpm_bvar_grad_cu(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, struct gkyl_array *grad_b
)
{
  assert(false);
}
