#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_harm_vol_priv.h>

void
gkyl_dg_pkpm_harm_vol_setup(struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
  int cdim = inp->conf_basis->ndim;
  assert(inp->conf_basis->poly_order == 1 && inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->qbym = inp->qbym;
  up->use_gpu = inp->use_gpu;
  assert(inp->pos_map);
  up->pos_map = 0;
  up->jacob_pos = inp->pos_map->jacob_pos;
  up->nop = PKPM_HO_NOP(cdim);
  // the largest |value| over the cell of the degree <= 2 modes k = k0 + 3 k1 + 9 k2, products of
  // the orthonormal Legendre polynomials, whose largest values are at the cell edges:
  // |Phat_0| = 1/sqrt 2, |Phat_1(1)| = sqrt(3/2), |Phat_2(1)| = sqrt(5/2)
  const double m1[3] = {0.7071067811865475, 1.224744871391589, 1.5811388300841898};
  up->n2 = 1;
  for (int d = 0; d < cdim; ++d) {
    up->n2 *= 3;
  }
  for (int k = 0; k < up->n2; ++k) {
    up->bmax[k] = 1.0;
    for (int d = 0, kk = k; d < cdim; ++d, kk /= 3) {
      up->bmax[k] *= m1[kk % 3];
    }
  }
  up->wmax = fmax(fabs(inp->phase_grid->lower[cdim]), fabs(inp->phase_grid->upper[cdim]));
  // the source pairs of the rate estimate
  up->nsrc = 0;
  for (int e = 0; e < PKPM_HO_NENT; ++e) {
    if (pkpm_ho_ent[e][2] != 0) {
      continue;
    }
    int p = 0;
    while (p < up->nsrc && !(up->src_o[p] == pkpm_ho_ent[e][0] && up->src_i[p] == pkpm_ho_ent[e][1])) {
      ++p;
    }
    if (p == up->nsrc) {
      assert(up->nsrc < PKPM_HM_MAXSRC);
      up->src_o[p] = pkpm_ho_ent[e][0];
      up->src_i[p] = pkpm_ho_ent[e][1];
      up->src_ent[p][0] = up->src_ent[p][1] = up->src_ent[p][2] = -1;
      up->nsrc += 1;
    }
    up->src_ent[p][pkpm_ho_ent[e][3]] = e;
  }
}

struct gkyl_dg_pkpm_harm_vol *
gkyl_dg_pkpm_harm_vol_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_harm_vol_cu_dev_inew(inp);
  }
#endif
  struct gkyl_dg_pkpm_harm_vol *up = gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_harm_vol_setup(up, inp);
  up->build_kern = ser_pkpm_harm_build_kernels[up->cdim - 1];
  up->mat_kern = ser_pkpm_harm_mat_kernels[up->cdim - 1];
  up->vol_kern = ser_pkpm_harm_vol_kernels[up->cdim - 1];
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

int
gkyl_dg_pkpm_harm_vol_num_coef(const struct gkyl_dg_pkpm_harm_vol *up)
{
  return up->nop * up->n2;
}

int
gkyl_dg_pkpm_harm_vol_num_wnode(const struct gkyl_dg_pkpm_harm_vol *up)
{
  int nn = 1;
  for (int d = 0; d < up->cdim; ++d) {
    nn *= 5;
  }
  return PKPM_HS_W_END * nn;
}

void
gkyl_dg_pkpm_harm_vol_coef(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, struct gkyl_array *coef,
  struct gkyl_array *wnode
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(coef)) {
    gkyl_dg_pkpm_harm_vol_coef_cu(
      up, conf_range, bvar, grad_b, u, grad_u, xi, grad_xi, grad_ln_n, theta, dtB_over_magB, em, coef,
      wnode
    );
    return;
  }
#endif
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, conf_range);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(conf_range, citer.idx);
    up->build_kern(up->qbym, gkyl_array_cfetch(bvar, cidx), gkyl_array_cfetch(grad_b, cidx),
      gkyl_array_cfetch(u, cidx), gkyl_array_cfetch(grad_u, cidx), gkyl_array_cfetch(xi, cidx),
      gkyl_array_cfetch(grad_xi, cidx), gkyl_array_cfetch(grad_ln_n, cidx), gkyl_array_cfetch(theta, cidx),
      gkyl_array_cfetch(dtB_over_magB, cidx), gkyl_array_cfetch(em, cidx), gkyl_array_fetch(coef, cidx),
      gkyl_array_fetch(wnode, cidx));
  }
}

void
gkyl_dg_pkpm_harm_vol_advance(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    gkyl_dg_pkpm_harm_vol_advance_cu(up, conf_range, phase_range, coef, fin, cflrate, rhs);
    return;
  }
#endif
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM];
  double mat[PKPM_HM_MAT];
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, conf_range);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(conf_range, citer.idx);
    // the matrices of the stored operator of the configuration cell and the source rate of its column
    const double *ho = gkyl_array_cfetch(coef, cidx);
    up->mat_kern(up->phase_grid.dx, gkyl_array_cfetch(up->jacob_pos, cidx), ho, mat);
    double rate = pkpm_hm_src_bound(up, ho);
    gkyl_copy_int_arr(cdim, citer.idx, idx);
    for (int v = phase_range->lower[cdim]; v <= phase_range->upper[cdim]; ++v) {
      idx[cdim] = v;
      long pidx = gkyl_range_idx(phase_range, idx);
      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
      up->vol_kern(xc, up->phase_grid.dx, mat, gkyl_array_cfetch(fin, pidx), gkyl_array_fetch(rhs, pidx));
      double *cfl = gkyl_array_fetch(cflrate, pidx);
      cfl[0] += rate;
    }
  }
}

void
gkyl_dg_pkpm_harm_vol_release(struct gkyl_dg_pkpm_harm_vol *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_harm_vol *
gkyl_dg_pkpm_harm_vol_cu_dev_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_pkpm_harm_vol_coef_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, struct gkyl_array *coef,
  struct gkyl_array *wnode
)
{
  assert(false);
}

void
gkyl_dg_pkpm_harm_vol_advance_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
  assert(false);
}

#endif
