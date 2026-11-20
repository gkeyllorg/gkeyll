#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_mom_diffusion.h>
#include <gkyl_dg_pkpm_mom_diffusion_priv.h>

void
gkyl_dg_pkpm_mom_diffusion_setup(
  struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_dg_pkpm_mom_diffusion_inp *inp
)
{
  assert(inp->conf_basis->poly_order == 1 && inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
  assert(inp->order == 2 || inp->order == 4);
  assert(inp->pos_map);
  up->conf_grid = *inp->conf_grid;
  up->cdim = inp->conf_basis->ndim;
  up->order = inp->order;
  up->use_gpu = inp->use_gpu;
  up->num_dir = 0;
  for (int d = 0; d < up->cdim; ++d) {
    if (inp->D[d] != 0.0) {
      up->dirs[up->num_dir] = d;
      up->D[up->num_dir] = inp->D[d];
      up->num_dir += 1;
    }
  }
  up->pos_map = 0;
  up->jacob_pos = inp->pos_map->jacob_pos;
}

struct gkyl_dg_pkpm_mom_diffusion *
gkyl_dg_pkpm_mom_diffusion_inew(const struct gkyl_dg_pkpm_mom_diffusion_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_mom_diffusion_cu_dev_inew(inp);
  }
#endif
  struct gkyl_dg_pkpm_mom_diffusion *up = gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_mom_diffusion_setup(up, inp);
  for (int k = 0; k < up->num_dir; ++k) {
    up->kern[k] = pkpm_md_choose_kern(up->order, up->cdim, up->dirs[k]);
  }
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

void
gkyl_dg_pkpm_mom_diffusion_advance(
  const struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    gkyl_dg_pkpm_mom_diffusion_advance_cu(up, conf_range, fluid, cflrate, rhs);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);
    double *cfl = gkyl_array_fetch(cflrate, loc);
    cfl[0] += pkpm_md_cell(up, conf_range, iter.idx, fluid, gkyl_array_fetch(rhs, loc));
  }
}

void
gkyl_dg_pkpm_mom_diffusion_release(struct gkyl_dg_pkpm_mom_diffusion *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_mom_diffusion *
gkyl_dg_pkpm_mom_diffusion_cu_dev_inew(const struct gkyl_dg_pkpm_mom_diffusion_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_pkpm_mom_diffusion_advance_cu(
  const struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
  assert(false);
}

#endif
