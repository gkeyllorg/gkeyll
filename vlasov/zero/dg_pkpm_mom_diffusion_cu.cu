/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_mom_diffusion.h>
#include <gkyl_dg_pkpm_mom_diffusion_priv.h>
}

// One thread per configuration cell.
__global__ static void
gkyl_dg_pkpm_mom_diffusion_advance_cu_kernel(
  const struct gkyl_dg_pkpm_mom_diffusion *up, struct gkyl_range conf_range, const struct gkyl_array *fluid,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    double *cfl = (double *)gkyl_array_fetch(cflrate, loc);
    cfl[0] += pkpm_md_cell(up, &conf_range, idx, fluid, (double *)gkyl_array_fetch(rhs, loc));
  }
}

void
gkyl_dg_pkpm_mom_diffusion_advance_cu(
  const struct gkyl_dg_pkpm_mom_diffusion *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
  gkyl_dg_pkpm_mom_diffusion_advance_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, fluid->on_dev, cflrate->on_dev, rhs->on_dev
  );
}

// The device kernel pointers of the generated kernels.
__global__ static void
dg_pkpm_mom_diffusion_set_cu_dev_ptrs(struct gkyl_dg_pkpm_mom_diffusion *up)
{
  for (int k = 0; k < up->num_dir; ++k) {
    up->kern[k] = pkpm_md_choose_kern(up->order, up->cdim, up->dirs[k]);
  }
}

struct gkyl_dg_pkpm_mom_diffusion *
gkyl_dg_pkpm_mom_diffusion_cu_dev_inew(const struct gkyl_dg_pkpm_mom_diffusion_inp *inp)
{
  struct gkyl_dg_pkpm_mom_diffusion *up =
    (struct gkyl_dg_pkpm_mom_diffusion *)gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_mom_diffusion_setup(up, inp);
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  // The device copy needs the device pointer of the Jacobian array; the host object keeps the
  // acquired map for its own lifetime management.
  up->jacob_pos = inp->pos_map->jacob_pos->on_dev;
  struct gkyl_dg_pkpm_mom_diffusion *up_cu =
    (struct gkyl_dg_pkpm_mom_diffusion *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);
  dg_pkpm_mom_diffusion_set_cu_dev_ptrs<<<1, 1>>>(up_cu);
  up->on_dev = up_cu;
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  return up;
}
