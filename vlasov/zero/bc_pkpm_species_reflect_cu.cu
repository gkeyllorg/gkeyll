/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_bc_pkpm_species_reflect.h>
#include <gkyl_bc_pkpm_species_reflect_priv.h>
}

// One thread per ghost phase-space cell of the wall.
__global__ static void
gkyl_bc_pkpm_species_reflect_advance_cu_kernel(
  const struct gkyl_bc_pkpm_species_reflect *up, struct gkyl_range ghost_r,
  const struct gkyl_array *bvar, struct gkyl_array *f
)
{
  int gidx[GKYL_MAX_DIM], sidx[GKYL_MAX_DIM];
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < ghost_r.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&ghost_r, linc, gidx);
    bool flip_w = bc_pkpm_species_reflect_skin(up, gidx, bvar, sidx);
    bc_pkpm_species_reflect_cell(
      up, flip_w, (const double *)gkyl_array_cfetch(f, gkyl_range_idx(&up->phase_range_ext, sidx)),
      (double *)gkyl_array_fetch(f, gkyl_range_idx(&up->phase_range_ext, gidx))
    );
  }
}

void
gkyl_bc_pkpm_species_reflect_advance_cu(
  const struct gkyl_bc_pkpm_species_reflect *up, const struct gkyl_array *bvar, struct gkyl_array *f
)
{
  int nblocks = up->ghost_r.nblocks, nthreads = up->ghost_r.nthreads;
  gkyl_bc_pkpm_species_reflect_advance_cu_kernel<<<nblocks, nthreads>>>(
    up->on_dev, up->ghost_r, bvar->on_dev, f->on_dev
  );
}

struct gkyl_bc_pkpm_species_reflect *
gkyl_bc_pkpm_species_reflect_cu_dev_inew(const struct gkyl_bc_pkpm_species_reflect_inp *inp)
{
  assert(bc_pkpm_species_reflect_symmetric_vel(inp->phase_grid, inp->cdim));
  struct gkyl_bc_pkpm_species_reflect *up =
    (struct gkyl_bc_pkpm_species_reflect *)gkyl_malloc(sizeof(*up));
  up->dir = inp->dir;
  up->edge = inp->edge;
  up->cdim = inp->cdim;
  up->num_basis = inp->phase_basis->num_basis;
  up->conf_num_basis = inp->conf_num_basis;
  up->num_blocks = inp->num_blocks;
  up->basis = inp->phase_basis_on_dev;
  up->skin_r = *inp->skin_r;
  up->ghost_r = *inp->ghost_r;
  up->phase_range_ext = *inp->phase_range_ext;
  up->conf_range_ext = *inp->conf_range_ext;
  up->use_gpu = inp->use_gpu;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  struct gkyl_bc_pkpm_species_reflect *up_cu =
    (struct gkyl_bc_pkpm_species_reflect *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;
  return up;
}
