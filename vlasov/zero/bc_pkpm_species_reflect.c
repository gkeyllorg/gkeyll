#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_bc_pkpm_species_reflect.h>
#include <gkyl_bc_pkpm_species_reflect_priv.h>

struct gkyl_bc_pkpm_species_reflect *
gkyl_bc_pkpm_species_reflect_inew(const struct gkyl_bc_pkpm_species_reflect_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_bc_pkpm_species_reflect_cu_dev_inew(inp);
  }
#endif
  assert(bc_pkpm_species_reflect_symmetric_vel(inp->phase_grid, inp->cdim));
  struct gkyl_bc_pkpm_species_reflect *up = gkyl_malloc(sizeof(*up));
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
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

void
gkyl_bc_pkpm_species_reflect_advance(
  const struct gkyl_bc_pkpm_species_reflect *up, const struct gkyl_array *bvar, struct gkyl_array *f
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(f)) {
    gkyl_bc_pkpm_species_reflect_advance_cu(up, bvar, f);
    return;
  }
#endif
  int sidx[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    bool flip_w = bc_pkpm_species_reflect_skin(up, iter.idx, bvar, sidx);
    bc_pkpm_species_reflect_cell(
      up, flip_w, gkyl_array_cfetch(f, gkyl_range_idx(&up->phase_range_ext, sidx)),
      gkyl_array_fetch(f, gkyl_range_idx(&up->phase_range_ext, iter.idx))
    );
  }
}

void
gkyl_bc_pkpm_species_reflect_release(struct gkyl_bc_pkpm_species_reflect *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_bc_pkpm_species_reflect *
gkyl_bc_pkpm_species_reflect_cu_dev_inew(const struct gkyl_bc_pkpm_species_reflect_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_bc_pkpm_species_reflect_advance_cu(
  const struct gkyl_bc_pkpm_species_reflect *up, const struct gkyl_array *bvar, struct gkyl_array *f
)
{
  assert(false);
}

#endif
