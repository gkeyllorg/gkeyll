#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion_gen.h>
#include <gkyl_dg_updater_diffusion_gen.h>
#include <gkyl_dg_updater_diffusion_gen_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_diffusion_gen_acquire_eqn(const struct gkyl_dg_updater_diffusion_gen *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_diffusion_gen*
gkyl_dg_updater_diffusion_gen_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_range *diff_range, bool use_gpu)
{
  struct gkyl_dg_updater_diffusion_gen *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_diffusion_gen));

  up->use_gpu = use_gpu;
  up->dgeqn = gkyl_dg_diffusion_gen_new(basis, diff_range, up->use_gpu);

  int ndim = basis->ndim;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2*GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    up_dirs[d] = d;
    zero_flux_flags[d] = zero_flux_flags[d+ndim] = 0;
  }

  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, ndim, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_diffusion_gen_advance(struct gkyl_dg_updater_diffusion_gen *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // Set arrays needed and call the specific advance method required
  gkyl_diffusion_gen_set_auxfields(up->dgeqn, (struct gkyl_dg_diffusion_gen_auxfields) { .Dij = coeff });
#ifdef GKYL_HAVE_CUDA
//    if (up->use_gpu)
//      // hyper_dg_gen_stencil NOT YET IMPLEMENTED ON DEVICE
//      gkyl_hyper_dg_gen_stencil_advance_cu(up->hyperdg, update_rng, fIn, cflrate, rhs);
//    else
//      gkyl_hyper_dg_gen_stencil_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
    assert(false);
#else
    gkyl_hyper_dg_gen_stencil_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
#endif
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_diffusion_gen_tm
gkyl_dg_updater_diffusion_gen_get_tm(const struct gkyl_dg_updater_diffusion_gen *up)
{
  return (struct gkyl_dg_updater_diffusion_gen_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_diffusion_gen_release(struct gkyl_dg_updater_diffusion_gen *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
