#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion_fluid.h>
#include <gkyl_dg_updater_diffusion_fluid.h>
#include <gkyl_dg_updater_diffusion_fluid_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_diffusion_fluid_acquire_eqn(const struct gkyl_dg_updater_diffusion_fluid *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_diffusion_fluid*
gkyl_dg_updater_diffusion_fluid_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, bool is_diff_const, int num_equations,
  const bool *diff_in_dir, int diff_order, const struct gkyl_range *diff_range,
  const bool *is_zero_flux_dir, bool use_gpu)
{
  struct gkyl_dg_updater_diffusion_fluid *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_diffusion_fluid));

  int ndim = basis->ndim;
  up->use_gpu = use_gpu;
  bool is_dir_diffusive[GKYL_MAX_CDIM];
  for (int d=0; d<ndim; d++) is_dir_diffusive[d] = diff_in_dir==NULL? true : diff_in_dir[d];

  up->dgeqn = gkyl_dg_diffusion_fluid_new(basis, is_diff_const, num_equations, is_dir_diffusive,
                                          diff_order, diff_range, up->use_gpu);

  int num_up_dirs = 0;
  for (int d=0; d<ndim; d++) num_up_dirs += is_dir_diffusive[d]? 1 : 0;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2*GKYL_MAX_DIM];
  int linc = 0;
  for (int d=0; d<ndim; ++d) {
    if (is_dir_diffusive[d]) up_dirs[linc] = d;
    linc += 1;

    zero_flux_flags[d] = zero_flux_flags[d+ndim] = is_zero_flux_dir[d]? 1 : 0;
  }

  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_diffusion_fluid_advance(struct gkyl_dg_updater_diffusion_fluid *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // Set arrays needed and call the specific advance method required
  gkyl_dg_diffusion_fluid_set_auxfields(up->dgeqn, (struct gkyl_dg_diffusion_fluid_auxfields) { .D = coeff });
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_diffusion_fluid_tm
gkyl_dg_updater_diffusion_fluid_get_tm(const struct gkyl_dg_updater_diffusion_fluid *up)
{
  return (struct gkyl_dg_updater_diffusion_fluid_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_diffusion_fluid_release(struct gkyl_dg_updater_diffusion_fluid *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
