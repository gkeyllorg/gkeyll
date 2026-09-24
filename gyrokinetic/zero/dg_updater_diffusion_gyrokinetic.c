#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion_gyrokinetic.h>
#include <gkyl_dg_updater_diffusion_gyrokinetic.h>
#include <gkyl_dg_updater_diffusion_gyrokinetic_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_diffusion_gyrokinetic_acquire_eqn(const struct gkyl_dg_updater_diffusion_gyrokinetic *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_diffusion_gyrokinetic*
gkyl_dg_updater_diffusion_gyrokinetic_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, bool is_diff_const, 
  const bool *diff_in_dir, int diff_order, const struct gkyl_range *diff_range,
  const bool *is_zero_flux_bc,
  const struct gkyl_array *coeff, const struct gkyl_array *jacobgeo_inv, bool use_gpu)
{
  struct gkyl_dg_updater_diffusion_gyrokinetic *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_diffusion_gyrokinetic));

  int pdim = basis->ndim;
  int cdim = cbasis->ndim;
  up->use_gpu = use_gpu;
  bool is_dir_diffusive[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) is_dir_diffusive[d] = diff_in_dir==NULL? true : diff_in_dir[d];

  up->dgeqn = gkyl_dg_diffusion_gyrokinetic_new(basis, cbasis, is_diff_const, is_dir_diffusive,
    diff_order, diff_range, up->use_gpu);

  gkyl_dg_diffusion_gyrokinetic_set_auxfields(up->dgeqn, (struct gkyl_dg_diffusion_gyrokinetic_auxfields) {
    .D = coeff, .jacobgeo_inv = jacobgeo_inv });

  int num_up_dirs = 0;
  for (int d=0; d<cdim; d++) num_up_dirs += is_dir_diffusive[d]? 1 : 0;
  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2*GKYL_MAX_DIM];
  int linc = 0;
  for (int d=0; d<cdim; ++d) {
    if (is_dir_diffusive[d]) up_dirs[linc] = d;
    linc += 1;

    zero_flux_flags[d] = is_zero_flux_bc[d]? 1 : 0;
    zero_flux_flags[d+pdim] = is_zero_flux_bc[d+pdim]? 1 : 0;
  }

  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_diffusion_gyrokinetic_advance(struct gkyl_dg_updater_diffusion_gyrokinetic *up,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_diffusion_gyrokinetic_tm
gkyl_dg_updater_diffusion_gyrokinetic_get_tm(const struct gkyl_dg_updater_diffusion_gyrokinetic *up)
{
  return (struct gkyl_dg_updater_diffusion_gyrokinetic_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_diffusion_gyrokinetic_release(struct gkyl_dg_updater_diffusion_gyrokinetic *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
