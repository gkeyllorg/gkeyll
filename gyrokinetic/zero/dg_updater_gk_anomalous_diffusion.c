#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_gk_anomalous_diffusion.h>
#include <gkyl_dg_updater_gk_anomalous_diffusion.h>
#include <gkyl_dg_updater_gk_anomalous_diffusion_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_gk_anomalous_diffusion_acquire_eqn(const struct gkyl_dg_updater_gk_anomalous_diffusion *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_gk_anomalous_diffusion*
gkyl_dg_updater_gk_anomalous_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  enum gkyl_gyrokinetic_bc_type bc_x_lower, enum gkyl_gyrokinetic_bc_type bc_x_upper,
  const struct gkyl_array *nu, const struct gkyl_array *jacobgeo_inv, bool use_gpu)
{
  struct gkyl_dg_updater_gk_anomalous_diffusion *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_gk_anomalous_diffusion));

  int pdim = basis->ndim;
  int cdim = cbasis->ndim;
  up->use_gpu = use_gpu;

  // Choice of boundary_surf and boundary_diag kernels:
  //   boundary_surf: zero_flux or local
  //   boundary_diag: local or recovery
  // MF 2025/09/10: as of now these options are meant for (here
  // N/A means not applicable):
  //            bound_surf  bound_diag  hyper_dg-zero_flux
  // SKIP:       N/A         recovery    no
  // ABSORB:     N/A         N/A         no
  // PERIODIC:   N/A         N/A         no
  // FIXED_FUNC: N/A         N/A         no
  // ZERO_FLUX:  zero_flux   N/A         yes
  // ELSE:       local       local       yes

  up->dgeqn = gkyl_gk_anomalous_diffusion_new(basis, cbasis,
    conf_range, bc_x_lower, bc_x_upper, up->use_gpu);

  gkyl_gk_anomalous_diffusion_set_auxfields(up->dgeqn,
    (struct gkyl_gk_anomalous_diffusion_auxfields) {.nu = nu, .jacobgeo_inv = jacobgeo_inv });

  int num_up_dirs = 1;
  int up_dirs[GKYL_MAX_DIM];
  up_dirs[0] = 0;

  // Determine if hyper_dg should apply boundary_surf kernels.
  int use_boundary_surf[2*GKYL_MAX_DIM] = {0};
  if ( !((bc_x_lower == GKYL_BC_GK_SKIP) ||
         (bc_x_lower == GKYL_BC_GK_SPECIES_ABSORB) ||
         (bc_x_lower == GKYL_BC_GK_SPECIES_PERIODIC) ||
         (bc_x_lower == GKYL_BC_GK_SPECIES_FIXED_FUNC)) )
    use_boundary_surf[0]      = 1;
  if ( !((bc_x_upper == GKYL_BC_GK_SKIP) ||
         (bc_x_upper == GKYL_BC_GK_SPECIES_ABSORB) ||
         (bc_x_upper == GKYL_BC_GK_SPECIES_PERIODIC) ||
         (bc_x_upper == GKYL_BC_GK_SPECIES_FIXED_FUNC)) )
    use_boundary_surf[0+pdim] = 1;

  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, num_up_dirs, up_dirs, use_boundary_surf, 1, up->use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_gk_anomalous_diffusion_advance(struct gkyl_dg_updater_gk_anomalous_diffusion *up,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_gk_anomalous_diffusion_tm
gkyl_dg_updater_gk_anomalous_diffusion_get_tm(const struct gkyl_dg_updater_gk_anomalous_diffusion *up)
{
  return (struct gkyl_dg_updater_gk_anomalous_diffusion_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_gk_anomalous_diffusion_release(struct gkyl_dg_updater_gk_anomalous_diffusion *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
