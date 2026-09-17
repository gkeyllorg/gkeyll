#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_fpo_vlasov_diff.h>
#include <gkyl_dg_fpo_vlasov_drag.h>
#include <gkyl_dg_updater_fpo_vlasov.h>
#include <gkyl_dg_updater_collisions_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_updater_collisions*
gkyl_dg_updater_fpo_vlasov_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *phase_range, bool use_gpu)
{
  struct gkyl_dg_updater_collisions *up = gkyl_malloc(sizeof(gkyl_dg_updater_collisions));
  up->use_gpu = use_gpu;
  up->coll_drag = gkyl_dg_fpo_vlasov_drag_new(pbasis, phase_range, up->use_gpu);
  up->coll_diff = gkyl_dg_fpo_vlasov_diff_new(pbasis, phase_range, up->use_gpu);

  int pdim = pbasis->ndim;
  int vdim = 3;
  int cdim = pdim - vdim;
  int num_up_dirs = vdim;
  int up_dirs[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<vdim; ++d)
    up_dirs[d] = d + pbasis->ndim - vdim;

  int zero_flux_flags[2*GKYL_MAX_DIM] = { 0 };
  for (int d=cdim; d<pdim; ++d)
    zero_flux_flags[d] = zero_flux_flags[d+pdim] = 1;
  
  up->diff = gkyl_hyper_dg_new(grid, pbasis, up->coll_diff, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);
  up->drag = gkyl_hyper_dg_new(grid, pbasis, up->coll_drag, num_up_dirs, up_dirs, zero_flux_flags, 1, up->use_gpu);

  up->diff_tm = 0.0; up->drag_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_fpo_vlasov_advance(struct gkyl_dg_updater_collisions *fpo,
  long offsets[36],
  const struct gkyl_range *update_rng,
  const struct gkyl_array *drag_coeff, 
  const struct gkyl_array *drag_coeff_surf,
  const struct gkyl_array *sgn_drag_coeff_surf,
  const struct gkyl_array *const_sgn_drag_coeff_surf,
  const struct gkyl_array *diff_coeff, 
  const struct gkyl_array *diff_coeff_surf,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  // Set arrays needed
  gkyl_fpo_vlasov_drag_set_auxfields(fpo->coll_drag,
    (struct gkyl_dg_fpo_vlasov_drag_auxfields) { 
      .drag_coeff = drag_coeff, 
      .drag_coeff_surf = drag_coeff_surf,
      .sgn_drag_coeff_surf = sgn_drag_coeff_surf,
      .const_sgn_drag_coeff_surf = const_sgn_drag_coeff_surf
  });

  gkyl_fpo_vlasov_diff_set_auxfields(fpo->coll_diff,
    (struct gkyl_dg_fpo_vlasov_diff_auxfields) { 
      .diff_coeff = diff_coeff,
      .diff_coeff_surf = diff_coeff_surf
  });

  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(fpo->drag, update_rng, fIn, cflrate, rhs);
  fpo->drag_tm += gkyl_time_diff_now_sec(wst);

  // Fokker-Planck diffusion requires generalized hyper dg operator due to 
  // off diagonal terms in diffusion tensor and mixed partial derivatives
  wst = gkyl_wall_clock();
  gkyl_hyper_dg_gen_stencil_advance(fpo->diff, offsets, update_rng, fIn, cflrate, rhs);
  fpo->diff_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_fpo_vlasov_tm
gkyl_dg_updater_fpo_vlasov_get_tm(const struct gkyl_dg_updater_collisions *coll)
{
  return (struct gkyl_dg_updater_fpo_vlasov_tm) {
    .diff_tm = coll->diff_tm,
    .drag_tm = coll->drag_tm
  };
}

void
gkyl_dg_updater_fpo_vlasov_release(struct gkyl_dg_updater_collisions* coll)
{
  gkyl_dg_eqn_release(coll->coll_diff);
  gkyl_dg_eqn_release(coll->coll_drag);
  gkyl_hyper_dg_release(coll->drag);
  gkyl_hyper_dg_release(coll->diff);
  gkyl_free(coll);
}
