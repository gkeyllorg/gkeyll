#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_pkpm.h>
#include <gkyl_dg_updater_pkpm_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

gkyl_dg_updater_pkpm*
gkyl_dg_updater_pkpm_new(const struct gkyl_rect_grid *conf_grid, const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const bool *is_zero_flux_dir, 
  struct gkyl_dg_vlasov_pkpm_auxfields *vlasov_pkpm_inp, struct gkyl_dg_euler_pkpm_auxfields *euler_pkpm_inp, 
  bool use_gpu)
{
  gkyl_dg_updater_pkpm *up = gkyl_malloc(sizeof(gkyl_dg_updater_pkpm));
  up->use_gpu = use_gpu;

  up->eqn_vlasov = gkyl_dg_vlasov_pkpm_new(conf_basis, phase_basis, conf_range, phase_range, up->use_gpu);
  gkyl_vlasov_pkpm_set_auxfields(up->eqn_vlasov, *vlasov_pkpm_inp);

  up->eqn_fluid = gkyl_dg_euler_pkpm_new(conf_basis, conf_range, up->use_gpu);
  gkyl_euler_pkpm_set_auxfields(up->eqn_fluid, *euler_pkpm_inp);

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim;
  int vdim = pdim-cdim;
  int up_dirs_conf[GKYL_MAX_DIM], zero_flux_flags_conf[2*GKYL_MAX_DIM];
  int up_dirs_phase[GKYL_MAX_DIM], zero_flux_flags_phase[2*GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) {
    up_dirs_conf[d] = d;
    up_dirs_phase[d] = d;
    zero_flux_flags_conf[d] = zero_flux_flags_conf[d+cdim] = is_zero_flux_dir[d] ? 1 : 0;
    zero_flux_flags_phase[d] = zero_flux_flags_phase[d+pdim] = is_zero_flux_dir[d] ? 1 : 0;
  }
  int num_up_dirs_conf = cdim;

  for (int d=cdim; d<pdim; ++d) {
    up_dirs_phase[d] = d;
    zero_flux_flags_phase[d] = zero_flux_flags_phase[d+pdim] = 1; // zero-flux BCs in vel-space
  }
  int num_up_dirs_phase = pdim;

  up->up_vlasov = gkyl_hyper_dg_new(phase_grid, phase_basis, up->eqn_vlasov, num_up_dirs_phase, up_dirs_phase, zero_flux_flags_phase, 1, up->use_gpu);
  up->up_fluid = gkyl_hyper_dg_new(conf_grid, conf_basis, up->eqn_fluid, num_up_dirs_conf, up_dirs_conf, zero_flux_flags_conf, 1, up->use_gpu);

  up->vlasov_tm = 0.0;
  up->fluid_tm = 0.0;
  
  return up;
}

void
gkyl_dg_updater_pkpm_advance(gkyl_dg_updater_pkpm *pkpm,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng, 
  const struct gkyl_array* GKYL_RESTRICT fIn, const struct gkyl_array* GKYL_RESTRICT fluidIn, 
  struct gkyl_array* GKYL_RESTRICT cflrate_f, struct gkyl_array* GKYL_RESTRICT cflrate_fluid, 
  struct gkyl_array* GKYL_RESTRICT rhs_f, struct gkyl_array* GKYL_RESTRICT rhs_fluid)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(pkpm->up_vlasov, update_phase_rng, fIn, cflrate_f, rhs_f);
  pkpm->vlasov_tm += gkyl_time_diff_now_sec(wst);

  wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(pkpm->up_fluid, update_conf_rng, fluidIn, cflrate_fluid, rhs_fluid);
  pkpm->fluid_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_pkpm_tm
gkyl_dg_updater_pkpm_get_tm(const gkyl_dg_updater_pkpm *pkpm)
{
  return (struct gkyl_dg_updater_pkpm_tm) {
    .vlasov_tm = pkpm->vlasov_tm,
    .fluid_tm = pkpm->fluid_tm,
  };
}

void
gkyl_dg_updater_pkpm_release(gkyl_dg_updater_pkpm* pkpm)
{
  gkyl_dg_eqn_release(pkpm->eqn_vlasov);
  gkyl_dg_eqn_release(pkpm->eqn_fluid);
  gkyl_hyper_dg_release(pkpm->up_vlasov);
  gkyl_hyper_dg_release(pkpm->up_fluid);
  gkyl_free(pkpm);
}
