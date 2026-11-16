#include <gkyl_vlasov_priv.h>

// Take a time-step with a first-order operator split: the SSP-RK3 method for
// everything but the fluid-EM coupling, then the implicit fluid-EM coupling
// over the time-step actually taken.
struct gkyl_update_status
vlasov_update_op_split(gkyl_vlasov_app *app, double dt0)
{
  struct gkyl_update_status st = vlasov_update_ssp_rk3(app, dt0);

  struct timespec wst = gkyl_wall_clock();
  vm_fluid_em_coupling_update(app, app->fl_em, app->tcurr, st.dt_actual);
  app->stat.fl_em_tm += gkyl_time_diff_now_sec(wst);

  return st;
}
