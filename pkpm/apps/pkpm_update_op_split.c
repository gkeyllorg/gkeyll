#include <gkyl_pkpm_priv.h>

// Take time-step using the SSP-RK3 method for the hyperbolic components
// Then, we use the actual timestep taken with the SSP-RK3 method to update
// fluid-EM coupling implicitly.
struct gkyl_update_status
pkpm_update_op_split(gkyl_pkpm_app* app, double dt0)
{
  struct gkyl_update_status st = pkpm_update_explicit_ssp_rk3(app, dt0);

  pkpm_fluid_em_coupling_update(app, app->pkpm_em, app->tcurr, st.dt_actual);

  return st;
}