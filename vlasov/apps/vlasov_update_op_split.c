#include <gkyl_vlasov_priv.h>

// Take time-step using the SSP-RK3 method for the hyperbolic components
// Then, we use the actual timestep taken with the SSP-RK3 method to update
// fluid-EM coupling and/or BGK collisions implicitly.
struct gkyl_update_status
vlasov_update_op_split(gkyl_vlasov_app* app, double dt0)
{
  struct gkyl_update_status st = vlasov_update_ssp_rk3(app,dt0);

  // Take the implicit timestep for BGK collisions
  if (app->has_implicit_coll_scheme) {
    vlasov_update_implicit_coll(app, st.dt_actual);
  }

  // Take the implicit timestep for fluid-EM coupling
  if (app->has_fluid_em_coupling) {
    vm_fluid_em_coupling_update(app, app->fl_em, app->tcurr, st.dt_actual);
  }

  return st;
}