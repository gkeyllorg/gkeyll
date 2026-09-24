#include <gkyl_vlasov_priv.h>

// Take time-step using the SSP-RK3 method for the hyperbolic components, then
// use the actual timestep taken to advance the implicit half (BGK collisions
// and/or fluid-EM coupling) in one pass.
struct gkyl_update_status
vlasov_update_op_split(gkyl_vlasov_app* app, double dt0)
{
  struct gkyl_update_status st = vlasov_update_ssp_rk3(app, dt0);

  vlasov_update_implicit(app, st.dt_actual);

  return st;
}