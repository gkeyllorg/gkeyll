#include <gkyl_gyrokinetic_priv.h>

// Advance explicit terms with SSP-RK3, then use its accepted timestep for
// implicit physical collisions and orbit relaxation.
struct gkyl_update_status
gyrokinetic_update_op_split(gkyl_gyrokinetic_app *app, double dt0)
{
  // Prepare frozen projectors before changing any species state.
  for (int i = 0; i < app->num_species; ++i) {
    if (!gk_species_orbit_average_prepare(app, &app->species[i])) {
      return (struct gkyl_update_status){.success = false};
    }
  }
  struct gkyl_update_status st = gyrokinetic_update_ssp_rk3(app, dt0);
  if (st.success) {
    if (app->has_implicit_coll_scheme) {
      gyrokinetic_update_implicit_coll(app, st.dt_actual);
    }
    gyrokinetic_update_orbit_average(app, st.dt_actual);
  }
  return st;
}
