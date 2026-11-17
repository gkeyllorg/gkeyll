#include <gkyl_vlasov_priv.h>

// Staggered operator split of the locally implicit fluid-Maxwell source update S and the SSP-RK3
// step T of the hyperbolic terms. Over a run the sequence is
//   S(h/2) T S(h) T S(h) T ... T S(h/2),
// the Strang split with the two source half steps that meet between consecutive transport steps
// taken as ONE implicit-midpoint step. Two half steps rotate the plasma oscillation by
// 4 atan(omega_pe h/4) per step, which reaches pi at omega_pe h = 4 and makes the curl coupling (the
// frozen c^2 curl B current in S, the B drift with frozen E in T) resonate; one step rotates it by
// 2 atan(omega_pe h/2) < pi. So between updates the source-updated state (fluids, E) lags tcurr by
// half the last step (E at half steps, B at whole steps).
//
// Each update applies the sources over the interval between the midpoint of the previous transport
// step and the midpoint of this one, (dt_prev + dt)/2. The step T takes is only known afterwards (it
// is bounded by the CFL condition), so the source step before T assumes the requested step dt0, capped
// by the step the previous update suggested, and a source step after T covers the difference. T takes
// dt0 unless its CFL bound is more than 1% smaller, so the difference is zero when the step is steady
// or growing, and negative when T took less (the implicit midpoint rule is time-reversible); a positive
// correction would add rotation to a step that may already be close to pi. The first update of a run
// is a Lie start (no source step before T). The cap by the suggested step also makes a restart (which
// restores both steps) continue as the uninterrupted run. Time-dependent applied fields, currents and
// accelerations are evaluated at the midpoint of each source interval (second order for one step, and
// the same time for a step and its inverse).
//
// The transport step drifts B with the E of the source step over its three stages, so a limited
// field (limit_em) is limited after the source steps of an update too: E enters the transport step
// as the stages of an explicit update leave it; the positivity limiter of the fluids acts there as
// well. The synchronization for output and the inverse step of a restart are the bare source step,
// so that one inverts the other.
struct gkyl_update_status
vlasov_update_op_split(gkyl_vlasov_app *app, double dt0)
{
  vlasov_op_split_unsync(app);

  if (app->op_split_dt_suggested > 0.0) {
    dt0 = fmin(dt0, app->op_split_dt_suggested);
  }
  double dt_prev = app->op_split_dt_prev;
  double dt_plan = dt_prev > 0.0 ? dt0 : 0.0;

  // sources from the midpoint of the previous transport step to the planned midpoint of this one
  double tau = 0.5 * (dt_prev + dt_plan);
  if (tau > 0.0) {
    struct timespec wst = gkyl_wall_clock();
    vm_fluid_em_coupling_update(app, app->fl_em, app->tcurr - 0.5 * dt_prev + 0.5 * tau, tau);
    vm_fluid_em_coupling_limit(app, app->fl_em);
    app->stat.fl_em_tm += gkyl_time_diff_now_sec(wst);
  }

  struct gkyl_update_status st = vlasov_update_ssp_rk3(app, dt0);

  // the actual midpoint differs from the planned one when the step changed
  double rem = 0.5 * (st.dt_actual - dt_plan);
  if (rem != 0.0) {
    struct timespec wst = gkyl_wall_clock();
    vm_fluid_em_coupling_update(app, app->fl_em, app->tcurr + 0.5 * dt_plan + 0.5 * rem, rem);
    vm_fluid_em_coupling_limit(app, app->fl_em);
    app->stat.fl_em_tm += gkyl_time_diff_now_sec(wst);
  }

  app->op_split_dt_prev = st.dt_actual;
  app->op_split_dt_suggested = st.dt_suggested;

  return st;
}

void
vlasov_op_split_sync(gkyl_vlasov_app *app)
{
  if (!app->has_fluid_em_coupling || app->op_split_synced || !(app->op_split_dt_prev > 0.0)) {
    return;
  }
  struct vm_fluid_em_coupling *fl_em = app->fl_em;
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    gkyl_array_copy(fl_em->fluid_stag[i], fl_em->species[i]->fluid->fluid);
  }
  gkyl_array_copy(fl_em->em_stag, app->field->em);
  // tcurr is the end of the last transport step; the sources stand at its midpoint (applied fields
  // at the midpoint of the half step, as in the inverse step of a restart)
  double dt_half = 0.5 * app->op_split_dt_prev;
  vm_fluid_em_coupling_update(app, fl_em, app->tcurr - 0.5 * dt_half, dt_half);
  app->op_split_synced = true;
}

void
vlasov_op_split_unsync(gkyl_vlasov_app *app)
{
  if (!app->op_split_synced) {
    return;
  }
  struct vm_fluid_em_coupling *fl_em = app->fl_em;
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    gkyl_array_copy(fl_em->species[i]->fluid->fluid, fl_em->fluid_stag[i]);
  }
  gkyl_array_copy(app->field->em, fl_em->em_stag);
  app->op_split_synced = false;
}

void
vlasov_op_split_restart(gkyl_vlasov_app *app, int frame)
{
  if (!app->has_fluid_em_coupling) {
    return;
  }
  // The frames hold the synchronized state of the staggered operator split; the uninterrupted run
  // held the sources half its last step behind. Step them back (the implicit midpoint rule is
  // time-reversible: this inverts the synchronization step) and restore the steps, so the next
  // update continues as the uninterrupted run. Frames without the steps restart as a new run (Lie
  // start).
  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz + 1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);
  struct vlasov_output_meta meta = vlasov_meta_from_file(fileNm);

  app->op_split_synced = false;
  app->op_split_dt_prev = meta.dt_prev;
  app->op_split_dt_suggested = meta.dt_suggested;
  if (meta.dt_prev > 0.0) {
    double dt_half = 0.5 * meta.dt_prev;
    vm_fluid_em_coupling_update(app, app->fl_em, app->tcurr - 0.5 * dt_half, -dt_half);
  }
}
