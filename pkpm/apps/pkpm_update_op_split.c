#include <gkyl_pkpm_priv.h>

// Staggered operator split of the locally implicit fluid-Maxwell source update S and the SSP-RK3
// step T of the hyperbolic terms (plan 17.34). Over a run the sequence is
//   S(h/2) T S(h) T S(h) T ... T S(h/2),
// the Strang split with the two source half steps that meet between consecutive transport steps
// taken as ONE implicit-midpoint step. Two half steps rotate the plasma oscillation by
// 4 atan(omega_pe h/4) per step, which reaches pi at omega_pe h = 4 and makes the curl coupling (the
// frozen c^2 curl B current in S, the B drift with frozen E in T) resonate; one step rotates it by
// 2 atan(omega_pe h/2) < pi. So between updates the source-updated state (fluid, E) lags tcurr by half
// the last step (E at half steps, B at whole steps).
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
struct gkyl_update_status
pkpm_update_op_split(gkyl_pkpm_app *app, double dt0)
{
  pkpm_op_split_unsync(app);

  if (app->op_split_dt_suggested > 0.0) {
    dt0 = fmin(dt0, app->op_split_dt_suggested);
  }
  double dt_prev = app->op_split_dt_prev;
  double dt_plan = dt_prev > 0.0 ? dt0 : 0.0;

  // sources from the midpoint of the previous transport step to the planned midpoint of this one
  double tau = 0.5 * (dt_prev + dt_plan);
  if (tau > 0.0) {
    pkpm_fluid_em_coupling_update(app, app->pkpm_em, app->tcurr - 0.5 * dt_prev + 0.5 * tau, tau);
  }

  struct gkyl_update_status st = pkpm_update_explicit_ssp_rk3(app, dt0);

  // the actual midpoint differs from the planned one when the step changed
  double rem = 0.5 * (st.dt_actual - dt_plan);
  if (rem != 0.0) {
    pkpm_fluid_em_coupling_update(app, app->pkpm_em, app->tcurr + 0.5 * dt_plan + 0.5 * rem, rem);
  }

  app->op_split_dt_prev = st.dt_actual;
  app->op_split_dt_suggested = st.dt_suggested;

  return st;
}

void
pkpm_op_split_sync(gkyl_pkpm_app *app)
{
  if (app->use_explicit_source || app->op_split_synced || !(app->op_split_dt_prev > 0.0)) {
    return;
  }
  for (int i = 0; i < app->num_species; ++i) {
    gkyl_array_copy(app->op_split_fluid[i], app->species[i].fluid);
  }
  gkyl_array_copy(app->op_split_em, app->field->em);
  // the positivity limiter of the source step is not part of the evolution: keep it out of both
  // records (the statistics of the dynvec and the per-cell counts of the frames)
  long pos_count[GKYL_MAX_SPECIES][PKPM_POS_NUM_STAT];
  for (int i = 0; i < app->num_species; ++i) {
    for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
      pos_count[i][c] = app->species[i].pos_count_since[c];
    }
    gkyl_array_copy(app->op_split_pos_count[i], app->species[i].pkpm_pos_count);
  }
  // tcurr is the end of the last transport step; the sources stand at its midpoint (applied fields
  // at the midpoint of the half step, as in the inverse step of a restart)
  double dt_half = 0.5 * app->op_split_dt_prev;
  pkpm_fluid_em_coupling_update(app, app->pkpm_em, app->tcurr - 0.5 * dt_half, dt_half);
  for (int i = 0; i < app->num_species; ++i) {
    for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
      app->species[i].pos_count_since[c] = pos_count[i][c];
    }
    gkyl_array_copy(app->species[i].pkpm_pos_count, app->op_split_pos_count[i]);
  }
  app->op_split_synced = true;
}

void
pkpm_op_split_unsync(gkyl_pkpm_app *app)
{
  if (!app->op_split_synced) {
    return;
  }
  for (int i = 0; i < app->num_species; ++i) {
    gkyl_array_copy(app->species[i].fluid, app->op_split_fluid[i]);
  }
  gkyl_array_copy(app->field->em, app->op_split_em);
  app->op_split_synced = false;
}
