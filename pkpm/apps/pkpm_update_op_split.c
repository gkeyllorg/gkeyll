#include <gkyl_pkpm_priv.h>

// Staggered operator split of the locally implicit fluid-Maxwell source update S and the SSP-RK3
// step T of the hyperbolic terms. Over a run the sequence is
//   T S(h/2) S(h) T S(h) T ... S(h) T S(h/2):
// the first update has no source step before T and applies S(h/2) after it, every later update
// applies one implicit-midpoint step S(h) before T, and the last S(h/2) is the synchronization
// of the output (pkpm_op_split_sync). Between updates the source-updated state (fluid, E) lags
// tcurr by half the last step (E at half steps, B at whole steps).
//
// Each update applies the sources over the interval between the midpoint of the previous
// transport step and the midpoint of this one, (dt_prev + dt)/2. The source step before T assumes
// the requested step dt0, capped by the step the previous update suggested (a restart restores
// both steps from the frame). A source step after T covers the difference between the planned
// and the actual midpoint: half the step in the first update, and a negative step when T took
// less than planned (the implicit midpoint rule is time-reversible). Time-dependent applied
// fields, currents and accelerations are evaluated at the midpoint of each source interval.
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
  if (app->use_explicit_sources || app->op_split_synced || !(app->op_split_dt_prev > 0.0)) {
    return;
  }
  for (int i = 0; i < app->num_species; ++i) {
    gkyl_array_copy(app->op_split_fluid[i], app->species[i].fluid->fluid);
  }
  gkyl_array_copy(app->op_split_em, app->field->em);
  // the positivity limiter of the source step is not part of the evolution: keep it out of both
  // records (the statistics of the dynvec and the per-cell counts of the frames)
  long pos_count[GKYL_MAX_SPECIES][PKPM_POS_NUM_STAT];
  for (int i = 0; i < app->num_species; ++i) {
    for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
      pos_count[i][c] = app->species[i].fluid->pos_count_since[c];
    }
    gkyl_array_copy(app->op_split_pos_count[i], app->species[i].fluid->pkpm_pos_count);
  }
  // tcurr is the end of the last transport step; the sources stand at its midpoint (applied fields
  // at the midpoint of the half step)
  double dt_half = 0.5 * app->op_split_dt_prev;
  pkpm_fluid_em_coupling_update(app, app->pkpm_em, app->tcurr - 0.5 * dt_half, dt_half);
  for (int i = 0; i < app->num_species; ++i) {
    for (int c = 0; c < PKPM_POS_NUM_STAT; ++c) {
      app->species[i].fluid->pos_count_since[c] = pos_count[i][c];
    }
    gkyl_array_copy(app->species[i].fluid->pkpm_pos_count, app->op_split_pos_count[i]);
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
    gkyl_array_copy(app->species[i].fluid->fluid, app->op_split_fluid[i]);
  }
  gkyl_array_copy(app->field->em, app->op_split_em);
  app->op_split_synced = false;
}
