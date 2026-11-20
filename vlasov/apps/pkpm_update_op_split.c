#include <gkyl_pkpm_priv.h>

// Staggered operator split of the implicit fluid-Maxwell source step S and the SSP-RK3 transport
// step T. Each update is the Strang step S(h/2) T(h) S(h/2) of the step h it takes, with the two
// source half steps between consecutive transport steps taken as one step; between updates the
// fluids and E lag tcurr by half the last step. The first update of a run takes a trial forward
// Euler step for h. If the CFL bound of a stage of T is more than 1% below the planned step, the
// fluids and the field are restored, the source step is repeated over the smaller interval and T
// is attempted again. Applied fields, currents and accelerations are evaluated at the midpoint of
// each source interval.
//
// The positivity limiter of the fluids is applied after the source step of an update, and after
// that of the synchronization for output, so that the fluid that is written is limited as the
// fluid that is evolved. Each frame is written with the staggered fluids and field (the files
// _restart_<frame>), from which a restart continues the run.

// Step the CFL condition allows from the present state, at most dt0: a forward Euler step into the
// stage arrays of the SSP-RK3 step, which overwrites them.
static double
op_split_trial_dt(gkyl_pkpm_app *app, double dt0)
{
  int ns = app->num_species;
  struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];
  const struct gkyl_array *fluidin[ns];
  struct gkyl_array *fluidout[ns];
  struct gkyl_update_status st = {.success = true};

  for (int i = 0; i < ns; ++i) {
    fin[i] = app->species[i].kinetic->f;
    fluidin[i] = app->species[i].fluid->fluid;
    fout[i] = app->species[i].kinetic->f1;
    fluidout[i] = app->species[i].fluid->fluid1;
  }
  pkpm_forward_euler(
    app, 1, app->tcurr, dt0, fin, fluidin, app->field->em, fout, fluidout, app->field->em1, &st
  );
  return st.dt_actual;
}

struct gkyl_update_status
pkpm_update_op_split(gkyl_pkpm_app *app, double dt0)
{
  pkpm_op_split_unsync(app);

  if (app->op_split_dt_suggested > 0.0) {
    dt0 = fmin(dt0, app->op_split_dt_suggested);
  }
  double dt_prev = app->op_split_dt_prev;
  if (!(dt_prev > 0.0)) {
    // first update of a run: no previous update suggested a step
    dt0 = op_split_trial_dt(app, dt0);
  }

  // the staggered state at the start of the update, for a transport step that takes less than
  // planned
  for (int i = 0; i < app->num_species; ++i) {
    gkyl_array_copy(app->op_split_fluid[i], app->species[i].fluid->fluid);
  }
  gkyl_array_copy(app->op_split_em, app->field->em);

  struct gkyl_update_status st;
  while (true) {
    // sources from the midpoint of the previous transport step (the start of a run) to the
    // planned midpoint of this one
    double tau = 0.5 * (dt_prev + dt0);
    struct timespec wst = gkyl_wall_clock();
    pkpm_fluid_em_coupling_update(app, app->pkpm_em, app->tcurr - 0.5 * dt_prev + 0.5 * tau, tau);
    pkpm_fluid_em_coupling_limit(app, app->pkpm_em);
    app->stat.pkpm_em_tm += gkyl_time_diff_now_sec(wst);

    st = pkpm_update_explicit_ssp_rk3_attempt(app, dt0);
    if (!(st.dt_actual < dt0)) {
      break;
    }
    // the transport step allows less than planned: repeat the source step for the smaller one
    for (int i = 0; i < app->num_species; ++i) {
      gkyl_array_copy(app->species[i].fluid->fluid, app->op_split_fluid[i]);
    }
    gkyl_array_copy(app->field->em, app->op_split_em);
    dt0 = st.dt_actual;
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
  // the limiter of the synchronization is not part of the evolution: keep it out of both records
  // (the statistics of the dynvec and the per-cell counts of the frames)
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
  pkpm_fluid_em_coupling_limit(app, app->pkpm_em);
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
