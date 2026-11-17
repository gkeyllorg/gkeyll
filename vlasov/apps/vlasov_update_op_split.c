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
// step and the midpoint of this one, (dt_prev + dt)/2. The step T can take is only known once its
// stages are evaluated (it is bounded by the CFL condition), so the source step before T assumes
// the requested step dt0, capped by the step the previous update suggested. T takes dt0 unless the
// CFL bound of one of its stages is more than 1% smaller; then nothing of the attempt is kept: the
// fluids and the field are restored to their state before the source step, the source step is
// repeated over the smaller interval (with the same currents of the kinetic species, which T has
// not moved), and T is attempted again. The first update of a run has no previous step to suggest
// one: a trial forward Euler step from the initial state (into the stage arrays, which T
// overwrites) gives the step the CFL condition allows, and the run starts with the source half
// step S(h/2) over it. Every update is thus the Strang step S(h/2) T(h) S(h/2) of the step h it
// took, which is second order for any sequence of steps, and after every update the
// source-updated state lags tcurr by half of that step, the half step of the synchronization for
// output. The cap by the suggested step also makes a restart (which restores both steps) continue
// as the uninterrupted run. Time-dependent applied fields, currents and accelerations are
// evaluated at the midpoint of each source interval (second order for one step, and the same time
// for a step and its inverse).
//
// Kinetic species (hybrid runs): the source step takes the current of each kinetic species at its
// start as a known current, and the transport step moves the species with the E of the source
// step, which it does not change. This is a leapfrog of the kinetic currents and E (second order,
// the energy error is second order without drift), stable for omega_p dt < 2 at a fixed step, with
// omega_p the plasma frequency of the kinetic species: the implicit fluids do not respond to the
// mode that alternates in sign from step to step. The update therefore bounds the step by
// cfl_frac / omega_p (vm_fluid_em_coupling_dt). The margin of two at cfl_frac = 1 is not derived: a
// changing step perturbs the leapfrog, and a run at omega_p dt = 1.6 whose step the CFL condition
// kept changing went unstable.
//
// The transport step drifts B with the E of the source step over its three stages, so a limited
// field (limit_em) is limited after the source steps of an update too: E enters the transport step
// as the stages of an explicit update leave it; the positivity limiter of the fluids acts there as
// well. The synchronization for output and the inverse step of a restart are the bare source step,
// so that one inverts the other.
// Step the CFL condition allows from the present state, at most dt0: a forward Euler step into the
// stage arrays of the SSP-RK3 step, which overwrites them.
static double
op_split_trial_dt(gkyl_vlasov_app *app, double dt0)
{
  int num_species = app->num_species;
  struct gkyl_array *fin[num_species];
  struct gkyl_array *fout[num_species];
  const struct gkyl_array *fluidin[num_species];
  struct gkyl_array *fluidout[num_species];
  struct gkyl_update_status st = {.success = true};

  vlasov_species_gather_rk_state(app, VM_RK_F, VM_RK_F1, fin, fout, fluidin, fluidout);
  vlasov_forward_euler(
    app, 1, app->tcurr, dt0, fin, fluidin, app->field->em, fout, fluidout, app->field->em1, &st
  );
  return st.dt_actual;
}

struct gkyl_update_status
vlasov_update_op_split(gkyl_vlasov_app *app, double dt0)
{
  vlasov_op_split_unsync(app);

  if (app->op_split_dt_suggested > 0.0) {
    dt0 = fmin(dt0, app->op_split_dt_suggested);
  }
  // the leapfrog of the kinetic currents resolves their plasma frequency
  double dt_src = vm_fluid_em_coupling_dt(app, app->fl_em);
  dt0 = fmin(dt0, dt_src);
  double dt_prev = app->op_split_dt_prev;
  if (!(dt_prev > 0.0)) {
    // first update of a run: no previous update suggested a step
    dt0 = op_split_trial_dt(app, dt0);
  }

  // the staggered state at the start of the update, for a transport step that takes less than
  // planned
  struct vm_fluid_em_coupling *fl_em = app->fl_em;
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    gkyl_array_copy(fl_em->fluid_stag[i], fl_em->species[i]->fluid->fluid);
  }
  gkyl_array_copy(fl_em->em_stag, app->field->em);

  struct gkyl_update_status st;
  while (true) {
    // sources from the midpoint of the previous transport step (the start of a run) to the
    // planned midpoint of this one
    double tau = 0.5 * (dt_prev + dt0);
    struct timespec wst = gkyl_wall_clock();
    vm_fluid_em_coupling_update(app, fl_em, app->tcurr - 0.5 * dt_prev + 0.5 * tau, tau);
    vm_fluid_em_coupling_limit(app, fl_em);
    app->stat.fl_em_tm += gkyl_time_diff_now_sec(wst);

    st = vlasov_update_ssp_rk3_attempt(app, dt0);
    if (!(st.dt_actual < dt0)) {
      break;
    }
    // the transport step allows less than planned: repeat the source step for the smaller one
    for (int i = 0; i < fl_em->num_fluid; ++i) {
      gkyl_array_copy(fl_em->species[i]->fluid->fluid, fl_em->fluid_stag[i]);
    }
    gkyl_array_copy(app->field->em, fl_em->em_stag);
    dt0 = st.dt_actual;
  }

  st.dt_suggested = fmin(st.dt_suggested, dt_src);
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
  // update continues as the uninterrupted run. Frames without the steps restart as a new run.
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
