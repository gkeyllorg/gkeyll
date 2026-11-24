#include <gkyl_vlasov_priv.h>

// Staggered operator split of the implicit fluid-Maxwell source step S and the SSP-RK3 transport
// step T. Each update is the Strang step S(h/2) T(h) S(h/2) of the step h it takes, with the two
// source half steps between consecutive transport steps taken as one step; between updates the
// fluids and E lag tcurr by half the last step. The first update of a run takes a trial forward
// Euler step for h. If the CFL bound of a stage of T is more than 1% below the planned step, the
// fluids and the field are restored, the source step is repeated over the smaller interval and T
// is attempted again. Applied fields, currents and accelerations are evaluated at the midpoint of
// each source interval.
//
// Kinetic species: the source step takes the current of each kinetic species at its start as a
// known current, and the transport step moves the species with the E of the source step. This
// leapfrog of the kinetic currents and E is stable for omega_p dt < 2, with omega_p the plasma
// frequency of the kinetic species, so the step is bounded by cfl_frac / omega_p
// (vm_fluid_em_coupling_dt).
//
// The positivity limiter of the fluids and the field limiter (limit_em) are applied after the
// source step of an update, and after that of the synchronization for output, so that the fluids
// and the field that are written are limited as those that are evolved. The species record the
// limiter of the updates that complete, not that of a repeated source step or of the
// synchronization. Each frame is written
// with the staggered fluids and field (the files _restart_<frame>), from which a restart
// continues the run.

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
  // the trial step is not part of the evolution
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_discard_step(app, &app->species[i]);
  }
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
    // what the species recorded in the source step and in the stages goes with them
    for (int i = 0; i < app->num_species; ++i) {
      vlasov_species_discard_step(app, &app->species[i]);
    }
    app->stat.nfl_em_retake += 1;
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
  // at the midpoint of the half step)
  double dt_half = 0.5 * app->op_split_dt_prev;
  struct timespec wst = gkyl_wall_clock();
  vm_fluid_em_coupling_update(app, fl_em, app->tcurr - 0.5 * dt_half, dt_half);
  vm_fluid_em_coupling_limit(app, fl_em);
  app->stat.fl_em_tm += gkyl_time_diff_now_sec(wst);
  // the limiter of the synchronization is not part of the evolution: the species keep it out of
  // their records
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    vlasov_species_discard_step(app, fl_em->species[i]);
  }
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

// File of the staggered fluid state of a species for a frame (the fluid files of a PKPM species
// are named <species>_fluid, its distributions take the name of the species).
static cstr
op_split_fluid_file(const gkyl_vlasov_app *app, const struct vlasov_species *sp, int frame)
{
  const char *fmt = sp->type == GKYL_SPECIES_PKPM ? "%s-%s_fluid_restart_%d.gkyl" :
                                                    "%s-%s_restart_%d.gkyl";
  return cstr_from_fmt(fmt, app->name, sp->name, frame);
}

// Write an array of the staggered state as the stepper holds it, with the steps of the split.
static void
op_split_write_array(
  gkyl_vlasov_app *app, const struct gkyl_array *arr, struct gkyl_array *arr_host, double tm,
  int frame, const char *fname
)
{
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
    .dt_prev = app->op_split_dt_prev,
    .dt_suggested = app->op_split_dt_suggested,
  });
  if (app->use_gpu) {
    // Copy data from device to host before writing it out.
    gkyl_array_copy(arr_host, arr);
    arr = arr_host;
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, arr, fname);
  vlasov_array_meta_release(mt);
}

void
vlasov_op_split_write_field(gkyl_vlasov_app *app, double tm, int frame)
{
  if (!app->has_fluid_em_coupling) {
    return;
  }
  // Before the first update nothing lags, and the field is the staggered one.
  const struct gkyl_array *em = app->op_split_synced ? app->fl_em->em_stag : app->field->em;
  cstr fileNm = cstr_from_fmt("%s-field_restart_%d.gkyl", app->name, frame);
  op_split_write_array(app, em, app->field->em_host, tm, frame, fileNm.str);
  cstr_drop(&fileNm);
}

void
vlasov_op_split_write_species(
  gkyl_vlasov_app *app, const struct vlasov_species *sp, double tm, int frame
)
{
  if (!app->has_fluid_em_coupling) {
    return;
  }
  struct vm_fluid_em_coupling *fl_em = app->fl_em;
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    if (fl_em->species[i] == sp) {
      struct vm_fluid_species *fs = sp->fluid;
      const struct gkyl_array *fluid = app->op_split_synced ? fl_em->fluid_stag[i] : fs->fluid;
      cstr fileNm = op_split_fluid_file(app, sp, frame);
      op_split_write_array(app, fluid, fs->fluid_host, tm, frame, fileNm.str);
      cstr_drop(&fileNm);
    }
  }
}

// Read an array of the staggered state of a frame; the status carries the steps of the split.
static struct gkyl_app_restart_status
op_split_read_array(
  gkyl_vlasov_app *app, struct gkyl_array *arr, struct gkyl_array *arr_host, const char *fname
)
{
  struct gkyl_app_restart_status rstat = vlasov_header_from_file(app, fname);
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, arr_host, fname);
    if (app->use_gpu) {
      gkyl_array_copy(arr, arr_host);
    }
  }
  return rstat;
}

struct gkyl_app_restart_status
vlasov_op_split_restart(gkyl_vlasov_app *app, int frame)
{
  struct gkyl_app_restart_status rstat = {.io_status = GKYL_ARRAY_RIO_SUCCESS};
  if (!app->has_fluid_em_coupling) {
    return rstat;
  }
  // The frames hold the synchronized state; the staggered fluids and field replace it, with the
  // steps of the split, so that the next update continues the run.
  struct vm_fluid_em_coupling *fl_em = app->fl_em;
  cstr fileNm = cstr_from_fmt("%s-field_restart_%d.gkyl", app->name, frame);
  struct vlasov_output_meta meta = vlasov_meta_from_file(fileNm.str);
  rstat = op_split_read_array(app, app->field->em, app->field->em_host, fileNm.str);
  cstr_drop(&fileNm);
  if (rstat.io_status != GKYL_ARRAY_RIO_SUCCESS) {
    return rstat;
  }
  vlasov_field_apply_bc(app, app->field->em);

  for (int i = 0; i < fl_em->num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    cstr fileNm_fluid = op_split_fluid_file(app, fl_em->species[i], frame);
    struct gkyl_app_restart_status rstat_fluid =
      op_split_read_array(app, fs->fluid, fs->fluid_host, fileNm_fluid.str);
    cstr_drop(&fileNm_fluid);
    if (rstat_fluid.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      return rstat_fluid;
    }
    vm_fluid_species_apply_bc(app, fs, fs->fluid);
  }

  app->op_split_synced = false;
  app->op_split_dt_prev = meta.dt_prev;
  app->op_split_dt_suggested = meta.dt_suggested;

  return rstat;
}
