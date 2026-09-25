#include <gkyl_vlasov_priv.h>

// Vlasov field dispatch: vlasov_field_new selects the concrete field type, and
// the wrappers below forward to the function pointers set by its constructor.

// --- Null field (GKYL_FIELD_NULL) --------------------------------------------
// Used when no field is present. A zeroed em is allocated (em1, emnew and
// em_host alias it) so the app holds a valid field and calls its methods
// uniformly; every method is a no-op.

static double
null_field_update(
  gkyl_vlasov_app *app, double tcurr, const struct gkyl_array *fin[], const struct gkyl_array *emin,
  struct gkyl_array *emout
)
{
  return DBL_MAX; // no field force -> no CFL constraint
}

static void
null_field_calc(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
}

static void
null_field_apply_ic(
  gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *fin[], double t0
)
{
}

static void
null_field_calc_energy(
  gkyl_vlasov_app *app, double tm, struct vm_field *field, const struct gkyl_array *fin[]
)
{
}

static void
null_field_write(gkyl_vlasov_app *app, double tm, int frame, const struct gkyl_array *fin[])
{
}

static void
null_field_write_energy(gkyl_vlasov_app *app)
{
}

static struct gkyl_app_restart_status
null_field_from_file(gkyl_vlasov_app *app, struct vm_field *field, const char *fname)
{
  return (struct gkyl_app_restart_status){
    .io_status = GKYL_ARRAY_RIO_SUCCESS,
    .frame = 0,
    .stime = 0.0,
  };
}

static void
null_field_release(const gkyl_vlasov_app *app, struct vm_field *f)
{
  gkyl_array_release(f->em_host);
  gkyl_array_release(f->emnew);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->em);
  gkyl_free(f);
}

static struct vm_field *
null_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  struct vm_field *f = gkyl_malloc(sizeof(struct vm_field));
  f->info = vm->field;
  f->field_id = GKYL_FIELD_NULL;

  // em1, emnew and em_host alias em; nothing is stepped or written.
  f->em = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->em, 0.0);
  f->em1 = gkyl_array_acquire(f->em);
  f->emnew = gkyl_array_acquire(f->em);
  f->em_host = gkyl_array_acquire(f->em);

  // No external fields/potentials/currents, nothing time-dependent.
  f->has_ext_em = f->ext_em_evolve = false;
  f->has_ext_pot = f->ext_pot_evolve = false;
  f->has_app_current = f->app_current_evolve = false;
  f->is_first_energy_write_call = true;

  f->update_func = null_field_update;
  f->combine_func = vp_field_combine; // no-op (shared with Vlasov-Poisson)
  f->copy_range_func = vp_field_copy_range; // no-op
  f->apply_ic_func = null_field_apply_ic;
  f->apply_bc_func = vp_field_apply_bc; // no-op
  f->limiter_func = vp_field_limiter; // no-op
  f->complete_update_func = vp_field_complete_update; // no-op
  f->calc_ext_em_func = null_field_calc;
  f->calc_app_current_func = null_field_calc;
  f->calc_ext_pot_func = null_field_calc;
  f->calc_energy_func = null_field_calc_energy;
  f->write_func = null_field_write;
  f->write_energy_func = null_field_write_energy;
  f->from_file_func = null_field_from_file;
  f->release_func = null_field_release;

  return f;
}

struct vm_field *
vlasov_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  // A field object is always created: the null field when no field is present.
  if (vm->skip_field) {
    return null_field_new(vm, app);
  }
  return vm->is_electrostatic ? vp_field_new(vm, app) : vm_field_new(vm, app);
}

double
vlasov_field_update(
  gkyl_vlasov_app *app, double tcurr, const struct gkyl_array *fin[], const struct gkyl_array *emin,
  struct gkyl_array *emout
)
{
  return app->field->update_func(app, tcurr, fin, emin, emout);
}

void
vlasov_field_combine(
  gkyl_vlasov_app *app, struct gkyl_array *out, double c1, const struct gkyl_array *arr1, double c2,
  const struct gkyl_array *arr2
)
{
  app->field->combine_func(app, app->field, out, c1, arr1, c2, arr2);
}

void
vlasov_field_copy_range(gkyl_vlasov_app *app, struct gkyl_array *out, const struct gkyl_array *inp)
{
  app->field->copy_range_func(app, app->field, out, inp);
}

void
vlasov_field_apply_ic(gkyl_vlasov_app *app, const struct gkyl_array *fin[], double t0)
{
  app->field->apply_ic_func(app, app->field, fin, t0);
}

void
vlasov_field_apply_bc(gkyl_vlasov_app *app, struct gkyl_array *em)
{
  app->field->apply_bc_func(app, app->field, em);
}

void
vlasov_field_limiter(gkyl_vlasov_app *app, struct gkyl_array *em)
{
  app->field->limiter_func(app, app->field, em);
}

void
vlasov_field_complete_update(
  gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  app->field->complete_update_func(app, dt, fin, fluidin, emin, emout);
}

void
vlasov_field_calc_ext_em(gkyl_vlasov_app *app, double tm)
{
  app->field->calc_ext_em_func(app, app->field, tm);
}

void
vlasov_field_calc_app_current(gkyl_vlasov_app *app, double tm)
{
  app->field->calc_app_current_func(app, app->field, tm);
}

void
vlasov_field_calc_ext_pot(gkyl_vlasov_app *app, double tm)
{
  app->field->calc_ext_pot_func(app, app->field, tm);
}

void
vlasov_field_calc_energy(gkyl_vlasov_app *app, double tm, const struct gkyl_array *fin[])
{
  app->field->calc_energy_func(app, tm, app->field, fin);
}

void
vlasov_field_write(gkyl_vlasov_app *app, double tm, int frame, const struct gkyl_array *fin[])
{
  app->field->write_func(app, tm, frame, fin);
}

void
vlasov_field_write_energy(gkyl_vlasov_app *app)
{
  app->field->write_energy_func(app);
}

struct gkyl_app_restart_status
vlasov_field_from_file(gkyl_vlasov_app *app, const char *fname)
{
  return app->field->from_file_func(app, app->field, fname);
}

struct gkyl_app_restart_status
vlasov_field_read_from_frame(gkyl_vlasov_app *app, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-field_%d.gkyl", app->name, frame);
  struct gkyl_app_restart_status rstat = vlasov_field_from_file(app, fileNm.str);
  cstr_drop(&fileNm);

  app->field->is_first_energy_write_call = false; // append to existing diagnostic

  // The app clock resumes at the frame time. Every species restart sets it
  // from the species file (all field types); the field sets it here too when
  // it has a file of its own (Maxwell, GR Maxwell), so a field-only restart is
  // consistent as well, including a rewind to frame zero. The Poisson and null
  // fields return a placeholder status and leave the clock to the species.
  bool has_field_file = app->field->field_id == GKYL_FIELD_E_B ||
                        app->field->field_id == GKYL_FIELD_GR_D_B;
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS && has_field_file) {
    app->tcurr = rstat.stime;
  }

  return rstat;
}

void
vlasov_field_release(gkyl_vlasov_app *app)
{
  app->field->release_func(app, app->field);
}
