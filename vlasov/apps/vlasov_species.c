#include <gkyl_vlasov_priv.h>

#include <string.h>

// Vlasov species dispatch over the species container (struct vlasov_species).
// Per-aspect operations apply to whichever aspects (kinetic, fluid) a species owns
// (a PKPM species owns both); the staging phases, the explicit field coupling,
// the completion of a step and the diagnostics, output, restart and release of
// the species as a whole go through the function pointers set by the
// constructors (those of a PKPM species in vm_pkpm_species.c). RK-state arrays
// (fin, fluidin, fout, fluidout) are indexed by species and are NULL where a
// species lacks that aspect.

// --- Methods assigned by the constructors -----------------------------------

// Self moments of a fluid species: the node values of its state for its collision partners.
static void
fluid_calc_self_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  sp->coll_cross->calc_prim(app, sp, fluidin);
}

// Wall losses of a kinetic species: the integrated boundary fluxes of its right-hand side.
static double
kinetic_wall_loss_rate(gkyl_vlasov_app *app, const struct vlasov_species *sp)
{
  return vm_species_bflux_loss_rate(app, &sp->kinetic->bflux);
}

// A fluid species has no wall-loss diagnostic.
static double
fluid_wall_loss_rate(gkyl_vlasov_app *app, const struct vlasov_species *sp)
{
  return 0.0;
}

// Nothing to keep at the end of a time step (kinetic and fluid species).
static void
species_no_complete_step(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
}

// Empty finish-stage method (kinetic and fluid species).
static void
species_no_finish_stage(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, double dt, struct gkyl_array *fout,
  struct gkyl_array *fluidout
)
{
}

// No explicit field coupling: the null field, and fluid species that carry no
// current or whose EM coupling is the implicit operator split (see
// vm_fluid_em_coupling.c).
static void
species_no_field_coupling(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
}

// Current of a momentum-carrying fluid coupled explicitly to the Maxwell field.
static void
fluid_accumulate_current(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
  vm_fluid_species_accumulate_current(app, sp->fluid, fluidin, target);
  double qbymeps = sp->charge / (sp->mass * app->field->info.epsilon0);
  vm_field_accumulate_ghost_current(app, app->field, qbymeps, sp->fluid->mom, target);
}

static void
kinetic_calc_self_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  vm_species_lbo_moms(app, sp->kinetic, &sp->kinetic->lbo, fin);
  vm_species_bgk_moms(app, sp->kinetic, &sp->kinetic->bgk, fin);
  vm_species_source_bgk_moms(app, sp->kinetic, &sp->kinetic->src_bgk, fin);
}

static void
kinetic_calc_coupled_vars(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  vm_species_lbo_cross_moms(app, sp->kinetic, &sp->kinetic->lbo, fin);
}

static void
fluid_calc_coupled_vars(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  // Primitive variables; run after all species' self moments.
  vm_fluid_species_prim_vars(app, sp->fluid, fluidin);
}

// Explicit coupling to Vlasov-Maxwell: accumulate this species' current onto
// the EM RHS ('target' is emout), J_s = q_s * m1i, scaled by -1/epsilon0.
static void
kinetic_accumulate_current(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
  struct vm_species *s = sp->kinetic;
  double qbyeps = sp->charge / app->field->info.epsilon0;

  vm_species_moment_calc(&s->m1i, s->local, app->local, fin);
  gkyl_array_accumulate_range(target, -qbyeps, s->m1i.marr, &app->local);
  vm_field_accumulate_ghost_current(app, app->field, qbyeps, s->m1i.marr, target);
}

// Explicit coupling of a triad species to GR-Maxwell: deposit the current
// through the GR kernel, q/eps0*(rho*beta - alpha*e^i_a*Jhat^a), which needs
// the charge density as well as the current.
static void
kinetic_accumulate_current_gr(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
  struct vm_species *s = sp->kinetic;
  double qbyeps = sp->charge / app->field->info.epsilon0;

  vm_species_moment_calc(&s->m0, s->local, app->local, fin);
  vm_species_moment_calc(&s->m1i, s->local, app->local, fin);
  gkyl_dg_gr_maxwell_current_deposition_advance(
    s->collisionless.calc_current_dep, &app->local, qbyeps, app->field->geom->lapse,
    app->field->geom->shift, app->field->geom->vierb_con, s->m0.marr, s->m1i.marr, target
  );
  vm_field_accumulate_ghost_current(app, app->field, qbyeps, s->m1i.marr, target);
}

// Explicit coupling to Vlasov-Poisson: deposit this species' number density
// onto the source of each potential the field solves for, weighted by the
// species' coupling to that potential (q for the electrostatic potential,
// alpha_g*m for the gravitational one). The sources live on the field's
// potentials, so 'target' is unused.
static void
kinetic_accumulate_poisson_sources(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
  struct vm_species *s = sp->kinetic;

  vm_species_moment_calc(&s->m0, s->local, app->local, fin);
  for (int p = 0; p < app->field->num_pots; ++p) {
    struct vp_potential *pot = &app->field->pots[p];
    double weight = pot->src_weight(pot, sp->charge, sp->mass);
    gkyl_array_accumulate_range(pot->rho, weight, s->m0.marr, &app->local);
  }
}

// Integrated diagnostic moments of the aspect of a kinetic or of a fluid species.
static void
species_calc_integrated_mom_aspects(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm)
{
  if (sp->kinetic) {
    vm_species_calc_integrated_mom(app, sp->kinetic, tm);
  }
  if (sp->fluid) {
    vm_fluid_species_calc_integrated_mom(app, sp->fluid, tm);
  }
}

// Write the evolved state of the aspect of a kinetic or of a fluid species.
static void
species_write_aspects(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame)
{
  if (sp->kinetic) {
    vm_species_write(app, sp->kinetic, tm, frame);
  }
  if (sp->fluid) {
    vm_fluid_species_write(app, sp->fluid, tm, frame);
  }
}

// Write the diagnostic moments of a kinetic species (a fluid species has none).
static void
species_write_mom_aspects(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame)
{
  if (sp->kinetic) {
    vm_species_write_mom(app, sp->kinetic, tm, frame);
  }
}

// Append the integrated diagnostic moments of the aspect of a kinetic or of a fluid species.
static void
species_write_integrated_mom_aspects(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->kinetic) {
    vm_species_write_integrated_mom(app, sp->kinetic);
  }
  if (sp->fluid) {
    vm_fluid_species_write_integrated_mom(app, sp->fluid);
  }
}

// Restart the aspect of a kinetic or of a fluid species from its frame file.
static struct gkyl_app_restart_status
species_read_from_frame_aspects(gkyl_vlasov_app *app, struct vlasov_species *sp, int frame)
{
  vlasov_species_apply_ic(app, sp, 0.0);

  struct gkyl_app_restart_status rstat = {.io_status = GKYL_ARRAY_RIO_SUCCESS};
  if (sp->kinetic) {
    rstat = vm_species_read_from_frame(app, sp->kinetic, frame);
  }
  if (sp->fluid && rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat = vm_fluid_species_read_from_frame(app, sp->fluid, frame);
  }

  // The app clock resumes at the frame time: the stepper evaluates
  // time-dependent forcing and sources at stage times taken from it.
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    app->tcurr = rstat.stime;
  }

  return rstat;
}

// Release the aspect of a kinetic or of a fluid species.
static void
species_release_aspects(const gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->kinetic) {
    vm_species_release(app, sp->kinetic);
    gkyl_free(sp->kinetic);
  }
  if (sp->fluid) {
    vm_fluid_species_release(app, sp->fluid);
    gkyl_free(sp->fluid);
  }
}

// --- Constructors ------------------------------------------------------------

// Construct a species from its input: validate the declared type against the
// input blocks and dispatch to the typed constructor.
void
vlasov_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
)
{
  switch (inp->type) {
    case GKYL_SPECIES_VLASOV:
      assert(inp->kinetic.cells[0] > 0); // kinetic species require a velocity grid
      assert(!inp->fluid.equation); // ... and must not carry a fluid equation object
      vlasov_kinetic_species_init(app, inp, sp);
      break;
    case GKYL_SPECIES_FLUID:
      assert(inp->fluid.equation); // fluid species require an equation object
      assert(inp->kinetic.cells[0] == 0); // ... and must not declare a velocity grid
      vlasov_fluid_species_init(app, inp, sp);
      break;
    case GKYL_SPECIES_PKPM:
      // PKPM species require the grid of the normalized velocity and no equation object (their
      // fluid flux is a moment of the kinetic flux).
      assert(inp->kinetic.cells[0] > 0);
      assert(!inp->fluid.equation);
      vlasov_pkpm_species_init(app, inp, sp);
      break;
    default:
      assert(false); // unsupported species type
      break;
  }
}

void
vlasov_kinetic_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
)
{
  *sp = (struct vlasov_species){};
  sp->type = GKYL_SPECIES_VLASOV;

  strcpy(sp->name, inp->name);
  sp->charge = inp->charge;
  sp->mass = inp->mass;

  sp->kinetic = gkyl_malloc(sizeof(struct vm_species));
  *sp->kinetic = (struct vm_species){};
  sp->kinetic->info = inp->kinetic;
  strcpy(sp->kinetic->name, inp->name);
  sp->kinetic->charge = inp->charge;
  sp->kinetic->mass = inp->mass;

  sp->coll_cross = 0;
  sp->calc_self_moms_func = kinetic_calc_self_moms;
  sp->calc_coupled_vars_func = kinetic_calc_coupled_vars;
  sp->finish_stage_func = species_no_finish_stage;
  sp->complete_step_func = species_no_complete_step;
  sp->calc_integrated_mom_func = species_calc_integrated_mom_aspects;
  sp->write_func = species_write_aspects;
  sp->write_mom_func = species_write_mom_aspects;
  sp->write_integrated_mom_func = species_write_integrated_mom_aspects;
  sp->read_from_frame_func = species_read_from_frame_aspects;
  sp->wall_loss_rate_func = kinetic_wall_loss_rate;
  sp->release_func = species_release_aspects;
  // Explicit field coupling by field type: Maxwell fields take the current
  // density (through the GR kernel for a triad species coupled to GR-Maxwell),
  // Poisson fields the weighted number density of each potential, the null
  // field nothing.
  enum gkyl_field_id field_id = app->field->field_id;
  switch (field_id) {
    case GKYL_FIELD_NULL:
      sp->accumulate_field_coupling_func = species_no_field_coupling;
      break;
    case GKYL_FIELD_PHI:
      sp->accumulate_field_coupling_func = kinetic_accumulate_poisson_sources;
      break;
    default:
      sp->accumulate_field_coupling_func =
        vm_species_has_gr_em_triad_coupling(app, field_id, inp->kinetic.model_id) ?
          kinetic_accumulate_current_gr :
          kinetic_accumulate_current;
      break;
  }
}

void
vlasov_fluid_species_init(
  struct gkyl_vlasov_app *app, const struct gkyl_vlasov_species *inp, struct vlasov_species *sp
)
{
  *sp = (struct vlasov_species){};
  sp->type = GKYL_SPECIES_FLUID;

  strcpy(sp->name, inp->name);
  sp->charge = inp->charge;
  sp->mass = inp->mass;

  sp->fluid = gkyl_malloc(sizeof(struct vm_fluid_species));
  *sp->fluid = (struct vm_fluid_species){};
  sp->fluid->info = inp->fluid;
  strcpy(sp->fluid->name, inp->name);
  sp->fluid->charge = inp->charge;
  sp->fluid->mass = inp->mass;

  sp->coll_cross = &sp->fluid->coll_cross;
  sp->calc_self_moms_func = fluid_calc_self_moms;
  sp->calc_coupled_vars_func = fluid_calc_coupled_vars;
  sp->finish_stage_func = species_no_finish_stage;
  sp->complete_step_func = species_no_complete_step;
  sp->calc_integrated_mom_func = species_calc_integrated_mom_aspects;
  sp->write_func = species_write_aspects;
  sp->write_mom_func = species_write_mom_aspects;
  sp->write_integrated_mom_func = species_write_integrated_mom_aspects;
  sp->read_from_frame_func = species_read_from_frame_aspects;
  sp->wall_loss_rate_func = fluid_wall_loss_rate;
  sp->release_func = species_release_aspects;
  // Field coupling of the stage: the current of a momentum-carrying (Euler or
  // isothermal Euler) fluid when the Maxwell field is coupled explicitly; nothing otherwise.
  bool has_explicit_current = app->use_explicit_sources && app->field->field_id == GKYL_FIELD_E_B &&
                              vm_fluid_em_coupling_supported(sp);
  sp->accumulate_field_coupling_func = has_explicit_current ? fluid_accumulate_current :
                                                              species_no_field_coupling;
}

// --- Wrappers over the constructor-assigned methods --------------------------

// Moments of the species' own stage state (LBO moments and boundary corrections,
// explicit BGK, explicit temperature source; staging phase 1).
void
vlasov_species_calc_self_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  sp->calc_self_moms_func(app, sp, stage, fin, fluidin);
}

// Coupled variables (staging phase 2): LBO cross-collision moments or fluid
// primitive variables. Run after self moments for all species.
void
vlasov_species_calc_coupled_vars(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin
)
{
  sp->calc_coupled_vars_func(app, sp, fin, fluidin);
}

// Accumulate this species' explicit source contribution onto the field's
// target array (Maxwell: current onto emout; Poisson: the weighted number
// density onto the source of each potential, target unused; no-op for the null
// field and for implicitly-coupled fluid species).
void
vlasov_species_accumulate_field_coupling(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *target
)
{
  sp->accumulate_field_coupling_func(app, sp, fin, fluidin, target);
}

// --- Per-aspect operations ---------------------------------------------------

// Compute the time-dependent applied acceleration of each present aspect at a
// stage time.
void
vlasov_species_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr)
{
  if (sp->kinetic) {
    vm_species_calc_app_accel(app, sp->kinetic, tcurr);
  }
  if (sp->fluid) {
    vm_fluid_species_calc_stage_app_accel(app, sp->fluid, tcurr);
  }
}

// Compute the RHS of each present aspect; returns the smallest of their stable
// time-steps.
double
vlasov_species_rhs(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, const struct gkyl_array *emin, struct gkyl_array *fout,
  struct gkyl_array *fluidout
)
{
  double dt = DBL_MAX;
  if (sp->kinetic) {
    dt = fmin(dt, vm_species_rhs(app, sp->kinetic, fin, emin, fout));
  }
  if (sp->fluid) {
    dt = fmin(dt, vm_fluid_species_rhs(app, sp->fluid, fluidin, emin, fluidout));
  }
  return dt;
}

// Compute the moments that rescale adaptive sources (kinetic aspect only). Run for
// every species before any source is accumulated: an adaptive source reads its
// partners' moments.
void
vlasov_species_calc_source_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin
)
{
  if (sp->kinetic) {
    vm_species_calc_source_moms(app, sp->kinetic, fin);
  }
}

// Accumulate the source onto the RHS for each present aspect. Run after every
// species' RHS: a boundary-flux source reads the boundary fluxes its source
// species computed there.
void
vlasov_species_source_rhs(
  gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr, struct gkyl_array *fout,
  struct gkyl_array *fluidout
)
{
  if (sp->kinetic) {
    vm_species_rhs_source(app, sp->kinetic, tcurr, fout);
  }
  if (sp->fluid) {
    vm_fluid_species_rhs_source(app, sp->fluid, fluidout);
  }
}

// Implicit update, phase 1: moments of the input distribution (kinetic aspect only;
// a no-op for species without implicit BGK collisions or temperature source).
void
vlasov_species_calc_implicit_moms(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin, double dt
)
{
  if (sp->kinetic) {
    vm_species_calc_implicit_moms(app, sp->kinetic, fin, dt);
  }
}

// Implicit update, phase 2: the backward-Euler step, fout = fin + dt*(implicit
// RHS). Run after phase 1 for all species (cross collisions read the partners'
// moments).
void
vlasov_species_rhs_implicit(
  gkyl_vlasov_app *app, struct vlasov_species *sp, const struct gkyl_array *fin,
  struct gkyl_array *fout, double dt
)
{
  if (sp->kinetic) {
    vm_species_rhs_implicit(app, sp->kinetic, fin, fout, dt);
  }
}

// Implicit update, phase 3: apply BCs and copy the update into the stage state.
void
vlasov_species_finish_implicit_update(
  gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fout,
  struct gkyl_array *fstage, double tcurr
)
{
  if (sp->kinetic) {
    vm_species_apply_bc(app, sp->kinetic, fout, tcurr);
    vm_species_copy_range(sp->kinetic, fstage, fout, &sp->kinetic->local_ext);
  }
}

// Forward-Euler accumulate: out = dt*out + in, for each present aspect.
void
vlasov_species_step_f(
  struct vlasov_species *sp, double dt, const struct gkyl_array *fin,
  const struct gkyl_array *fluidin, struct gkyl_array *fout, struct gkyl_array *fluidout
)
{
  if (sp->kinetic) {
    vm_species_step_f(sp->kinetic, fout, dt, fin);
  }
  if (sp->fluid) {
    vm_fluid_species_step_f(sp->fluid, fluidout, dt, fluidin);
  }
}

// Finish-stage method of the species, between step_f and the BCs (empty for
// kinetic and fluid species).
void
vlasov_species_finish_stage(
  gkyl_vlasov_app *app, struct vlasov_species *sp, int stage, double dt, struct gkyl_array *fout,
  struct gkyl_array *fluidout
)
{
  sp->finish_stage_func(app, sp, stage, dt, fout, fluidout);
}

// Complete a time step of the species, once its last stage has been accepted
// (empty for kinetic and fluid species).
void
vlasov_species_complete_step(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  sp->complete_step_func(app, sp);
}

// RK buffer of the kinetic (fluid) aspect selected by the stepper.
static struct gkyl_array *
kinetic_rk_buf(struct vm_species *d, enum vm_rk_buf b)
{
  return b == VM_RK_F ? d->f : (b == VM_RK_F1 ? d->f1 : d->fnew);
}

static struct gkyl_array *
fluid_rk_buf(struct vm_fluid_species *f, enum vm_rk_buf b)
{
  return b == VM_RK_F ? f->fluid : (b == VM_RK_F1 ? f->fluid1 : f->fluidnew);
}

// Combine RK buffers, out = c1*b1 + c2*b2, for each present aspect.
void
vlasov_species_combine(
  gkyl_vlasov_app *app, struct vlasov_species *sp, enum vm_rk_buf out, double c1, enum vm_rk_buf b1,
  double c2, enum vm_rk_buf b2
)
{
  if (sp->kinetic) {
    struct vm_species *d = sp->kinetic;
    vm_species_combine(
      d, kinetic_rk_buf(d, out), c1, kinetic_rk_buf(d, b1), c2, kinetic_rk_buf(d, b2), &d->local_ext
    );
  }
  if (sp->fluid) {
    struct vm_fluid_species *f = sp->fluid;
    vm_fluid_species_combine(
      f, fluid_rk_buf(f, out), c1, fluid_rk_buf(f, b1), c2, fluid_rk_buf(f, b2), &app->local_ext
    );
  }
}

// Copy one RK buffer into another, out = inp, for each present aspect.
void
vlasov_species_copy_range(
  gkyl_vlasov_app *app, struct vlasov_species *sp, enum vm_rk_buf out, enum vm_rk_buf inp
)
{
  if (sp->kinetic) {
    struct vm_species *d = sp->kinetic;
    vm_species_copy_range(d, kinetic_rk_buf(d, out), kinetic_rk_buf(d, inp), &d->local_ext);
  }
  if (sp->fluid) {
    struct vm_fluid_species *f = sp->fluid;
    vm_fluid_species_copy_range(f, fluid_rk_buf(f, out), fluid_rk_buf(f, inp), &app->local_ext);
  }
}

void
vlasov_species_gather_rk_state(
  gkyl_vlasov_app *app, enum vm_rk_buf in, enum vm_rk_buf out, struct gkyl_array *fin[],
  struct gkyl_array *fout[], const struct gkyl_array *fluidin[], struct gkyl_array *fluidout[]
)
{
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    fin[i] = sp->kinetic ? kinetic_rk_buf(sp->kinetic, in) : 0;
    fout[i] = sp->kinetic ? kinetic_rk_buf(sp->kinetic, out) : 0;
    fluidin[i] = sp->fluid ? fluid_rk_buf(sp->fluid, in) : 0;
    fluidout[i] = sp->fluid ? fluid_rk_buf(sp->fluid, out) : 0;
  }
}

// Apply boundary conditions to each present aspect.
void
vlasov_species_apply_bc(
  gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *f, struct gkyl_array *fluid,
  double tcurr
)
{
  if (sp->kinetic) {
    vm_species_apply_bc(app, sp->kinetic, f, tcurr);
  }
  if (sp->fluid) {
    vm_fluid_species_apply_bc(app, sp->fluid, fluid);
  }
}

// Apply the limiters to the stepped fluid state (fluid aspect only).
void
vlasov_species_limiter(gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fluid)
{
  if (sp->fluid) {
    vm_fluid_species_limiter(app, sp->fluid, fluid);
  }
}

// --- Initial conditions, diagnostics, I/O, restart, and release --------------

// Project initial conditions for each present aspect.
void
vlasov_species_apply_ic(gkyl_vlasov_app *app, struct vlasov_species *sp, double t0)
{
  app->tcurr = t0;
  if (sp->kinetic) {
    struct timespec wtm = gkyl_wall_clock();
    vm_species_apply_ic(app, sp->kinetic, t0);
    app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
  }
  if (sp->fluid) {
    struct timespec wtm = gkyl_wall_clock();
    vm_fluid_species_apply_ic(app, sp->fluid, t0);
    app->stat.init_fluid_species_tm += gkyl_time_diff_now_sec(wtm);
  }
}

// Compute integrated diagnostic moments for each present aspect.
void
vlasov_species_calc_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm)
{
  sp->calc_integrated_mom_func(app, sp, tm);
}

// Compute the integrated L2 norm of the distribution (kinetic aspect only).
void
vlasov_species_calc_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm)
{
  if (sp->kinetic) {
    vm_species_calc_L2(app, sp->kinetic, tm);
  }
}

// Write the evolved state (and attendant diagnostics: sources, emission
// spectra, fluid primitive variables) of each present aspect for this frame.
void
vlasov_species_write(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame)
{
  sp->write_func(app, sp, tm, frame);
}

// Write diagnostic moments (kinetic aspect only).
void
vlasov_species_write_mom(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame)
{
  sp->write_mom_func(app, sp, tm, frame);
}

// Append integrated diagnostic moments for each present aspect.
void
vlasov_species_write_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  sp->write_integrated_mom_func(app, sp);
}

// Append the integrated L2 norm of the distribution (kinetic aspect only).
void
vlasov_species_write_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->kinetic) {
    vm_species_write_L2(app, sp->kinetic);
  }
}

// Append the LTE-correction iteration status (kinetic aspect only).
void
vlasov_species_write_lte_corr_status(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->kinetic) {
    vm_species_lte_write_max_corr_status(app, sp->kinetic);
  }
}

// Read the evolved state of the species from the named file through the method
// of its aspect, which also rebuilds what a restart does not carry (the file of a
// PKPM species holds its distributions; its fluid state is read with a frame).
struct gkyl_app_restart_status
vlasov_species_from_file(gkyl_vlasov_app *app, struct vlasov_species *sp, const char *fname)
{
  struct gkyl_app_restart_status rstat = {.io_status = GKYL_ARRAY_RIO_SUCCESS};
  if (sp->kinetic) {
    rstat = vm_species_from_file(app, sp->kinetic, fname);
  }
  if (sp->type == GKYL_SPECIES_FLUID) {
    rstat = vm_fluid_species_from_file(app, sp->fluid, fname);
  }
  return rstat;
}

// Restart a species from a frame. Every aspect is seeded from the initial
// conditions first, so that fixed-function BC buffers are filled before the
// read overwrites the interior; each aspect then reads its frame file and marks
// its diagnostic dynvectors to append, and the app clock is set to the frame time.
struct gkyl_app_restart_status
vlasov_species_read_from_frame(gkyl_vlasov_app *app, struct vlasov_species *sp, int frame)
{
  return sp->read_from_frame_func(app, sp, frame);
}

void
vlasov_species_gather_state(
  gkyl_vlasov_app *app, struct gkyl_array *distf[], struct gkyl_array *fluid[]
)
{
  int num_species = app->num_species;
  for (int i = 0; i < num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (distf) {
      distf[i] = sp->kinetic ? sp->kinetic->f : 0;
    }
    if (fluid) {
      fluid[i] = sp->fluid ? sp->fluid->fluid : 0;
    }
  }
}

void
vlasov_species_gather_dist(gkyl_vlasov_app *app, const struct gkyl_array *fin[])
{
  int num_species = app->num_species;
  for (int i = 0; i < num_species; ++i) {
    fin[i] = app->species[i].kinetic ? app->species[i].kinetic->f : 0;
  }
}

void
vlasov_species_init_kinetic_aspects(struct gkyl_vm *vm, gkyl_vlasov_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].type == GKYL_SPECIES_VLASOV) {
      vm_species_init(vm, app, app->species[i].kinetic);
    }
  }
}

void
vlasov_species_link_kinetic_aspects(gkyl_vlasov_app *app)
{
  // Emission walls need their impact species' boundary-flux objects.
  for (int i = 0; i < app->num_species; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (app->species[i].type != GKYL_SPECIES_VLASOV) {
      continue;
    }
    if (vms->emit_lo) {
      vm_species_emission_cross_init(app, vms, &vms->bc_emission_lo);
    }
    if (vms->emit_up) {
      vm_species_emission_cross_init(app, vms, &vms->bc_emission_up);
    }
  }

  // Cross collisions need the partners' collision objects.
  for (int i = 0; i < app->num_species; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (app->species[i].type != GKYL_SPECIES_VLASOV) {
      continue;
    }
    vm_species_lbo_cross_init(app, vms, &vms->lbo);
    vm_species_bgk_cross_init(app, vms, &vms->bgk);
  }

  // Sources may create a boundary-flux updater on their source species.
  for (int i = 0; i < app->num_species; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (app->species[i].type == GKYL_SPECIES_VLASOV && vms->source_id) {
      vm_species_source_init(app, vms, &vms->src);
    }
  }
}

void
vlasov_species_init_fluid_aspects(struct gkyl_vm *vm, gkyl_vlasov_app *app)
{
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].type == GKYL_SPECIES_FLUID) {
      vm_fluid_species_init(vm, app, app->species[i].fluid);
    }
  }

  for (int i = 0; i < app->num_species; ++i) {
    struct vm_fluid_species *vmf = app->species[i].fluid;
    if (app->species[i].type == GKYL_SPECIES_FLUID && vmf->source_id) {
      vm_fluid_species_source_init(app, vmf, &vmf->src);
    }
  }

  // Collisions refer to other species: after the containers of all species exist.
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].type == GKYL_SPECIES_FLUID) {
      vm_fluid_species_coll_cross_init(app, &app->species[i]);
    }
  }
}

void
vlasov_species_init_pkpm(struct gkyl_vm *vm, gkyl_vlasov_app *app)
{
  app->num_pkpm_species = 0;
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].type == GKYL_SPECIES_PKPM) {
      app->pkpm_species[app->num_pkpm_species++] = &app->species[i];
    }
  }

  for (int i = 0; i < app->num_pkpm_species; ++i) {
    vm_pkpm_species_init(vm, app, app->pkpm_species[i]);
  }

  // Sources refer to other species (the boundary-flux rate): after every species exists.
  for (int i = 0; i < app->num_pkpm_species; ++i) {
    struct vlasov_species *sp = app->pkpm_species[i];
    vm_pkpm_species_source_init(app, sp, &sp->kinetic->pkpm.src);
  }

  // Cross-species terms need the collision objects of the partners.
  for (int i = 0; i < app->num_pkpm_species; ++i) {
    vm_pkpm_species_coll_cross_init(app, app->pkpm_species[i]);
  }
}

double
vlasov_species_wall_loss_rate(gkyl_vlasov_app *app, const struct vlasov_species *sp)
{
  return sp->wall_loss_rate_func(app, sp);
}

// Release the aspects of the species and the container's allocations.
void
vlasov_species_release(const gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  sp->release_func(app, sp);
}
