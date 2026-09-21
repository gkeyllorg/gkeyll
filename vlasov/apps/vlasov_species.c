#include <gkyl_vlasov_priv.h>

#include <string.h>

// Vlasov species dispatch over the species container (struct vlasov_species).
// Per-aspect operations apply to whichever aspects (kinetic, fluid) a species owns;
// the staging phases and the explicit field coupling go through the function
// pointers set by the constructors. RK-state arrays (fin, fluidin, fout,
// fluidout) are indexed over the overall species count and are NULL where a
// species lacks that aspect.

// --- Methods assigned by the constructors -----------------------------------

// No-op moments phase for species types with no work in that phase.
static void
species_no_calc_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin) { }

// No explicit field coupling: the null field, and fluid species (whose EM
// coupling is the implicit operator split, see vm_fluid_em_coupling.c).
static void
species_no_field_coupling(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target) { }

static void
kinetic_calc_self_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  vm_species_lbo_moms(app, sp->kinetic, &sp->kinetic->lbo, fin);
  vm_species_bgk_moms(app, sp->kinetic, &sp->kinetic->bgk, fin);
}

static void
kinetic_calc_coupled_vars(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  vm_species_lbo_cross_moms(app, sp->kinetic, &sp->kinetic->lbo, fin);
}

static void
fluid_calc_coupled_vars(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  // Primitive variables; run after all species' self moments.
  vm_fluid_species_prim_vars(app, sp->fluid, fluidin);
}

// Ghost-current correction to dE/dt = -J in 1x: replace the current by its
// global average (cell average of J/(eps0*nx) summed over the domain).
static void
kinetic_accumulate_ghost_current(gkyl_vlasov_app *app, struct vm_species *s,
  double qbyeps, struct gkyl_array *target)
{
  if (app->field->use_ghost_current) {
    double avals_ghost_current[1], avals_ghost_current_global[1];
    // First set the scalar ghost current array to the cell average
    // current/(epsilon0*nx) where nx is the number of x cells.
    gkyl_array_set_range(app->field->ghost_current, qbyeps/app->grid.cells[0], s->m1i.marr, &app->local);
    // Integrate the current over the whole domain to find the globally averaged ghost current.
    if (app->use_gpu) {
      gkyl_array_reduce_range(app->field->red_ghost_current, app->field->ghost_current, GKYL_SUM, &app->local);
      gkyl_cu_memcpy(avals_ghost_current, app->field->red_ghost_current, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      gkyl_array_reduce_range(avals_ghost_current, app->field->ghost_current, GKYL_SUM, &app->local);
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, avals_ghost_current, avals_ghost_current_global);
    // Set the scalar ghost current array to the global average current and accumulate to the electric field.
    gkyl_array_clear(app->field->ghost_current, avals_ghost_current_global[0]);
    gkyl_array_accumulate_range(target, 1.0, app->field->ghost_current, &app->local);
  }
}

// Explicit coupling to Vlasov-Maxwell: accumulate this species' current onto
// the EM RHS ('target' is emout), J_s = q_s * m1i, scaled by -1/epsilon0.
static void
kinetic_accumulate_current(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target)
{
  struct vm_species *s = sp->kinetic;
  double qbyeps = sp->charge/app->field->info.epsilon0;

  vm_species_moment_calc(&s->m1i, s->local, app->local, fin);
  gkyl_array_accumulate_range(target, -qbyeps, s->m1i.marr, &app->local);
  kinetic_accumulate_ghost_current(app, s, qbyeps, target);
}

// Explicit coupling of a triad species to GR-Maxwell: deposit the current
// through the GR kernel, q/eps0*(rho*beta - alpha*e^i_a*Jhat^a), which needs
// the charge density as well as the current.
static void
kinetic_accumulate_current_gr(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target)
{
  struct vm_species *s = sp->kinetic;
  double qbyeps = sp->charge/app->field->info.epsilon0;

  vm_species_moment_calc(&s->m0, s->local, app->local, fin);
  vm_species_moment_calc(&s->m1i, s->local, app->local, fin);
  gkyl_dg_gr_maxwell_current_deposition_advance(s->collisionless.calc_current_dep,
    &app->local, qbyeps, app->field->geom->lapse, app->field->geom->shift,
    app->field->geom->vierb_con, s->m0.marr, s->m1i.marr, target);
  kinetic_accumulate_ghost_current(app, s, qbyeps, target);
}

// Explicit coupling to Vlasov-Poisson: accumulate this species' charge density
// onto the Poisson source ('target' is the field's rho_c), rho_s = q_s * m0.
static void
kinetic_accumulate_charge_dens(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *target)
{
  struct vm_species *s = sp->kinetic;

  vm_species_moment_calc(&s->m0, s->local, app->local, fin);
  gkyl_array_accumulate_range(target, sp->charge, s->m0.marr, &app->local);
}

// --- Constructors ------------------------------------------------------------

// Construct a species from its input: validate the declared type against the
// input blocks and dispatch to the typed constructor.
void
vlasov_species_init(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_species *inp, struct vlasov_species *sp)
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
    default:
      assert(false); // GKYL_SPECIES_PKPM is reserved and not supported here
      break;
  }
}

void
vlasov_kinetic_species_init(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_species *inp, struct vlasov_species *sp)
{
  *sp = (struct vlasov_species) { };
  sp->type = GKYL_SPECIES_VLASOV;

  strcpy(sp->name, inp->name);
  sp->charge = inp->charge;
  sp->mass = inp->mass;

  sp->kinetic = gkyl_malloc(sizeof(struct vm_species));
  *sp->kinetic = (struct vm_species) { };
  sp->kinetic->info = inp->kinetic;
  strcpy(sp->kinetic->name, inp->name);
  sp->kinetic->charge = inp->charge;
  sp->kinetic->mass = inp->mass;

  sp->calc_self_moms_func = kinetic_calc_self_moms;
  sp->calc_coupled_vars_func = kinetic_calc_coupled_vars;
  // Explicit field coupling by field type: Maxwell fields take the current
  // density (through the GR kernel for a triad species coupled to GR-Maxwell),
  // Poisson fields the charge density, the null field nothing.
  enum gkyl_field_id field_id = app->field->field_id;
  switch (field_id) {
    case GKYL_FIELD_NULL:
      sp->accumulate_field_coupling_func = species_no_field_coupling;
      break;
    case GKYL_FIELD_PHI:
      sp->accumulate_field_coupling_func = kinetic_accumulate_charge_dens;
      break;
    default:
      sp->accumulate_field_coupling_func =
        vm_species_has_gr_em_triad_coupling(app, field_id, inp->kinetic.model_id) ?
        kinetic_accumulate_current_gr : kinetic_accumulate_current;
      break;
  }
}

void
vlasov_fluid_species_init(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_species *inp, struct vlasov_species *sp)
{
  *sp = (struct vlasov_species) { };
  sp->type = GKYL_SPECIES_FLUID;

  strcpy(sp->name, inp->name);
  sp->charge = inp->charge;
  sp->mass = inp->mass;

  sp->fluid = gkyl_malloc(sizeof(struct vm_fluid_species));
  *sp->fluid = (struct vm_fluid_species) { };
  sp->fluid->info = inp->fluid;
  strcpy(sp->fluid->name, inp->name);
  sp->fluid->charge = inp->charge;
  sp->fluid->mass = inp->mass;

  sp->calc_self_moms_func = species_no_calc_moms;
  sp->calc_coupled_vars_func = fluid_calc_coupled_vars;
  sp->accumulate_field_coupling_func = species_no_field_coupling;
}

// --- Wrappers over the constructor-assigned methods --------------------------

// Self-collision moments and boundary corrections (staging phase 1).
void
vlasov_species_calc_self_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  sp->calc_self_moms_func(app, sp, fin);
}

// Coupled variables (staging phase 2): LBO cross-collision moments or fluid
// primitive variables. Run after self moments for all species.
void
vlasov_species_calc_coupled_vars(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  sp->calc_coupled_vars_func(app, sp, fin, fluidin);
}

// Accumulate this species' explicit source contribution onto the field's
// target array (Maxwell: current onto emout; Poisson: charge density onto
// rho_c; no-op for the null field and for implicitly-coupled fluid species).
void
vlasov_species_accumulate_field_coupling(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, struct gkyl_array *target)
{
  sp->accumulate_field_coupling_func(app, sp, fin, fluidin, target);
}

// --- Per-aspect operations ---------------------------------------------------

// Compute time-dependent applied acceleration (kinetic aspect only).
void
vlasov_species_calc_app_accel(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr)
{
  if (sp->kinetic && sp->kinetic->collisionless.app_accel_evolve)
    vm_species_collisionless_app_accel(app, &sp->kinetic->collisionless, tcurr);
}

// Compute the species RHS, returning the maximum stable time-step across
// aspects. Cross-aspect coupling is staged in the moments phases, so the aspect
// RHSs are independent here.
double
vlasov_species_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, const struct gkyl_array *emin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  double dt = DBL_MAX;
  if (sp->kinetic)  dt = fmin(dt, vm_species_rhs(app, sp->kinetic, fin, emin, fout));
  if (sp->fluid) dt = fmin(dt, vm_fluid_species_rhs(app, sp->fluid, fluidin, emin, fluidout));
  return dt;
}

// Pre-compute moments needed to rescale adaptive sources (kinetic aspect only).
// Done in its own pass so the species RHS (and any boundary fluxes) are complete.
void
vlasov_species_calc_source_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  if (sp->kinetic && sp->kinetic->source_id)
    vm_species_source_adapt_moms(app, sp->kinetic, &sp->kinetic->src, fin);
}

// Accumulate the source onto the RHS for each present aspect. Takes the full
// RK-state arrays since the source may read/write across species (e.g.
// boundary-flux sources).
void
vlasov_species_source_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fluidout[])
{
  if (sp->kinetic && sp->kinetic->source_id) {
    if (sp->kinetic->src.evolve_source)
      vm_species_source_calc(app, sp->kinetic, &sp->kinetic->src, tcurr);
    vm_species_source_adapt(app, sp->kinetic, &sp->kinetic->src);
    vm_species_source_rhs(app, sp->kinetic, &sp->kinetic->src, fin, fout);
  }
  if (sp->fluid && sp->fluid->source_id)
    vm_fluid_species_source_rhs(app, sp->fluid, &sp->fluid->src, fluidin, fluidout);
}

// Implicit collision update, phase 1: moments of the input distribution
// (kinetic aspect only; a no-op for species without implicit BGK collisions).
void
vlasov_species_calc_implicit_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  if (sp->kinetic)
    vm_species_bgk_moms_implicit(app, sp->kinetic, &sp->kinetic->bgk, fin);
}

// Implicit collision update, phase 2: the implicit RHS. Run after phase 1 for
// all species (the RHS also computes cross moments).
void
vlasov_species_rhs_implicit(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *fout, double dt)
{
  if (sp->kinetic)
    vm_species_rhs_implicit(app, sp->kinetic, fin, fout, dt);
}

// Implicit collision update, phase 3: BCs and copy-back into the solution.
void
vlasov_species_finish_implicit_update(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *fout, double tcurr)
{
  if (sp->kinetic) {
    vm_species_apply_bc(app, sp->kinetic, fout, tcurr);
    gkyl_array_copy_range(sp->kinetic->f, fout, &sp->kinetic->local_ext);
  }
}

// Forward-Euler accumulate: out = dt*out + in, for each present aspect.
void
vlasov_species_step_f(struct vlasov_species *sp, double dt,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  if (sp->kinetic)  vm_species_step_f(sp->kinetic, fout, dt, fin);
  if (sp->fluid) vm_fluid_species_step_f(sp->fluid, fluidout, dt, fluidin);
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
vlasov_species_combine(gkyl_vlasov_app *app, struct vlasov_species *sp, enum vm_rk_buf out,
  double c1, enum vm_rk_buf b1, double c2, enum vm_rk_buf b2)
{
  if (sp->kinetic) {
    struct vm_species *d = sp->kinetic;
    vm_species_combine(d, kinetic_rk_buf(d, out), c1, kinetic_rk_buf(d, b1),
      c2, kinetic_rk_buf(d, b2), &d->local_ext);
  }
  if (sp->fluid) {
    struct vm_fluid_species *f = sp->fluid;
    vm_fluid_species_combine(f, fluid_rk_buf(f, out), c1, fluid_rk_buf(f, b1),
      c2, fluid_rk_buf(f, b2), &app->local_ext);
  }
}

// Copy one RK buffer into another, out = inp, for each present aspect.
void
vlasov_species_copy_range(gkyl_vlasov_app *app, struct vlasov_species *sp,
  enum vm_rk_buf out, enum vm_rk_buf inp)
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
vlasov_species_gather_rk_state(gkyl_vlasov_app *app, enum vm_rk_buf in, enum vm_rk_buf out,
  const struct gkyl_array *fin[], struct gkyl_array *fout[],
  const struct gkyl_array *fluidin[], struct gkyl_array *fluidout[])
{
  for (int i=0; i<app->num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    fin[i]      = sp->kinetic ? kinetic_rk_buf(sp->kinetic, in)  : 0;
    fout[i]     = sp->kinetic ? kinetic_rk_buf(sp->kinetic, out) : 0;
    fluidin[i]  = sp->fluid   ? fluid_rk_buf(sp->fluid, in)      : 0;
    fluidout[i] = sp->fluid   ? fluid_rk_buf(sp->fluid, out)     : 0;
  }
}

// Apply boundary conditions to each present aspect.
void
vlasov_species_apply_bc(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *f, struct gkyl_array *fluid, double tcurr)
{
  if (sp->kinetic)  vm_species_apply_bc(app, sp->kinetic, f, tcurr);
  if (sp->fluid) vm_fluid_species_apply_bc(app, sp->fluid, fluid);
}

// Apply the post-update slope limiter (fluid aspect only).
void
vlasov_species_limiter(gkyl_vlasov_app *app, struct vlasov_species *sp, struct gkyl_array *fluid)
{
  if (sp->fluid) vm_fluid_species_limiter(app, sp->fluid, fluid);
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
  if (sp->kinetic) {
    vm_species_calc_integrated_mom(app, sp->kinetic, tm);
    if (sp->kinetic->src.write_source)
      vm_species_source_calc_integrated_mom(app, sp->kinetic, &sp->kinetic->src, tm);
  }
  if (sp->fluid)
    vm_fluid_species_calc_integrated_mom(app, sp->fluid, tm);
}

// Compute the integrated L2 norm of the distribution (kinetic aspect only).
void
vlasov_species_calc_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm)
{
  if (sp->kinetic)
    vm_species_calc_L2(app, sp->kinetic, tm);
}

// Write the evolved state (and attendant diagnostics: sources, emission
// spectra, fluid primitive variables) of each present aspect for this frame.
void
vlasov_species_write(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame)
{
  if (sp->kinetic) {
    struct vm_species *vms = sp->kinetic;
    vm_species_write(app, vms, tm, frame);

    if (vms->src.write_source)
      vm_species_source_write(app, vms, &vms->src, tm, frame);

    struct gkyl_msgpack_data *mt = vlasov_array_meta_new( (struct vlasov_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = vms->basis.id
      }
    );
    if (vms->emit_lo)
      vm_species_emission_write(app, vms, &vms->bc_emission_lo, mt, frame);
    if (vms->emit_up)
      vm_species_emission_write(app, vms, &vms->bc_emission_up, mt, frame);
    vlasov_array_meta_release(mt);
  }
  if (sp->fluid)
    vm_fluid_species_write(app, sp->fluid, tm, frame);
}

// Write diagnostic moments (kinetic aspect only).
void
vlasov_species_write_mom(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame)
{
  if (sp->kinetic) {
    vm_species_write_mom(app, sp->kinetic, tm, frame);
    if (sp->kinetic->src.write_source)
      vm_species_source_write_mom(app, sp->kinetic, &sp->kinetic->src, tm, frame);
  }
}

// Append integrated diagnostic moments for each present aspect.
void
vlasov_species_write_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->kinetic) {
    vm_species_write_integrated_mom(app, sp->kinetic);
    if (sp->kinetic->src.write_source)
      vm_species_source_write_integrated_mom(app, sp->kinetic, &sp->kinetic->src);
  }
  if (sp->fluid)
    vm_fluid_species_write_integrated_mom(app, sp->fluid);
}

// Append the integrated L2 norm of the distribution (kinetic aspect only).
void
vlasov_species_write_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->kinetic)
    vm_species_write_L2(app, sp->kinetic);
}

// Append the LTE-correction iteration status (kinetic aspect only).
void
vlasov_species_write_lte_corr_status(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->kinetic)
    vm_species_lte_write_max_corr_status(app, sp->kinetic);
}

// Read each present aspect's evolved state from the named file and rebuild what
// a restart does not carry: the velocity- and configuration-space Jacobian
// rescale and boundary fluxes (kinetic), BCs, sources, and static applied
// accelerations.
struct gkyl_app_restart_status
vlasov_species_from_file(gkyl_vlasov_app *app, struct vlasov_species *sp, const char *fname)
{
  struct gkyl_app_restart_status rstat = vlasov_header_from_file(app, fname);

  if (sp->kinetic) {
    struct vm_species *vms = sp->kinetic;
    if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
      rstat.io_status =
        gkyl_comm_array_read(vms->comm, &vms->grid, &vms->local, vms->f_host, fname);
      if (app->use_gpu) {
        gkyl_array_copy(vms->f, vms->f_host);
      }
      if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
        // Rescale distribution function by velocity-space Jacobian if present
        // since output distribution function does not include velocity-space Jacobian.
        // Need to do this before applying boundary conditions since we only know f on
        // the local range for the rescaling.
        gkyl_vlasov_velocity_map_rescale_jacobvel(vms->vel_map, &app->basis, &vms->basis,
          &vms->local, vms->f, vms->f_no_J);
        gkyl_array_copy(vms->f, vms->f_no_J);
        // Also rescale by the configuration-space Jacobian since the output
        // distribution does not include it either (stored f is J_x*J_v*f).
        gkyl_vlasov_position_map_rescale_jacobpos(vms->pos_map, &vms->basis,
          &vms->local, vms->f, vms->f_no_J);
        gkyl_array_copy(vms->f, vms->f_no_J);

        if (vms->calc_bflux) {
          vm_species_bflux_rhs(app, vms, &vms->bflux, vms->f, vms->f);
        }
        vm_species_apply_bc(app, vms, vms->f, rstat.stime);
        if (vms->source_id) {
          vm_species_source_calc(app, vms, &vms->src, 0.0);
        }
      }
    }
    vm_species_collisionless_app_accel(app, &vms->collisionless, rstat.stime);
  }

  if (sp->fluid) {
    struct vm_fluid_species *vm_fs = sp->fluid;
    if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
      rstat.io_status =
        gkyl_comm_array_read(app->comm, &app->grid, &app->local, vm_fs->fluid_host, fname);
      if (app->use_gpu) {
        gkyl_array_copy(vm_fs->fluid, vm_fs->fluid_host);
      }
      if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
        vm_fluid_species_apply_bc(app, vm_fs, vm_fs->fluid);
        if (vm_fs->source_id) {
          vm_fluid_species_source_calc(app, vm_fs, 0.0);
        }
      }
    }
    vm_fluid_species_calc_app_accel(app, vm_fs, rstat.stime);
  }

  return rstat;
}

// Restart a species from a frame. The kinetic aspect is seeded from the initial
// conditions first, so that fixed-function BC buffers are filled before the
// read overwrites the interior; the diagnostic dynvectors are then marked to
// append.
struct gkyl_app_restart_status
vlasov_species_read_from_frame(gkyl_vlasov_app *app, struct vlasov_species *sp, int frame)
{
  if (sp->kinetic)
    vlasov_species_apply_ic(app, sp, 0.0);

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, sp->name, frame);
  struct gkyl_app_restart_status rstat = vlasov_species_from_file(app, sp, fileNm.str);
  cstr_drop(&fileNm);

  if (sp->kinetic) {
    sp->kinetic->is_first_integ_write_call = false; // append to existing diagnostic
    sp->kinetic->is_first_integ_L2_write_call = false; // append to existing diagnostic
  }
  if (sp->fluid) {
    sp->fluid->is_first_integ_write_call = false; // append to existing diagnostic
  }

  return rstat;
}

void
vlasov_species_gather_state(gkyl_vlasov_app *app, struct gkyl_array *distf[],
  struct gkyl_array *fluid[])
{
  int num_species = app->num_species;
  for (int i=0; i<num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (distf) distf[i] = sp->kinetic ? sp->kinetic->f : 0;
    if (fluid) fluid[i] = sp->fluid ? sp->fluid->fluid : 0;
  }
}

void
vlasov_species_gather_dist(gkyl_vlasov_app *app, const struct gkyl_array *fin[])
{
  int num_species = app->num_species;
  for (int i=0; i<num_species; ++i)
    fin[i] = app->species[i].kinetic ? app->species[i].kinetic->f : 0;
}

void
vlasov_species_init_kinetic_aspects(struct gkyl_vm *vm, gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i)
    if (app->species[i].kinetic)
      vm_species_init(vm, app, app->species[i].kinetic);
}

void
vlasov_species_link_kinetic_aspects(gkyl_vlasov_app *app)
{
  // Emission walls need their impact species' boundary-flux objects.
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (!vms) continue;
    if (vms->emit_lo)
      vm_species_emission_cross_init(app, vms, &vms->bc_emission_lo);
    if (vms->emit_up)
      vm_species_emission_cross_init(app, vms, &vms->bc_emission_up);
  }

  // Cross collisions need the partners' collision objects.
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (!vms) continue;
    vm_species_lbo_cross_init(app, vms, &vms->lbo);
    vm_species_bgk_cross_init(app, vms, &vms->bgk);
  }

  // Sources may create a boundary-flux updater on their source species.
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (vms && vms->source_id)
      vm_species_source_init(app, vms, &vms->src);
  }
}

void
vlasov_species_init_fluid_aspects(struct gkyl_vm *vm, gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i)
    if (app->species[i].fluid)
      vm_fluid_species_init(vm, app, app->species[i].fluid);

  for (int i=0; i<app->num_species; ++i) {
    struct vm_fluid_species *vmf = app->species[i].fluid;
    if (vmf && vmf->source_id)
      vm_fluid_species_source_init(app, vmf, &vmf->src);
  }
}

// Release each present aspect and the container's allocations.
void
vlasov_species_release(const gkyl_vlasov_app *app, struct vlasov_species *sp)
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
