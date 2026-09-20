#include <gkyl_vlasov_priv.h>

#include <string.h>

// Vlasov species dispatch over the species container (struct vlasov_species).
// Per-aspect operations apply to whichever aspects (dist, fluid) a species owns;
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
  vm_species_lbo_moms(app, sp->dist, &sp->dist->lbo, fin);
  vm_species_bgk_moms(app, sp->dist, &sp->dist->bgk, fin);
}

static void
kinetic_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  vm_species_lbo_cross_moms(app, sp->dist, &sp->dist->lbo, fin);
}

static void
fluid_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  // Primitive moments; run in the cross-moms phase so they follow all self moments.
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
  struct vm_species *s = sp->dist;
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
  struct vm_species *s = sp->dist;
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
  struct vm_species *s = sp->dist;

  vm_species_moment_calc(&s->m0, s->local, app->local, fin);
  gkyl_array_accumulate_range(target, sp->charge, s->m0.marr, &app->local);
}

// --- Constructors ------------------------------------------------------------

// Construct a species from its input: validate the declared type against the
// input blocks and dispatch to the typed constructor.
void
vlasov_species_new(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_species *inp, struct vlasov_species *sp)
{
  switch (inp->type) {
    case GKYL_SPECIES_VLASOV:
      assert(inp->kinetic.cells[0] > 0); // kinetic species require a velocity grid
      assert(!inp->fluid.equation); // ... and must not carry a fluid equation object
      vlasov_kinetic_species_new(app, inp, sp);
      break;
    case GKYL_SPECIES_FLUID:
      assert(inp->fluid.equation); // fluid species require an equation object
      assert(inp->kinetic.cells[0] == 0); // ... and must not declare a velocity grid
      vlasov_fluid_species_new(app, inp, sp);
      break;
    default:
      assert(false); // GKYL_SPECIES_PKPM is reserved and not supported here
      break;
  }
}

void
vlasov_kinetic_species_new(struct gkyl_vlasov_app *app,
  const struct gkyl_vlasov_species *inp, struct vlasov_species *sp)
{
  *sp = (struct vlasov_species) { };
  sp->type = GKYL_SPECIES_VLASOV;

  strcpy(sp->name, inp->name);
  sp->charge = inp->charge;
  sp->mass = inp->mass;

  sp->dist = gkyl_malloc(sizeof(struct vm_species));
  *sp->dist = (struct vm_species) { };
  sp->dist->info = inp->kinetic;
  strcpy(sp->dist->name, inp->name);
  sp->dist->charge = inp->charge;
  sp->dist->mass = inp->mass;

  sp->calc_self_moms_func = kinetic_calc_self_moms;
  sp->calc_cross_moms_func = kinetic_calc_cross_moms;
  // Explicit field coupling by field type: Maxwell fields take the current
  // density (through the GR kernel for a triad species coupled to GR-Maxwell),
  // Poisson fields the charge density, the null field nothing.
  enum gkyl_field_id field_id = app->field->field_id;
  switch (field_id) {
    case GKYL_FIELD_NULL:
      sp->accumulate_field_coupling_func = species_no_field_coupling;
      break;
    case GKYL_FIELD_PHI:
    case GKYL_FIELD_PHI_EXT_POTENTIALS:
    case GKYL_FIELD_PHI_EXT_FIELDS:
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
vlasov_fluid_species_new(struct gkyl_vlasov_app *app,
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
  sp->calc_cross_moms_func = fluid_calc_cross_moms;
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

// Cross-collision moments / fluid primitive variables (staging phase 2).
// Run after self moments for all species, so cross moments see them.
void
vlasov_species_calc_cross_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin)
{
  sp->calc_cross_moms_func(app, sp, fin, fluidin);
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
  if (sp->dist && sp->dist->collisionless.app_accel_evolve)
    vm_species_collisionless_app_accel(app, &sp->dist->collisionless, tcurr);
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
  if (sp->dist)  dt = fmin(dt, vm_species_rhs(app, sp->dist, fin, emin, fout));
  if (sp->fluid) dt = fmin(dt, vm_fluid_species_rhs(app, sp->fluid, fluidin, emin, fluidout));
  return dt;
}

// Pre-compute moments needed to rescale adaptive sources (kinetic aspect only).
// Done in its own pass so the species RHS (and any boundary fluxes) are complete.
void
vlasov_species_calc_source_moms(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin)
{
  if (sp->dist && sp->dist->source_id)
    vm_species_source_adapt_moms(app, sp->dist, &sp->dist->src, fin);
}

// Accumulate the source onto the RHS for each present aspect. Takes the full
// RK-state arrays since the source may read/write across species (e.g.
// boundary-flux sources).
void
vlasov_species_source_rhs(gkyl_vlasov_app *app, struct vlasov_species *sp, double tcurr,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fluidout[])
{
  if (sp->dist && sp->dist->source_id) {
    if (sp->dist->src.evolve_source)
      vm_species_source_calc(app, sp->dist, &sp->dist->src, tcurr);
    vm_species_source_adapt(app, sp->dist, &sp->dist->src);
    vm_species_source_rhs(app, sp->dist, &sp->dist->src, fin, fout);
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
  if (sp->dist)
    vm_species_bgk_moms_implicit(app, sp->dist, &sp->dist->bgk, fin);
}

// Implicit collision update, phase 2: the implicit RHS. Run after phase 1 for
// all species (the RHS also computes cross moments).
void
vlasov_species_rhs_implicit(gkyl_vlasov_app *app, struct vlasov_species *sp,
  const struct gkyl_array *fin, struct gkyl_array *fout, double dt)
{
  if (sp->dist)
    vm_species_rhs_implicit(app, sp->dist, fin, fout, dt);
}

// Implicit collision update, phase 3: BCs and copy-back into the solution.
void
vlasov_species_finish_implicit_update(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *fout, double tcurr)
{
  if (sp->dist) {
    vm_species_apply_bc(app, sp->dist, fout, tcurr);
    gkyl_array_copy_range(sp->dist->f, fout, &sp->dist->local_ext);
  }
}

// Forward-Euler accumulate: out = dt*out + in, for each present aspect.
void
vlasov_species_step_f(struct vlasov_species *sp, double dt,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin,
  struct gkyl_array *fout, struct gkyl_array *fluidout)
{
  if (sp->dist)  vm_species_step_f(sp->dist, fout, dt, fin);
  if (sp->fluid) vm_fluid_species_step_f(sp->fluid, fluidout, dt, fluidin);
}

// Combine RK stages into the first-stage buffer: f1 = c1*f + c2*fnew, for each
// present aspect (the RK3 stepper's combine pattern).
void
vlasov_species_combine(gkyl_vlasov_app *app, struct vlasov_species *sp, double c1, double c2)
{
  if (sp->dist) {
    struct vm_species *d = sp->dist;
    vm_species_combine(d, d->f1, c1, d->f, c2, d->fnew, &d->local_ext);
  }
  if (sp->fluid) {
    struct vm_fluid_species *f = sp->fluid;
    vm_fluid_species_combine(f, f->fluid1, c1, f->fluid, c2, f->fluidnew, &app->local_ext);
  }
}

// Copy the first-stage buffer back into the solution: f = f1, for each aspect.
void
vlasov_species_copy_range(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->dist) {
    struct vm_species *d = sp->dist;
    vm_species_copy_range(d, d->f, d->f1, &d->local_ext);
  }
  if (sp->fluid) {
    struct vm_fluid_species *f = sp->fluid;
    vm_fluid_species_copy_range(f, f->fluid, f->fluid1, &app->local_ext);
  }
}

// Apply boundary conditions to each present aspect.
void
vlasov_species_apply_bc(gkyl_vlasov_app *app, struct vlasov_species *sp,
  struct gkyl_array *f, struct gkyl_array *fluid, double tcurr)
{
  if (sp->dist)  vm_species_apply_bc(app, sp->dist, f, tcurr);
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
  if (sp->dist) {
    struct timespec wtm = gkyl_wall_clock();
    vm_species_apply_ic(app, sp->dist, t0);
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
  if (sp->dist) {
    vm_species_calc_integrated_mom(app, sp->dist, tm);
    if (sp->dist->src.write_source)
      vm_species_source_calc_integrated_mom(app, sp->dist, &sp->dist->src, tm);
  }
  if (sp->fluid)
    vm_fluid_species_calc_integrated_mom(app, sp->fluid, tm);
}

// Compute the integrated L2 norm of the distribution (kinetic aspect only).
void
vlasov_species_calc_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm)
{
  if (sp->dist)
    vm_species_calc_L2(app, sp->dist, tm);
}

// Write the evolved state (and attendant diagnostics: sources, emission
// spectra, fluid primitive variables) of each present aspect for this frame.
void
vlasov_species_write(gkyl_vlasov_app *app, struct vlasov_species *sp, double tm, int frame)
{
  if (sp->dist) {
    struct vm_species *vms = sp->dist;
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
  if (sp->dist) {
    vm_species_write_mom(app, sp->dist, tm, frame);
    if (sp->dist->src.write_source)
      vm_species_source_write_mom(app, sp->dist, &sp->dist->src, tm, frame);
  }
}

// Append integrated diagnostic moments for each present aspect.
void
vlasov_species_write_integrated_mom(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->dist) {
    vm_species_write_integrated_mom(app, sp->dist);
    if (sp->dist->src.write_source)
      vm_species_source_write_integrated_mom(app, sp->dist, &sp->dist->src);
  }
  if (sp->fluid)
    vm_fluid_species_write_integrated_mom(app, sp->fluid);
}

// Append the integrated L2 norm of the distribution (kinetic aspect only).
void
vlasov_species_write_integrated_L2_f(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->dist)
    vm_species_write_L2(app, sp->dist);
}

// Append the LTE-correction iteration status (kinetic aspect only).
void
vlasov_species_write_lte_corr_status(gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->dist)
    vm_species_lte_write_max_corr_status(app, sp->dist);
}

// Read each present aspect's evolved state from the named file and rebuild what
// a restart does not carry: the velocity- and configuration-space Jacobian
// rescale and boundary fluxes (kinetic), BCs, sources, and static applied
// accelerations.
struct gkyl_app_restart_status
vlasov_species_from_file(gkyl_vlasov_app *app, struct vlasov_species *sp, const char *fname)
{
  struct gkyl_app_restart_status rstat = vlasov_header_from_file(app, fname);

  if (sp->dist) {
    struct vm_species *vms = sp->dist;
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
  if (sp->dist)
    vlasov_species_apply_ic(app, sp, 0.0);

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, sp->name, frame);
  struct gkyl_app_restart_status rstat = vlasov_species_from_file(app, sp, fileNm.str);
  cstr_drop(&fileNm);

  if (sp->dist) {
    sp->dist->is_first_integ_write_call = false; // append to existing diagnostic
    sp->dist->is_first_integ_L2_write_call = false; // append to existing diagnostic
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
  int num_species = app->num_species + app->num_fluid_species;
  for (int i=0; i<num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (distf) distf[i] = sp->dist ? sp->dist->f : 0;
    if (fluid) fluid[i] = sp->fluid ? sp->fluid->fluid : 0;
  }
}

void
vlasov_species_gather_dist(gkyl_vlasov_app *app, const struct gkyl_array *fin[])
{
  int num_species = app->num_species + app->num_fluid_species;
  for (int i=0; i<num_species; ++i)
    fin[i] = app->species[i].dist ? app->species[i].dist->f : 0;
}

// Release each present aspect and the aspect allocations themselves.
void
vlasov_species_release(const gkyl_vlasov_app *app, struct vlasov_species *sp)
{
  if (sp->dist) {
    vm_species_release(app, sp->dist);
    gkyl_free(sp->dist);
  }
  if (sp->fluid) {
    vm_fluid_species_release(app, sp->fluid);
    gkyl_free(sp->fluid);
  }
}
