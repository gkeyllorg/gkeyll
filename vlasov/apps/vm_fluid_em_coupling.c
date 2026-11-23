#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

bool
vm_fluid_em_coupling_supported(const struct vlasov_species *sp)
{
  // The fluid aspects that carry a momentum: that of a PKPM species, and the Euler and
  // isothermal Euler equations of a fluid species.
  if (sp->type == GKYL_SPECIES_PKPM) {
    return true;
  }
  return sp->type == GKYL_SPECIES_FLUID && (sp->fluid->info.equation->type == GKYL_EQN_EULER ||
                                            sp->fluid->info.equation->type == GKYL_EQN_ISO_EULER);
}

static double
fluid_em_coupling_dt_none(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  return DBL_MAX;
}

static double
fluid_em_coupling_dt_kinetic(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  gkyl_array_clear(fl_em->omega_p_sq, 0.0);
  for (int i = 0; i < fl_em->num_kinetic; ++i) {
    struct vm_species *vms = fl_em->kinetic[i]->kinetic;
    vm_species_accumulate_omega_p_sq(app, vms, vms->f, fl_em->m0_avg, fl_em->omega_p_sq);
  }
  gkyl_array_reduce_range(fl_em->omega_p_sq_max, fl_em->omega_p_sq, GKYL_MAX, &app->local);

  double omega_p_sq_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_p_sq_ho, fl_em->omega_p_sq_max, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_p_sq_ho[0] = fl_em->omega_p_sq_max[0];
  }
  double omega_p_sq_local = omega_p_sq_ho[0], omega_p_sq_global;
  gkyl_comm_allreduce_host(
    app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &omega_p_sq_local, &omega_p_sq_global
  );
  return omega_p_sq_global > 0.0 ? app->cfl / sqrt(omega_p_sq_global) : DBL_MAX;
}

// The solve of a source step. It advances the momentum of the fluids and the electric field (with
// the resistive layer, if any), and the energy of a fluid that has one by the change of the bulk
// kinetic energy (the work of the field): the density is unchanged, and so is the internal energy
// at the nodes of the solve. It acts on the physical fluids and field (without the position-map
// Jacobian).
static void
fluid_em_coupling_advance(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double dt,
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES], struct gkyl_array *fluids[GKYL_MAX_SPECIES]
)
{
  gkyl_dg_calc_fluid_em_coupling_advance(
    fl_em->slvr, dt, app_accels, app->field->ext_em, fl_em->curl_current, fl_em->sigma, fluids,
    fl_em->em_solve
  );
}

// Domain mean of E_x of the field of the solve. On a J-weighted field the mean is taken of the
// field with the position-map Jacobian, whose cell averages sum to the integral over the
// physical domain; the evolved field is used as scratch for it (the update rescales into it at
// the end).
static double
fluid_em_coupling_ex_mean(struct gkyl_vlasov_app *app)
{
  struct vm_field *field = app->field;
  if (field->weight_by_pos_jacob) {
    gkyl_vlasov_position_map_rescale_jacobpos_conf(
      app->pos_map, &app->local, field->em_no_J, field->em
    );
  }
  return vm_field_ghost_current_mean(app, field, 1.0, field->em);
}

// The solve with the ghost current (1x periodic): a uniform current along x whose strength keeps
// the average of E_x at its value at the start of the step, so that the average current of the
// fluids over the step drives no uniform field (the kinetic species remove their own average
// when they deposit their current). The solve is linear in the current at a fixed initial state:
// two trial solves, without and with a unit ghost current, give the strength, and a third takes
// the step (the energy of an Euler fluid is not linear in the current).
static void
fluid_em_coupling_advance_ghost(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double dt,
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES], struct gkyl_array *fluids[GKYL_MAX_SPECIES]
)
{
  struct vm_field *field = app->field;
  int num_fluid = fl_em->num_fluid;

  for (int i = 0; i < num_fluid; ++i) {
    gkyl_array_copy(fl_em->fluid_ghost[i], fluids[i]);
  }
  gkyl_array_copy(fl_em->em_ghost, fl_em->em_solve);
  double ex_mean_start = vm_field_ghost_current_mean(app, field, 1.0, field->em);

  double ghost_current[3] = {0.0, 1.0, 0.0}; // the trial strengths; the last is the solution
  double ex_mean_trial[2];
  double ghost_current_prev = 0.0;
  for (int trial = 0; trial < 3; ++trial) {
    if (trial == 2) {
      ghost_current[2] = (ex_mean_start - ex_mean_trial[0]) / (ex_mean_trial[1] - ex_mean_trial[0]);
    }
    gkyl_array_shiftc(fl_em->curl_current, ghost_current[trial] - ghost_current_prev, 0);
    ghost_current_prev = ghost_current[trial];
    for (int i = 0; i < num_fluid; ++i) {
      gkyl_array_copy(fluids[i], fl_em->fluid_ghost[i]);
    }
    gkyl_array_copy(fl_em->em_solve, fl_em->em_ghost);
    fluid_em_coupling_advance(app, fl_em, dt, app_accels, fluids);
    if (trial < 2) {
      ex_mean_trial[trial] = fluid_em_coupling_ex_mean(app);
    }
  }
}

// initialize fluid-EM coupling object
struct vm_fluid_em_coupling *
vm_fluid_em_coupling_new(struct gkyl_vlasov_app *app)
{
  struct vm_fluid_em_coupling *fl_em = gkyl_malloc(sizeof(struct vm_fluid_em_coupling));

  // Gather the species whose fluid aspect carries a momentum (fluid and PKPM species) and the
  // kinetic species in declaration order.
  int num_species = app->num_species;
  bool has_energy[GKYL_MAX_SPECIES]; // the fluid state carries an energy (not isothermal Euler)
  fl_em->num_fluid = 0;
  fl_em->num_kinetic = 0;
  for (int i = 0; i < num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (vm_fluid_em_coupling_supported(sp)) {
      fl_em->species[fl_em->num_fluid] = sp;
      fl_em->qbym[fl_em->num_fluid] = sp->charge / sp->mass;
      has_energy[fl_em->num_fluid] = sp->type == GKYL_SPECIES_PKPM ||
                                     sp->fluid->eqn_type == GKYL_EQN_EULER;
      fl_em->num_fluid += 1;
    }
    if (sp->type == GKYL_SPECIES_VLASOV) {
      fl_em->kinetic[fl_em->num_kinetic] = sp;
      fl_em->num_kinetic += 1;
    }
  }

  // The source update is solved at the Gauss-Legendre nodes of the p = 1 basis.
  if (app->poly_order != 1) {
    gkyl_exit("vm_fluid_em_coupling: fluid-EM coupling requires polynomial order 1.");
  }

  // Initialize solver
  fl_em->slvr = gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
    .conf_basis = &app->basis,
    .mem_range = &app->local,
    .num_species = fl_em->num_fluid,
    .qbym = fl_em->qbym,
    .has_energy = has_energy,
    .epsilon0 = app->field->info.epsilon0,
    .is_field_static = app->field->info.is_static,
    .use_gpu = app->use_gpu,
  });

  int num_basis = app->basis.num_basis;
  long num_cells = app->local_ext.volume;
  fl_em->curl_rhs = mkarr(app->use_gpu, 8 * num_basis, num_cells);
  fl_em->curl_current = mkarr(app->use_gpu, 3 * num_basis, num_cells);
  gkyl_array_clear(fl_em->curl_current, 0.0);
  // The solve works with the physical field: em_no_J of a field that carries the position-map
  // Jacobian (divided out before the solve and restored after it), the evolved field otherwise.
  fl_em->em_solve = app->field->em;
  fl_em->curl_rhs_no_J = 0;
  if (app->field->weight_by_pos_jacob) {
    fl_em->em_solve = app->field->em_no_J;
    fl_em->curl_rhs_no_J = mkarr(app->use_gpu, 8 * num_basis, num_cells);
  }
  // The resistive layer damps E in the solve; a static field has none.
  fl_em->sigma = app->field->has_sigma && !app->field->info.is_static ? app->field->sigma : 0;
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    const struct gkyl_array *fluid = fl_em->species[i]->fluid->fluid;
    fl_em->fluid_stag[i] = mkarr(app->use_gpu, fluid->ncomp, fluid->size);
  }
  fl_em->em_stag = mkarr(app->use_gpu, app->field->em->ncomp, app->field->em->size);

  // Ghost current: the fluids' average current over a source step is removed as well (a static
  // field takes no current).
  fl_em->advance_func = fluid_em_coupling_advance;
  fl_em->em_ghost = 0;
  if (app->field->use_ghost_current && !app->field->info.is_static) {
    fl_em->advance_func = fluid_em_coupling_advance_ghost;
    for (int i = 0; i < fl_em->num_fluid; ++i) {
      const struct gkyl_array *fluid = fl_em->species[i]->fluid->fluid;
      fl_em->fluid_ghost[i] = mkarr(app->use_gpu, fluid->ncomp, fluid->size);
    }
    fl_em->em_ghost = mkarr(app->use_gpu, app->field->em->ncomp, app->field->em->size);
  }

  // A static field takes no current, so its kinetic species set no step.
  fl_em->omega_p_sq = 0;
  fl_em->m0_avg = 0;
  fl_em->omega_p_sq_max = 0;
  fl_em->src_dt_func = fluid_em_coupling_dt_none;
  if (fl_em->num_kinetic > 0 && !app->field->info.is_static) {
    fl_em->omega_p_sq = mkarr(app->use_gpu, 1, num_cells);
    fl_em->m0_avg = mkarr(app->use_gpu, 1, num_cells);
    fl_em->omega_p_sq_max = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) :
                                           gkyl_malloc(sizeof(double));
    fl_em->src_dt_func = fluid_em_coupling_dt_kinetic;
  }

  return fl_em;
}

void
vm_fluid_em_coupling_update(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double tcurr, double dt
)
{
  int num_fluid = fl_em->num_fluid;

  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];

  // The solve acts on the physical state: the fluid of a PKPM species and a field with PKPM
  // species carry the position-map Jacobian, which is divided out before and restored after
  // (identity map: copies).
  for (int i = 0; i < num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    fluids[i] = fs->fluid;
    if (fl_em->species[i]->type == GKYL_SPECIES_PKPM) {
      gkyl_vlasov_position_map_divide_jacobpos_conf(
        app->pos_map, &app->local, fs->fluid->ncomp, fs->fluid, fs->pkpm.fluid_no_J
      );
      fluids[i] = fs->pkpm.fluid_no_J;
    }
    app_accels[i] = fs->app_accel;
    // Compute applied accelerations if present and time-dependent.
    // Note: applied accelerations use proj_on_basis
    // so does copy to GPU every call if app->use_gpu = true.
    if (fs->app_accel_evolve) {
      vm_fluid_species_calc_app_accel(app, fs, tcurr);
    }
  }
  // Compute external EM field or applied currents if present and time-dependent
  // (this object only exists when there are fluid species and a field).
  // Note: external EM field and applied currents use proj_on_basis
  // so does copy to GPU every call if app->use_gpu = true.
  if (app->field->ext_em_evolve) {
    vlasov_field_calc_ext_em(app, tcurr);
  }
  if (app->field->app_current_evolve) {
    vlasov_field_calc_app_current(app, tcurr);
  }

  // The currents of the solve besides the fluids': the applied current, the effective current
  // of the curl and the current of the kinetic species. A static field takes none of them
  // (curl_current stays zero), so its E does not change.
  struct vm_field *field = app->field;
  if (field->weight_by_pos_jacob) {
    gkyl_vlasov_position_map_divide_jacobpos_conf(
      app->pos_map, &app->local, field->em->ncomp, field->em, field->em_no_J
    );
  }
  if (!field->info.is_static) {
    // Stormer-Verlet in the curl terms: E advances only here, with c^2 curl B of the field at
    // the start of the step (B is not changed by the sources) entering as the effective current
    // -eps0 c^2 curl B; the transport step then drifts B with this E.
    gkyl_array_clear(fl_em->curl_rhs, 0.0);
    gkyl_hyper_dg_advance(field->slvr, &app->local, field->em, field->cflrate, fl_em->curl_rhs);
    // Leapfrog in the kinetic currents: each kinetic species deposits -q/eps0 M1 of its
    // distribution at the start of the step on the E rows, as it does on the stage right-hand
    // side of an explicit update; the transport step then moves it with the E of this step.
    for (int i = 0; i < fl_em->num_kinetic; ++i) {
      struct vlasov_species *sp = fl_em->kinetic[i];
      vlasov_species_accumulate_field_coupling(app, sp, sp->kinetic->f, 0, fl_em->curl_rhs);
    }
    // the physical currents of the solve
    const struct gkyl_array *curl_rhs = fl_em->curl_rhs, *app_current = field->app_current;
    if (field->weight_by_pos_jacob) {
      gkyl_vlasov_position_map_divide_jacobpos_conf(
        app->pos_map, &app->local, fl_em->curl_rhs->ncomp, fl_em->curl_rhs, fl_em->curl_rhs_no_J
      );
      curl_rhs = fl_em->curl_rhs_no_J;
      app_current = field->app_current_no_J;
    }
    gkyl_array_set(fl_em->curl_current, 1.0, app_current);
    gkyl_array_accumulate_offset(fl_em->curl_current, -field->info.epsilon0, curl_rhs, 0);
  }

  fl_em->advance_func(app, fl_em, dt, app_accels, fluids);

  // Restore the Jacobian and apply BCs to the updated fluids and field.
  for (int i = 0; i < num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    if (fl_em->species[i]->type == GKYL_SPECIES_PKPM) {
      gkyl_vlasov_position_map_rescale_jacobpos_conf(
        app->pos_map, &app->local, fs->pkpm.fluid_no_J, fs->fluid
      );
    }
    vm_fluid_species_apply_bc(app, fs, fs->fluid);
  }
  if (field->weight_by_pos_jacob) {
    gkyl_vlasov_position_map_rescale_jacobpos_conf(
      app->pos_map, &app->local, field->em_no_J, field->em
    );
  }
  vlasov_field_apply_bc(app, app->field->em);
}

double
vm_fluid_em_coupling_dt(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  return fl_em->src_dt_func(app, fl_em);
}

void
vm_fluid_em_coupling_limit(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  // The positivity limiter keeps the internal energy positive at the vertices of the fluids, and
  // a limited field (limit_em) is limited: E enters the transport step as the stages of an
  // explicit update leave it.
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    vm_fluid_species_apply_positivity(app, fs, fs->fluid);
    vm_fluid_species_apply_bc(app, fs, fs->fluid);
  }
  vlasov_field_limiter(app, app->field->em);
}

void
vm_fluid_em_coupling_release(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  gkyl_dg_calc_fluid_em_coupling_release(fl_em->slvr);
  gkyl_array_release(fl_em->curl_rhs);
  if (fl_em->curl_rhs_no_J) {
    gkyl_array_release(fl_em->curl_rhs_no_J);
  }
  gkyl_array_release(fl_em->curl_current);
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    gkyl_array_release(fl_em->fluid_stag[i]);
  }
  gkyl_array_release(fl_em->em_stag);
  if (fl_em->em_ghost) {
    for (int i = 0; i < fl_em->num_fluid; ++i) {
      gkyl_array_release(fl_em->fluid_ghost[i]);
    }
    gkyl_array_release(fl_em->em_ghost);
  }
  if (fl_em->omega_p_sq) {
    gkyl_array_release(fl_em->omega_p_sq);
    gkyl_array_release(fl_em->m0_avg);
    if (app->use_gpu) {
      gkyl_cu_free(fl_em->omega_p_sq_max);
    } else {
      gkyl_free(fl_em->omega_p_sq_max);
    }
  }
  gkyl_free(fl_em);
}
