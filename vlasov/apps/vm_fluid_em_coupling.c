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

// The step bound of the densities that couple explicitly to the field, cfl / omega with omega^2
// the largest sum over the rate list of the weight times the cell average of the physical
// density at the start of the update. Maxwell split: the current of a kinetic species at the
// start of a source step is a known current of the step, so the step resolves its plasma
// frequency. Poisson field: every species couples explicitly to the potentials (a fluid aspect
// through the kick, a kinetic species in the stages of the transport step), a leapfrog stable
// for omega dt < 2 at a fixed step (derived for the plasma oscillation of a cold plasma at every
// angle to B and every Omega dt; the Jeans rate of self-gravity is resolved for accuracy).
static double
fluid_em_coupling_dt_rate(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  gkyl_array_clear(fl_em->omega_sq, 0.0);
  const struct vlasov_rate_list *rl = &fl_em->rates;
  for (int k = 0; k < rl->num; ++k) {
    struct vlasov_species *sp = &app->species[rl->idx[k]];
    const struct gkyl_array *fin = sp->kinetic ? sp->kinetic->f : 0;
    const struct gkyl_array *fluidin = sp->fluid ? sp->fluid->fluid : 0;
    vlasov_species_accumulate_omega_sq(
      app, sp, fin, fluidin, rl->weight[k], fl_em->dens_avg, fl_em->omega_sq
    );
  }
  gkyl_array_reduce_range(fl_em->omega_sq_max, fl_em->omega_sq, GKYL_MAX, &app->local);

  double omega_sq_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_sq_ho, fl_em->omega_sq_max, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_sq_ho[0] = fl_em->omega_sq_max[0];
  }
  double omega_sq_local = omega_sq_ho[0], omega_sq_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &omega_sq_local, &omega_sq_global);
  return omega_sq_global > 0.0 ? app->cfl / sqrt(omega_sq_global) : DBL_MAX;
}

// The arrays of the step bound and its function, once the rate lists are set.
static void
fluid_em_coupling_rate_init(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  fl_em->omega_sq = 0;
  fl_em->dens_avg = 0;
  fl_em->omega_sq_max = 0;
  fl_em->src_dt_func = fluid_em_coupling_dt_none;
  if (fl_em->rates.num > 0) {
    fl_em->omega_sq = mkarr(app->use_gpu, 1, app->local_ext.volume);
    fl_em->dens_avg = mkarr(app->use_gpu, 1, app->local_ext.volume);
    fl_em->omega_sq_max = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) :
                                         gkyl_malloc(sizeof(double));
    fl_em->src_dt_func = fluid_em_coupling_dt_rate;
  }
}

// The solve of a source step. It advances the momentum of the fluids and the electric field (with
// the resistive layer, if any), and the energy of a fluid that has one by the change of the bulk
// kinetic energy (the work of the field): the density is unchanged, and so is the internal energy
// at the nodes of the solve. It acts on the physical fluids and field (without the position-map
// Jacobian).
static void
fluid_em_coupling_advance(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double dt,
  struct gkyl_array *fluids[GKYL_MAX_SPECIES]
)
{
  const struct gkyl_array *accels[GKYL_MAX_SPECIES];
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    accels[i] = fl_em->accel_zero;
  }
  gkyl_dg_calc_fluid_em_coupling_advance(
    fl_em->slvr, dt, accels, app->field->ext_em, fl_em->curl_current, fl_em->sigma, fluids,
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
  struct gkyl_array *fluids[GKYL_MAX_SPECIES]
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
    fluid_em_coupling_advance(app, fl_em, dt, fluids);
    if (trial < 2) {
      ex_mean_trial[trial] = fluid_em_coupling_ex_mean(app);
    }
  }
}

// The field part of the staggered state: the Maxwell field (E advances in the source step), none
// for the Poisson field (its potentials follow the densities, which the source step leaves alone).
static void
fluid_em_coupling_save_field_maxwell(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  gkyl_array_copy(fl_em->em_stag, app->field->em);
}

static void
fluid_em_coupling_restore_field_maxwell(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em
)
{
  gkyl_array_copy(app->field->em, fl_em->em_stag);
}

static void
fluid_em_coupling_field_none(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
}

// The physical fluid states of the solve: the fluid of a PKPM species carries the position-map
// Jacobian, which is divided out (identity map: a copy) and restored by
// fluid_em_coupling_restore_fluids.
static void
fluid_em_coupling_physical_fluids(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em,
  struct gkyl_array *fluids[GKYL_MAX_SPECIES]
)
{
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    fluids[i] = fs->fluid;
    if (fl_em->species[i]->type == GKYL_SPECIES_PKPM) {
      gkyl_vlasov_position_map_divide_jacobpos_conf(
        app->pos_map, &app->local, fs->fluid->ncomp, fs->fluid, fs->pkpm.fluid_no_J
      );
      fluids[i] = fs->pkpm.fluid_no_J;
    }
  }
}

static void
fluid_em_coupling_restore_fluids(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    if (fl_em->species[i]->type == GKYL_SPECIES_PKPM) {
      gkyl_vlasov_position_map_rescale_jacobpos_conf(
        app->pos_map, &app->local, fs->pkpm.fluid_no_J, fs->fluid
      );
    }
    vm_fluid_species_apply_bc(app, fs, fs->fluid);
  }
}

// The source step of the Maxwell field.
static void
fluid_em_coupling_update_maxwell(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double tcurr, double dt
)
{
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];

  // The solve acts on the physical state: the fluid of a PKPM species and a field with PKPM
  // species carry the position-map Jacobian, which is divided out before and restored after
  // (identity map: copies).
  fluid_em_coupling_physical_fluids(app, fl_em, fluids);

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

  fl_em->advance_func(app, fl_em, dt, fluids);

  // Restore the Jacobian and apply BCs to the updated fluids and field.
  fluid_em_coupling_restore_fluids(app, fl_em);
  if (field->weight_by_pos_jacob) {
    gkyl_vlasov_position_map_rescale_jacobpos_conf(
      app->pos_map, &app->local, field->em_no_J, field->em
    );
  }
  vlasov_field_apply_bc(app, app->field->em);
}

// The kick of the Poisson field: the fluids take the forces of the potentials of the present
// densities, which the kick does not change (so the potentials stand still over it), with the
// external field: the solve with a static field, whose E is -grad(phi + phi_ext), and the
// gravitational acceleration of every fluid. It advances the momenta, implicit in the external
// magnetic field, and the energies by the work of the forces; the electric field is unchanged.
static void
fluid_em_coupling_update_poisson(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double tcurr, double dt
)
{
  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  fluid_em_coupling_physical_fluids(app, fl_em, fluids);

  // Time-dependent external fields and potentials at the midpoint of the source interval.
  if (app->field->ext_em_evolve) {
    vlasov_field_calc_ext_em(app, tcurr);
  }
  if (app->field->ext_pot_evolve) {
    vlasov_field_calc_ext_pot(app, tcurr);
  }
  vp_field_calc_kick(app, app->field);

  gkyl_dg_calc_fluid_em_coupling_advance(
    fl_em->slvr, dt, fl_em->accel_kick, app->field->ext_em, fl_em->curl_current, 0, fluids,
    app->field->em_pot
  );

  fluid_em_coupling_restore_fluids(app, fl_em);
}

// Gather the species whose fluid aspect carries a momentum (fluid and PKPM species) and the
// kinetic species in declaration order, and whether each fluid state carries an energy.
static void
fluid_em_coupling_gather(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, bool has_energy[GKYL_MAX_SPECIES]
)
{
  fl_em->num_fluid = 0;
  fl_em->num_kinetic = 0;
  for (int i = 0; i < app->num_species; ++i) {
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
}

// The coupling to the Poisson field: the kick of the staggered split (see vlasov_update_op_split).
static void
fluid_em_coupling_init_poisson(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em,
  const bool has_energy[GKYL_MAX_SPECIES]
)
{
  // The solve with a static field: the E rows keep E, so their scale is immaterial (a run with
  // self-gravity alone has no permittivity).
  double epsilon0 = app->field->info.epsilon0 > 0.0 ? app->field->info.epsilon0 : 1.0;
  fl_em->slvr = gkyl_dg_calc_fluid_em_coupling_inew(&(struct gkyl_dg_calc_fluid_em_coupling_inp){
    .conf_basis = &app->basis,
    .mem_range = &app->local,
    .num_species = fl_em->num_fluid,
    .qbym = fl_em->qbym,
    .has_energy = has_energy,
    .epsilon0 = epsilon0,
    .is_field_static = true,
    .use_gpu = app->use_gpu,
  });

  int num_basis = app->basis.num_basis;
  long num_cells = app->local_ext.volume;
  // No currents: the applied current of the solve is zero.
  fl_em->curl_rhs = 0;
  fl_em->curl_rhs_no_J = 0;
  fl_em->curl_current = mkarr(app->use_gpu, 3 * num_basis, num_cells);
  gkyl_array_clear(fl_em->curl_current, 0.0);
  fl_em->em_solve = 0;
  fl_em->sigma = 0;
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    const struct gkyl_array *fluid = fl_em->species[i]->fluid->fluid;
    fl_em->fluid_stag[i] = mkarr(app->use_gpu, fluid->ncomp, fluid->size);
  }
  fl_em->em_stag = 0;
  fl_em->accel_zero = mkarr(app->use_gpu, 3 * num_basis, num_cells);
  gkyl_array_clear(fl_em->accel_zero, 0.0);
  // Every fluid takes the gravitational acceleration of the field (zero without self-gravity).
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    fl_em->accel_kick[i] = app->field->accel_grav;
  }
  fl_em->advance_func = 0;
  fl_em->em_ghost = 0;

  // The species whose density sets a rate of the potentials, in declaration order.
  struct vlasov_rate_list *rl = &fl_em->rates;
  rl->num = 0;
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    double rate_sq = 0.0;
    for (int p = 0; p < app->field->num_pots && sp->mass > 0.0; ++p) {
      const struct vp_potential *pot = &app->field->pots[p];
      rate_sq += fabs(pot->rate_sq_weight(pot, sp->charge, sp->mass));
    }
    if (!(rate_sq > 0.0)) {
      continue;
    }
    if (sp->type == GKYL_SPECIES_VLASOV) {
      rl->idx[rl->num] = i;
      rl->weight[rl->num] = rate_sq;
      rl->num += 1;
    } else if (vm_fluid_em_coupling_supported(sp)) {
      rl->idx[rl->num] = i;
      rl->weight[rl->num] = rate_sq / sp->mass; // per unit mass density
      rl->num += 1;
    }
  }
  fluid_em_coupling_rate_init(app, fl_em);

  fl_em->update_func = fluid_em_coupling_update_poisson;
  fl_em->save_field_func = fluid_em_coupling_field_none;
  fl_em->restore_field_func = fluid_em_coupling_field_none;
}

// initialize fluid-EM coupling object
struct vm_fluid_em_coupling *
vm_fluid_em_coupling_new(struct gkyl_vlasov_app *app)
{
  struct vm_fluid_em_coupling *fl_em = gkyl_malloc(sizeof(struct vm_fluid_em_coupling));

  bool has_energy[GKYL_MAX_SPECIES]; // the fluid state carries an energy (not isothermal Euler)
  fluid_em_coupling_gather(app, fl_em, has_energy);

  // The source update is solved at the Gauss-Legendre nodes of the p = 1 basis.
  if (app->poly_order != 1) {
    gkyl_exit("vm_fluid_em_coupling: fluid-EM coupling requires polynomial order 1.");
  }
  if (app->field->field_id == GKYL_FIELD_PHI) {
    fluid_em_coupling_init_poisson(app, fl_em, has_energy);
    return fl_em;
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
  fl_em->accel_zero = mkarr(app->use_gpu, 3 * num_basis, num_cells);
  gkyl_array_clear(fl_em->accel_zero, 0.0);

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

  // The plasma frequency of the charged kinetic species sets the step (the fluids are implicit in
  // the source step); a static field takes no current, so they set none.
  struct vlasov_rate_list *rl = &fl_em->rates;
  rl->num = 0;
  for (int i = 0; i < app->num_species && !app->field->info.is_static; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (sp->type == GKYL_SPECIES_VLASOV && sp->charge != 0.0) {
      rl->idx[rl->num] = i;
      rl->weight[rl->num] = sp->charge * sp->charge / (app->field->info.epsilon0 * sp->mass);
      rl->num += 1;
    }
  }
  fluid_em_coupling_rate_init(app, fl_em);

  fl_em->update_func = fluid_em_coupling_update_maxwell;
  fl_em->save_field_func = fluid_em_coupling_save_field_maxwell;
  fl_em->restore_field_func = fluid_em_coupling_restore_field_maxwell;

  return fl_em;
}

void
vm_fluid_em_coupling_update(
  struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em, double tcurr, double dt
)
{
  fl_em->update_func(app, fl_em, tcurr, dt);
}

void
vm_fluid_em_coupling_save_field(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  fl_em->save_field_func(app, fl_em);
}

void
vm_fluid_em_coupling_restore_field(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  fl_em->restore_field_func(app, fl_em);
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
  if (fl_em->curl_rhs) {
    gkyl_array_release(fl_em->curl_rhs);
  }
  if (fl_em->curl_rhs_no_J) {
    gkyl_array_release(fl_em->curl_rhs_no_J);
  }
  gkyl_array_release(fl_em->curl_current);
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    gkyl_array_release(fl_em->fluid_stag[i]);
  }
  if (fl_em->em_stag) {
    gkyl_array_release(fl_em->em_stag);
  }
  gkyl_array_release(fl_em->accel_zero);
  if (fl_em->em_ghost) {
    for (int i = 0; i < fl_em->num_fluid; ++i) {
      gkyl_array_release(fl_em->fluid_ghost[i]);
    }
    gkyl_array_release(fl_em->em_ghost);
  }
  if (fl_em->omega_sq) {
    gkyl_array_release(fl_em->omega_sq);
    gkyl_array_release(fl_em->dens_avg);
    if (app->use_gpu) {
      gkyl_cu_free(fl_em->omega_sq_max);
    } else {
      gkyl_free(fl_em->omega_sq_max);
    }
  }
  gkyl_free(fl_em);
}
