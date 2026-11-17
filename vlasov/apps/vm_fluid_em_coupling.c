#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

bool
vm_fluid_em_coupling_supported(const struct vm_fluid_species *f)
{
  // Euler only: isothermal Euler has the momentum layout the coupling expects
  // but no fluid DG solver yet, and is rejected at fluid initialization.
  return f->eqn_type == GKYL_EQN_EULER;
}

// initialize fluid-EM coupling object
struct vm_fluid_em_coupling *
vm_fluid_em_coupling_new(struct gkyl_vlasov_app *app)
{
  struct vm_fluid_em_coupling *fl_em = gkyl_malloc(sizeof(struct vm_fluid_em_coupling));

  // Gather the coupled (momentum-carrying) fluid species in declaration order.
  int num_species = app->num_species;
  fl_em->num_fluid = 0;
  for (int i = 0; i < num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (sp->fluid && vm_fluid_em_coupling_supported(sp->fluid)) {
      fl_em->species[fl_em->num_fluid] = sp;
      fl_em->qbym[fl_em->num_fluid] = sp->charge / sp->mass;
      fl_em->num_fluid += 1;
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
    .epsilon0 = app->field->info.epsilon0,
    .is_field_static = app->field->info.is_static,
    .use_gpu = app->use_gpu,
  });

  int num_basis = app->basis.num_basis;
  long num_cells = app->local_ext.volume;
  fl_em->curl_rhs = mkarr(app->use_gpu, 8 * num_basis, num_cells);
  fl_em->curl_current = mkarr(app->use_gpu, 3 * num_basis, num_cells);
  gkyl_array_clear(fl_em->curl_current, 0.0);
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    const struct gkyl_array *fluid = fl_em->species[i]->fluid->fluid;
    fl_em->fluid_stag[i] = mkarr(app->use_gpu, fluid->ncomp, fluid->size);
  }
  fl_em->em_stag = mkarr(app->use_gpu, app->field->em->ncomp, app->field->em->size);

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

  for (int i = 0; i < num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    fluids[i] = fs->fluid;
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

  // The currents of the solve besides the species': the applied current and the effective
  // current of the curl. A static field takes none of them (curl_current stays zero), so its E
  // does not change.
  struct vm_field *field = app->field;
  if (!field->info.is_static) {
    // Stormer-Verlet in the curl terms: E advances only here, with c^2 curl B of the field at
    // the start of the step (B is not changed by the sources) entering as the effective current
    // -eps0 c^2 curl B; the transport step then drifts B with this E.
    gkyl_array_clear(fl_em->curl_rhs, 0.0);
    gkyl_hyper_dg_advance(field->slvr, &app->local, field->em, field->cflrate, fl_em->curl_rhs);
    gkyl_array_set(fl_em->curl_current, 1.0, field->app_current);
    gkyl_array_accumulate_offset(fl_em->curl_current, -field->info.epsilon0, fl_em->curl_rhs, 0);
  }

  gkyl_dg_calc_fluid_em_coupling_advance(
    fl_em->slvr, dt, app_accels, field->ext_em, fl_em->curl_current, fluids, field->em
  );

  // The solve advances the momentum and the electric field, and the fluid energy by the change
  // of the bulk kinetic energy (the work of the field): the density is unchanged, and so is the
  // internal energy at the nodes of the solve.
  for (int i = 0; i < num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    vm_fluid_species_apply_bc(app, fs, fluids[i]);
  }
  vlasov_field_apply_bc(app, app->field->em);
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
  gkyl_array_release(fl_em->curl_current);
  for (int i = 0; i < fl_em->num_fluid; ++i) {
    gkyl_array_release(fl_em->fluid_stag[i]);
  }
  gkyl_array_release(fl_em->em_stag);
  gkyl_free(fl_em);
}
