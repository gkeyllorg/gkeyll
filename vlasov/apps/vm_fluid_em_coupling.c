#include <assert.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

// initialize fluid-EM coupling object
struct vm_fluid_em_coupling*
vm_fluid_em_coupling_init(struct gkyl_vlasov_app *app)
{
  struct vm_fluid_em_coupling *fl_em = gkyl_malloc(sizeof(struct vm_fluid_em_coupling));

  // Gather the fluid-bearing species in declaration order.
  int num_species = app->num_species + app->num_fluid_species;
  fl_em->num_fluid = 0;
  for (int i=0; i<num_species; ++i) {
    struct vlasov_species *sp = &app->species[i];
    if (sp->fluid) {
      fl_em->species[fl_em->num_fluid] = sp;
      fl_em->qbym[fl_em->num_fluid] = sp->charge/sp->mass;
      fl_em->num_fluid += 1;
    }
  }

  // Initialize solver
  fl_em->slvr = gkyl_dg_calc_fluid_em_coupling_new(&app->basis, &app->local,
    fl_em->num_fluid, fl_em->qbym, app->field->info.epsilon0, app->use_gpu);

  return fl_em;
}

void
vm_fluid_em_coupling_update(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em,
  double tcurr, double dt)
{
  int num_fluid = fl_em->num_fluid;

  struct gkyl_array *fluids[GKYL_MAX_SPECIES];
  const struct gkyl_array *app_accels[GKYL_MAX_SPECIES];

  for (int i=0; i<num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;
    fluids[i] = fs->fluid;

    if (fs->eqn_type == GKYL_EQN_EULER) {
      gkyl_dg_calc_fluid_vars_advance(fs->calc_fluid_vars,
        fluids[i], fs->cell_avg_prim, fs->u, fs->u_surf);
      gkyl_dg_calc_fluid_vars_ke(fs->calc_fluid_vars, &app->local,
        fluids[i], fs->u, fs->ke_old);
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
  // (this object only exists when a dynamic Maxwell field is present).
  // Note: external EM field and applied currents use proj_on_basis
  // so does copy to GPU every call if app->use_gpu = true.
  if (app->field->ext_em_evolve) {
    vlasov_field_calc_ext_em(app, tcurr);
  }
  if (app->field->app_current_evolve) {
    vlasov_field_calc_app_current(app, tcurr);
  }

  gkyl_dg_calc_fluid_em_coupling_advance(fl_em->slvr, dt,
    app_accels, app->field->ext_em, app->field->app_current,
    fluids, app->field->em);

  for (int i=0; i<num_fluid; ++i) {
    struct vm_fluid_species *fs = fl_em->species[i]->fluid;

    // Compute the updated energy from the old and new kinetic energies
    if (fs->eqn_type == GKYL_EQN_EULER) {
      gkyl_dg_calc_fluid_vars_advance(fs->calc_fluid_vars,
        fluids[i], fs->cell_avg_prim, fs->u, fs->u_surf);
      gkyl_dg_calc_fluid_vars_ke(fs->calc_fluid_vars, &app->local,
        fluids[i], fs->u, fs->ke_new);
      gkyl_dg_calc_fluid_em_coupling_energy(fl_em->slvr,
        fs->ke_old, fs->ke_new, fluids[i]);
    }
    vm_fluid_species_apply_bc(app, fs, fluids[i]);
  }
  vlasov_field_apply_bc(app, app->field->em);
}

void
vm_fluid_em_coupling_release(struct gkyl_vlasov_app *app, struct vm_fluid_em_coupling *fl_em)
{
  gkyl_dg_calc_fluid_em_coupling_release(fl_em->slvr);
  gkyl_free(fl_em);
}
