#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static void
gk_field_rhs_phi_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Solve the Poisson equation in 1x with the parallel FEM projection.
//  gk_field_fem_projection_par(app, field, field->rho_c, field->phi_smooth);
  // divide rho_c by epsilon
  gkyl_dg_div_op_range(app->species[0].m0.mem_geo, &app->basis, 0, field->phi_smooth, 0, field->rho_c, 0,
                       field->epsilon, &app->local);

  // Apply periodic BCs to phi.
  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, app->periodic_dirs, field->phi_smooth);
}

static void
gk_field_fem_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->phi_smooth);
  gkyl_array_release(f->phi_fem);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->apardot_fem);
  }

  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_array_release(f->epsilon);
  
  gkyl_fem_parproj_release(f->fem_parproj);

  gkyl_array_integrate_release(f->calc_em_energy);
}

void
gk_field_fem_new_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  f->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  
  // Allocate arrays for electrostatic potential.
  f->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Allocate electromagnetic arrays if needed.
  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    f->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->apardot_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }
  
  // Allocate phi_host for I/O.
  f->phi_host = f->phi_smooth;
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  f->rhs_phi_func = gk_field_rhs_phi_1x;

  // Allocate array for the polarization weight times geometric coefficients.
  f->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);

  double polarization_weight = 0.0;
  double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
  // Linearized polarization density
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
  }
  // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
  gkyl_array_copy(f->epsilon, app->gk_geom->geo_int.jacobgeo);
  gkyl_array_scale(f->epsilon, polarization_weight);
  gkyl_array_scale(f->epsilon, f->info.kperpSq);

  double es_energy_fac_1d_adiabatic = 0.0;
  if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_adiabatic;
    // Add the contribution from adiabatic electrons
    double n_s0 = f->info.electron_density;
    double q_s = f->info.electron_charge;
    double T_s = f->info.electron_temp;
    double quasineut_contr = q_s*n_s0*q_s/T_s;
    es_energy_fac_1d_adiabatic = 0.5*quasineut_contr;
    
    struct gkyl_array *epsilon_adiab = mkarr(app->use_gpu, f->epsilon->ncomp, f->epsilon->size);
    gkyl_array_copy(epsilon_adiab, app->gk_geom->geo_int.jacobgeo);
    gkyl_array_scale(epsilon_adiab, quasineut_contr);
    gkyl_array_accumulate(f->epsilon, 1., epsilon_adiab);
    gkyl_array_release(epsilon_adiab);
  } else {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_poisson;
  }

  // Gather epsilon for (global) smoothing in z.
  struct gkyl_array *epsilon_global = mkarr(app->use_gpu, f->epsilon->ncomp, app->global_ext.volume);
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, f->epsilon, epsilon_global);

  // Potential smoothing (in z) updater
  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d=0; d<app->num_periodic_dir; ++d)
    if (app->periodic_dirs[d] == app->cdim-1) {
      fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;
    }

  f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc, 0, epsilon_global, 0, app->use_gpu);

  f->es_energy_fac_1d = 0.5*polarization_weight*f->info.kperpSq + es_energy_fac_1d_adiabatic;

  f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
    1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);

  // Create operator needed for FLR effects.
  f->use_flr = false;
  f->invert_flr = gk_field_invert_flr_none;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    if (s->info.flr.type)
      f->use_flr = f->use_flr || s->info.flr.type;
  }

  gkyl_array_release(epsilon_global);

  f->release_func = gk_field_fem_release_1x;
}
