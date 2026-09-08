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
  gk_field_fem_projection_par(app, field, field->rho_c, field->phi_smooth);
}

static void 
gk_field_ohm_solve_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out){
  struct timespec wst = gkyl_wall_clock();

  // Compute the LHS factor (= k_perp^2/mu_0 + sum_s q_s^2 n_s/m_s)
  // gkyl_array_set(field->dApartdtSlvr_lhs_factor, 1.0, field->dApartdtSlvr_kSq);
  gkyl_array_clear(field->dApartdtSlvr_lhs_factor, 0.0);
  gkyl_array_accumulate_range(field->dApartdtSlvr_lhs_factor,-1.0, field->dApartdtSlvr_kSq, &app->local);
  gkyl_array_accumulate_range(field->dApartdtSlvr_lhs_factor, 1.0, field->lapWeightAmpere, &app->local);

  // Weak division method dApar/dt = sum_s q_s int dv vpar d/dt(F_s) / ( (k_perp^2/mu_0 + sum_s q_s^2/m_s int dv F_s) )
  gkyl_dg_div_op_range(field->div_mem, &app->basis, 0, out, 0, field->currentDensdot,
    0, field->dApartdtSlvr_lhs_factor, &app->local);
  
  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
  
}

static void 
gk_field_ampere_solve_1x_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out){
  struct timespec wst = gkyl_wall_clock();

  // Weak division method Apar = sum_s q_s int dv vpar F_s / (k_perp^2/mu_0)
  gkyl_dg_div_op_range(field->div_mem, &app->basis, 0, out, 0, field->currentDens,
    0, field->lapWeightAmpere, &app->local);

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_field_ampere_solve_1x_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *out){
  // Do nothing.
}

static void
gk_field_em_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  gk_field_accumulate_current_dens_dot(app, field, rhs_in);
  gk_field_accumulate_ohms_kSq(app, field, f_in);
  gk_field_ohm_solve_1x(app, field, field->apardot);
}

static void
gk_field_em_rhs_none(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  // Do nothing.
}

static void
gk_field_fem_projection_par_none(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Do nothing.
}

static void
gk_field_fem_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->phi_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->apar);
  gkyl_array_release(f->apar_curr);
  gkyl_array_release(f->apar1);
  gkyl_array_release(f->aparnew);
  gkyl_array_release(f->apardot);
  gkyl_array_release(f->amperesol);

  if (f->is_em) {
    gkyl_array_release(f->currentDens);
    gkyl_array_release(f->currentDensdot);
    gkyl_array_release(f->currentDens_global);
    gkyl_array_release(f->dApartdtSlvr_lhs_factor);
    gkyl_array_release(f->lapWeightAmpere);
    gkyl_array_release(f->dApartdtSlvr_kSq);
    gkyl_dg_bin_op_mem_release(f->div_mem);
    if (app->use_gpu) {
      gkyl_array_release(f->apar_host);
      gkyl_array_release(f->apardot_host);
      gkyl_array_release(f->amperesol_host);
    }
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
  f->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar_curr = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar1 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->aparnew = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->amperesol = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Allocate additional arrays for EM case.
  if (f->is_em) {
    
    f->apar_host = f->apar;
    f->apardot_host = f->apardot;
    f->amperesol_host = f->amperesol;
    if (app->use_gpu) {
      f->apar_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      f->apardot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      f->amperesol_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    }

    f->currentDens = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->currentDensdot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->lapWeightAmpere = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
    f->dApartdtSlvr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    
    f->currentDens_global = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    f->dApartdtSlvr_lhs_factor = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->div_mem = app->use_gpu? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis)
      : gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    
    // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
    gkyl_array_copy(f->lapWeightAmpere, app->gk_geom->geo_int.jacobgeo);
    gkyl_array_scale(f->lapWeightAmpere, f->info.kperpSq/f->info.mu0);

    f->fem_parproj_ampere_bc = GKYL_FEM_PARPROJ_NONE;
    for (int d=0; d<app->num_periodic_dir; ++d)
      if (app->periodic_dirs[d] == app->cdim-1) f->fem_parproj_ampere_bc = GKYL_FEM_PARPROJ_PERIODIC;

    assert(f->info.mu0 > 0.0);
    f->apar_energy_fac_1d = 0.5/f->info.mu0 * f->info.kperpSq;
    f->apar_energy_fac_1d *= -1.0; // the -1 is to be consistent with ES energy.
  }
  
  // Allocate phi_host for I/O.
  f->phi_host = f->phi_smooth;
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  // Set function pointers.
  f->rhs_phi_func = gk_field_rhs_phi_1x;
  if (f->is_em) {
    f->em_rhs_func = f->info.is_apar_static ? gk_field_em_rhs_none : gk_field_em_rhs_enabled;
    f->ampere_solve = gk_field_ampere_solve_1x_enabled;
    f->fem_projection_par_apar_func = gk_field_fem_projection_par_none;
  } else {
    f->em_rhs_func = gk_field_em_rhs_none;
    f->ampere_solve = gk_field_ampere_solve_1x_none;
  }

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
