#include <gkyl_alloc.h>
#include <gkyl_app_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_flash_priv.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_wv_maxwell.h>

// Initialize field object
struct flash_field*
flash_field_new(struct gkyl_flash_inp *inp, struct gkyl_flash_app *app)
{
  struct flash_field *f = gkyl_malloc(sizeof(struct flash_field));
  int cdim = app->cdim;

  f->info = inp->field;

  f->mu0 = inp->field.mu0;
  f->epsilon0 = inp->field.epsilon0;

  // Electromagnetic field variables (Ex, Ey, Ez, Bx, By, Bz, phi, psi).
  f->em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  f->tot_em = mkarr(app->use_gpu, 8*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  
  // Magnetic field unit vector and integer array for determining if |B|^2 uses only cell averages for weak division.
  f->bvar = mkarr(app->use_gpu, 9*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  f->cell_avg_magB2 = mk_int_arr(app->use_gpu, 1, app->local_kinetic_ext.volume);

  // Surface expansion of magnetic field unit vector and tensor.
  int Ncomp_surf = 2*cdim*4;
  int Nbasis_surf = app->confBasis.num_basis/(app->confBasis.poly_order + 1); // *only valid for tensor bases for cdim > 1*
  f->bvar_surf = mkarr(app->use_gpu, Ncomp_surf*Nbasis_surf, app->local_kinetic_ext.volume);

  // Volume expansion of div(b)
  f->div_b = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_kinetic_ext.volume);

  // Surface expansion of max b penalization for streaming in PKPM system max(|b_i_l|, |b_i_r|)
  f->max_b = mkarr(app->use_gpu, 2*cdim*Nbasis_surf, app->local_kinetic_ext.volume);

   f->ext_em = mkarr(app->use_gpu, 6*app->confBasis.num_basis, app->local_kinetic_ext.volume);
  // Check if external EM fields are specified
  f->has_ext_em = false;
  if (f->info.ext_em_init) {
    f->has_ext_em = true;
    int num_ret_vals = 6;
    f->proj_ext_em = gkyl_proj_on_basis_new(&app->grid_kinetic, &app->confBasis, 
      app->poly_order+1, num_ret_vals, f->info.ext_em_init, f->info.ext_em_ctx);
  }

  // Check if user specified limiter_fac to adjust diffusion in slope limiter.
  // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_em_vars.h)
  double limiter_fac = f->info.limiter_fac == 0 ? 0.0 : f->info.limiter_fac;

  // Geometry object for Maxwell updater to compute div(b)
  f->geom = gkyl_wave_geom_new(&app->grid_kinetic, &app->global_kinetic_ext, 0, 0, false);

  // Wave equation updater for Maxwell's equations to compute div(b)
  double e_fact = 0.0, b_fact = 0.0;
  double c = 1/sqrt(f->epsilon0*f->mu0);
  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, e_fact, b_fact, app->use_gpu);
  f->calc_bvar = gkyl_dg_calc_em_vars_new(&app->grid_kinetic, &app->confBasis, &app->local_kinetic_ext,
    maxwell, f->geom, limiter_fac, 0, app->use_gpu);
  gkyl_wv_eqn_release(maxwell);

  return f;
}

void
flash_field_apply_ic(struct gkyl_flash_app *app, struct flash_field *field, double t0)
{
  // Project 
  int num_ret_vals = 8;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid_kinetic, &app->confBasis, 
    app->poly_order+1, num_ret_vals, field->info.init, field->info.ctx);

  gkyl_proj_on_basis_advance(proj, t0, &app->local_kinetic_ext, field->em);
  gkyl_proj_on_basis_release(proj);

  // Project external electromagnetic fields onto basis.
  if (field->has_ext_em) {
    gkyl_proj_on_basis_advance(field->proj_ext_em, t0, &app->local_kinetic_ext, field->ext_em);
  }
}

void
flash_field_calc_bvar(gkyl_flash_app *app, struct flash_field *field, struct gkyl_array *em)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(field->tot_em, 0.0);
  gkyl_array_set(field->tot_em, 1.0, em); // em as input arg here to allow for computing div(b) at intermediate steps.

  // Check for external fields and accumulate them onto tot_em
  if (field->has_ext_em) {
    gkyl_array_accumulate_range(field->tot_em, 1.0, field->ext_em, &app->local_kinetic_ext);
  }

  // Compute magnetic field unit vector and unit tensor.
  gkyl_dg_calc_em_vars_advance(field->calc_bvar, field->tot_em, 
    field->cell_avg_magB2, field->bvar, field->bvar_surf);

  // Compute div(b) and max_b = max(|b_i_l|, |b_i_r|)
  gkyl_array_clear(field->div_b, 0.0); // Incremented in each dimension, so clear beforehand
  gkyl_dg_calc_em_vars_div_b(field->calc_bvar, &app->local_kinetic, 
    field->bvar_surf, field->bvar,
    field->max_b, field->div_b);

  app->stat.field_vars_tm += gkyl_time_diff_now_sec(wst);
}

void
flash_field_release(struct flash_field *f)
{
  gkyl_array_release(f->em);
  gkyl_array_release(f->tot_em);
  gkyl_array_release(f->ext_em);
  gkyl_array_release(f->bvar);
  gkyl_array_release(f->bvar_surf);
  gkyl_array_release(f->cell_avg_magB2);
  gkyl_array_release(f->div_b);
  gkyl_array_release(f->max_b);

  if (f->has_ext_em)
    gkyl_proj_on_basis_release(f->proj_ext_em);

  gkyl_wave_geom_release(f->geom);

  gkyl_dg_calc_em_vars_release(f->calc_bvar);

  gkyl_free(f);
}