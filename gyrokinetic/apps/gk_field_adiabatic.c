#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>

#include <assert.h>
#include <time.h>

// Adiabatic species response in the field equation. Each adiabatic species s enters through
// its background charge density q_s n_s0 J and K_s = (q_s^2 n_s0/T_s)*J, computed from its
// t=0 moments (see gk_species_scaling_adiabatic_init); the field holds their sums rho_bg and
// K. The dimensionality is selected once in gk_field_adiabatic_new, which sets the methods
// used afterwards:
// - 1x: there is no flux-surface average and K is added to the weight of the parallel
//   projection that solves the field equation.
// - 2x/3x: the field equation
//     -div(eps grad_perp phi) + K (phi - <phi>) = rho,
//   with <phi> = int J phi dy dz / int J dy dz, is solved by the fem_poisson_perp_adiabatic
//   updater (Woodbury identity), and the flux-surface average (FSA) and its extension are
//   computed by flux_surf_avg. With Pade FLR effects the solve is for Phi_0 (phi = A Phi_0,
//   A = 1 - rho^2*nabla_perp^2) and the response K (phi - <phi>) becomes
//     [K - K_0 rho^2 nabla_perp^2] (Phi_0 - <Phi_0>),
//   with K_0 the J-weighted average of K (see gk_field_adiabatic_k_avg), i.e. the updater's
//   Laplacian term with factor K_0 rho^2/eps_pol.
// Only adiabatic electrons are implemented in 2x/3x; methods specific to them are named
// gk_field_adiab_elc_*.

// Add the background charge density of the adiabatic species to rho_c.
static void
gk_field_adiabatic_accumulate_rho_c_bg(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  gkyl_array_accumulate_range(field->rho_c, 1.0, field->gk_adiab.rho_bg, &app->local);
}

// Integral of a conf-space array over the global domain.
static double
gk_field_adiabatic_integral(struct gkyl_gyrokinetic_app *app, const struct gkyl_array *arr)
{
  double *integ = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  double *integ_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  struct gkyl_array_integrate *int_op = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu
  );
  gkyl_array_integrate_advance(int_op, arr, 1.0, NULL, &app->local, NULL, integ);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, integ, integ_red);

  double integ_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(integ_ho, integ_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(integ);
    gkyl_cu_free(integ_red);
  } else {
    integ_ho[0] = integ_red[0];
    gkyl_free(integ);
    gkyl_free(integ_red);
  }
  gkyl_array_integrate_release(int_op);
  return integ_ho[0];
}

// J-weighted average of the adiabatic response coefficient K = sum_s q_s^2 n_s0/T_s over the
// global domain (exact for uniform profiles).
static double
gk_field_adiabatic_k_avg(struct gkyl_gyrokinetic_app *app, const struct gk_field *gkf)
{
  return gk_field_adiabatic_integral(app, gkf->gk_adiab.kJ) /
         gk_field_adiabatic_integral(app, app->gk_geom->geo_int.jacobgeo);
}

/** 1x: K is added to the weight of the parallel projection. **/

static void
gk_field_adiabatic_init_1x(gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkyl_array_accumulate(gkf->epsilon, 1.0, gkf->gk_adiab.kJ);
  gk_field_fem_parproj_new(app, gkf, gkf->epsilon);

  // The 1x energy uses a scalar factor, the J-weighted average of K.
  gkf->es_energy_fac_1d += 0.5 * gk_field_adiabatic_k_avg(app, gkf);
}

static void
gk_field_adiabatic_dphi_1x(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  gkyl_array_copy(out, phi);
}

// The parallel projection is released with the 1x field solver.
static void
gk_field_adiabatic_release_1x(const gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
}

/** 2x/3x: adiabatic electrons with a flux-surface average. **/

// Smooth a right-hand side along z like the charge density, before each perpendicular solve.
static void
gk_field_adiab_elc_fem_parproj_op(struct gkyl_array *rhs, void *ctx)
{
  gkyl_gyrokinetic_app *app = ctx;
  app->field->fem_projection_par_rho_func(app, app->field, rhs, rhs);
}

static void
gk_field_adiab_elc_init_2x3x(gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  // 1D (x) basis and ranges of the flux-surface averages.
  gkyl_cart_modal_serendip(&ad->basis_x, 1, app->poly_order);
  gkyl_range_init(&ad->local_x_ext, 1, &app->local_ext.lower[0], &app->local_ext.upper[0]);
  gkyl_sub_range_init(&ad->local_x, &ad->local_x_ext, &app->local.lower[0], &app->local.upper[0]);
  ad->phi_fsa = mkarr(app->use_gpu, ad->basis_x.num_basis, ad->local_x_ext.volume);
  ad->phi_buff = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Flux-surface average (1D in x) and its extension, reduced over the ranks decomposing z.
  ad->fsa = gkyl_flux_surf_avg_inew(&(struct gkyl_flux_surf_avg_inp){
    .grid = &app->grid,
    .basis = app->basis,
    .basis_x = ad->basis_x,
    .local = &app->local,
    .local_ext = &app->local_ext,
    .local_x = &ad->local_x,
    .local_x_ext = &ad->local_x_ext,
    .jacobgeo = app->gk_geom->geo_int.jacobgeo,
    .is_deflated = true,
    .comm = app->comm,
    .use_gpu = app->use_gpu,
  });

  // Energy diagnostic, (1/2) int K (phi-<phi>)^2.
  ad->calc_energy = gkyl_array_integrate_new(
    &app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED, app->use_gpu
  );
  ad->energy_red = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));

  // Factor of the Laplacian term in the response with FLR effects, K_0 rho^2/eps_pol, with the
  // constant K_0 the J-weighted average of K (from the t=0 moments of the adiabatic species).
  double lap_fac = 0.0;
  if (gkf->use_flr) {
    double k_avg = gk_field_adiabatic_k_avg(app, gkf);
    assert(k_avg > 0.0);
    lap_fac = k_avg * gkf->flr_local_fac;
    gkyl_gyrokinetic_app_cout(
      app, stdout, "Adiabatic response FLR term: K_0 = %.6e, Laplacian factor = %.6e.\n", k_avg,
      lap_fac
    );
  }

  // Field solver with the adiabatic response, smoothing every right-hand side along z.
  struct timespec wst = gkyl_wall_clock();
  ad->solver = gkyl_fem_poisson_perp_adiabatic_inew(&(struct gkyl_fem_poisson_perp_adiabatic_inp){
    .grid = &app->grid,
    .basis = app->basis,
    .local = &app->local,
    .local_ext = &app->local_ext,
    .basis_x = ad->basis_x,
    .local_x = &ad->local_x,
    .local_x_ext = &ad->local_x_ext,
    .fsa = ad->fsa,
    .bcs = &gkf->poisson_bcs,
    .bias_line_list = gkf->info.bias_line_list,
    .epsilon = gkf->epsilon,
    .kJ = ad->kJ,
    .lap_fac = lap_fac,
    .fem_parproj_op = gk_field_adiab_elc_fem_parproj_op,
    .fem_parproj_ctx = app,
    .use_gpu = app->use_gpu,
  });

  gkyl_gyrokinetic_app_cout(
    app, stdout, "Adiabatic response with FSA: m = %ld, cond_1(I-G) = %.4e, setup %.3e sec.\n",
    ad->local_x.volume * ad->basis_x.num_basis, gkyl_fem_poisson_perp_adiabatic_cond(ad->solver),
    gkyl_time_diff_now_sec(wst)
  );
}

// Solve the field equation with the adiabatic response, apply the FLR correction and smooth
// phi along z.
static void
gk_field_adiab_elc_rhs_phi_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Solve -div(eps grad_perp phi) + K (phi - <phi>) = rho (for Phi_0 with FLR effects).
  gkyl_fem_poisson_perp_adiabatic_advance(field->gk_adiab.solver, field->rho_c, field->phi_smooth);

  field->invert_flr(app, field, field->phi_smooth);

  // Smooth the potential along z.
  field->fem_projection_par_phi_func(app, field, field->phi_smooth, field->phi_smooth);
}

// Compute out = phi - <phi>, with <phi> the flux-surface average of phi.
static void
gk_field_adiab_elc_dphi_2x3x(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  const struct gk_field_adiabatic *ad = &gkf->gk_adiab;
  gkyl_flux_surf_avg_advance(ad->fsa, phi, ad->phi_fsa);
  gkyl_flux_surf_avg_inflate(ad->fsa, ad->phi_fsa, out);
  gkyl_array_scale_range(out, -1.0, &app->local);
  gkyl_array_accumulate_range(out, 1.0, phi, &app->local);
}

// Add factor*(1/2) int K (phi-<phi>)^2 over the local range to out (device memory if use_gpu).
static void
gk_field_adiab_elc_energy_2x3x(
  gkyl_gyrokinetic_app *app, const struct gk_field *gkf, double factor, double *out
)
{
  const struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  gk_field_adiab_elc_dphi_2x3x(app, gkf, gkf->phi_smooth, ad->phi_buff);
  gkyl_array_integrate_advance(
    ad->calc_energy, ad->phi_buff, 0.5 * factor, ad->kJ, &app->local, &app->local, ad->energy_red
  );

  if (app->use_gpu) {
    double en[2];
    gkyl_cu_memcpy(&en[0], out, sizeof(double), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_memcpy(&en[1], ad->energy_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
    en[0] += en[1];
    gkyl_cu_memcpy(out, &en[0], sizeof(double), GKYL_CU_MEMCPY_H2D);
  } else {
    out[0] += ad->energy_red[0];
  }
}

static void
gk_field_adiab_elc_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  gkyl_fem_poisson_perp_adiabatic_release(ad->solver);
  gkyl_flux_surf_avg_release(ad->fsa);
  gkyl_array_release(ad->phi_fsa);
  gkyl_array_release(ad->phi_buff);
  gkyl_array_integrate_release(ad->calc_energy);
  if (app->use_gpu) {
    gkyl_cu_free(ad->energy_red);
  } else {
    gkyl_free(ad->energy_red);
  }
}

void
gk_field_adiabatic_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  assert(gkf->gkfield_id == GKYL_GK_FIELD_ES);
  // With FLR effects phi = A Phi_0 must be applied with the FEM operator (the local term
  // assumes no adiabatic response).
  assert(!gkf->use_flr || gkf->info.flr.use_fem_operator);

  ad->rho_bg = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  ad->kJ = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  gkf->accumulate_rhoc_bg_func = gk_field_adiabatic_accumulate_rho_c_bg;
  if (app->cdim > 1) {
    // Adiabatic electrons respond to phi - <phi>, only valid on closed field lines.
    assert(!app->gk_geom->has_LCFS);
    assert(!gkf->info.bias_line_list || gkf->info.bias_line_list->num_bias_line == 0);
    for (int k = 0; k < ad->num_species; k++) {
      assert(app->species[ad->species_idx[k]].info.charge < 0.0); // Adiabatic ions not implemented.
    }
    ad->init_func = gk_field_adiab_elc_init_2x3x;
    ad->dphi_func = gk_field_adiab_elc_dphi_2x3x;
    ad->release_func = gk_field_adiab_elc_release_2x3x;
    gkf->rhs_phi_func = gk_field_adiab_elc_rhs_phi_2x3x;
    gkf->energy_adiab_func = gk_field_adiab_elc_energy_2x3x;
  } else {
    ad->init_func = gk_field_adiabatic_init_1x;
    ad->dphi_func = gk_field_adiabatic_dphi_1x;
    ad->release_func = gk_field_adiabatic_release_1x;
  }
}

void
gk_field_adiabatic_init(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  struct gk_field_adiabatic *ad = &gkf->gk_adiab;

  // Response coefficients of each adiabatic species from its t=0 moments, summed into K and
  // the background charge density.
  gkyl_array_clear(ad->rho_bg, 0.0);
  gkyl_array_clear(ad->kJ, 0.0);
  for (int k = 0; k < ad->num_species; k++) {
    struct gk_species *gks = &app->species[ad->species_idx[k]];
    gk_species_scaling_adiabatic_init(app, gks, &gks->sca);
    gkyl_array_accumulate(ad->rho_bg, 1.0, gks->sca.rho_bg);
    gkyl_array_accumulate(ad->kJ, 1.0, gks->sca.kJ);
  }

  ad->init_func(app, gkf);
}

void
gk_field_adiabatic_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf)
{
  gkf->gk_adiab.release_func(app, gkf);
  gkyl_array_release(gkf->gk_adiab.rho_bg);
  gkyl_array_release(gkf->gk_adiab.kJ);
}
