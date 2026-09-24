#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_neut_species_collisionless_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_neut_species_collisionless_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_dg_updater_vlasov_advance(gkcls->vlasov_slvr, &species->local, 
    fin, species->cflrate, rhs);

  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_collisionless_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, double tm, int frame)
{
}

static void
gk_neut_species_collisionless_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_neut_species_collisionless_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, 
  struct gk_collisionless *gkcls)
{
  gkcls->collisionless_id = gkns->info.collisionless.type;
  gkcls->write_diagnostics = gkns->info.collisionless.write_diagnostics;

  gkcls->write_diags_func_neut = gk_neut_species_collisionless_write_diags_disabled;
  gkcls->rhs_func_neut = gk_neut_species_collisionless_rhs_disabled;

  if (gkcls->collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL) {

    int cdim = app->cdim, vdim = gkns->info.vdim;
    int pdim = cdim+vdim;

    // Determine which directions are zero-flux. By default
    // we do not have zero-flux boundary conditions in any direction.
    bool is_zero_flux[2*GKYL_MAX_DIM] = {false};
    for (int dir=0; dir<app->cdim; ++dir) {
      if (gkns->lower_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        is_zero_flux[dir] = true;
      if (gkns->upper_bc[dir].type == GKYL_BC_GK_SPECIES_ZERO_FLUX)
        is_zero_flux[dir+pdim] = true;
    }

    // Need to figure out size of alpha_surf and sgn_alpha_surf by finding size of surface basis set 
    struct gkyl_basis surf_basis, surf_quad_basis;
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
    gkyl_cart_modal_tensor(&surf_quad_basis, pdim-1, app->poly_order);
  
    int alpha_surf_sz = (cdim+vdim)*surf_basis.num_basis; 
    int sgn_alpha_surf_sz = (cdim+vdim)*surf_quad_basis.num_basis; // sign(alpha) is store at quadrature points
  
    // Allocate arrays to store fields:
    // 1. alpha_surf (surface phase space velocity)
    // 2. sgn_alpha_surf (sign(alpha_surf) at quadrature points)
    // 3. const_sgn_alpha (boolean for if sign(alpha_surf) is a constant, either +1 or -1)
    gkcls->alpha_surf = mkarr(app->use_gpu, alpha_surf_sz, gkns->local_ext.volume);
    gkcls->sgn_alpha_surf = mkarr(app->use_gpu, sgn_alpha_surf_sz, gkns->local_ext.volume);
    gkcls->const_sgn_alpha = mk_int_arr(app->use_gpu, cdim+vdim, gkns->local_ext.volume);
  
    // Pre-compute alpha_surf, sgn_alpha_surf, const_sgn_alpha, and cot_vec since they are time-independent
    struct gkyl_dg_calc_canonical_pb_vars *calc_vars = gkyl_dg_calc_canonical_pb_vars_new(&gkns->grid, 
      &app->basis, &gkns->basis, app->use_gpu);
    gkyl_dg_calc_canonical_pb_vars_alpha_surf(calc_vars, &app->local, &gkns->local, &gkns->local_ext, gkns->hamil,
      gkcls->alpha_surf, gkcls->sgn_alpha_surf, gkcls->const_sgn_alpha);
    gkyl_dg_calc_canonical_pb_vars_release(calc_vars);
  
    struct gkyl_dg_canonical_pb_auxfields aux_inp = {.hamil = gkns->hamil, .alpha_surf = gkcls->alpha_surf, 
      .sgn_alpha_surf = gkcls->sgn_alpha_surf, .const_sgn_alpha = gkcls->const_sgn_alpha};
  
    gkcls->vlasov_slvr = gkyl_dg_updater_vlasov_new(&gkns->grid, &app->basis, &gkns->basis, 
      &app->local, &gkns->local_vel, &gkns->local, is_zero_flux, gkns->model_id, gkns->field_id, &aux_inp, app->use_gpu);

    // Methods chosen at runtime.
    gkcls->rhs_func_neut = gk_neut_species_collisionless_rhs_enabled;
    if (gkcls->write_diagnostics) {
      gkcls->write_diags_func_neut = gk_neut_species_collisionless_write_diags_enabled;
    }
  }
}

void
gk_neut_species_collisionless_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_collisionless *gkcls, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkcls->rhs_func_neut(app, species, gkcls, fin, rhs);
}

void
gk_neut_species_collisionless_write_diags(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_collisionless *gkcls, double tm, int frame)
{
  gkcls->write_diags_func_neut(app, gkns, gkcls, tm, frame);
}

void
gk_neut_species_collisionless_release(const struct gkyl_gyrokinetic_app *app, const struct gk_collisionless *gkcls)
{
  if (gkcls->collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL) {
    gkyl_array_release(gkcls->alpha_surf);
    gkyl_array_release(gkcls->sgn_alpha_surf);
    gkyl_array_release(gkcls->const_sgn_alpha);

    gkyl_dg_updater_vlasov_release(gkcls->vlasov_slvr);

    if (gkcls->write_diagnostics) {
    }
  }
}
