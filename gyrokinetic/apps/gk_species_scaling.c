#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_scaling_cross_moms_disabled(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  // Do nothing.
}

static void
gk_species_scaling_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // Do nothing.
}

static void
gk_species_scaling_apply_fixed_fraction(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  struct gk_species *gks_ref = &app->species[sca->ref_species_idx];

  // Reciprocal of the density of this species.
  struct gkyl_array *m0_inv = gks->lte.moms.marr;
  gk_species_moment_calc(&gks->m0, gks->local, app->local, fin);
  gkyl_dg_inv_op_range(&app->basis, 0, m0_inv, 0, gks->m0.marr, &app->local);

  // Density of the reference species.
  gk_species_moment_calc(&gks_ref->m0, gks_ref->local, app->local, gks_ref->f);

  // Ratio of the densities.
  struct gkyl_array *m0ratio = gks->m0.marr;
  gkyl_dg_mul_op_range(&app->basis, 0, m0ratio, 0, m0_inv, 0, gks_ref->m0.marr, &app->local);

  // Multiply this species by ratio of densities times fixed_fraction.
  gkyl_array_scale_range(m0ratio, sca->fixed_fraction, &app->local);

  gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, fin, m0ratio, fin, &app->local, &gks->local);
}

static void
copy_lower_z_ghost_to_all_z_conf(gkyl_gyrokinetic_app *app, struct gkyl_array *arr_out,
  struct gkyl_array *arr_in, int comp_in, struct gkyl_array *buff)
{
  // Copy the component 'comp_in' of the lower ghost cell of 'arr_in' into the
  // interior cells of 'arr_out' for configuration-space ranges and arrays.
  int cdim = app->cdim;
  int idx_par = cdim-1;

  // Copy lower ghost of input field into buffer.
  struct gkyl_range* lower_ghost = &app->local_lower_ghost[idx_par];
  gkyl_array_set_offset_range(buff, 1.0, arr_in, comp_in, lower_ghost); 

  // Create a range for each z cell, and copy from buffer into that range.
  int num_cells_z = app->local.upper[idx_par]-app->local.lower[idx_par]+1;
  int lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM];
  for (int d=0; d<idx_par; ++d) {
    lower[d] = app->local.lower[d];
    upper[d] = app->local.upper[d];
  }
  struct gkyl_range rng_curr;
  for (int i=0; i<num_cells_z; ++i) {
    lower[idx_par] = app->local.lower[idx_par]+i;
    upper[idx_par] = lower[idx_par];
    gkyl_sub_range_init(&rng_curr, &app->local_ext, lower, upper);
    gkyl_array_copy_range_to_range(arr_out, buff, &rng_curr, lower_ghost); 
  }
}

static void
gk_species_scaling_apply_boltzmann(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  struct gk_species *gks_ref = &app->species[sca->ref_species_idx];

  // Compute density and temperature.
  gk_species_moment_calc(&gks->lte.moms, gks->local, app->local, fin);
  gkyl_dg_div_op_range(gks->lte.moms.mem_geo, &app->basis, 0, gks->lte.moms.marr,
    0, gks->lte.moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local);

  // Boltzmann density = n_sheath * exp(-q * (phi-phi_sheath)/T ).
  // Copy the sheath potential from the lower ghost cell to the interior cells.
  int idx_par = app->cdim-1;
  int off = 2*idx_par;
  copy_lower_z_ghost_to_all_z_conf(app, sca->sheath_val, app->field->sheath_vals[off], 1*app->basis.num_basis, sca->buffer_conf);

  // Compute ( phi-phi_sheath)/(T/m) ).
  gkyl_array_copy_range(sca->buffer_conf, app->field->phi_smooth, &app->local);
  gkyl_array_accumulate_range(sca->buffer_conf, -1.0, sca->sheath_val, &app->local);
  gkyl_dg_div_op_range(gks->lte.moms.mem_geo, &app->basis, 0, sca->buffer_conf,
    0, sca->buffer_conf, 2, gks->lte.moms.marr, &app->local);

  // Compute exp(-q * (phi-phi_sheath)/T ).
  double qDm = gks->info.charge/gks->info.mass;
  struct gkyl_array *m0_boltz = gks->m0.marr;
  gkyl_proj_exp_on_basis_advance(sca->proj_exp, &app->local, 1.0, -qDm, sca->buffer_conf, m0_boltz);

  // Copy the sheath density from the lower ghost cell to the interior cells.
  copy_lower_z_ghost_to_all_z_conf(app, sca->sheath_val, app->field->sheath_vals[off], 0*app->basis.num_basis, sca->buffer_conf);

  // Compute n_sheath * exp(-q * (phi-phi_sheath)/T ).
  gkyl_dg_mul_op_range(&app->basis, 0, sca->buffer_conf, 0, sca->sheath_val, 0, m0_boltz, &app->local);

  // Reciprocal of the density.
  struct gkyl_array *m0_inv = gks->m0.marr;
  gkyl_dg_inv_op_range(&app->basis, 0, m0_inv, 0, gks->lte.moms.marr, &app->local);

  // Ratio of the densities.
  struct gkyl_array *m0ratio = gks->m0.marr;
  gkyl_dg_mul_op_range(&app->basis, 0, m0ratio, 0, m0_inv, 0, sca->buffer_conf, &app->local);

  // Multiply this species by ratio of densities times fixed_fraction.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, fin, m0ratio, fin, &app->local, &gks->local);
}

static void
gk_species_scaling_apply_disabled(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  // Do nothing.
}

static void
gk_species_scaling_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  // React diagnostics usually written from gk_species.
  // In the case of static gk_species, write_diagnostics flag
  // can be used to check reaction rates from gk_species.
  struct timespec wst = gkyl_wall_clock();
//  // Compute reaction rate
//  const struct gkyl_array *fin[app->num_species];
//  const struct gkyl_array *fin_neut[app->num_species];
//  for (int i=0; i<app->num_species; ++i) {
//    fin[i] = app->species[i].f;
//  }
//  for (int i=0; i<app->num_species; ++i) {
//    fin_neut[i] = app->species[i].f;
//  }
//  gk_species_scaling_cross_moms(app, gks, gkr, fin, fin_neut);
  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  
  struct timespec wtm = gkyl_wall_clock();
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gk_species_scaling_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  // Do nothing
}

void 
gk_species_scaling_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_scaling *sca)
{
  struct gkyl_gyrokinetic_scaling_inp *sca_inp = &gks->info.scaling;

  sca->type = 0;
  sca->cross_moms_func = gk_species_scaling_cross_moms_disabled;
  sca->rhs_func = gk_species_scaling_rhs_disabled;
  sca->apply_func = gk_species_scaling_apply_disabled;
  sca->write_func = gk_species_scaling_write_disabled;

  if (sca_inp->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
    sca->type = sca_inp->type;
    sca->write_diagnostics = sca_inp->write_diagnostics;

    sca->fixed_fraction = sca_inp->fixed_fraction;

    sca->apply_func = gk_species_scaling_apply_fixed_fraction;
    if (sca->write_diagnostics)
      sca->write_func = gk_species_scaling_write_enabled;
    else
      sca->write_func = gk_species_scaling_write_disabled;
  }
  else if (sca_inp->type == GKYL_GK_SPECIES_SCALING_BOLTZMANN) {

    sca->type = sca_inp->type;
    sca->sheath_val = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    sca->buffer_conf = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    sca->proj_exp = gkyl_proj_exp_on_basis_new(&app->basis, app->basis.poly_order+1, app->use_gpu);

    sca->apply_func = gk_species_scaling_apply_boltzmann;
    if (sca->write_diagnostics)
      sca->write_func = gk_species_scaling_write_enabled;
    else
      sca->write_func = gk_species_scaling_write_disabled;
  }
}

void 
gk_species_scaling_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
    struct gkyl_gyrokinetic_scaling_inp *sca_inp = &gks->info.scaling;

    // Fetch index of species for indexing arrays.
    sca->ref_species_idx = gk_find_species_idx(app, sca_inp->ref_species_name);
  }
}

void 
gk_species_scaling_apply_ic_cross(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
  }
}

void
gk_species_scaling_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  sca->cross_moms_func(app, gks, sca, fin, fin_neut);
}

void
gk_species_scaling_rhs(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  sca->rhs_func(app, gks, sca, fin, rhs);
}

void
gk_species_scaling_apply(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_scaling *sca, struct gkyl_array *fin, struct gkyl_array **bflux[])
{
  sca->apply_func(app, gks, sca, fin, bflux);
}

void
gk_species_scaling_write(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_scaling *sca, int ridx, double tm, int frame)
{
  sca->write_func(app, gks, sca, ridx, tm, frame);
}

void 
gk_species_scaling_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_scaling *sca)
{
  if (sca->type == GKYL_GK_SPECIES_SCALING_FIXED_FRACTION) {
  } 
  else if (sca->type == GKYL_GK_SPECIES_SCALING_BOLTZMANN) {
    gkyl_proj_exp_on_basis_release(sca->proj_exp);
    gkyl_array_release(sca->sheath_val);
    gkyl_array_release(sca->buffer_conf);
  }
}
