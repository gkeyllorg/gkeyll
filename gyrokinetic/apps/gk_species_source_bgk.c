#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
proj_on_basis_c2p_phase_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gk_proj_on_basis_c2p_func_ctx *c2p_ctx = ctx;
  int cdim = c2p_ctx->cdim; // Assumes update range is a phase range.
  gkyl_position_map_eval_mc2nu(c2p_ctx->pos_map, xcomp, xphys);
  gkyl_velocity_map_eval_c2p(c2p_ctx->vel_map, &xcomp[cdim], &xphys[cdim]);
}

static double
gk_species_source_bgk_volume_integrate(gkyl_gyrokinetic_app *app, struct gk_source_bgk *src, const struct gkyl_array *arrin)
{
  // Compute the volume integral of arrin.
  gkyl_array_integrate_advance(src->vol_integ_op, arrin, 1.0,
    0, &app->local, 0, src->volint_local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, src->volint_local, src->volint_global);
  double volint_global = 0.0;
  if (app->use_gpu)
    gkyl_cu_memcpy(&volint_global, src->volint_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&volint_global, src->volint_global, sizeof(double));
  return volint_global;
}

static void
gk_species_source_bgk_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_species_source_bgk_rhs_feq_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  // compute -nu * (f - feq) and add to the rhs.
  // Use f_lte to store Jrate*feq (coming from the projection done at initialization)
  gkyl_array_clear(species->lte.f_lte, 0.0);
  gkyl_array_accumulate(species->lte.f_lte, 1.0, src->Jrate_df);
  // Compute Jrate_df
  gkyl_array_clear(src->Jrate_df, 0.0);
  gkyl_bgk_collisions_advance(src->bgk_op, &app->local, &species->local, 
    src->rate, species->lte.f_lte, fin, src->implicit_step, src->dt_implicit, src->Jrate_df, species->cflrate);
  // Add to the rhs.
  gkyl_array_accumulate(rhs, 1.0, src->Jrate_df);
  // Compute integrated diagnostics if needed.
  src->update_integrated_diags_rhs_func(app, species, src, 0.0);
  // Reset Jrate_df to Jrate*feq.
  gkyl_array_clear(src->Jrate_df, 0.0);
  gkyl_array_accumulate(src->Jrate_df, 1.0, species->lte.f_lte);
}

static void
gk_species_source_bgk_rhs_heating_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{

  struct timespec wst = gkyl_wall_clock();

  // Compute Maxwellian moments (n, u_par, T/m).
  gk_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, species->lte.moms.marr, 
    0, species->lte.moms.marr, 0, app->gk_geom->geo_int.jacobgeo, &app->local);  

  // Volume integrate Jrate times the thermal M2.
  gkyl_dg_mul_op_range(&app->basis, 0, src->Jrate_mom, 0, src->Jrate, 0, species->lte.moms.marr, &app->local);
  gkyl_dg_mul_op_range(&app->basis, 0, src->Jrate_mom, 0, src->Jrate_mom, 2, species->lte.moms.marr, &app->local);
  double Jrate_M2thermal_int = GKYL_MAX2(0.0, gk_species_source_bgk_volume_integrate(app, src, src->Jrate_mom));

  // Volume integrate Jrate times the vtsq_shape time M0.
  gkyl_dg_mul_op_range(&app->basis, 0, src->Jrate_mom, 0, src->Jrate_vtsq_shape, 0, species->lte.moms.marr, &app->local);
  double Jrate_vtsq_shape_M0_int = GKYL_MAX2(0.0, gk_species_source_bgk_volume_integrate(app, src, src->Jrate_mom));

  // Thermal speed squared of the Maxwellian.
  src->vtsq_amplitude = (src->norm_power + Jrate_M2thermal_int)/Jrate_vtsq_shape_M0_int;
  gkyl_array_set_offset_range(species->lte.moms.marr, src->vtsq_amplitude, src->vtsq_shape, 2*app->basis.num_basis, &app->local);

  // Compute the Maxwellian (this overwrite f_lte).
  gk_species_lte_from_moms(app, species, &species->lte, species->lte.moms.marr);

  // Multiply the Maxwellian by Jrate and keep it in lte.f_lte
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, species->lte.f_lte, 
    src->Jrate, species->lte.f_lte, &app->local, &species->local);

  // Assemble the BGK-like term and add it to rhs.
  gkyl_array_clear(src->Jrate_df, 0.0);
  gkyl_bgk_collisions_advance(src->bgk_op, &app->local, &species->local, 
    src->rate, species->lte.f_lte, fin, src->implicit_step, src->dt_implicit, src->Jrate_df, species->cflrate);
  gkyl_array_accumulate(rhs, 1.0, src->Jrate_df);

  app->stat.species_source_bgk_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_source_bgk_rhs_accumulate_maxwellian(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *out)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute Maxwellian moments (n, u_par, T/m).
  gk_species_moment_calc(&species->lte.moms, species->local, app->local, fin);
  gk_species_moment_diag_jacobgeo_div(app, &species->lte.moms, species->lte.moms.marr, species->lte.moms.marr);

  // Compute M0,M1,M2
  gk_species_moment_calc(&src->correct_mom_op, species->local, app->local, fin);
  gk_species_moment_diag_jacobgeo_div(app, &src->correct_mom_op, src->correct_mom_op.marr, src->correct_mom_op.marr);

  // Set a minimum on the density
  //gkyl_array_set_offset(src->Jrate_cap, src->damping_factor, species->lte.moms.marr, 0*app->basis.num_basis);

  // Divide M0dot by the rate and add M0
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, src->Jrate_mom, 0, src->M0dot, 0, src->rate, &app->local);  
  gkyl_array_accumulate_offset(src->Jrate_mom, 1.0, species->lte.moms.marr, 0*app->basis.num_basis);
  // Set the density
  gkyl_array_set_offset(species->lte.moms.marr, 1.0, src->Jrate_mom, 0*app->basis.num_basis);

  // Now do momentum
  // Divide M1dot by the rate, add on M1, divide by density to get upar
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, src->Jrate_mom, 0, src->M1dot, 0, src->rate, &app->local);  
  gkyl_array_accumulate_offset(src->Jrate_mom, 1.0, src->correct_mom_op.marr, 1*app->basis.num_basis);
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, src->Jrate_mom, 0, src->Jrate_mom, 0, species->lte.moms.marr, &app->local);  
  // Set the LTE moments for projection and project
  gkyl_array_set_offset(species->lte.moms.marr, 1.0, src->Jrate_mom, 1*app->basis.num_basis);

  // Now Do Energy
  // Set a minimum on vtsq so it doesn't go negative
  gkyl_array_set_offset(src->Jrate_cap, src->damping_factor, species->lte.moms.marr, 2*app->basis.num_basis);
  // T/m = M2dot/nu + M2 - n_s upar_s^2
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, src->Jrate_mom, 0, src->M2dot, 0, src->rate, &app->local);  
  gkyl_array_accumulate_offset(src->Jrate_mom, 1.0, src->correct_mom_op.marr, 2*app->basis.num_basis);

  gkyl_dg_mul_op_range(&app->basis, 1, src->correct_mom_op.marr, 1, species->lte.moms.marr, 1, species->lte.moms.marr, &app->local);
  gkyl_dg_mul_op_range(&app->basis, 1, src->correct_mom_op.marr, 1, src->correct_mom_op.marr, 0, species->lte.moms.marr, &app->local);

  gkyl_array_accumulate_offset(src->Jrate_mom,-1.0, src->correct_mom_op.marr, 1*app->basis.num_basis);

  gkyl_array_scale(src->Jrate_mom, 1.0/3.0);
  gkyl_dg_div_op_range(species->lte.moms.mem_geo, &app->basis, 0, src->Jrate_mom, 0, src->Jrate_mom, 0, species->lte.moms.marr, &app->local);  
  // Apply the cap so we don't drive the temperature negative
  gkyl_array_max_by_cell_per_cell_avg_range(src->Jrate_mom, src->Jrate_cap, &app->local);
  // Set the temperature
  gkyl_array_set_offset(species->lte.moms.marr, 1.0, src->Jrate_mom, 2*app->basis.num_basis);

  // Project with LTE moments
  gk_species_lte_from_moms(app, species, &species->lte, species->lte.moms.marr);

  // Add this Maxwellian onto the average of Maxwellians.
  gkyl_array_accumulate_range(out, 1.0, species->lte.f_lte, &species->local);

  app->stat.species_source_bgk_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_source_bgk_rhs_external_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkyl_array_clear(src->Jrate_df, 0.0);

  gk_species_source_bgk_rhs_accumulate_maxwellian(app, species, src, fin, src->Jrate_df);

  // Multiply the Maxwellian by Jrate.
  gkyl_dg_mul_conf_phase_op_range(&app->basis, &species->basis, species->lte.f_lte, 
    src->Jrate, src->Jrate_df, &app->local, &species->local);
  // Assemble the BGK-like term and add it to Jrate_feq to conserve -nu*(f-feq) term for use in diagnostics.
  gkyl_array_clear(src->Jrate_df, 0.0);
  gkyl_bgk_collisions_advance(src->bgk_op, &app->local, &species->local, src->rate, species->lte.f_lte,
    fin, src->implicit_step, src->dt_implicit, src->Jrate_df, species->cflrate);
  gkyl_array_accumulate(rhs, 1.0, src->Jrate_df);
}

static void
gk_species_source_bgk_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, double tm, int frame)
{
}

static void
gk_species_source_bgk_write_diags_heating_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // Write the Maxwellian square thermal speed amplitude.
  gkyl_dynvec_append(src->vtsq_amp_diag, tm, &src->vtsq_amplitude);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    const char *fmt = "%s-%s_source_bgk_vtsq_amplitude.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name);
    
    if (src->is_first_diag_dynvec_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Squared thermal speed amplitude." }
      };
      int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(src->vtsq_amp_diag, fileNm, mt);
      src->is_first_diag_dynvec_write_call = false;
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(src->vtsq_amp_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(src->vtsq_amp_diag);
  app->stat.n_diag_io += 1;
  
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_source_bgk_write_diags_external_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, double tm, int frame)
{
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_conf_len, gks->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_conf_len, gks->io_meta_conf, "frame", frame);
  struct gkyl_msgpack_map_elem desc_bgk_moms[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING,
      .cval = "BGK source particle (M0), momentum (M1) or kinetic energy (M2) source/sink rate." }
  };
  int io_meta_len[] = {gks->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_conf, app->gk_geom->io_meta_basic, desc_bgk_moms};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  if (app->use_gpu) {
    gkyl_array_copy(src->M0dot_host, src->M0dot);
    gkyl_array_copy(src->M1dot_host, src->M1dot);
    gkyl_array_copy(src->M2dot_host, src->M2dot);
  }

  cstr fileNm = cstr_from_fmt("%s-%s_source_bgk_M0dot_%d.gkyl", app->name, gks->info.name, frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, src->M0dot_host, fileNm.str);
  cstr_drop(&fileNm);
  fileNm = cstr_from_fmt("%s-%s_source_bgk_M1dot_%d.gkyl", app->name, gks->info.name, frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, src->M1dot_host, fileNm.str);
  cstr_drop(&fileNm);
  fileNm = cstr_from_fmt("%s-%s_source_bgk_M2dot_%d.gkyl", app->name, gks->info.name, frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, src->M2dot_host, fileNm.str);
  cstr_drop(&fileNm);

  gkyl_msgpack_data_release(mt); 
}

static void
gk_species_source_bgk_update_integrated_diags_disabled(gkyl_gyrokinetic_app* app,
  struct gk_species *gks, struct gk_source_bgk *src, double tm)
{
  // Do nothing.
}

static void
gk_species_source_bgk_update_integrated_diags_enabled(gkyl_gyrokinetic_app* app,
  struct gk_species *gks, struct gk_source_bgk *src, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  int num_mom = src->integ_mom_op.num_mom;
  double avals_global[num_mom];
  
  // Compute integrated moments of source term.
  gk_species_moment_calc(&src->integ_mom_op, gks->local, app->local, src->Jrate_df); 
  app->stat.n_mom += 1;

  // Reduce (sum) over whole domain, append to diagnostics.
  gkyl_array_reduce_range(src->red_integ_diag, src->integ_mom_op.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
    src->red_integ_diag, src->red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, src->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, src->red_integ_diag_global, sizeof(double[num_mom]));
  }
  for (int i=0; i<num_mom; ++i) {
    src->int_mom_global[i] = avals_global[i];
  }
}

static void
gk_species_source_bgk_calc_integrated_diags_disabled(gkyl_gyrokinetic_app* app,
  struct gk_species *gks, struct gk_source_bgk *src, double tm)
{
  // Do nothing.
}

static void
gk_species_source_bgk_calc_integrated_diags_enabled(gkyl_gyrokinetic_app* app,
  struct gk_species *gks, struct gk_source_bgk *src, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  src->update_integrated_diags_func(app, gks, src, tm);
  gkyl_dynvec_append(src->integ_diag, tm, src->int_mom_global);
  
  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

static void
gk_species_source_bgk_write_integrated_diags_disabled(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_source_bgk *src)
{
  // Empty.
}

static void
gk_species_source_bgk_write_integrated_diags_enabled(gkyl_gyrokinetic_app *app,
  struct gk_species *gks, struct gk_source_bgk *src)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s_source_bgk_integrated_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, gkyl_distribution_moments_strs[GKYL_F_MOMENT_M0M1M2]);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, gkyl_distribution_moments_strs[GKYL_F_MOMENT_M0M1M2]);

    if (src->is_first_diag_dynvec_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Volume integrated moment of the BGK source." }
      };
      int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(src->integ_diag, fileNm, mt);
      src->is_first_diag_dynvec_write_call = false;
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(src->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(src->integ_diag);
  app->stat.n_diag_io += 1;
  
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_source_bgk_write_array(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, int frame, double stime, char* file_suffix, char* description,
  struct gkyl_msgpack_map_elem *iom, int iom_len, struct gkyl_rect_grid grid, struct gkyl_range local,
  struct gkyl_array *arrout)
{
  // Write out a conf-space or a phase-space array.

  // Package metadata.
  gkyl_msgpack_map_elem_set_double(iom_len, iom, "time", stime);
  gkyl_msgpack_map_elem_set_uint(iom_len, iom, "frame", frame);
  struct gkyl_msgpack_map_elem io_meta_f[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = description},
  };
  int io_meta_f_len = sizeof(io_meta_f)/sizeof(io_meta_f[0]);
  int io_meta_len[] = {iom_len, app->gk_geom->io_meta_basic_len, io_meta_f_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {iom, app->gk_geom->io_meta_basic, io_meta_f};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Construct the file handles for collision frequency and primitive moments.
  const char *fmt = "%s-%s_%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, file_suffix, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, file_suffix, frame);

  struct gkyl_array *arr_ho;
  if (app->use_gpu) {  
    arr_ho = mkarr(false, arrout->ncomp, arrout->size);
    // Copy data from device to host before writing it out.
    gkyl_array_copy(arr_ho, arrout);
  }
  else {
    arr_ho = gkyl_array_acquire(arrout);
  }

  gkyl_comm_array_write(app->comm, &grid, &local, mt, arr_ho, fileNm);
  gkyl_msgpack_data_release(mt); 
  gkyl_array_release(arr_ho);
}

void 
gk_species_source_bgk_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_source_bgk *src)
{
  src->source_bgk_id = src->source_bgk_id ? src->source_bgk_id : gks->info.source_bgk.source_bgk_id;
  src->write_diagnostics = src->write_diagnostics ? src->write_diagnostics : gks->info.source_bgk.write_diagnostics;

  src->write_diags_func = gk_species_source_bgk_write_diags_disabled;
  src->rhs_func = gk_species_source_bgk_rhs_disabled;
  src->update_integrated_diags_rhs_func = gk_species_source_bgk_update_integrated_diags_disabled;
  src->update_integrated_diags_func = gk_species_source_bgk_update_integrated_diags_disabled;
  src->calc_integrated_diags_func = gk_species_source_bgk_calc_integrated_diags_disabled;
  src->write_integrated_diags_func = gk_species_source_bgk_write_integrated_diags_disabled;

  if (src->source_bgk_id){
    // Collision rate
    src->rate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    // Jacobian times the collision rate.
    src->Jrate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    // Jrate times the equilibrium distribution function.
    src->Jrate_df = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
    // BGK operator.
    src->bgk_op = gkyl_bgk_collisions_new(&app->basis, &gks->basis, app->use_gpu);
    src->implicit_step = false;
    src->dt_implicit = 1e9;
    // Diagnostics.
    if (src->write_diagnostics) {
      src->update_integrated_diags_func = gk_species_source_bgk_update_integrated_diags_enabled;
      src->calc_integrated_diags_func = gk_species_source_bgk_calc_integrated_diags_enabled;
      src->write_integrated_diags_func = gk_species_source_bgk_write_integrated_diags_enabled;
      gk_species_moment_init(app, gks, &src->integ_mom_op, GKYL_F_MOMENT_M0M1M2, true);
      if (app->use_gpu) {
        src->red_integ_diag = gkyl_cu_malloc(sizeof(double[src->integ_mom_op.num_mom]));
        src->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[src->integ_mom_op.num_mom]));
      } 
      else {
        src->red_integ_diag = gkyl_malloc(sizeof(double[src->integ_mom_op.num_mom]));
        src->red_integ_diag_global = gkyl_malloc(sizeof(double[src->integ_mom_op.num_mom]));
      }
      src->int_mom_global = gkyl_malloc(sizeof(double[src->integ_mom_op.num_mom]));
      src->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, src->integ_mom_op.num_mom);
      src->is_first_diag_dynvec_write_call = true;
    }

    if (src->source_bgk_id == GKYL_SOURCE_BGK_STATIC) {
      // source_bgk rate.
      struct gkyl_array *rate_host = app->use_gpu? mkarr(false, src->rate->ncomp, src->rate->size)
                                                : gkyl_array_acquire(src->rate);
      gkyl_proj_on_basis *proj_rate = gkyl_proj_on_basis_new(&app->grid, &app->basis,
        app->poly_order+1, 1, gks->info.source_bgk.rate_profile, gks->info.source_bgk.rate_profile_ctx);

      gkyl_proj_on_basis_advance(proj_rate, 0.0, &app->local, rate_host);
      gkyl_array_copy(src->rate, rate_host);
      gkyl_proj_on_basis_release(proj_rate);
      gkyl_array_release(rate_host);
      // Multiply the rate by the conf-space Jacobian.
      gkyl_dg_mul_op_range(&app->basis, 0, src->Jrate, 0, app->gk_geom->geo_int.jacobgeo, 0, src->rate, &app->local);

      // Project the user provided function (we use here Jrate_df as a temporary array).
      struct gkyl_array *Jrate_fmax_host = app->use_gpu? mkarr(false, src->Jrate_df->ncomp, src->Jrate_df->size)
                                                      : gkyl_array_acquire(src->Jrate_df);
      struct gk_proj_on_basis_c2p_func_ctx proj_feq_shape_c2p_ctx = {
        .cdim = app->cdim,
        .vdim = gks->local_vel.ndim,
        .vel_map = gks->vel_map,
        .pos_map = app->position_map,
      };
      gkyl_proj_on_basis *proj_feq_shape = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp) {
          .grid = &gks->grid,
          .basis = &gks->basis,
          .qtype = GKYL_GAUSS_QUAD,
          .num_quad = gks->basis.poly_order+1,
          .num_ret_vals = 1,
          .eval = gks->info.source_bgk.feq_shape,
          .ctx = gks->info.source_bgk.feq_shape_ctx,
          .c2p_func = proj_on_basis_c2p_phase_func,
          .c2p_func_ctx = &proj_feq_shape_c2p_ctx,
        }
      );
      gkyl_proj_on_basis_advance(proj_feq_shape, 0.0, &gks->local, Jrate_fmax_host);
      gkyl_array_copy(src->Jrate_df, Jrate_fmax_host);
      gkyl_proj_on_basis_release(proj_feq_shape);
      gkyl_array_release(Jrate_fmax_host);

      // Multiply feq by the same Jacobians as the distribution function.
      // (Similar to gk_species_projection_calc.)
      gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, src->Jrate_df,
        app->gk_geom->geo_corn.bmag, src->Jrate_df, &app->local, &gks->local);
      gkyl_array_scale_by_cell(src->Jrate_df, gks->vel_map->jacobvel);
      gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, src->Jrate_df,
        app->gk_geom->geo_int.jacobgeo, src->Jrate_df, &app->local, &gks->local);

      // Methods chosen at runtime.
      src->rhs_func = gk_species_source_bgk_rhs_feq_enabled;
      if (src->write_diagnostics) {
        // Replace the update int. diags func. with the one that is called within the RHS.
        src->update_integrated_diags_func = gk_species_source_bgk_update_integrated_diags_disabled;
        src->update_integrated_diags_rhs_func = gk_species_source_bgk_update_integrated_diags_enabled;
        gk_species_source_bgk_write_array(app, gks, src, 0, 0.0, "source_bgk_rate", "BGK source relaxation rate.", 
          gks->io_meta_conf, gks->io_meta_conf_len, app->grid, app->local, src->rate);
        gk_species_source_bgk_write_array(app, gks, src, 0, 0.0, "source_bgk_feq", "BGK source equilibrium function.", 
          gks->io_meta_phase, gks->io_meta_phase_len, gks->grid, gks->local, src->Jrate_df);
      }

      // Multiply feq by the rate to get Jrate_df = nu*feq.
      gkyl_dg_mul_conf_phase_op_range(&app->basis, &gks->basis, src->Jrate_df,
        src->rate, src->Jrate_df, &app->local, &gks->local);
    }

    if (src->source_bgk_id == GKYL_SOURCE_BGK_HEATING) {
      int vdim_phys = gks->info.vdim == 1? 1 : 3;
      src->norm_power = 2.0*gks->info.source_bgk.power/(vdim_phys*gks->info.mass);

      // source_bgk rate.
      struct gkyl_array *rate_host = app->use_gpu? mkarr(false, src->rate->ncomp, src->rate->size)
                                                : gkyl_array_acquire(src->rate);
      gkyl_proj_on_basis *proj_rate = gkyl_proj_on_basis_new(&app->grid, &app->basis,
        app->poly_order+1, 1, gks->info.source_bgk.rate_profile, gks->info.source_bgk.rate_profile_ctx);

      gkyl_proj_on_basis_advance(proj_rate, 0.0, &app->local, rate_host);
      gkyl_array_copy(src->rate, rate_host);
      gkyl_proj_on_basis_release(proj_rate);
      gkyl_array_release(rate_host);
      // Multiply the rate by the conf-space Jacobian.
      gkyl_dg_mul_op_range(&app->basis, 0, src->Jrate, 0, app->gk_geom->geo_int.jacobgeo, 0, src->rate, &app->local);

      // source_bgk rate.
      src->vtsq_shape = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      struct gkyl_array *vtsq_shape_host = app->use_gpu? mkarr(false, src->vtsq_shape->ncomp, src->vtsq_shape->size)
                                                      : gkyl_array_acquire(src->vtsq_shape);
      gkyl_proj_on_basis *proj_vtsq_shape = gkyl_proj_on_basis_new(&app->grid, &app->basis,
        app->poly_order+1, 1, gks->info.source_bgk.temp_shape, gks->info.source_bgk.temp_shape_ctx);
      gkyl_proj_on_basis_advance(proj_vtsq_shape, 0.0, &app->local, vtsq_shape_host);
      gkyl_array_copy(src->vtsq_shape, vtsq_shape_host);
      gkyl_proj_on_basis_release(proj_vtsq_shape);
      gkyl_array_release(vtsq_shape_host);

      // Multiply Jrate by the shape of v_t^2.
      src->Jrate_vtsq_shape = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      gkyl_dg_mul_op_range(&app->basis, 0, src->Jrate_vtsq_shape, 0, src->Jrate, 0, src->vtsq_shape, &app->local);

      // Rate times a velocity moment.
      src->Jrate_mom = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      // Volume integrator.
      src->vol_integ_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
      if (app->use_gpu) {
        src->volint_local = gkyl_cu_malloc(sizeof(double));
        src->volint_global = gkyl_cu_malloc(sizeof(double));
      } 
      else {
        src->volint_local = gkyl_malloc(sizeof(double));
        src->volint_global = gkyl_malloc(sizeof(double));
      }

      if (src->write_diagnostics) {
        src->vtsq_amp_diag = gkyl_dynvec_new(GKYL_DOUBLE, 1);
        // Write out the source_bgk rate and vtsq shape.
        gk_species_source_bgk_write_array(app, gks, src, 0, 0.0, "source_bgk_rate", "BGK source relaxation rate.", 
          gks->io_meta_conf, gks->io_meta_conf_len, app->grid, app->local, src->rate);
        gk_species_source_bgk_write_array(app, gks, src, 0, 0.0, "source_bgk_temp_shape", "BGK source temperature shape.", 
          gks->io_meta_conf, gks->io_meta_conf_len, app->grid, app->local, src->vtsq_shape);
      }

      // Methods chosen at runtime.
      src->rhs_func = gk_species_source_bgk_rhs_heating_enabled;
      if (src->write_diagnostics) {
        src->write_diags_func = gk_species_source_bgk_write_diags_heating_enabled;
      }
    }

    if (src->source_bgk_id == GKYL_SOURCE_BGK_EXTERNAL) {
      // source_bgk rate.
      gkyl_array_shiftc(src->rate, pow(sqrt(2.0),app->cdim)/src->injection_time, 0); // Sets rate = 1/injection_time

      // Multiply the rate by the conf-space Jacobian.
      gkyl_dg_mul_op_range(&app->basis, 0, src->Jrate, 0, app->gk_geom->geo_int.jacobgeo, 0, src->rate, &app->local);

      // External source rates
      src->M0dot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      src->M1dot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      src->M2dot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      src->M0dot_host =  src->M0dot;
      src->M1dot_host =  src->M1dot;
      src->M2dot_host =  src->M2dot;

      if(app->use_gpu) {
        src->M0dot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
        src->M1dot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
        src->M2dot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      }

      // Rate times a velocity moment.
      src->Jrate_mom = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      src->Jrate_cap = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

      // Correction moments of the operator.
      gk_species_moment_init(app, gks, &src->correct_mom_op, GKYL_F_MOMENT_M0M1M2, false);

      // Methods chosen at runtime.
      src->rhs_func = gk_species_source_bgk_rhs_external_enabled;
      if (src->write_diagnostics) {
        src->write_diags_func = gk_species_source_bgk_write_diags_external_enabled;
      }
    }
  }
}

void
gk_species_source_bgk_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_source_bgk *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  src->rhs_func(app, species, src, fin, rhs);
}

void
gk_species_source_bgk_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, double tm, int frame)
{
  src->write_diags_func(app, gks, src, tm, frame);
}

void
gk_species_source_bgk_calc_integrated_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_source_bgk *src, double tm)
{
  src->calc_integrated_diags_func(app, gks, src, tm);
}

void
gk_species_source_bgk_write_integrated_diags(gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_source_bgk *src)
{
  src->write_integrated_diags_func(app, gks, src);
}

void
gk_species_source_bgk_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source_bgk *src)
{
  if (src->source_bgk_id) {
    gkyl_array_release(src->rate);
    gkyl_array_release(src->Jrate);
    gkyl_array_release(src->Jrate_df);
    gkyl_bgk_collisions_release(src->bgk_op);

    if (src->write_diagnostics) {
      gk_species_moment_release(app, &src->integ_mom_op); 
      gkyl_dynvec_release(src->integ_diag);
      if (app->use_gpu) {
        gkyl_cu_free(src->red_integ_diag);
        gkyl_cu_free(src->red_integ_diag_global);
      }
      else {
        gkyl_free(src->red_integ_diag);
        gkyl_free(src->red_integ_diag_global);
      }
      gkyl_free(src->int_mom_global);
    }

    if (src->source_bgk_id == GKYL_SOURCE_BGK_HEATING) {
      gkyl_array_release(src->vtsq_shape);
      gkyl_array_release(src->Jrate_vtsq_shape);
      gkyl_array_release(src->Jrate_mom);

      gkyl_array_integrate_release(src->vol_integ_op);
      if (app->use_gpu) {
        gkyl_cu_free(src->volint_local);
        gkyl_cu_free(src->volint_global);
      }
      else {
        gkyl_free(src->volint_local);
        gkyl_free(src->volint_global);
      }

      if (src->write_diagnostics) {
        gkyl_dynvec_release(src->vtsq_amp_diag);
      }
    } 
    else if (src->source_bgk_id == GKYL_SOURCE_BGK_EXTERNAL) {
      gkyl_array_release(src->Jrate_mom);
      gkyl_array_release(src->Jrate_cap);

      gkyl_array_release(src->M0dot);
      gkyl_array_release(src->M1dot);
      gkyl_array_release(src->M2dot);
      if (app->use_gpu) {
        gkyl_array_release(src->M0dot_host);
        gkyl_array_release(src->M1dot_host);
        gkyl_array_release(src->M2dot_host);
      }

      gk_species_moment_release(app, &src->correct_mom_op); 
    }
  }
}
