#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gkns_pos_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, double tm, int frame)
{
  // Empty.
}

static void
gkns_pos_deltaf_moms_calc(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos)
{
  gk_neut_species_moment_calc(&pos->moms, gkns->local, app->local, pos->fbuffer_ptr);

  // Rescale moment by inverse of Jacobian if needed.
  gk_species_moment_diag_jacobgeo_div(app, &pos->moms, pos->moms.marr, pos->moms.marr);
}

static void
gkns_pos_deltaf_moms_clear(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos)
{
  gkyl_array_clear(pos->moms.marr, 0.0);
  pos->deltaf_moms_func_neut = gkns_pos_deltaf_moms_calc;
}

static void
gkns_pos_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, double tm, int frame)
{
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gkns->io_meta_conf_len, gkns->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gkns->io_meta_conf_len, gkns->io_meta_conf, "frame", frame);
  struct gkyl_msgpack_map_elem desc[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING,
      .cval = "M0M1M2 moments of the change in the distribution by the positivity shift." }
  };
  int io_meta_len[] = {gkns->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_conf, app->gk_geom->io_meta_basic, desc};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  struct timespec wst = gkyl_wall_clock();
  // We placed the change in f in fbuffer_ptr.
  pos->deltaf_moms_func_neut(app, gkns, pos);
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

  if (app->use_gpu)
    gkyl_array_copy(pos->moms.marr_host, pos->moms.marr);

  const char *fmt = "%s-%s_positivity_%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name,
    gkyl_distribution_moments_strs[GKYL_F_MOMENT_M0M1M2], frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name,
    gkyl_distribution_moments_strs[GKYL_F_MOMENT_M0M1M2], frame);
  
  struct timespec wtm = gkyl_wall_clock();
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
    pos->moms.marr_host, fileNm);
  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.n_neut_diag_io += 1;
  
  gkyl_msgpack_data_release(mt); 

  app->stat.n_neut_diag += 1;
}

static void
gkns_pos_calc_integrated_diags_disabled(gkyl_gyrokinetic_app* app,
  struct gk_neut_species *gkns, struct gk_positivity *pos, double tm)
{
  // Empty.
}

static void
gkns_pos_deltaf_integ_moms_calc(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos)
{
  gk_neut_species_moment_calc(&pos->integ_moms, gkns->local, app->local, pos->fbuffer_ptr); 
}

static void
gkns_pos_deltaf_integ_moms_clear(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos)
{
  gkyl_array_clear(pos->integ_moms.marr, 0.0);
  pos->deltaf_integ_moms_func_neut = gkns_pos_deltaf_integ_moms_calc;
}

static void
gkns_pos_calc_integrated_diags_enabled(gkyl_gyrokinetic_app* app,
  struct gk_neut_species *gkns, struct gk_positivity *pos, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  int num_mom = pos->integ_moms.num_mom;
  double avals_global[num_mom];
  
  // The change in f after the positivity operation is stored in fbuffer_ptr.
  pos->deltaf_integ_moms_func_neut(app, gkns, pos);
  app->stat.n_mom += 1;

  // Reduce (sum) over whole domain, append to diagnostics.
  gkyl_array_reduce_range(pos->red_integ_diag, pos->integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
    pos->red_integ_diag, pos->red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, pos->red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, pos->red_integ_diag_global, sizeof(double[num_mom]));
  }
  gkyl_dynvec_append(pos->integ_diag, tm, avals_global);
  
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_neut_diag += 1;
}

static void
gkns_pos_write_integrated_diags_disabled(gkyl_gyrokinetic_app *app,
  struct gk_neut_species *gkns, struct gk_positivity *pos)
{
  // Empty.
}

static void
gkns_pos_write_integrated_diags_enabled(gkyl_gyrokinetic_app *app,
  struct gk_neut_species *gkns, struct gk_positivity *pos)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s_positivity_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, "integrated_moms");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, "integrated_moms");

    if (pos->is_first_integ_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Volume integrated moments of change in distribution due to positivity shift." }
      };
      int io_meta_len[] = {gkns->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(pos->integ_diag, fileNm, mt);
      pos->is_first_integ_write_call = false;
      pos->integ_diag_file_exists = true;
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(pos->integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(pos->integ_diag);
  app->stat.n_neut_diag_io += 1;

  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gkns_pos_apply_disabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  struct gk_positivity *pos, struct gkyl_array *fbuffer, struct gkyl_array *fout)
{
  // Do nothing.
}

void
gkns_pos_apply_enabled(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, struct gkyl_array *fbuffer, struct gkyl_array *fout)
{
  struct timespec wtm = gkyl_wall_clock();
  // Copy f so we can calculate the moments of delta f later. 
  pos->fbuffer_ptr = fbuffer;
  gkyl_array_set(pos->fbuffer_ptr, -1.0, fout);

  // Shift each species.
  gkyl_positivity_shift_vlasov_advance(pos->shift_op_vlasov, &app->local, &gkns->local,
    fout, gkns->m0.marr, pos->delta_m0);

  app->stat.species_pos_shift_tm += gkyl_time_diff_now_sec(wtm);
}

void
gk_neut_species_positivity_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_positivity *pos)
{
  pos->type = gkns->info.positivity.type;
  pos->write_diagnostics = gkns->info.positivity.write_diagnostics;
  pos->quasineut_rescale = false;

  pos->apply_func_neut = gkns_pos_apply_disabled;
  pos->write_diags_func_neut = gkns_pos_write_diags_disabled;
  pos->calc_integrated_diags_func_neut = gkns_pos_calc_integrated_diags_disabled;
  pos->write_integrated_diags_func_neut = gkns_pos_write_integrated_diags_disabled;

  if (pos->type) {

    pos->delta_m0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    // Positivity shift updater.
    pos->shift_op_vlasov = gkyl_positivity_shift_vlasov_new(app->basis, gkns->basis,
      gkns->grid, &app->local_ext, app->use_gpu);

    // Methods chosen at runtime.
    pos->apply_func_neut = gkns_pos_apply_enabled;

    if (pos->write_diagnostics) {
      // Allocate data for diagnostic moments.
      gk_neut_species_moment_init(app, gkns, &pos->moms, GKYL_F_MOMENT_M0M1M2, false);

      // Integrated moments of delta f.
      gk_neut_species_moment_init(app, gkns, &pos->integ_moms, GKYL_F_MOMENT_M0M1M2, true);

      if (app->use_gpu) {
        pos->red_integ_diag = gkyl_cu_malloc(sizeof(double[pos->integ_moms.num_mom]));
        pos->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[pos->integ_moms.num_mom]));
      } 
      else {
        pos->red_integ_diag = gkyl_malloc(sizeof(double[pos->integ_moms.num_mom]));
        pos->red_integ_diag_global = gkyl_malloc(sizeof(double[pos->integ_moms.num_mom]));
      }
      pos->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, pos->integ_moms.num_mom);
      pos->is_first_integ_write_call = true;

      // Methods chosen at runtime.
      pos->deltaf_moms_func_neut = gkns_pos_deltaf_moms_clear;
      pos->deltaf_integ_moms_func_neut = gkns_pos_deltaf_integ_moms_clear;
      pos->write_diags_func_neut = gkns_pos_write_diags_enabled;
      pos->calc_integrated_diags_func_neut = gkns_pos_calc_integrated_diags_enabled;
      pos->write_integrated_diags_func_neut = gkns_pos_write_integrated_diags_enabled;
    }

  }
}

void
gk_neut_species_positivity_apply(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, struct gkyl_array *fbuffer, struct gkyl_array *fout)
{
  pos->apply_func_neut(app, gkns, pos, fbuffer, fout);
}

void
gk_neut_species_positivity_write_diags(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, double tm, int frame)
{
  pos->write_diags_func_neut(app, gkns, pos, tm, frame);
}

void
gk_neut_species_positivity_calc_integrated_diags(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns,
  struct gk_positivity *pos, double tm)
{
  pos->calc_integrated_diags_func_neut(app, gkns, pos, tm);
}

void
gk_neut_species_positivity_write_integrated_diags(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gk_positivity *pos)
{
  pos->write_integrated_diags_func_neut(app, gkns, pos);
}

void
gk_neut_species_positivity_release(const struct gkyl_gyrokinetic_app *app, const struct gk_positivity *pos)
{
  if (pos->type) {

    gkyl_array_release(pos->delta_m0);
    gkyl_positivity_shift_vlasov_release(pos->shift_op_vlasov);

    if (pos->write_diagnostics) {
      gk_neut_species_moment_release(app, &pos->moms);
      gk_neut_species_moment_release(app, &pos->integ_moms); 
      gkyl_dynvec_release(pos->integ_diag);
      if (app->use_gpu) {
        gkyl_cu_free(pos->red_integ_diag);
        gkyl_cu_free(pos->red_integ_diag_global);
      }
      else {
        gkyl_free(pos->red_integ_diag);
        gkyl_free(pos->red_integ_diag_global);
      }
    }

  }
}

void
gk_neut_species_positivity_reset(gkyl_gyrokinetic_app* app, double tm,
  struct gk_neut_species *gkns, struct gk_positivity *pos, struct gkyl_gyrokinetic_positivity pos_inp)
{
  gk_neut_species_positivity_release(app, pos);

  gkns->info.positivity = pos_inp;
  gk_neut_species_positivity_init(app, gkns, pos);

  if (pos_inp.type) {
    if (pos->integ_diag_file_exists) {
      pos->is_first_integ_write_call = false;
      pos->deltaf_integ_moms_func_neut = gkns_pos_deltaf_integ_moms_calc;
    }
  }
}
