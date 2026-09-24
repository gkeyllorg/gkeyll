#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_anomalous_diff_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

static void
gk_species_anomalous_diff_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_dg_updater_gk_anomalous_diffusion_advance(gkad->slvr, &species->local, 
    fin, species->cflrate, rhs);

  app->stat.species_diffusion_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_anomalous_diff_write_diags_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, double tm, int frame)
{
}

static void
gk_species_anomalous_diff_write_diags_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_anomalous_diff_write_conf_array(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, int frame, double stime, char* file_suffix, char* description,
  struct gkyl_array *arrout, struct gkyl_array *arrout_host)
{
  // Write out a conf-space array.
  
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_conf_len, gks->io_meta_conf, "time", stime);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_conf_len, gks->io_meta_conf, "frame", frame);
  struct gkyl_msgpack_map_elem desc[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = description }
  };
  int io_meta_len[] = {gks->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_conf, app->gk_geom->io_meta_basic, desc};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Construct the file handles for collision frequency and primitive moments.
  const char *fmt = "%s-%s_%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, file_suffix, frame);
  char fileNm[sz+1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, file_suffix, frame);

  struct gkyl_array *arr_ho;
  if (app->use_gpu) {  
    if (arrout_host)
      arr_ho = gkyl_array_acquire(arrout_host);
    else {
      arr_ho = mkarr(false, arrout->ncomp, arrout->size);
    }
    // Copy data from device to host before writing it out.
    gkyl_array_copy(arr_ho, arrout);
  }
  else {
    if (arrout_host)
      arr_ho = gkyl_array_acquire(arrout_host);
    else
      arr_ho = gkyl_array_acquire(arrout);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, arr_ho, fileNm);
  gkyl_msgpack_data_release(mt); 
  gkyl_array_release(arr_ho);
}

static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}

void 
gk_species_anomalous_diff_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks, 
  struct gk_anomalous_diff *gkad)
{
  gkad->anom_diff_id = gks->info.anomalous_diffusion.anomalous_diff_id;
  gkad->write_diagnostics = gks->info.anomalous_diffusion.write_diagnostics;

  gkad->write_diags_func = gk_species_anomalous_diff_write_diags_disabled;
  gkad->rhs_func = gk_species_anomalous_diff_rhs_disabled;

  if (gkad->anom_diff_id) {

    // Project the diffusivity.
    gkad->diffD = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *diffD_ho = app->use_gpu? mkarr(false, gkad->diffD->ncomp, gkad->diffD->size)
                                              : gkyl_array_acquire(gkad->diffD);
    struct gkyl_eval_on_nodes *diffDproj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .num_ret_vals = 1,
      .eval = gks->info.anomalous_diffusion.D_profile,
      .ctx = gks->info.anomalous_diffusion.D_profile_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
    });
    gkyl_eval_on_nodes_advance(diffDproj, 0.0, &app->local, diffD_ho);
    gkyl_array_copy(gkad->diffD, diffD_ho);

    gkyl_eval_on_nodes_release(diffDproj);
    gkyl_array_release(diffD_ho);

    // Multiply diffD by g^xx*jacobgeo.
    gkyl_dg_mul_op(&app->basis, 0, gkad->diffD, 0, app->gk_geom->geo_int.gxxj, 0, gkad->diffD);

    // Sync diffusivity.
    int num_periodic_dir = app->num_periodic_dir;
    gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
      num_periodic_dir, app->periodic_dirs, gkad->diffD); 
    // ABSORB and FIXED_FUNC BCs need to fill the ghost cell. MF 2025/11/03: The
    // only sensible way to fill the ghost cell is to use the value at the
    // boundary.
    for (int b=0; b<2; ++b) {
      if ((b == 0 && ((gks->lower_bc[0].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) || (gks->lower_bc[0].type == GKYL_BC_GK_SPECIES_ABSORB))) || 
          (b == 1 && ((gks->upper_bc[0].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) || (gks->upper_bc[0].type == GKYL_BC_GK_SPECIES_ABSORB))) ) {
        int dir = 0;
        enum gkyl_edge_loc edge = b==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
        struct gkyl_range *skin_r = b==0? &app->local_lower_skin[dir] : &app->local_upper_skin[dir];
        struct gkyl_range *ghost_r = b==0? &app->local_lower_ghost[dir] : &app->local_upper_ghost[dir];

        long vol = skin_r->volume;
        long buff_sz = 1;
        buff_sz = buff_sz > vol ? buff_sz : vol;
        struct gkyl_array *bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

        struct gkyl_bc_basic_gyrokinetic *gfss_bc_op = gkyl_bc_basic_gyrokinetic_new(dir, edge,
          GKYL_BC_GK_FIELD_BOUNDARY_VALUE, app->basis_on_dev, skin_r, ghost_r, app->basis.num_basis, app->cdim, app->use_gpu);
        gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, bc_buffer, gkad->diffD);
        gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, bc_buffer, app->gk_geom->geo_int.jacobgeo_inv);

        gkyl_bc_basic_gyrokinetic_release(gfss_bc_op);
        gkyl_array_release(bc_buffer);
      }
    }
    gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, gkad->diffD);

    // Create solver.
    gkad->slvr = gkyl_dg_updater_gk_anomalous_diffusion_new(&gks->grid, &gks->basis, &app->basis,
      &app->local, gks->lower_bc[0].type, gks->upper_bc[0].type,
      gkad->diffD, app->gk_geom->geo_int.jacobgeo_inv, app->use_gpu);

    if (gkad->write_diagnostics) {
      // Write out the diffusivity.
      gk_anomalous_diff_write_conf_array(app, gks, gkad, 0, 0.0, "anom_diff", "Anomalous diffusivity.", gkad->diffD, 0);
    }

    // Methods chosen at runtime.
    gkad->rhs_func = gk_species_anomalous_diff_rhs_enabled;
    if (gkad->write_diagnostics) {
      gkad->write_diags_func = gk_species_anomalous_diff_write_diags_enabled;
    }
  }
}

void
gk_species_anomalous_diff_rhs(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_anomalous_diff *gkad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  gkad->rhs_func(app, species, gkad, fin, rhs);
}

void
gk_species_anomalous_diff_write_diags(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_anomalous_diff *gkad, double tm, int frame)
{
  gkad->write_diags_func(app, gks, gkad, tm, frame);
}

void
gk_species_anomalous_diff_release(const struct gkyl_gyrokinetic_app *app, const struct gk_anomalous_diff *gkad)
{
  if (gkad->anom_diff_id) {

    gkyl_array_release(gkad->diffD);
    gkyl_dg_updater_gk_anomalous_diffusion_release(gkad->slvr);

//    if (gkad->write_diagnostics) {
//    }
  }
}
