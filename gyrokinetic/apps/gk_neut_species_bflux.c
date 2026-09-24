#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_translate_dim.h>

static int
gk_neut_species_bflux_idx(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge) {
  // Given a direction 'dir' and an edge 'edge' return the boundary index.
  for (int b=0; b<bflux->num_boundaries; ++b) {
    if (dir == bflux->boundaries_dir[b] && edge == bflux->boundaries_edge[b])
      return b;
  }
  return -1;
}

static void
gk_neut_species_bflux_clear_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    for (int m=0; m<bflux->num_calc_moms; m++)
      gkyl_array_clear_range(fin[b*bflux->num_calc_moms+m], val, bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_clear_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  // Do nothing.
}

void
gk_neut_species_bflux_clear(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  bflux->bflux_clear_func(app, bflux, fin, val);
}

static void
gk_neut_species_bflux_scale_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    for (int m=0; m<bflux->num_calc_moms; m++)
      gkyl_array_scale_range(fin[b*bflux->num_calc_moms+m], val, bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_scale_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  // Do nothing.
}

void
gk_neut_species_bflux_scale(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fin, double val)
{
  bflux->bflux_scale_func(app, bflux, fin, val);
}

static void
gk_neut_species_bflux_step_f_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double dt, struct gkyl_array **fin)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    for (int m=0; m<bflux->num_calc_moms; m++)
      gkyl_array_accumulate_range(gkyl_array_scale_range(fout[b*bflux->num_calc_moms+m], dt, bflux->boundaries_conf_ghost[b]),
        1.0, fin[b*bflux->num_calc_moms+m], bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_step_f_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double dt, struct gkyl_array **fin)
{
  // Do nothing.
}

void
gk_neut_species_bflux_step_f(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double dt, struct gkyl_array **fin)
{
  bflux->bflux_step_f_func(app, bflux, fout, dt, fin);
}

static void
gk_neut_species_bflux_set_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac, struct gkyl_array **fin)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    for (int m=0; m<bflux->num_calc_moms; m++)
      gkyl_array_set_range(fout[b*bflux->num_calc_moms+m],
        fac, fin[b*bflux->num_calc_moms+m], bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_set_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac, struct gkyl_array **fin)
{
}

void
gk_neut_species_bflux_set(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac, struct gkyl_array **fin)
{
  bflux->bflux_set_func(app, bflux, fout, fac, fin);
}

static void
gk_neut_species_bflux_accumulate_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac, struct gkyl_array **fin)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    for (int m=0; m<bflux->num_calc_moms; m++)
      gkyl_array_accumulate_range(fout[b*bflux->num_calc_moms+m],
        fac, fin[b*bflux->num_calc_moms+m], bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_accumulate_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac, struct gkyl_array **fin)
{
}

void
gk_neut_species_bflux_accumulate(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac, struct gkyl_array **fin)
{
  bflux->bflux_accumulate_func(app, bflux, fout, fac, fin);
}

static void
gk_neut_species_bflux_combine_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    for (int m=0; m<bflux->num_calc_moms; m++)
      gkyl_array_accumulate_range(gkyl_array_set_range(fout[b*bflux->num_calc_moms+m],
        fac1, fin1[b*bflux->num_calc_moms+m], bflux->boundaries_conf_ghost[b]),
        fac2, fin2[b*bflux->num_calc_moms+m], bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_combine_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2)
{
}

void
gk_neut_species_bflux_combine(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, double fac1, struct gkyl_array **fin1, double fac2, struct gkyl_array **fin2)
{
  bflux->bflux_combine_func(app, bflux, fout, fac1, fin1, fac2, fin2);
}

static void
gk_neut_species_bflux_copy_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin)
{
  for (int b=0; b<bflux->num_boundaries; ++b)
    for (int m=0; m<bflux->num_calc_moms; m++)
      gkyl_array_copy_range(fout[b*bflux->num_calc_moms+m], fin[b*bflux->num_calc_moms+m], bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_copy_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin)
{
}

void
gk_neut_species_bflux_copy(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  struct gkyl_array **fout, struct gkyl_array **fin)
{
  bflux->bflux_copy_func(app, bflux, fout, fin);
}

void
gk_neut_species_bflux_rhs_calc(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  for (int b=0; b<bflux->num_boundaries; ++b) {
    // Ghost cells of the rhs array are filled with the bflux. This is overwritten
    // by the boundary conditions, but it is used before that happens.
    gkyl_array_clear_range(rhs, 0.0, bflux->boundaries_phase_ghost[b]);
    gkyl_boundary_flux_advance(bflux->flux_slvr[b], fin, rhs);
    gkyl_array_copy_range_to_range(bflux->flux[b], rhs, &bflux->boundaries_phase_ghost_nosub[b], bflux->boundaries_phase_ghost[b]);
  }
  app->stat.neut_species_bflux_calc_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_bflux_rhs_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
}

void
gk_neut_species_bflux_rhs(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  bflux->bflux_rhs_func(app, bflux, fin, rhs);
}

static void
gk_neut_species_bflux_calc_moms_enabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  struct timespec wst = gkyl_wall_clock();
  // Compute moments of boundary fluxes.
  for (int b=0; b<bflux->num_boundaries; ++b) {
    for (int m=0; m<bflux->num_calc_moms; m++) {
      gk_neut_species_moment_calc(&bflux->moms_op[m], *bflux->boundaries_phase_ghost[b], *bflux->boundaries_conf_ghost[b], rhs);

      gkyl_array_copy_range(bflux_moms[b*bflux->num_calc_moms+m], bflux->moms_op[m].marr, bflux->boundaries_conf_ghost[b]);
    }
  }
  app->stat.neut_species_bflux_moms_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_bflux_calc_moms_disabled(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  // Do nothing.
}

void
gk_neut_species_bflux_calc_moms(gkyl_gyrokinetic_app *app, struct gk_boundary_fluxes *bflux,
  const struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  bflux->bflux_calc_moms_func(app, bflux, rhs, bflux_moms);
}

static void
gk_neut_species_bflux_get_flux_mom_enabled(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  enum gkyl_distribution_moments mom_type, struct gkyl_array **bflux_moms, struct gkyl_array *out, const struct gkyl_range *out_rng)
{
  int b = gk_neut_species_bflux_idx(bflux, dir, edge);
  int mom_idx = -1;
  for (int m=0; m<bflux->num_calc_moms; m++) {
    if (bflux->calc_mom_names[m] == mom_type) {
      mom_idx = m;
      break;
    }
  }
  gkyl_array_copy_range_to_range(out, bflux_moms[b*bflux->num_calc_moms+mom_idx], out_rng, bflux->boundaries_conf_ghost[b]);
}

static void
gk_neut_species_bflux_get_flux_mom_disabled(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  enum gkyl_distribution_moments mom_type, struct gkyl_array **bflux_moms, struct gkyl_array *out, const struct gkyl_range *out_rng)
{
  // Do nothing.
}

void
gk_neut_species_bflux_get_flux_mom(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  enum gkyl_distribution_moments mom_type, struct gkyl_array **bflux_moms, struct gkyl_array *out, const struct gkyl_range *out_rng)
{
  bflux->bflux_get_flux_mom_func(bflux, dir, edge, mom_type, bflux_moms, out, out_rng);
}

static void
gk_neut_species_bflux_get_flux_enabled(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  struct gkyl_array *out, const struct gkyl_range *out_rng)
{
  int b = gk_neut_species_bflux_idx(bflux, dir, edge);
  gkyl_array_copy_range_to_range(out, bflux->flux[b], out_rng, bflux->boundaries_phase_ghost[b]);
}

static void
gk_neut_species_bflux_get_flux_disabled(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  struct gkyl_array *out, const struct gkyl_range *out_rng)
{
}

void
gk_neut_species_bflux_get_flux(struct gk_boundary_fluxes *bflux, int dir, enum gkyl_edge_loc edge,
  struct gkyl_array *out, const struct gkyl_range *out_rng)
{
  bflux->bflux_get_flux_func(bflux, dir, edge, out, out_rng);
}

void
gk_neut_species_bflux_calc_integrated_mom_enabled(gkyl_gyrokinetic_app* app,
  void *spec_in, struct gk_boundary_fluxes *bflux, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  const struct gk_neut_species *gkns = spec_in;

  int num_diag_int_mom = gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments;
  for (int m=0; m<num_diag_int_mom; m++) {
    int int_mom_idx = bflux->diag_int_mom_idx[m];
    int num_mom_comp = bflux->moms_op[int_mom_idx].num_mom; 
    double avals_global[num_mom_comp];

    for (int b=0; b<bflux->num_boundaries; ++b) {
      // Integrated moment of the boundary flux.
      int dir = bflux->boundaries_dir[b];
      gkyl_array_integrate_advance(bflux->integ_op[m], bflux->f[b*bflux->num_calc_moms+int_mom_idx], 1.0, 0,
        bflux->boundaries_conf_ghost[b], 0, bflux->int_moms_local);

      gkyl_comm_allreduce(app->comm_plane[dir], GKYL_DOUBLE, GKYL_SUM, num_mom_comp, 
        bflux->int_moms_local, bflux->int_moms_global);
      if (app->use_gpu) {
        gkyl_cu_memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom_comp]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom_comp]));
      }

      gkyl_dynvec_append(bflux->intmom[b*num_diag_int_mom+m], tm, avals_global);
      
    } 
  } 
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_bflux_calc_voltime_integrated_mom_enabled(gkyl_gyrokinetic_app* app,
  void *spec_in, struct gk_boundary_fluxes *bflux, double tm)
{
  const struct gk_neut_species *gkns = spec_in;

  int num_diag_int_mom = gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments;
  for (int m=0; m<num_diag_int_mom; m++) {
    int int_mom_idx = bflux->diag_int_mom_idx[m];
    int num_mom_comp = bflux->moms_op[int_mom_idx].num_mom; 
    double avals_global[num_mom_comp];
    double *intmom_cumm_buff = bflux->intmom_cumm_buff[m];

    for (int b=0; b<bflux->num_boundaries; ++b) {
      // Integrated moment of the boundary flux.
      int dir = bflux->boundaries_dir[b];
      gkyl_array_integrate_advance(bflux->integ_op[m], bflux->f[b*bflux->num_calc_moms+int_mom_idx], 1., 0,
        bflux->boundaries_conf_ghost[b], 0, bflux->int_moms_local);

      gkyl_comm_allreduce(app->comm_plane[dir], GKYL_DOUBLE, GKYL_SUM, num_mom_comp, 
        bflux->int_moms_local, bflux->int_moms_global);
      if (app->use_gpu) {
        gkyl_cu_memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom_comp]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        memcpy(avals_global, bflux->int_moms_global, sizeof(double[num_mom_comp]));
      }

      for (int k=0; k<num_mom_comp; k++)
        intmom_cumm_buff[b*num_mom_comp+k] += avals_global[k];
    }
  }
}

static void
gk_neut_species_bflux_calc_voltime_integrated_mom_disabled(gkyl_gyrokinetic_app* app,
  void *species, struct gk_boundary_fluxes *bflux, double tm)
{
}

void
gk_neut_species_bflux_calc_voltime_integrated_mom(gkyl_gyrokinetic_app* app,
  void *species, struct gk_boundary_fluxes *bflux, double tm)
{
  bflux->bflux_calc_voltime_int_mom_func(app, species, bflux, tm);
}

static void
gk_neut_species_bflux_append_integrated_mom(gkyl_gyrokinetic_app* app,
  void *species, struct gk_boundary_fluxes *bflux, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  // Append the time integrated moment of the boundary flux.
  const struct gk_neut_species *gkns = species;

  int num_diag_int_mom = gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments;
  for (int m=0; m<num_diag_int_mom; m++) {
    int int_mom_idx = bflux->diag_int_mom_idx[m];
    int num_mom_comp = bflux->moms_op[int_mom_idx].num_mom;
    double *intmom_cumm_buff = bflux->intmom_cumm_buff[m];
    for (int b=0; b<bflux->num_boundaries; ++b)
      gkyl_dynvec_append(bflux->intmom[b*num_diag_int_mom+m], tm, &intmom_cumm_buff[b*num_mom_comp]);
  }
  app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_neut_species_bflux_calc_integrated_mom_disabled(gkyl_gyrokinetic_app* app,
  void *species, struct gk_boundary_fluxes *bflux, double tm)
{
}

void
gk_neut_species_bflux_calc_integrated_mom(gkyl_gyrokinetic_app* app,
  void *species, struct gk_boundary_fluxes *bflux, double tm)
{
  bflux->bflux_calc_integrated_mom_func(app, species, bflux, tm);
}
  
static void
gk_neut_species_bflux_write_integrated_mom_enabled(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux)
{
  struct timespec wtm = gkyl_wall_clock();
  const struct gk_neut_species *gkns = species;

  int rank, comm_size;
  gkyl_comm_get_rank(app->comm, &rank);
  gkyl_comm_get_size(app->comm, &comm_size);

  const char *vars[] = {"x","y","z"};
  const char *edge[] = {"lower","upper"};

  int num_diag_int_mom = gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments;
  for (int b=0; b<bflux->num_boundaries; ++b) {
    int dir = bflux->boundaries_dir[b];
    int edi = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? 0 : 1;

    if ((edi == 0 && rank == 0) || (edi == 1 && rank == comm_size-1)) {
      for (int m=0; m<num_diag_int_mom; m++) {
        // Write integrated moments of the boundary fluxes.
        const char *fmt = "%s-%s_bflux_%s%s_integrated_%s.gkyl";
        const char *mom_name = gkyl_distribution_moments_strs[gkns->info.boundary_flux_diagnostics.integrated_diag_moments[m]];

        int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, vars[dir], edge[edi], mom_name);
        char fileNm[sz+1]; // ensures no buffer overflow
        snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, vars[dir], edge[edi], mom_name);

        struct timespec wtm = gkyl_wall_clock();
        if (bflux->is_first_intmom_write_call[b]) {
          struct gkyl_msgpack_map_elem io_meta_phi[] = {
            { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Volume integrated moments of the boundary flux." }
          };
          int io_meta_len[] = {gkns->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
          const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi};
          struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

          gkyl_dynvec_write_wmeta(bflux->intmom[b*num_diag_int_mom+m], fileNm, mt);
          gkyl_msgpack_data_release(mt);
        }
        else {
          gkyl_dynvec_awrite(bflux->intmom[b*num_diag_int_mom+m], fileNm);
        }

        app->stat.n_diag_io += 1;

        gkyl_dynvec_clear(bflux->intmom[b*num_diag_int_mom+m]);
      
        if (bflux->is_first_intmom_write_call[b])
          bflux->is_first_intmom_write_call[b] = false;
      }
    }
  }
  app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gk_neut_species_bflux_write_integrated_mom_disabled(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux)
{
}

void
gk_neut_species_bflux_write_integrated_mom(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux)
{
  bflux->bflux_write_integrated_mom_func(app, species, bflux);
}

static void
gk_neut_species_bflux_write_mom_enabled(gkyl_gyrokinetic_app* app, void *species,
  struct gk_boundary_fluxes *bflux, double tm, int frame)
{
  const struct gk_neut_species *gkns = species;

  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gkns->io_meta_conf_len, gkns->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gkns->io_meta_conf_len, gkns->io_meta_conf, "frame", frame);
  struct gkyl_msgpack_map_elem desc[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Velocity-space moment of the boundary flux." }
  };
  int io_meta_len[] = {gkns->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gkns->io_meta_conf, app->gk_geom->io_meta_basic, desc};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  int rank, comm_size;
  gkyl_comm_get_rank(app->comm, &rank);
  gkyl_comm_get_size(app->comm, &comm_size);

  const char *vars[] = {"x","y","z"};
  const char *edge[] = {"lower","upper"};

  int num_diag_mom = gkns->info.boundary_flux_diagnostics.num_diag_moments;

  for (int b=0; b<bflux->num_boundaries; ++b) {
    int dir = bflux->boundaries_dir[b];
    int edi = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? 0 : 1;

    if ((edi == 0 && rank == 0) || (edi == 1 && rank == comm_size-1)) {
      for (int m=0; m<num_diag_mom; ++m) {
        struct timespec wst = gkyl_wall_clock();
        int mom_idx = bflux->diag_mom_idx[m];
        struct gkyl_array *mom_arr = bflux->f[b*bflux->num_calc_moms+mom_idx];

        // For now copy the moment to the skin ghost and write it out.
        gkyl_array_copy_range_to_range(mom_arr, mom_arr,
          bflux->boundaries_conf_skin[b], bflux->boundaries_conf_ghost[b]);

        // Rescale by dx/2 in the direction of the boundary to account for the
        // normalization in the boundary surf kernels.
        gkyl_array_scale_range(mom_arr, 0.5*app->grid.dx[dir], bflux->boundaries_conf_skin[b]);
        app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
          
        const char *fmt = "%s-%s_bflux_%s%s_%s_%d.gkyl";
        const char *mom_name = gkyl_distribution_moments_strs[gkns->info.boundary_flux_diagnostics.diag_moments[mom_idx]];
        int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, vars[dir], edge[edi], mom_name, frame);
        char fileNm[sz+1]; // ensures no buffer overflow
        snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, vars[dir], edge[edi], mom_name, frame);
        
        struct timespec wtm = gkyl_wall_clock();
        if (app->cdim > 1) {
          struct timespec wst = gkyl_wall_clock();
          // Project the moment down to lower dimensions.
          int num_mom_comp = bflux->moms_op[mom_idx].num_mom;
          gkyl_translate_dim_advance(bflux->transdim[b], bflux->boundaries_conf_skin_fullx[b], &bflux->surf_local[dir],
            mom_arr, num_mom_comp, bflux->mom_surf[b*num_diag_mom+m]);
          app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

          struct timespec wtm = gkyl_wall_clock();
          if (app->use_gpu)
            gkyl_array_copy(bflux->mom_surf_ho[b*num_diag_mom+m], bflux->mom_surf[b*num_diag_mom+m]);
  
          gkyl_comm_array_write(bflux->comm_surf[dir], &bflux->grid_surf[dir], &bflux->surf_local[dir], mt,
            bflux->mom_surf_ho[b*num_diag_mom+m], fileNm);
          app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
        }
        else {
          // Don't project down to 0D; the infrastructure doesn't make it easy to do so.
          struct timespec wtm = gkyl_wall_clock();
          if (app->use_gpu)
            gkyl_array_copy(bflux->moms_op[mom_idx].marr_host, mom_arr);
  
          gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
            bflux->moms_op[mom_idx].marr_host, fileNm);

          app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
        }
  
        app->stat.n_diag_io += 1;
      }
    }
  }

  gkyl_msgpack_data_release(mt);  
}

static void
gk_neut_species_bflux_write_mom_disabled(gkyl_gyrokinetic_app* app, void *species,
  struct gk_boundary_fluxes *bflux, double tm, int frame)
{
}

void
gk_neut_species_bflux_write_mom(gkyl_gyrokinetic_app* app, void *species,
  struct gk_boundary_fluxes *bflux, double tm, int frame)
{
  bflux->bflux_write_mom_func(app, species, bflux, tm, frame);
}

static int *
bflux_unionize_moms(int num_add_moms, enum gkyl_distribution_moments add_moms[BFLUX_MAX_MOM_NAMES],
  int *num_moms, enum gkyl_distribution_moments moms[BFLUX_MAX_MOM_NAMES])
{
  // Check if each of the num_add_moms moments in add_moms is included in the
  // moms list of num_moms moments. If it's not, include it and increment
  // num_moms. Return a list of the indices of each add_moms in the moms list.
  int *add_mom_idx = gkyl_malloc(num_add_moms*sizeof(int));
  int num_moms_base = num_moms[0];
  for (int i=0; i<num_add_moms; i++) {
    enum gkyl_distribution_moments mom_name = add_moms[i];
    bool included = false;
    for (int j=0; j<num_moms_base; j++) {
      if (mom_name == moms[j]) {
        included = true;
        add_mom_idx[i] = j;
        break;
      }
    }

    if (!included) {
      add_mom_idx[i] = num_moms[0];
      moms[num_moms[0]] = mom_name;
      num_moms[0] += 1;
    }
  }
  return add_mom_idx;
}

void 
gk_neut_species_bflux_init(struct gkyl_gyrokinetic_app *app, void *species,
  struct gk_boundary_fluxes *bflux, enum gkyl_species_bflux_type bflux_type,
  struct gkyl_phase_diagnostics_inp add_moms_inp)
{ 
  struct gk_neut_species *gkns = species;

  bflux->allocated_solver = false;
  bflux->allocated_moms = false;
  bflux->allocated_diags = false;

  // Set function pointers to empty functions.
  bflux->bflux_rhs_func = gk_neut_species_bflux_rhs_disabled;
  bflux->bflux_calc_moms_func = gk_neut_species_bflux_calc_moms_disabled;
  bflux->bflux_get_flux_func = gk_neut_species_bflux_get_flux_disabled;
  bflux->bflux_get_flux_mom_func = gk_neut_species_bflux_get_flux_mom_disabled;
  bflux->bflux_clear_func = gk_neut_species_bflux_clear_disabled;
  bflux->bflux_scale_func = gk_neut_species_bflux_scale_disabled;
  bflux->bflux_step_f_func = gk_neut_species_bflux_step_f_disabled;
  bflux->bflux_set_func = gk_neut_species_bflux_set_disabled;
  bflux->bflux_accumulate_func = gk_neut_species_bflux_accumulate_disabled;
  bflux->bflux_combine_func = gk_neut_species_bflux_combine_disabled;
  bflux->bflux_copy_func = gk_neut_species_bflux_copy_disabled;
  bflux->bflux_calc_integrated_mom_func = gk_neut_species_bflux_calc_integrated_mom_disabled;
  bflux->bflux_write_integrated_mom_func = gk_neut_species_bflux_write_integrated_mom_disabled;
  bflux->bflux_calc_voltime_int_mom_func = gk_neut_species_bflux_calc_voltime_integrated_mom_disabled;
  bflux->bflux_write_mom_func = gk_neut_species_bflux_write_mom_disabled;

  if (bflux_type != GK_SPECIES_BFLUX_NONE) {
    bflux->allocated_solver = true;

    // Set function pointer to compute bfluxes.
    bflux->bflux_rhs_func = gk_neut_species_bflux_rhs_calc; 
    bflux->bflux_get_flux_func = gk_neut_species_bflux_get_flux_enabled;

    // Identify the non-periodic, non-zero-flux boundaries to compute boundary fluxes at.
    bflux->num_boundaries = 0;
    for (int d=0; d<app->cdim; ++d) {
      for (int e=0; e<2; ++e) {
        if ( gkns->bc_is_np[d] &&
             ((e == 0 && gkns->lower_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX) ||
              (e == 1 && gkns->upper_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX)) ) {
          bflux->boundaries_dir[bflux->num_boundaries] = d;
          bflux->boundaries_edge[bflux->num_boundaries] = e==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
          bflux->boundaries_conf_skin[bflux->num_boundaries] = e==0? &app->local_lower_skin[d] : &app->local_upper_skin[d];
          bflux->boundaries_conf_ghost[bflux->num_boundaries] = e==0? &app->local_lower_ghost[d] : &app->local_upper_ghost[d];
          bflux->boundaries_phase_ghost[bflux->num_boundaries] = e==0? &gkns->local_lower_ghost[d] : &gkns->local_upper_ghost[d];
          bflux->num_boundaries++;
        }
      }
    }
  
    // Create an array of equation objects with terms that produce boundary fluxes.
    bflux->num_eqns = 0;
    if (gkns->collisionless.collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL)
      bflux->num_eqns += 1; // Collisionless terms.

    bflux->eqns = gkyl_malloc(bflux->num_eqns*sizeof(struct gkyl_dg_eqn *));

    int eqc = 0;
    if (gkns->collisionless.collisionless_id == GKYL_GK_COLLISIONLESS_NEUTRAL)
      bflux->eqns[eqc++] = gkyl_dg_updater_vlasov_acquire_eqn(gkns->collisionless.vlasov_slvr);
  
    // Allocate updater that computes boundary fluxes.
    for (int b=0; b<bflux->num_boundaries; ++b) {
      int dir = bflux->boundaries_dir[b];
      struct gkyl_range *skin_r = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &gkns->local_lower_skin[dir]
                                                                            : &gkns->local_upper_skin[dir];
      struct gkyl_range *ghost_r = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &gkns->local_lower_ghost[dir]
                                                                             : &gkns->local_upper_ghost[dir];
      bflux->flux_slvr[b] = gkyl_boundary_flux_new(dir, bflux->boundaries_edge[b], &gkns->grid,
        skin_r, ghost_r, bflux->num_eqns, bflux->eqns, app->use_gpu);
    }

    // Create a ghost range that the flux lives on, and allocate the array that stores the flux.
    int ndim = gkns->local.ndim;
    bflux->boundaries_phase_ghost_nosub = gkyl_malloc(bflux->num_boundaries*sizeof(struct gkyl_range));
    bflux->flux = gkyl_malloc(bflux->num_boundaries*sizeof(struct gkyl_array *));
    for (int b=0; b<bflux->num_boundaries; ++b) {
      int rlower[ndim], rupper[ndim];
      for (int d=0; d<ndim; d++) {
        rlower[d] = bflux->boundaries_phase_ghost[b]->lower[d];
        rupper[d] = bflux->boundaries_phase_ghost[b]->upper[d];
      }
      gkyl_range_init(&bflux->boundaries_phase_ghost_nosub[b], ndim, rlower, rupper);
      bflux->flux[b] = mkarr(app->use_gpu, gkns->basis.num_basis, bflux->boundaries_phase_ghost_nosub[b].volume);
    }
  }

  int num_diag_mom = gkns->info.boundary_flux_diagnostics.num_diag_moments;
  int num_diag_int_mom = gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments;
  
  if (bflux_type == GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS || bflux_type == GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS) {
    bflux->allocated_moms = true;

    // Set methods for time-stepping boundary fluxes needed for diagnostics.
    bflux->bflux_get_flux_mom_func = gk_neut_species_bflux_get_flux_mom_enabled;
    bflux->bflux_calc_moms_func = gk_neut_species_bflux_calc_moms_enabled;
    bflux->bflux_clear_func = gk_neut_species_bflux_clear_enabled;
    bflux->bflux_scale_func = gk_neut_species_bflux_scale_enabled;
    bflux->bflux_step_f_func = gk_neut_species_bflux_step_f_enabled;
    bflux->bflux_set_func = gk_neut_species_bflux_set_enabled;
    bflux->bflux_accumulate_func = gk_neut_species_bflux_accumulate_enabled;
    bflux->bflux_combine_func = gk_neut_species_bflux_combine_enabled;
    bflux->bflux_copy_func = gk_neut_species_bflux_copy_enabled;

    // Create a union of the diag_moms, int_diag_moms and add_moms lists. Also store
    // the index of each mom in this union list.
    int num_add_mom = add_moms_inp.num_diag_moments;
    assert(num_add_mom+num_diag_mom+num_diag_int_mom < BFLUX_MAX_MOM_NAMES+1);

    bflux->num_calc_moms = 0;
    if (num_diag_mom > 0)
      bflux->diag_mom_idx = bflux_unionize_moms(num_diag_mom,
        gkns->info.boundary_flux_diagnostics.diag_moments, &bflux->num_calc_moms, bflux->calc_mom_names);

    if (num_diag_int_mom > 0)
      bflux->diag_int_mom_idx = bflux_unionize_moms(num_diag_int_mom,
        gkns->info.boundary_flux_diagnostics.integrated_diag_moments, &bflux->num_calc_moms, bflux->calc_mom_names);

    if (num_add_mom > 0) {
      int *add_mom_idx = bflux_unionize_moms(num_add_mom,
        add_moms_inp.diag_moments, &bflux->num_calc_moms, bflux->calc_mom_names);
      gkyl_free(add_mom_idx);
    }

    // Create a moments app for each mom needed.
    bflux->moms_op = gkyl_malloc(sizeof(struct gk_species_moment[bflux->num_calc_moms]));
    bool need_hamil_ghost = false;
    for (int m=0; m<bflux->num_calc_moms; m++) {
      gk_neut_species_moment_init(app, gkns, &bflux->moms_op[m], bflux->calc_mom_names[m], false);

      need_hamil_ghost = (bflux->calc_mom_names[m] == GKYL_F_MOMENT_M1_FROM_H)
        || (bflux->calc_mom_names[m] == GKYL_F_MOMENT_ENERGY);
    }

    if (need_hamil_ghost) {
      // For moments that use the Hamiltonian, we must fill the ghost cell of H. Use an option
      // in bc_basic that fills the ghost cells keeping the value at the skin-ghost boundary the same.
      long buff_sz = 0;
      for (int b=0; b<bflux->num_boundaries; ++b) {
        int dir = bflux->boundaries_dir[b];
        struct gkyl_range *skin_r = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? &gkns->local_lower_skin[dir]
                                                                              : &gkns->local_upper_skin[dir];
    
        // MF 2025/09/29: The option `GKYL_BC_GK_SPECIES_REFLECT` here is a
        // just a place holder and almost certainly wrong. Currently REFLECT is
        // meant to reflect particles. Need something that mirrors the
        // Hamiltonian into the ghost cell.
        bflux->gfss_bc_op[b] = gkyl_bc_basic_gyrokinetic_new(bflux->boundaries_dir[b], bflux->boundaries_edge[b],
          GKYL_BC_GK_SPECIES_REFLECT, &gkns->basis, skin_r, bflux->boundaries_phase_ghost[b], 1, app->cdim, app->use_gpu);
        
        long vol = skin_r->volume;
        buff_sz = buff_sz > vol ? buff_sz : vol;
      }
      bflux->bc_buffer = mkarr(app->use_gpu, gkns->basis.num_basis, buff_sz);
    
      // Fill ghost cell of H.
      for (int b=0; b<bflux->num_boundaries; ++b)
        gkyl_bc_basic_gyrokinetic_advance(bflux->gfss_bc_op[b], bflux->bc_buffer, gkns->hamil);

      gkyl_array_release(bflux->bc_buffer);
      for (int b=0; b<bflux->num_boundaries; ++b)
        gkyl_bc_basic_gyrokinetic_release(bflux->gfss_bc_op[b]);
    }
  
    bflux->f = gkyl_malloc(bflux->num_boundaries*bflux->num_calc_moms*sizeof(struct gkyl_array *));
    bflux->f1 = gkyl_malloc(bflux->num_boundaries*bflux->num_calc_moms*sizeof(struct gkyl_array *));
    bflux->fnew = gkyl_malloc(bflux->num_boundaries*bflux->num_calc_moms*sizeof(struct gkyl_array *));
    for (int b=0; b<bflux->num_boundaries; ++b) {
      for (int m=0; m<bflux->num_calc_moms; m++) {
        // Allocate arrays storing moments of the boundary flux.
        int num_mom_comp = bflux->moms_op[m].num_mom;
        bflux->f[b*bflux->num_calc_moms+m] = mkarr(app->use_gpu, num_mom_comp*app->basis.num_basis, app->local_ext.volume);
        bflux->f1[b*bflux->num_calc_moms+m] = mkarr(app->use_gpu, num_mom_comp*app->basis.num_basis, app->local_ext.volume);
        bflux->fnew[b*bflux->num_calc_moms+m] = mkarr(app->use_gpu, num_mom_comp*app->basis.num_basis, app->local_ext.volume);
      }
    }
  }
    
  if (bflux_type == GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS) {
    assert(num_diag_mom > 0 || num_diag_int_mom > 0);

    bflux->allocated_diags = true;
  
    // Set methods for time-stepping boundary fluxes needed for diagnostics.
    if (gkns->info.boundary_flux_diagnostics.time_integrated) {
      bflux->bflux_calc_integrated_mom_func = gk_neut_species_bflux_append_integrated_mom;
      bflux->bflux_calc_voltime_int_mom_func = gk_neut_species_bflux_calc_voltime_integrated_mom_enabled;
    }
    else {
      bflux->bflux_calc_integrated_mom_func = gk_neut_species_bflux_calc_integrated_mom_enabled;
    }
    bflux->bflux_write_integrated_mom_func = gk_neut_species_bflux_write_integrated_mom_enabled;
    bflux->bflux_write_mom_func = gk_neut_species_bflux_write_mom_enabled;

    if (num_diag_mom > 0) {
      // Objects needed to output lower dimensional moments.
      if (app->cdim > 1) {
        struct gkyl_basis basis_conf_surf;
        switch (app->basis.b_type) {
          case GKYL_BASIS_MODAL_SERENDIPITY:
            gkyl_cart_modal_serendip(&basis_conf_surf, app->cdim-1, app->basis.poly_order);
            break;
          default:
            assert(false);
            break;
        }

        bflux->transdim = gkyl_malloc(bflux->num_boundaries*sizeof(struct gkyl_translate_dim *));
        bool diag_in_dir[GKYL_MAX_CDIM] = {0};
        for (int b=0; b<bflux->num_boundaries; ++b) {
          int dir = bflux->boundaries_dir[b];
          enum gkyl_edge_loc edge = bflux->boundaries_edge[b];
          // Updater that projects to lower dim.
          bflux->transdim[b] = gkyl_translate_dim_new(app->cdim, app->basis, app->cdim-1,
            basis_conf_surf, dir, edge==GKYL_LOWER_EDGE? GKYL_UPPER_EDGE : GKYL_LOWER_EDGE, app->use_gpu);

          if (!diag_in_dir[dir]) {

            // Create a communicator associated with a lower dimensional surface range.
            // Identify ranks on the same plane as this one.
            int num_ranks_surf = 0;
            int ranks_surf[app->decomp->ndecomp]; 
            for (int i=0; i<app->decomp->ndecomp; i++) {
              if (app->decomp->ranges[i].lower[dir] == app->local.lower[dir]) {
                ranks_surf[num_ranks_surf] = i;
                num_ranks_surf++;
              }
            }
            // Create a range tangentially global.
            int surf_dim = app->cdim-1;
            int lower_surf[surf_dim], upper_surf[surf_dim];
            int c = 0;
            for (int d=0; d<app->cdim; ++d) {
              if (d != dir) {
                lower_surf[c] = app->global.lower[d];
                upper_surf[c] = app->global.upper[d];
                c++;
              }
            }
            struct gkyl_range range_surf;
            gkyl_range_init(&range_surf, surf_dim, lower_surf, upper_surf);
        
            // Create decomp.
            int cuts_plane[GKYL_MAX_CDIM], cuts_surf[surf_dim];
            gkyl_rect_decomp_get_cuts(app->decomp_plane[dir], cuts_plane);
            c = 0;
            for (int d=0; d<app->cdim; ++d) {
              if (d != dir) {
                cuts_surf[c] = cuts_plane[d];
                c++;
              }
            }
            bflux->decomp_surf[dir] = gkyl_rect_decomp_new_from_cuts(surf_dim, cuts_surf, &range_surf);
        
            // Create a new communicator with ranks on surf.
            bool is_comm_valid;
            bflux->comm_surf[dir] = gkyl_comm_create_comm_from_ranks(app->comm, num_ranks_surf,
              ranks_surf, bflux->decomp_surf[dir], &is_comm_valid);
            assert(is_comm_valid);

            // Local and local extended surface range.
            int rank;
            gkyl_comm_get_rank(bflux->comm_surf[dir], &rank);
            int ghost[] = { 1, 1, 1 };
            gkyl_create_ranges(&bflux->decomp_surf[dir]->ranges[rank], ghost, &bflux->surf_local_ext[dir], &bflux->surf_local[dir]);

            // Create a surface grid.
            double grid_surf_lower[surf_dim], grid_surf_upper[surf_dim];
            int grid_surf_cells[surf_dim];
            c = 0;
            for (int d=0; d<app->cdim; ++d) {
              if (d != dir) {
                grid_surf_lower[c] = app->grid.lower[d];
                grid_surf_upper[c] = app->grid.upper[d];
                grid_surf_cells[c] = app->grid.cells[d];
                c++;
              }
            }
            gkyl_rect_grid_init(&bflux->grid_surf[dir], surf_dim, grid_surf_lower, grid_surf_upper, grid_surf_cells);

            diag_in_dir[dir] = true;
          }
        }

        // Allocate a lower dimensional array for each moment.
        bflux->mom_surf = gkyl_malloc(bflux->num_boundaries*num_diag_mom*sizeof(struct gkyl_array *));
        bflux->mom_surf_ho = gkyl_malloc(bflux->num_boundaries*num_diag_mom*sizeof(struct gkyl_array *));
        for (int b=0; b<bflux->num_boundaries; ++b) {
          int dir = bflux->boundaries_dir[b];
          for (int m=0; m<num_diag_mom; m++) {
            // Allocate arrays storing moments of the boundary flux.
            int mom_idx = bflux->diag_mom_idx[m];
            int num_mom_comp = bflux->moms_op[mom_idx].num_mom;
            bflux->mom_surf[b*num_diag_mom+m] = mkarr(app->use_gpu, num_mom_comp*basis_conf_surf.num_basis, bflux->surf_local_ext[dir].volume);
            bflux->mom_surf_ho[b*num_diag_mom+m] = app->use_gpu? 
              mkarr(false, num_mom_comp*basis_conf_surf.num_basis, bflux->surf_local_ext[dir].volume) :
              gkyl_array_acquire(bflux->mom_surf[b*num_diag_mom+m]);
          }
        }
      }
    }

    if (num_diag_int_mom > 0) {
      // Object to integrate moments of the bflux and dynvectors to store them.
      bflux->integ_op = gkyl_malloc(num_diag_int_mom*sizeof(struct gkyl_array_integrate *));
      bflux->intmom = gkyl_malloc(num_diag_int_mom*bflux->num_boundaries*sizeof(gkyl_dynvec));
      int num_mom_comp_max = 1;
      for (int m=0; m<num_diag_int_mom; m++) {
        int num_mom_comp = bflux->moms_op[bflux->diag_int_mom_idx[m]].num_mom;
        num_mom_comp_max = GKYL_MAX2(num_mom_comp_max, num_mom_comp);
        // Updater to compute the volume integral of the boundary flux moments.
        bflux->integ_op[m] = gkyl_array_integrate_new(&app->grid, &app->basis,
          num_mom_comp, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
        // Allocate a dynvector for each moment.
        for (int b=0; b<bflux->num_boundaries; ++b)
          bflux->intmom[b*bflux->num_calc_moms+m] = gkyl_dynvec_new(GKYL_DOUBLE, num_mom_comp);
      }
  
      if (app->use_gpu) {
        bflux->int_moms_local = gkyl_cu_malloc(num_mom_comp_max*sizeof(double));
        bflux->int_moms_global = gkyl_cu_malloc(num_mom_comp_max*sizeof(double));
      }
      else {
        bflux->int_moms_local = gkyl_malloc(num_mom_comp_max*sizeof(double));
        bflux->int_moms_global = gkyl_malloc(num_mom_comp_max*sizeof(double));
      }
    
      for (int b=0; b<bflux->num_boundaries; ++b)
        bflux->is_first_intmom_write_call[b] = true;
    
      if (gkns->info.boundary_flux_diagnostics.time_integrated) {
        // Cummulative integrated moments of boundary fluxes.
        bflux->intmom_cumm_buff = gkyl_malloc(num_diag_int_mom*sizeof(double *));
        for (int m=0; m<num_diag_int_mom; m++) {
          int num_mom_comp = bflux->moms_op[bflux->diag_int_mom_idx[m]].num_mom;
          bflux->intmom_cumm_buff[m] = gkyl_malloc(bflux->num_boundaries*num_mom_comp*sizeof(double));
          double *buff = bflux->intmom_cumm_buff[m];
          for (int b=0; b<bflux->num_boundaries; ++b) {
            for (int k=0; k<num_mom_comp; k++)
              buff[b*num_mom_comp+k] = 0.0;
          }
        }
      }

    }
  }

}

void
gk_neut_species_bflux_read_voltime_integrated_mom(gkyl_gyrokinetic_app *app,
  void *species, struct gk_boundary_fluxes *bflux)
{
  const struct gk_neut_species *gkns = species;

  if (gkns->info.boundary_flux_diagnostics.time_integrated) {
    int rank, comm_size;
    gkyl_comm_get_rank(app->comm, &rank);
    gkyl_comm_get_size(app->comm, &comm_size);

    const char *vars[] = {"x","y","z"};
    const char *edge[] = {"lower","upper"};

    int num_diag_int_mom = gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments;
    for (int b=0; b<bflux->num_boundaries; ++b) {
      int dir = bflux->boundaries_dir[b];
      int edi = bflux->boundaries_edge[b]==GKYL_LOWER_EDGE? 0 : 1;

      if ((edi == 0 && rank == 0) || (edi == 1 && rank == comm_size-1)) {
        for (int m=0; m<num_diag_int_mom; m++) {
          // Write integrated moments of the boundary fluxes.
          const char *fmt = "%s-%s_bflux_%s%s_integrated_%s.gkyl";
          const char *mom_name = gkyl_distribution_moments_strs[gkns->info.boundary_flux_diagnostics.integrated_diag_moments[m]];

          int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, vars[dir], edge[edi], mom_name);
          char fileNm[sz+1]; // ensures no buffer overflow
          snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, vars[dir], edge[edi], mom_name);

          // We didn't calculate int_mom at restart, so read the value from the previous sim
          // and append it. This is the best solution given the flow in present input files.
          bool res = gkyl_dynvec_read(bflux->intmom[b*num_diag_int_mom+m], fileNm);
          int num_mom_comp = gkyl_dynvec_ncomp(bflux->intmom[b*num_diag_int_mom+m]);

          double vals_prev[num_mom_comp];
          gkyl_dynvec_getlast(bflux->intmom[b*num_diag_int_mom+m], vals_prev);
          gkyl_dynvec_clear(bflux->intmom[b*num_diag_int_mom+m]);

          double *intmom_cumm_buff = bflux->intmom_cumm_buff[m];
          for (int k=0; k<num_mom_comp; k++)
            intmom_cumm_buff[b*num_mom_comp+k] += vals_prev[k];

        }
      }
    }
  }
}

void
gk_neut_species_bflux_release(const struct gkyl_gyrokinetic_app *app, const void *species,
  const struct gk_boundary_fluxes *bflux)
{
  const struct gk_neut_species *gkns = species;

  if (bflux->allocated_solver) {
    for (int b=0; b<bflux->num_boundaries; ++b) {
      gkyl_boundary_flux_release(bflux->flux_slvr[b]);
    }
    for (int i=0; i<bflux->num_eqns; i++) {
      gkyl_dg_eqn_release(bflux->eqns[i]);
    }
    gkyl_free(bflux->eqns);

    for (int b=0; b<bflux->num_boundaries; ++b) {
      gkyl_array_release(bflux->flux[b]);
    }
    gkyl_free(bflux->flux);
    gkyl_free(bflux->boundaries_phase_ghost_nosub);
  }

  int num_diag_mom = gkns->info.boundary_flux_diagnostics.num_diag_moments;
  int num_diag_int_mom = gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments;

  if (bflux->allocated_moms) {
    if (num_diag_mom > 0)
      gkyl_free(bflux->diag_mom_idx);

    if (num_diag_int_mom > 0)
      gkyl_free(bflux->diag_int_mom_idx);

    for (int m=0; m<bflux->num_calc_moms; m++)
      gk_neut_species_moment_release(app, &bflux->moms_op[m]); 

    gkyl_free(bflux->moms_op);
    gkyl_free(bflux->is_hamiltonian_mom);
    if (bflux->a_hamiltonian_mom) {
      gkyl_array_release(bflux->bc_buffer);
      for (int b=0; b<bflux->num_boundaries; ++b)
        gkyl_bc_basic_gyrokinetic_release(bflux->gfss_bc_op[b]);
    }
    for (int b=0; b<bflux->num_boundaries; ++b) {
      for (int m=0; m<bflux->num_calc_moms; m++) {
        gkyl_array_release(bflux->f[b*bflux->num_calc_moms+m]);
        gkyl_array_release(bflux->f1[b*bflux->num_calc_moms+m]);
        gkyl_array_release(bflux->fnew[b*bflux->num_calc_moms+m]);
      }
    }
    gkyl_free(bflux->f);
    gkyl_free(bflux->f1);
    gkyl_free(bflux->fnew);
  }

  if (bflux->allocated_diags) {
    if (num_diag_mom > 0) {
      if (app->cdim > 1) {
        // Objects needed to output lower dimensional moments.
        for (int b=0; b<bflux->num_boundaries; ++b)
          gkyl_translate_dim_release(bflux->transdim[b]);

        gkyl_free(bflux->transdim);

        bool diag_in_dir[GKYL_MAX_CDIM] = {0};
        for (int b=0; b<bflux->num_boundaries; ++b) {
          int dir = bflux->boundaries_dir[b];
          if (!diag_in_dir[dir]) {
            gkyl_rect_decomp_release(bflux->decomp_surf[dir]);
            gkyl_comm_release(bflux->comm_surf[dir]);
            diag_in_dir[dir] = true;
          }
        }

        for (int b=0; b<bflux->num_boundaries; ++b) {
          for (int m=0; m<num_diag_mom; m++) {
            gkyl_array_release(bflux->mom_surf[b*num_diag_mom+m]);
            gkyl_array_release(bflux->mom_surf_ho[b*num_diag_mom+m]);
          }
        }
        gkyl_free(bflux->mom_surf);
        gkyl_free(bflux->mom_surf_ho);
      }
    }

    if (num_diag_int_mom > 0) {
      for (int m=0; m<num_diag_int_mom; m++) {
        gkyl_array_integrate_release(bflux->integ_op[m]);
        for (int b=0; b<bflux->num_boundaries; ++b)
          gkyl_dynvec_release(bflux->intmom[b*bflux->num_calc_moms+m]);
      }
      gkyl_free(bflux->integ_op);
      gkyl_free(bflux->intmom);
      if (app->use_gpu) {
        gkyl_cu_free(bflux->int_moms_local);
        gkyl_cu_free(bflux->int_moms_global);
      }
      else {
        gkyl_free(bflux->int_moms_local);
        gkyl_free(bflux->int_moms_global);
      }
  
      if (gkns->info.boundary_flux_diagnostics.time_integrated) {
        for (int m=0; m<num_diag_int_mom; m++)
          gkyl_free(bflux->intmom_cumm_buff[m]);
        gkyl_free(bflux->intmom_cumm_buff);
      }
    }
  }
}
