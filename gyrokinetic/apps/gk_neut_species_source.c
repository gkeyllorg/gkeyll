#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_neut_species_source_init(
  struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, struct gk_source *src
)
{
  src->source_id = gkns->info.source.source_id;

  if (src->source_id) {
    int vdim = gkns->info.vdim;
    // we need to ensure source has same shape as distribution function
    src->source = mkarr(app->use_gpu, gkns->basis.num_basis, gkns->local_ext.volume);
    src->source_host = src->source;
    if (app->use_gpu) {
      src->source_host = mkarr(false, gkns->basis.num_basis, gkns->local_ext.volume);
    }

    src->evolve = gkns->info.source.evolve; // Whether the source is time dependent.

    src->num_sources = gkns->info.source.num_sources;
    if (src->source_id == GKYL_SOURCE_FROMFILE) {
      // A file-imported source shape is static: no time dependence.
      assert(!src->evolve);
      struct gkyl_array *fdo = mkarr(app->use_gpu, gkns->basis.num_basis, gkns->local_ext.volume);
      gkyl_array_clear(src->source, 0.0);
      for (int k = 0; k < gkns->info.source.num_sources; k++) {
        assert(!gkns->info.source.source_import[k].enforce_positivity
        ); // Not supported for sources yet.
        gk_neut_species_kinetic_read_distf_from_file(
          app, gkns, gkns->info.source.source_import[k], fdo
        );
        gkyl_array_accumulate(src->source, 1., fdo);
      }
      gkyl_array_release(fdo);
    } else {
      for (int k = 0; k < gkns->info.source.num_sources; k++) {
        gk_neut_species_projection_init(
          app, gkns, gkns->info.source.projection[k], &src->proj_source[k]
        );
      }
    }

    // Allocate data and updaters for diagnostic moments.
    src->num_diag_mom = gkns->info.num_diag_moments;
    gkns->src.moms = gkyl_malloc(sizeof(struct gk_species_moment[src->num_diag_mom]));
    for (int m = 0; m < src->num_diag_mom; ++m) {
      gk_neut_species_moment_init(app, gkns, &gkns->src.moms[m], gkns->info.diag_moments[m], false);
    }

    // Allocate data and updaters for integrated moments.
    enum gkyl_distribution_moments integ_mom_types[] = {GKYL_F_MOMENT_M0M1M2};
    gk_integ_moms_neut_init(app, gkns, &gkns->src.integ_moms, 1, integ_mom_types);
  }
}

void
gk_neut_species_source_calc(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, struct gk_source *src,
  struct gkyl_array *f_buffer, double tm
)
{
  // A file-imported source's shape was already read once in gk_neut_species_source_init.
  if (src->source_id && src->source_id != GKYL_SOURCE_FROMFILE) {
    for (int k = 0; k < gkns->info.source.num_sources; k++) {
      gk_neut_species_projection_calc(app, gkns, &src->proj_source[k], f_buffer, tm);
      gkyl_array_accumulate(src->source, 1., f_buffer);
    }
  }
}

// Compute rhs of the source
void
gk_neut_species_source_rhs(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, struct gk_source *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();
  if (src->source_id) {
    gkyl_array_accumulate(rhs, 1.0, src->source);
  }
  app->stat.neut_species_src_tm += gkyl_time_diff_now_sec(wst);
}

// Write functions
void
gk_neut_species_source_write(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, double tm, int frame
)
{
  if (gkns->src.source_id && (gkns->src.evolve || frame == 0)) {
    struct timespec wtm = gkyl_wall_clock();

    // Package metadata.
    gkyl_msgpack_map_elem_set_double(gkns->io_meta_phase_len, gkns->io_meta_phase, "time", tm);
    gkyl_msgpack_map_elem_set_uint(gkns->io_meta_phase_len, gkns->io_meta_phase, "frame", frame);
    struct gkyl_msgpack_map_elem desc_src[] = {
      {.key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Neutral species source."},
    };
    int io_meta_len[] = {gkns->io_meta_phase_len, app->gk_geom->io_meta_basic_len, 1};
    const struct gkyl_msgpack_map_elem *io_meta[] = {
      gkns->io_meta_phase, app->gk_geom->io_meta_basic, desc_src
    };
    struct gkyl_msgpack_data *mt =
      gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

    // Write out the source distribution function
    const char *fmt = "%s-%s_source_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, frame);
    char fileNm[sz + 1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, frame);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gkns->src.source_host, gkns->src.source);
    }

    gkyl_comm_array_write(gkns->comm, &gkns->grid, &gkns->local, mt, gkns->src.source_host, fileNm);
    app->stat.n_neut_io += 1;

    gkyl_msgpack_data_release(mt);
    app->stat.neut_species_io_tm += gkyl_time_diff_now_sec(wtm);
  }
}

void
gk_neut_species_source_write_mom(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, double tm, int frame
)
{
  if (gkns->src.source_id && (gkns->src.evolve || frame == 0)) {
    // Package metadata.
    gkyl_msgpack_map_elem_set_double(gkns->io_meta_conf_len, gkns->io_meta_conf, "time", tm);
    gkyl_msgpack_map_elem_set_uint(gkns->io_meta_conf_len, gkns->io_meta_conf, "frame", frame);
    struct gkyl_msgpack_map_elem desc_src_mom[] = {
      {
        .key = "Description",
        .elem_type = GKYL_MP_STRING,
        .cval = "Velocity-space moment of the neutral species source.",
      },
    };
    int io_meta_len[] = {gkns->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
    const struct gkyl_msgpack_map_elem *io_meta[] = {
      gkns->io_meta_conf, app->gk_geom->io_meta_basic, desc_src_mom
    };
    struct gkyl_msgpack_data *mt =
      gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

    for (int m = 0; m < gkns->info.num_diag_moments; ++m) {
      struct timespec wst = gkyl_wall_clock();
      gk_neut_species_moment_calc(
        app, &gkns->src.moms[m], gkns->local, app->local, gkns->src.source
      );
      app->stat.n_neut_mom += 1;

      // Rescale moment by inverse of Jacobian if needed.
      gk_species_moment_diag_jacobgeo_div(
        app, &gkns->src.moms[m], gkns->src.moms[m].marr, gkns->src.moms[m].marr
      );
      app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

      struct timespec wtm = gkyl_wall_clock();
      if (app->use_gpu) {
        gkyl_array_copy(gkns->src.moms[m].marr_host, gkns->src.moms[m].marr);
      }

      const char *fmt = "%s-%s_source_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(
        fmt, app->name, gkns->info.name,
        gkyl_distribution_moments_strs[gkns->info.source.diagnostics.diag_moments[m]], frame
      );
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(
        fileNm, sizeof fileNm, fmt, app->name, gkns->info.name,
        gkyl_distribution_moments_strs[gkns->info.source.diagnostics.diag_moments[m]], frame
      );

      gkyl_comm_array_write(
        app->comm, &app->grid, &app->local, mt, gkns->src.moms[m].marr_host, fileNm
      );
      app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
      app->stat.n_neut_diag_io += 1;
    }
    gkyl_msgpack_data_release(mt);

    app->stat.n_neut_diag += 1;
  }
}

void
gk_neut_species_source_calc_integrated_mom(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, double tm
)
{
  if (gkns->src.source_id && gkns->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    struct gk_integ_moms *im = &gkns->src.integ_moms;
    gk_neut_species_moment_calc(app, &im->ops[0], gkns->local, app->local, gkns->src.source);
    app->stat.n_neut_mom += 1;

    // Reduce to compute sum over whole domain, append to diagnostics.
    struct gkyl_array *arrs[] = {im->ops[0].marr};
    gk_integ_moms_reduce_append(app, im, arrs, tm);

    app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_neut_diag += 1;
  }
}

void
gk_neut_species_source_write_integrated_mom(gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns)
{
  if (gkns->src.source_id && gkns->src.evolve) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s_source_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name, "integrated_moms");
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name, "integrated_moms");

      if (gkns->src.integ_moms.is_first_write) {
        struct gkyl_msgpack_map_elem io_meta_phi[] = {
          {
            .key = "Description",
            .elem_type = GKYL_MP_STRING,
            .cval = "Volume integrated moments of the neutral source.",
          },
        };
        int io_meta_len[] = {gkns->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
        const struct gkyl_msgpack_map_elem *io_meta[] = {
          gkns->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi
        };
        struct gkyl_msgpack_data *mt =
          gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(gkns->src.integ_moms.diag[0], fileNm, mt);
        gkns->src.integ_moms.is_first_write = false;
        gkyl_msgpack_data_release(mt);
      } else {
        gkyl_dynvec_awrite(gkns->src.integ_moms.diag[0], fileNm);
      }
    }
    gkyl_dynvec_clear(gkns->src.integ_moms.diag[0]);

    app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_neut_diag_io += 1;
  }
}

// Release function
void
gk_neut_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  if (src->source_id) {
    gkyl_array_release(src->source);
    if (app->use_gpu) {
      gkyl_array_release(src->source_host);
    }
    if (src->source_id != GKYL_SOURCE_FROMFILE) {
      for (int k = 0; k < src->num_sources; k++) {
        gk_neut_species_projection_release(app, &src->proj_source[k]);
      }
    }

    // Release moment data.
    for (int i = 0; i < src->num_diag_mom; ++i) {
      gk_neut_species_moment_release(app, &src->moms[i]);
    }
    gkyl_free(src->moms);
    gk_integ_moms_release(app, &src->integ_moms);
  }
}
