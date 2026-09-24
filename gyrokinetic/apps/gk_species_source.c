#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void
gk_species_source_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // Empty.
}

void
gk_species_source_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();
  // Package metadata.
  gkyl_msgpack_map_elem_set_double(gks->io_meta_phase_len, gks->io_meta_phase, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_phase_len, gks->io_meta_phase, "frame", frame);
  struct gkyl_msgpack_map_elem io_meta_source[] = {
    { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Source." }
  };
  int io_meta_len[] = {gks->io_meta_phase_len, app->gk_geom->io_meta_basic_len, 1};
  const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_phase, app->gk_geom->io_meta_basic, io_meta_source};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

  // Write out the source distribution function
  const char *fmt = "%s-%s_source_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, frame);

  // Copy data from device to host before writing it out.
  if (app->use_gpu)
    gkyl_array_copy(gks->src.source_host, gks->src.source);

  gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->src.source_host, fileNm);
  app->stat.n_io += 1;

  gkyl_msgpack_data_release(mt); 
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_source_write_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_source_write_enabled(app, gks, tm, frame);
  gks->src.write_func = gk_species_source_write_disabled;
}

void
gk_species_source_write_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  // Empty.
}

void
gk_species_source_write_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  struct timespec wst = gkyl_wall_clock();

  // Update time/frame metadata (shared fields, set once before the loop).
  gkyl_msgpack_map_elem_set_double(gks->io_meta_conf_len, gks->io_meta_conf, "time", tm);
  gkyl_msgpack_map_elem_set_uint(gks->io_meta_conf_len, gks->io_meta_conf, "frame", frame);

  for (int m=0; m<gks->src.num_diag_mom; ++m) {
    struct timespec wst = gkyl_wall_clock();
    gk_species_moment_calc(&gks->src.moms[m], gks->local, app->local, gks->src.source);
    app->stat.n_mom += 1;

    // Rescale moment by inverse of Jacobian if needed.
    gk_species_moment_diag_jacobgeo_div(app, &gks->src.moms[m], gks->src.moms[m].marr, gks->src.moms[m].marr);
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

    struct timespec wtm = gkyl_wall_clock();
    if (app->use_gpu)
      gkyl_array_copy(gks->src.moms[m].marr_host, gks->src.moms[m].marr);

    struct gkyl_msgpack_map_elem io_meta_mom[] = {
      { .key = "Description", .elem_type = GKYL_MP_STRING,
        .cval = (char*)gkyl_distribution_moments_descriptions[gks->info.source.diagnostics.diag_moments[m]] }
    };
    int io_meta_len[] = {gks->io_meta_conf_len, app->gk_geom->io_meta_basic_len, 1};
    const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_conf, app->gk_geom->io_meta_basic, io_meta_mom};
    struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

    const char *fmt = "%s-%s_source_%s_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
      gkyl_distribution_moments_strs[gks->info.source.diagnostics.diag_moments[m]], frame);
    char fileNm[sz+1]; // Ensures no buffer overflow.
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
      gkyl_distribution_moments_strs[gks->info.source.diagnostics.diag_moments[m]], frame);

    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt,
      gks->src.moms[m].marr_host, fileNm);
    gkyl_msgpack_data_release(mt);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_diag_io += 1;
  }

  app->stat.n_diag += 1;
}

void
gk_species_source_write_mom_init_only(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gk_species_source_write_mom_enabled(app, gks, tm, frame);
  gks->src.write_mom_func = gk_species_source_write_mom_disabled;
}

void
gk_species_source_calc_integrated_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  // Empty.
}

void
gk_species_source_calc_integrated_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock();

  double tm_prev = gkyl_dynvec_getlast_tm(gks->src.integ_diag);

  int num_mom = gks->src.integ_moms.num_mom;
  double avals_global[num_mom];

  gk_species_moment_calc(&gks->src.integ_moms, gks->local, app->local, gks->src.source); 
  app->stat.n_mom += 1;

  // Reduce to compute sum over whole domain, append to diagnostics
  gkyl_array_reduce_range(gks->src.red_integ_diag, gks->src.integ_moms.marr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
    gks->src.red_integ_diag, gks->src.red_integ_diag_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(avals_global, gks->src.red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  }
  else {
    memcpy(avals_global, gks->src.red_integ_diag_global, sizeof(double[num_mom]));
  }

  if (gks->info.source.diagnostics.time_integrated) {
    // This assumes time-independent sources. For time dependent ones
    // step the source contributions in RK3 like we do for boundary fluxes.
    double avals_global_prev[num_mom];
    for (int k=0; k<num_mom; k++)
      avals_global_prev[k] = 0.0;
    gkyl_dynvec_getlast(gks->src.integ_diag, avals_global_prev);
  
    double tau = tm - tm_prev;
    for (int k=0; k<num_mom; k++)
      avals_global[k] = avals_global_prev[k] + tau*avals_global[k];
  }
  gkyl_dynvec_append(gks->src.integ_diag, tm, avals_global);

  if (gks->src.num_adapt_sources > 0) {

    // Adaptive sources diagnostic.
    double temp[gks->src.num_adapt_sources];
    double part[gks->src.num_adapt_sources];
    for (int k=0; k < gks->src.num_adapt_sources; ++k) {
      temp[k] = gks->src.adapt[k].temperature_curr;
      part[k] = gks->src.adapt[k].particle_src_curr;
    }
    gkyl_dynvec_append(gks->src.part_diag, tm, part);
    gkyl_dynvec_append(gks->src.temp_diag, tm, temp);
  }

  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

void
gk_species_source_write_integrated_mom_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  // empty
}

void
gk_species_source_write_integrated_mom_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  struct timespec wst = gkyl_wall_clock();

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // Write out integrated diagnostic moments.
    const char *fmt = "%s-%s_source_%s.gkyl";

    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
    char fileNm[sz+1]; // Ensures no buffer overflow.
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

    if (gks->src.is_first_integ_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Volume integrated moments of the source." }
      };
      int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(gks->src.integ_diag, fileNm, mt);
      gks->src.is_first_integ_write_call = false;
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(gks->src.integ_diag, fileNm);
    }
  }
  gkyl_dynvec_clear(gks->src.integ_diag);

  if (gks->src.num_adapt_sources > 0) {
    if (rank == 0) {
      // Write out particle diagnostics for adaptive sources.
      const char *fmt = "%s-%s_adapt_sources_%s.gkyl";

      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "particle");
      char fileNm[sz+1]; // Ensures no buffer overflow.
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "particle");

      if (gks->src.is_first_integ_write_call_adapt) {
        struct gkyl_msgpack_map_elem io_meta_phi[] = {
          { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Source particle injection rate." }
        };
        int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
        const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi};
        struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);
  
        gkyl_dynvec_write_wmeta(gks->src.part_diag, fileNm, mt);
        gkyl_msgpack_data_release(mt);
      }
      else {
        gkyl_dynvec_awrite(gks->src.part_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(gks->src.part_diag);

    if (rank == 0) {
      // Write out temperature diagnostics for adaptive sources.
      const char *fmt = "%s-%s_adapt_sources_%s.gkyl";

      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "temperature");
      char fileNm[sz+1]; // Ensures no buffer overflow.
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "temperature");

      if (gks->src.is_first_integ_write_call_adapt) {
        struct gkyl_msgpack_map_elem io_meta_phi[] = {
          { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Source temperature." }
        };
        int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 1};
        const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_phi};
        struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(gks->src.temp_diag, fileNm, mt);
        gkyl_msgpack_data_release(mt);
      }
      else {
        gkyl_dynvec_awrite(gks->src.temp_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(gks->src.temp_diag);
    gks->src.is_first_integ_write_call_adapt = false;
  }

  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag_io += 1;
}


void
gk_species_source_calc(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, struct gkyl_array *f_buffer, double tm)
{
  if (src->source_id) {
    gkyl_array_clear(src->source, 0.0);
    for (int k=0; k<s->info.source.num_sources; k++) {
      gk_species_projection_calc(app, s, &src->proj_source[k], f_buffer, tm);
      gkyl_array_accumulate(src->source, 1., f_buffer);
    }
  }
}

void
gk_species_source_adapt(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, struct gkyl_array *f_buffer, struct gkyl_array **bflux_moms[], double tm) 
{  
  struct timespec wst = gkyl_wall_clock();

  src->adapt_func(app, s, src, f_buffer, bflux_moms, tm);
  
  app->stat.species_src_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_species_source_adapt_disabled(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, struct gkyl_array *f_buffer, struct gkyl_array **bflux_moms[], double tm)
{
  // Empty.
}

static void
gk_species_source_adapt_enabled(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, struct gkyl_array *f_buffer, struct gkyl_array **bflux_moms[], double tm)
{
  for (int k=0; k < s->info.source.num_adapt_sources; ++k) {
    struct gk_adapt_source *adapt_src = &src->adapt[k];
    struct gk_species *s_adapt = adapt_src->adapt_species;
    struct gkyl_array **s_adapt_bflux_moms = bflux_moms[adapt_src->adapt_species_idx];

    // Accumulate energy and particle losses through the boundaries.
    double sum_particle_loss_local = 0.0;
    double sum_energy_loss_local = 0.0;
    for (int j=0; j < adapt_src->num_boundaries; ++j) {

      double integ_m0_local_j, integ_m2_local_j; // Integrated boundary flux moments summed over boundaries.

      gk_species_bflux_get_flux_mom(&s_adapt->bflux, adapt_src->dir[j], adapt_src->edge[j],
        GKYL_F_MOMENT_M0, s_adapt_bflux_moms, adapt_src->bflux_m0, &adapt_src->boundaries_conf_ghost[j]);
      gk_species_bflux_get_flux_mom(&s_adapt->bflux, adapt_src->dir[j], adapt_src->edge[j],
        GKYL_F_MOMENT_M2, s_adapt_bflux_moms, adapt_src->bflux_m2, &adapt_src->boundaries_conf_ghost[j]);

      // Integrate the boundary flux moments to get the total loss through the j-th boundary.
      gkyl_array_integrate_advance(adapt_src->integrate_op, adapt_src->bflux_m0, 1.0, 0,
        &adapt_src->boundaries_conf_ghost[j], 0, adapt_src->integ_m0);
      gkyl_array_integrate_advance(adapt_src->integrate_op, adapt_src->bflux_m2, 1.0, 0,
        &adapt_src->boundaries_conf_ghost[j], 0, adapt_src->integ_m2);
      if (app->use_gpu) {
        gkyl_cu_memcpy(&integ_m0_local_j, adapt_src->integ_m0, sizeof(double), GKYL_CU_MEMCPY_D2H);
        gkyl_cu_memcpy(&integ_m2_local_j, adapt_src->integ_m2, sizeof(double), GKYL_CU_MEMCPY_D2H);
      } else {
        memcpy(&integ_m0_local_j, adapt_src->integ_m0, sizeof(double));
        memcpy(&integ_m2_local_j, adapt_src->integ_m2, sizeof(double));
      }
      // Accumulate to the total losses.
      sum_particle_loss_local += integ_m0_local_j; // n
      sum_energy_loss_local += 0.5 * s_adapt->info.mass * integ_m2_local_j; // 1/2 * m * v^2
    }
    // Sum over all MPI processes.
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &sum_particle_loss_local, &adapt_src->particle_rate_loss);
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &sum_energy_loss_local, &adapt_src->energy_rate_loss);

    double particle_input = s->info.source.projection[k].total_num_particles;
    double energy_input = s->info.source.projection[k].total_kin_energy;

    // Particle and energy rate update.
    // balance = user target + loss
    double density_compensation = adapt_src->adapt_particle_fraction * adapt_src->particle_rate_loss;
    double particle_src_new = adapt_src->adapt_particle? 
      particle_input + density_compensation : particle_input;

    double energy_compensation = adapt_src->adapt_energy_fraction * adapt_src->energy_rate_loss;
    double energy_src_new = adapt_src->adapt_energy?
      energy_input + energy_compensation : energy_input;

    // Avoid negative particle source.
    // This is important to avoid division by zero in the temperature calculation.
    particle_src_new = fmax(particle_input, particle_src_new);
    
    // Compute the target temperature of the source following the rule:
    // T = 2/3 * Q/G (T: src temperature, Q: src energy rate, G: total particle rate)
    const double vdim_phys = s->info.vdim == 1? 1.0 : 3.0;
    double temperature_new = (2./vdim_phys) * energy_src_new/particle_src_new;

    // Impose the temperature to be within the limits.  
    temperature_new = fmin(temperature_new, s->info.source.projection[k].temp_max);
    temperature_new = fmax(temperature_new, s->info.source.projection[k].temp_min);

    // Update the density and temperature moments of the source
    gkyl_array_clear(src->proj_source[k].prim_moms, 0.0);
    gkyl_array_set_offset(src->proj_source[k].prim_moms, particle_src_new, src->proj_source[k].gaussian_profile, 0*app->basis.num_basis);
    // The parallel velocity is left to be 0
    double dg_norm = pow(sqrt(2.0), app->cdim);
    gkyl_array_shiftc(src->proj_source[k].prim_moms, dg_norm * temperature_new / s->info.mass, 2*app->basis.num_basis);

    // Refresh the current values of particle, energy and temperature (can be used for control).
    adapt_src->particle_src_curr = particle_src_new;
    adapt_src->energy_src_curr = energy_src_new;
    adapt_src->temperature_curr = temperature_new;
  }

  // Reproject the source
  gk_species_source_calc(app, s, &s->src, f_buffer, tm);
}

static void 
gk_species_source_adapt_after_first_step(gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src, struct gkyl_array *f_buffer, struct gkyl_array **bflux_moms[], double tm)
{
  // Re-point to the function that adapts the source, so it doesn't happen the
  // first time we call gk_species_source_adapt.
  src->adapt_func = gk_species_source_adapt_enabled;
}

void 
gk_species_source_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gk_source *src)
{
  src->source_id = s->info.source.source_id;

  // Default function pointers.
  src->write_func = gk_species_source_write_disabled;
  src->write_mom_func = gk_species_source_write_mom_disabled;
  src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_disabled;
  src->write_integrated_mom_func = gk_species_source_write_integrated_mom_disabled;
  src->adapt_func = gk_species_source_adapt_disabled;

  if (src->source_id) {
    // Allocate source array.
    src->source = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
    src->source_host = src->source;
    if (app->use_gpu) {
      src->source_host = mkarr(false, src->source->ncomp, src->source->size); 
    }

    src->evolve = s->info.source.evolve || s->info.source.num_adapt_sources > 0; // Whether the source is time dependent.

    src->num_sources = s->info.source.num_sources;
    for (int k=0; k<s->info.source.num_sources; k++)
      gk_species_projection_init(app, s, s->info.source.projection[k], &src->proj_source[k]);

    // Allocate data and updaters for diagnostic moments.
    src->num_diag_mom = s->info.source.diagnostics.num_diag_moments;
    if (src->num_diag_mom == 0) {
      src->num_diag_mom = s->info.num_diag_moments;
      for (int m=0; m<src->num_diag_mom; ++m)
        s->info.source.diagnostics.diag_moments[m] = s->info.diag_moments[m];
    }

    src->moms = gkyl_malloc(sizeof(struct gk_species_moment[src->num_diag_mom]));
    for (int m=0; m<src->num_diag_mom; ++m)
      gk_species_moment_init(app, s, &src->moms[m], s->info.source.diagnostics.diag_moments[m], false);

    // Allocate data and updaters for integrated moments.
    src->num_diag_int_mom = s->info.source.diagnostics.num_integrated_diag_moments;
    assert(src->num_diag_int_mom < 2); // 1 int moment allowed now.
    if (src->evolve || src->num_diag_int_mom > 0) {
      gk_species_moment_init(app, s, &src->integ_moms,
        src->num_diag_int_mom == 0? GKYL_F_MOMENT_M0M1M2PARM2PERP : s->info.source.diagnostics.integrated_diag_moments[0], true);
      int num_mom = src->integ_moms.num_mom;
      if (app->use_gpu) {
        src->red_integ_diag = gkyl_cu_malloc(sizeof(double[num_mom]));
        src->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[num_mom]));
      } 
      else {
        src->red_integ_diag = gkyl_malloc(sizeof(double[num_mom]));
        src->red_integ_diag_global = gkyl_malloc(sizeof(double[num_mom]));
      }
      // Allocate dynamic-vector to store all-reduced integrated moments.
      src->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, num_mom);
      src->is_first_integ_write_call = true;
    }
    
    // Set function pointers chosen at runtime.
    if (src->evolve) {
      src->write_func = gk_species_source_write_enabled;
      src->write_mom_func = gk_species_source_write_mom_enabled;
      src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_enabled;
      src->write_integrated_mom_func = gk_species_source_write_integrated_mom_enabled;
    }
    else {
      src->write_func = gk_species_source_write_init_only;
      src->write_mom_func = gk_species_source_write_mom_init_only;
      if (src->num_diag_int_mom > 0) {
        // User requested integrated diagnostics.
        src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_enabled;
        src->write_integrated_mom_func = gk_species_source_write_integrated_mom_enabled;
      } else {
        src->calc_integrated_mom_func = gk_species_source_calc_integrated_mom_disabled;
        src->write_integrated_mom_func = gk_species_source_write_integrated_mom_disabled;
      }
    }

    // Set up the adaptive source.
    src->num_adapt_sources = s->info.source.num_adapt_sources;
    assert(src->num_adapt_sources <= src->num_sources); // Adaptive source should be a subset of the sources.
    if (src->num_adapt_sources > 0){
      assert(s->info.vdim > 1); // MF 2025/10/24: hasn't been tested in 1v.
      src->adapt_func = gk_species_source_adapt_after_first_step;

      if (src->num_diag_int_mom > 0){
        // Allocate dynvecs to store the temperature and particle count diagnostics of the adaptive sources.
        src->temp_diag = gkyl_dynvec_new(GKYL_DOUBLE, src->num_adapt_sources);
        src->part_diag = gkyl_dynvec_new(GKYL_DOUBLE, src->num_adapt_sources);
        src->is_first_integ_write_call_adapt = true;
      }

      for (int k = 0; k < src->num_adapt_sources; ++k) {
        // Adaptive source must be a Maxwellian Gaussian projection.
        assert(src->proj_source[k].proj_id == GKYL_PROJ_MAXWELLIAN_GAUSSIAN);

        struct gk_adapt_source *adapt_src = &src->adapt[k];

        adapt_src->adapt_particle = s->info.source.adapt[k].adapt_particle;
        adapt_src->adapt_energy = s->info.source.adapt[k].adapt_energy;
        adapt_src->adapt_particle_fraction = s->info.source.adapt[k].has_adapt_particle_fraction ?
          s->info.source.adapt[k].adapt_particle_fraction : 1.0; // Default to full adaptation if not specified.
        adapt_src->adapt_energy_fraction = s->info.source.adapt[k].has_adapt_energy_fraction ?
          s->info.source.adapt[k].adapt_energy_fraction : 1.0; // Default to full adaptation if not specified.

        adapt_src->adapt_species = gk_find_species(app, s->info.source.adapt[k].adapt_to_species);
        assert(adapt_src->adapt_species != NULL); // Make sure the adaptive species is found.
        adapt_src->adapt_species_idx = gk_find_species_idx(app, adapt_src->adapt_species->info.name);

        adapt_src->particle_src_curr = s->info.source.projection[k].total_num_particles;
        adapt_src->energy_src_curr = s->info.source.projection[k].total_kin_energy;
        // The temperature computation makes sense only if we inject particles.
        const double vdim_phys = s->info.vdim == 1? 1.0 : 3.0;
        adapt_src->temperature_curr = s->info.source.projection[k].total_num_particles > 0?
          (2./vdim_phys) * adapt_src->energy_src_curr/adapt_src->particle_src_curr : s->info.source.projection[k].temp_min;

        gk_species_moment_init(app, adapt_src->adapt_species, &adapt_src->integ_threemoms, GKYL_F_MOMENT_M0M1M2, true);

        // Initialize the infrastructure to compute integrated moments of the boundary fluxes.
        adapt_src->bflux_m0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
        adapt_src->bflux_m2 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
        adapt_src->integrate_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);

        int num_mom = adapt_src->integ_threemoms.num_mom;
        if (app->use_gpu){
          adapt_src->red_integ_mom = gkyl_cu_malloc(sizeof(double[num_mom]));
          adapt_src->red_integ_mom_global = gkyl_cu_malloc(sizeof(double[num_mom]));
          adapt_src->integ_m0 = gkyl_cu_malloc(sizeof(double));
          adapt_src->integ_m2 = gkyl_cu_malloc(sizeof(double));
        }
        else {
          adapt_src->red_integ_mom = gkyl_malloc(sizeof(double[num_mom]));
          adapt_src->red_integ_mom_global = gkyl_malloc(sizeof(double[num_mom]));
          adapt_src->integ_m0 = gkyl_malloc(sizeof(double));
          adapt_src->integ_m2 = gkyl_malloc(sizeof(double));
        }

        adapt_src->num_boundaries = s->info.source.adapt[k].num_boundaries;
        bool is_dir_periodic[GKYL_MAX_CDIM] = {0};
        for (int j=0; j < app->num_periodic_dir; ++j) {
          is_dir_periodic[app->periodic_dirs[j]] = true;
        }
        for (int j=0; j < adapt_src->num_boundaries; ++j) {
          int dir = s->info.source.adapt[k].dir[j];
          int edge = s->info.source.adapt[k].edge[j];

          // Source adaptation on periodic, zero flux, or reflect boundary is not allowed.
          assert(is_dir_periodic[dir] == 0);
          if (edge == GKYL_LOWER_EDGE) {
            assert(s->lower_bc[dir].type != GKYL_BC_GK_SPECIES_ZERO_FLUX);
            assert(s->lower_bc[dir].type != GKYL_BC_GK_SPECIES_REFLECT);
          } else {                                     
            assert(s->upper_bc[dir].type != GKYL_BC_GK_SPECIES_ZERO_FLUX);
            assert(s->upper_bc[dir].type != GKYL_BC_GK_SPECIES_REFLECT);
          }

          // Default scenario: we set the ranges to the full range of the ghost cells.
          adapt_src->boundaries_phase_ghost[j] = edge == GKYL_LOWER_EDGE? s->local_lower_ghost[dir]
                                                                        : s->local_upper_ghost[dir];
          adapt_src->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE? app->local_lower_ghost[dir]
                                                                       : app->local_upper_ghost[dir];
          adapt_src->dir[j]  = dir;
          adapt_src->edge[j] = edge;

          // Specific scenario if we are in a inner wall limited case. We select only SOL range in parallel direction.
          if (dir == app->cdim-1 && app->gk_geom->has_LCFS) {
            adapt_src->boundaries_phase_ghost[j] = edge == GKYL_LOWER_EDGE? s->local_lower_ghost_par_sol
                                                                          : s->local_upper_ghost_par_sol;
            adapt_src->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE? app->local_lower_ghost_par_sol
                                                                         : app->local_upper_ghost_par_sol;
          }
        }
      }
    }
  }
}

void
gk_species_source_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_source *src, const struct gkyl_array *fin, struct gkyl_array *rhs)
{

  struct timespec wst = gkyl_wall_clock();
  if (src->source_id) {
    gkyl_array_accumulate(rhs, 1.0, src->source);
  }
  app->stat.species_src_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_source_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->src.write_func(app, gks, tm, frame);
}

void
gk_species_source_write_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  gks->src.write_mom_func(app, gks, tm, frame);
}

void
gk_species_source_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  gks->src.calc_integrated_mom_func(app, gks, tm);
}

void
gk_species_source_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  gks->src.write_integrated_mom_func(app, gks);
}

void
gk_species_source_release(const struct gkyl_gyrokinetic_app *app, const struct gk_source *src)
{
  if (src->source_id) {
    gkyl_array_release(src->source);
    if (app->use_gpu) {
      gkyl_array_release(src->source_host);
    }
    for (int k=0; k<src->num_sources; k++) {
      gk_species_projection_release(app, &src->proj_source[k]);
    }

    // Release moment data.
    for (int i=0; i<src->num_diag_mom; ++i) {
      gk_species_moment_release(app, &src->moms[i]);
    }
    gkyl_free(src->moms);

    if (src->evolve || src->num_diag_int_mom > 0) {
      gk_species_moment_release(app, &src->integ_moms); 
      if (app->use_gpu) {
        gkyl_cu_free(src->red_integ_diag);
        gkyl_cu_free(src->red_integ_diag_global);
      }
      else {
        gkyl_free(src->red_integ_diag);
        gkyl_free(src->red_integ_diag_global);
      }  
      gkyl_dynvec_release(src->integ_diag);
    }
    if (src->num_adapt_sources > 0) {
      if (src->num_diag_int_mom > 0){
        gkyl_dynvec_release(src->part_diag);
        gkyl_dynvec_release(src->temp_diag);
      }
      for (int k=0; k < src->num_adapt_sources; ++k) {
        const struct gk_adapt_source *adapt_src = &src->adapt[k];
        gk_species_moment_release(app, &adapt_src->integ_threemoms);
        gkyl_array_integrate_release(adapt_src->integrate_op);
        gkyl_array_release(adapt_src->bflux_m0);
        gkyl_array_release(adapt_src->bflux_m2);
        if (app->use_gpu) {
          gkyl_cu_free(adapt_src->red_integ_mom);
          gkyl_cu_free(adapt_src->red_integ_mom_global);
          gkyl_cu_free(adapt_src->integ_m0);
          gkyl_cu_free(adapt_src->integ_m2);
        }
        else {
          gkyl_free(adapt_src->red_integ_mom);
          gkyl_free(adapt_src->red_integ_mom_global);
          gkyl_free(adapt_src->integ_m0);
          gkyl_free(adapt_src->integ_m2);
        }
      }
    }
  }
}
