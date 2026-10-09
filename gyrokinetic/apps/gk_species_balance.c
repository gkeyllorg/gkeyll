#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

static void
gk_species_balance_step_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, double dt, bool is_old)
{
}

static void
gk_species_balance_step_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, double dt, bool is_old)
{
  // Store dt and the integrated M0 of f/dt at the beginning (is_old=true)
  // or end of the time step.
  bal->dt = dt;
  gk_species_moment_calc(&bal->f_mom, gks->local, app->local, gks->f);
  gkyl_array_set(is_old? bal->mom_old : bal->mom_new, 1.0/dt, bal->f_mom.marr);
}

static double
gk_species_balance_reduce(gkyl_gyrokinetic_app* app, struct gk_species_balance *bal,
  const struct gkyl_array *arr)
{
  // Sum a single component array over the local range and across ranks.
  double out;
  gkyl_array_reduce_range(bal->red_local, arr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, bal->red_local, bal->red_global);
  if (app->use_gpu)
    gkyl_cu_memcpy(&out, bal->red_global, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&out, bal->red_global, sizeof(double));
  return out;
}

static void
gk_species_balance_calc_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, double tm)
{
}

static void
gk_species_balance_calc_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, double tm)
{
  double row[bal->num_comp];
  for (int k=0; k<bal->num_comp; k++) row[k] = 0.0;

  // Particle number.
  gk_species_moment_calc(&bal->f_mom, gks->local, app->local, gks->f);
  double N = gk_species_balance_reduce(app, bal, bal->f_mom.marr);

  if (bal->dt > 0.0) {
    // Rates are only defined once a step has been taken.
    gkyl_array_set(bal->f_mom.marr, 1.0, bal->mom_new);
    gkyl_array_accumulate(bal->f_mom.marr, -1.0, bal->mom_old);
    double dNdt = gk_species_balance_reduce(app, bal, bal->f_mom.marr);

    double src = 0.0;
    if (bal->has_source) {
      gk_species_moment_calc(&bal->src_mom, gks->local, app->local, gks->src.source);
      src = gk_species_balance_reduce(app, bal, bal->src_mom.marr);
    }

    // Total flux through the boundaries.
    struct gk_boundary_fluxes *bflux = &gks->bflux;
    double bflux_tot_local = 0.0, bflux_tot = 0.0;
    for (int b=0; b<bal->num_boundaries; ++b) {
      int dir = bflux->boundaries_dir[b];
      enum gkyl_edge_loc edge = bflux->boundaries_edge[b];
      // Local ghost ranges also exist at inter-rank boundaries: only count the
      // contribution of ranks that touch the boundary.
      bool at_boundary = edge == GKYL_LOWER_EDGE?
        app->local.lower[dir] == app->global.lower[dir] :
        app->local.upper[dir] == app->global.upper[dir];
      if (at_boundary) {
        gk_species_bflux_get_flux_mom(bflux, dir, edge, bal->mom_type,
          bflux->f, bal->bflux_mom, bflux->boundaries_conf_ghost[b]);
        gkyl_array_integrate_advance(bal->integ_op, bal->bflux_mom, 1.0, 0,
          bflux->boundaries_conf_ghost[b], 0, bal->red_local);
        double flux_b;
        if (app->use_gpu)
          gkyl_cu_memcpy(&flux_b, bal->red_local, sizeof(double), GKYL_CU_MEMCPY_D2H);
        else
          memcpy(&flux_b, bal->red_local, sizeof(double));
        bflux_tot_local += flux_b;
      }
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &bflux_tot_local, &bflux_tot);

    // Particles gained/lost in the last step, also relative to N (round-off ~1e-16 if conservative).
    double mom_err = src - bflux_tot - dNdt;
    row[0] = dNdt;
    row[1] = src;
    row[2] = bflux_tot;
    row[3] = mom_err;
    row[4] = mom_err*bal->dt/N;
  }

  gkyl_dynvec_append(bal->diag, tm, row);
}

static void
gk_species_balance_write_disabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, bool is_first_write)
{
}

static void
gk_species_balance_write_enabled(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, bool is_first_write)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    // File named after the moment, e.g. <name>-<species>_balance_M0.gkyl.
    const char *mom_name = gkyl_distribution_moments_strs[bal->mom_type];
    const char *fmt = "%s-%s_balance_%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, mom_name);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, mom_name);

    if (is_first_write) {
      struct gkyl_msgpack_map_elem io_meta_bal[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING,
          .cval = "Species moment balance." },
        { .key = "Columns", .elem_type = GKYL_MP_STRING,
          .cval = "mom_dot, src_mom, bflux_tot_mom, mom_err, mom_err_norm" },
      };
      int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 2};
      const struct gkyl_msgpack_map_elem* io_meta[] = {gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_bal};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(bal->diag, fileNm, mt);
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(bal->diag, fileNm);
    }
  }
  gkyl_dynvec_clear(bal->diag);
  app->stat.n_diag_io += 1;
}

void
gk_species_balance_pre_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_species_balance *bal, enum gkyl_species_bflux_type *bflux_type,
  struct gkyl_phase_diagnostics_inp *add_bflux_moms_inp)
{
  // Set up what must exist before the bfluxes and the df/dt diagnostics are initialized.
  *bal = (struct gk_species_balance) { };
  bal->step_func = gk_species_balance_step_disabled;
  bal->calc_func = gk_species_balance_calc_disabled;
  bal->write_func = gk_species_balance_write_disabled;

  // Only the particle (M0) balance is available for now.
  assert(gks->info.num_balance_moments == 0 ||
    (gks->info.num_balance_moments == 1 && gks->info.balance_moments[0] == GKYL_F_MOMENT_M0));
  bal->enabled = gks->info.num_balance_moments > 0 && !gks->info.is_static;
  if (!bal->enabled)
    return;

  bal->mom_type = gks->info.balance_moments[0];

  // dN/dt and dt come from the df/dt diagnostics.
  gks->info.time_rate_diagnostics = true;

  // Step the balance moment of the boundary fluxes (same boundaries as gk_species_bflux_init).
  int num_bound = 0;
  for (int d=0; d<app->cdim; ++d) {
    if (gks->bc_is_np[d]) {
      if (gks->lower_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX) num_bound++;
      if (gks->upper_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX) num_bound++;
    }
  }
  if (num_bound > 0) {
    // Activate the step M0 moment calculation is needed.
    if (*bflux_type == GK_SPECIES_BFLUX_NONE || *bflux_type == GK_SPECIES_BFLUX_CALC_FLUX)
      *bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS;
    add_bflux_moms_inp->diag_moments[add_bflux_moms_inp->num_diag_moments++] = bal->mom_type;
  }
}

void
gk_species_balance_init(struct gkyl_gyrokinetic_app *app, struct gk_species *gks,
  struct gk_species_balance *bal)
{
  if (!bal->enabled)
    return;

  bal->num_boundaries = gks->bflux.num_boundaries;
  bal->num_comp = 5; // fdot, src, bflux_tot, mom_err, mom_err_norm=mom_err*dt/N.
  bal->diag = gkyl_dynvec_new(GKYL_DOUBLE, bal->num_comp);
  bal->dt = 0.0;

  gk_species_moment_init(app, gks, &bal->f_mom, bal->mom_type, true);
  bal->mom_old = mkarr(app->use_gpu, bal->f_mom.marr->ncomp, bal->f_mom.marr->size);
  bal->mom_new = mkarr(app->use_gpu, bal->f_mom.marr->ncomp, bal->f_mom.marr->size);

  bal->has_source = gks->info.source.source_id != 0;
  if (bal->has_source)
    gk_species_moment_init(app, gks, &bal->src_mom, bal->mom_type, true);

  if (bal->num_boundaries > 0) {
    bal->integ_op = gkyl_array_integrate_new(&app->grid, &app->basis, 1,
      GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu);
    bal->bflux_mom = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  if (app->use_gpu) {
    bal->red_local = gkyl_cu_malloc(sizeof(double));
    bal->red_global = gkyl_cu_malloc(sizeof(double));
  }
  else {
    bal->red_local = gkyl_malloc(sizeof(double));
    bal->red_global = gkyl_malloc(sizeof(double));
  }

  bal->step_func = gk_species_balance_step_enabled;
  bal->calc_func = gk_species_balance_calc_enabled;
  bal->write_func = gk_species_balance_write_enabled;
}

void
gk_species_balance_step(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, double dt, bool is_old)
{
  bal->step_func(app, gks, bal, dt, is_old);
}

void
gk_species_balance_calc(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, double tm)
{
  bal->calc_func(app, gks, bal, tm);
}

void
gk_species_balance_write(gkyl_gyrokinetic_app* app, struct gk_species *gks,
  struct gk_species_balance *bal, bool is_first_write)
{
  bal->write_func(app, gks, bal, is_first_write);
}

void
gk_species_balance_release(const struct gkyl_gyrokinetic_app *app,
  const struct gk_species_balance *bal)
{
  if (!bal->enabled)
    return;

  gkyl_dynvec_release(bal->diag);
  gk_species_moment_release(app, &bal->f_mom);
  gkyl_array_release(bal->mom_old);
  gkyl_array_release(bal->mom_new);
  if (bal->has_source)
    gk_species_moment_release(app, &bal->src_mom);
  if (bal->num_boundaries > 0) {
    gkyl_array_integrate_release(bal->integ_op);
    gkyl_array_release(bal->bflux_mom);
  }
  if (app->use_gpu) {
    gkyl_cu_free(bal->red_local);
    gkyl_cu_free(bal->red_global);
  }
  else {
    gkyl_free(bal->red_local);
    gkyl_free(bal->red_global);
  }
}
