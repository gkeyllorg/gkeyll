#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// Sum an array over the local range and across ranks, returning one component.
static double
gk_species_balance_reduce(
  gkyl_gyrokinetic_app *app, struct gk_species_balance *bal, const struct gkyl_array *arr, int comp
)
{
  double out;
  gkyl_array_reduce_range(bal->red_local, arr, GKYL_SUM, &app->local);
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, arr->ncomp, bal->red_local, bal->red_global);
  if (app->use_gpu) {
    gkyl_cu_memcpy(&out, bal->red_global + comp, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    out = bal->red_global[comp];
  }
  return out;
}

// Balance of a moment of f: d(mom)/dt = source - boundary fluxes, from the
// balance's own moments of f, of the source and of the stepped boundary-flux
// moments. For the energy of a single species the residual is the energy
// exchanged with the field; the sum over species closes with the field energy.
static void
gk_species_balance_calc_diag(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal, int d,
  double tm
)
{
  double bal_vals[bal->num_comp];
  for (int k = 0; k < bal->num_comp; k++) {
    bal_vals[k] = 0.0;
  }
  int comp = bal->mom_comp[d];

  // Balanced moment.
  gk_species_moment_calc(app, &bal->f_mom[d], gks->local, app->local, gks->f);
  double mom = gk_species_balance_reduce(app, bal, bal->f_mom[d].marr, comp);
  bal->mom[d] = mom;

  if (bal->dt > 0.0) {
    // Rates are only defined once a step has been taken.
    gkyl_array_set(bal->f_mom[d].marr, 1.0, bal->mom_new[d]);
    gkyl_array_accumulate(bal->f_mom[d].marr, -1.0, bal->mom_old[d]);
    double mom_dot = gk_species_balance_reduce(app, bal, bal->f_mom[d].marr, comp);

    double src = 0.0;
    if (bal->has_source) {
      gk_species_moment_calc(app, &bal->src_mom[d], gks->local, app->local, gks->src.source);
      src = gk_species_balance_reduce(app, bal, bal->src_mom[d].marr, comp);
    }
    // BGK source terms, combined over the RK stages of the last step.
    src += gk_species_balance_reduce(app, bal, bal->src_bgk[d], comp);

    // Total flux through the boundaries.
    struct gk_boundary_fluxes *bflux = &gks->bflux;
    double bflux_tot_local = 0.0, bflux_tot = 0.0;
    for (int b = 0; b < bal->num_boundaries; ++b) {
      int dir = bflux->boundaries_dir[b];
      enum gkyl_edge_loc edge = bflux->boundaries_edge[b];
      // Local ghost ranges also exist at inter-rank boundaries: only count the
      // contribution of ranks that touch the boundary.
      bool at_boundary = edge == GKYL_LOWER_EDGE ? app->local.lower[dir] == app->global.lower[dir] :
                                                   app->local.upper[dir] == app->global.upper[dir];
      if (at_boundary) {
        gk_species_bflux_get_flux_mom(
          bflux, dir, edge, bal->mom_types[d], bflux->f, bal->bflux_mom[d],
          bflux->boundaries_conf_ghost[b]
        );
        gkyl_array_integrate_advance(
          bal->integ_op[d], bal->bflux_mom[d], 1.0, 0, bflux->boundaries_conf_ghost[b], 0,
          bal->red_local
        );
        double flux_b;
        if (app->use_gpu) {
          gkyl_cu_memcpy(&flux_b, bal->red_local + comp, sizeof(double), GKYL_CU_MEMCPY_D2H);
        } else {
          flux_b = bal->red_local[comp];
        }
        bflux_tot_local += flux_b;
      }
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &bflux_tot_local, &bflux_tot);

    // Moment gained/lost in the last step, also relative to mom.
    double mom_err = src - bflux_tot - mom_dot;
    row[0] = mom_err;
    row[1] = mom_err * bal->dt / mom;
    row[2] = mom_dot;
    row[3] = src;
    row[4] = bflux_tot;
  }

  gkyl_dynvec_append(bal->diag[d], tm, row);
}

static void
gk_species_balance_rhs_begin_disabled(struct gk_species_balance *bal)
{
}

static void
gk_species_balance_rhs_begin_enabled(struct gk_species_balance *bal)
{
  // The first stage of a step may land in any slot (a failed stage restarts
  // the step): gk_species_balance_step moves it to the first slot.
  bal->rk_stage = GKYL_MIN2(bal->rk_stage + 1, 2);
  for (int d = 0; d < bal->num_diags; d++) {
    gkyl_array_clear(bal->src_bgk_stage[d][bal->rk_stage], 0.0);
  }
}

static void
gk_species_balance_source_bgk_rhs_disabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal,
  const struct gkyl_array *src_bgk_rhs
)
{
}

static void
gk_species_balance_source_bgk_rhs_enabled(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal,
  const struct gkyl_array *src_bgk_rhs
)
{
  for (int d = 0; d < bal->num_diags; d++) {
    gk_species_moment_calc(app, &bal->src_mom[d], gks->local, app->local, src_bgk_rhs);
    gkyl_array_accumulate(bal->src_bgk_stage[d][bal->rk_stage], 1.0, bal->src_mom[d].marr);
  }
}

void
gk_species_balance_pre_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal,
  enum gkyl_species_bflux_type *bflux_type, struct gkyl_phase_diagnostics_inp *add_bflux_moms_inp
)
{
  // Collect the balances requested by the species (written) and by the app
  // (only computed, for the sum over species), and request the stepping of
  // the boundary-flux moments they need.
  *bal = (struct gk_species_balance){};
  bal->rhs_begin_func = gk_species_balance_rhs_begin_disabled;
  bal->source_bgk_rhs_func = gk_species_balance_source_bgk_rhs_disabled;
  struct gkyl_gyrokinetic_species *info = &gks->info;
  for (int d = 0; d < info->num_balance_diags; d++) {
    bal->diag_types[d] = info->balance_diags[d];
    bal->write[d] = true;
    bal->num_diags++;
  }
  for (int d = 0; d < app->bal.num_diags; d++) {
    bool found = false;
    for (int k = 0; k < bal->num_diags; k++) {
      found = found || bal->diag_types[k] == app->bal.diag_types[d];
    }
    if (!found) {
      assert(bal->num_diags < 12);
      bal->diag_types[bal->num_diags] = app->bal.diag_types[d];
      bal->write[bal->num_diags] = false;
      bal->num_diags++;
    }
  }
  bal->enabled = bal->num_diags > 0 && !info->is_static;
  if (!bal->enabled) {
    bal->num_diags = 0;
    return;
  }

  for (int d = 0; d < bal->num_diags; d++) {
    switch (bal->diag_types[d]) {
      case GKYL_BALANCE_GK_PARTICLE:
        bal->mom_types[d] = GKYL_F_MOMENT_M0;
        bal->mom_comp[d] = 0;
        break;
      case GKYL_BALANCE_GK_ENERGY:
        bal->mom_types[d] = GKYL_F_MOMENT_HAMILTONIAN; // M0, m*M1 and H moments.
        bal->mom_comp[d] = 2;
        break;
      default:
        assert(false); // Balance diagnostic not available.
    }
  }

  bal->has_source = info->source.source_id != 0;

  // Boundaries with a boundary flux (same as in gk_species_bflux_init).
  bal->num_boundaries = 0;
  for (int d = 0; d < app->cdim; ++d) {
    if (gks->bc_is_np[d]) {
      if (gks->lower_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX) {
        bal->num_boundaries++;
      }
      if (gks->upper_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX) {
        bal->num_boundaries++;
      }
    }
  }
  if (bal->num_boundaries > 0) {
    // Step the balanced moments of the boundary fluxes (without writing them).
    if (*bflux_type == GK_SPECIES_BFLUX_NONE || *bflux_type == GK_SPECIES_BFLUX_CALC_FLUX) {
      *bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS;
    }
    for (int d = 0; d < bal->num_diags; d++) {
      bool found = false;
      for (int m = 0; m < add_bflux_moms_inp->num_diag_moments; m++) {
        found = found || add_bflux_moms_inp->diag_moments[m] == bal->mom_types[d];
      }
      if (!found) {
        add_bflux_moms_inp->diag_moments[add_bflux_moms_inp->num_diag_moments++] =
          bal->mom_types[d];
      }
    }
  }
}

void
gk_species_balance_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal
)
{
  if (!bal->enabled) {
    return;
  }

  bal->num_comp = 5; // mom_err, mom_err_norm=mom_err*dt/mom, mom_dot, src, bflux_tot.
  bal->dt = 0.0;
  int max_num_mom = 1;
  for (int d = 0; d < bal->num_diags; d++) {
    gk_species_moment_init(app, gks, &bal->f_mom[d], bal->mom_types[d], true);
    bal->mom_old[d] = mkarr(app->use_gpu, bal->f_mom[d].marr->ncomp, bal->f_mom[d].marr->size);
    bal->mom_new[d] = mkarr(app->use_gpu, bal->f_mom[d].marr->ncomp, bal->f_mom[d].marr->size);
    gk_species_moment_init(app, gks, &bal->src_mom[d], bal->mom_types[d], true);
    for (int k = 0; k < 3; k++) {
      bal->src_bgk_stage[d][k] =
        mkarr(app->use_gpu, bal->src_mom[d].marr->ncomp, bal->src_mom[d].marr->size);
    }
    bal->src_bgk[d] = mkarr(app->use_gpu, bal->src_mom[d].marr->ncomp, bal->src_mom[d].marr->size);
    gkyl_array_clear(bal->src_bgk[d], 0.0);
    if (bal->num_boundaries > 0) {
      bal->integ_op[d] = gkyl_array_integrate_new(
        &app->grid, &app->basis, bal->f_mom[d].num_mom, GKYL_ARRAY_INTEGRATE_OP_NONE, app->use_gpu
      );
      bal->bflux_mom[d] =
        mkarr(app->use_gpu, bal->f_mom[d].num_mom * app->basis.num_basis, app->local_ext.volume);
    }
    max_num_mom = GKYL_MAX2(max_num_mom, bal->f_mom[d].num_mom);
    bal->diag[d] = gkyl_dynvec_new(GKYL_DOUBLE, bal->num_comp);
    bal->is_first_write[d] = true;
  }
  bal->rk_stage = 0;
  bal->rhs_begin_func = gk_species_balance_rhs_begin_enabled;
  bal->source_bgk_rhs_func = gk_species_balance_source_bgk_rhs_enabled;

  if (app->use_gpu) {
    bal->red_local = gkyl_cu_malloc(max_num_mom * sizeof(double));
    bal->red_global = gkyl_cu_malloc(max_num_mom * sizeof(double));
  } else {
    bal->red_local = gkyl_malloc(max_num_mom * sizeof(double));
    bal->red_global = gkyl_malloc(max_num_mom * sizeof(double));
  }
}

void
gk_species_balance_step(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal, double dt,
  bool is_old
)
{
  // Store dt and the balanced moments of f/dt at the beginning (is_old=true)
  // or end of the time step.
  bal->dt = dt;
  for (int d = 0; d < bal->num_diags; d++) {
    gk_species_moment_calc(app, &bal->f_mom[d], gks->local, app->local, gks->f);
    gkyl_array_set(is_old ? bal->mom_old[d] : bal->mom_new[d], 1.0 / dt, bal->f_mom[d].marr);

    struct gkyl_array **stage = bal->src_bgk_stage[d];
    if (is_old) {
      // The last rhs evaluation is the first stage of this step.
      struct gkyl_array *first = stage[bal->rk_stage];
      stage[bal->rk_stage] = stage[0];
      stage[0] = first;
    } else {
      // Combine the BGK source moments of the stages as SSP-RK3 combines f.
      gkyl_array_set(bal->src_bgk[d], 1.0 / 6.0, stage[0]);
      gkyl_array_accumulate(bal->src_bgk[d], 1.0 / 6.0, stage[1]);
      gkyl_array_accumulate(bal->src_bgk[d], 2.0 / 3.0, stage[2]);
    }
  }
  if (is_old) {
    bal->rk_stage = 0;
  }
}

void
gk_species_balance_rhs_begin(struct gk_species_balance *bal)
{
  bal->rhs_begin_func(bal);
}

void
gk_species_balance_source_bgk_rhs(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal,
  const struct gkyl_array *src_bgk_rhs
)
{
  bal->source_bgk_rhs_func(app, gks, bal, src_bgk_rhs);
}

void
gk_species_balance_calc(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal, double tm
)
{
  for (int d = 0; d < bal->num_diags; d++) {
    gk_species_balance_calc_diag(app, gks, bal, d, tm);
  }
}

void
gk_species_balance_write(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal
)
{
  struct timespec wst = gkyl_wall_clock();
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  for (int d = 0; d < bal->num_diags; d++) {
    if (bal->write[d] && rank == 0) {
      // File named after the balance, e.g. <name>-<species>_balance_particle.gkyl.
      const char *fmt = "%s-%s_balance_%s.gkyl";
      const char *diag_name = gkyl_gyrokinetic_balance_diag_strs[bal->diag_types[d]];
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, diag_name);
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, diag_name);

      if (bal->is_first_write[d]) {
        struct gkyl_msgpack_map_elem io_meta_bal[] = {
          {
            .key = "Description",
            .elem_type = GKYL_MP_STRING,
            .cval = "Species balance: mom_err = src_mom - bflux_tot_mom - mom_dot, "
                    "mom_err_norm = mom_err*dt/mom.",
          },
          {
            .key = "Columns",
            .elem_type = GKYL_MP_STRING,
            .cval = "mom_err, mom_err_norm, mom_dot, src_mom, bflux_tot_mom",
          },
        };
        int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 2};
        const struct gkyl_msgpack_map_elem *io_meta[] = {
          gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_bal
        };
        struct gkyl_msgpack_data *mt =
          gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(bal->diag[d], fileNm, mt);
        gkyl_msgpack_data_release(mt);
        bal->is_first_write[d] = false;
      } else {
        gkyl_dynvec_awrite(bal->diag[d], fileNm);
      }
      app->stat.n_diag_io += 1;
    }
    gkyl_dynvec_clear(bal->diag[d]);
  }
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_balance_release(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_balance *bal
)
{
  if (!bal->enabled) {
    return;
  }

  for (int d = 0; d < bal->num_diags; d++) {
    gk_species_moment_release(app, &bal->f_mom[d]);
    gkyl_array_release(bal->mom_old[d]);
    gkyl_array_release(bal->mom_new[d]);
    gk_species_moment_release(app, &bal->src_mom[d]);
    for (int k = 0; k < 3; k++) {
      gkyl_array_release(bal->src_bgk_stage[d][k]);
    }
    gkyl_array_release(bal->src_bgk[d]);
    if (bal->num_boundaries > 0) {
      gkyl_array_integrate_release(bal->integ_op[d]);
      gkyl_array_release(bal->bflux_mom[d]);
    }
    gkyl_dynvec_release(bal->diag[d]);
  }
  if (app->use_gpu) {
    gkyl_cu_free(bal->red_local);
    gkyl_cu_free(bal->red_global);
  } else {
    gkyl_free(bal->red_local);
    gkyl_free(bal->red_global);
  }
}
