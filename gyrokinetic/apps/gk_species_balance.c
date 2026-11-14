#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// Component of the quantity balanced by a diagnostic in an integrated moment,
// -1 if the moment does not contain it.
static int
gk_species_balance_mom_comp(
  enum gkyl_gyrokinetic_balance_diag diag_type, enum gkyl_distribution_moments mom_type
)
{
  switch (diag_type) {
    case GKYL_BALANCE_GK_PARTICLE:
      // The density is the first component of these moments.
      return mom_type == GKYL_F_MOMENT_M0 || mom_type == GKYL_F_MOMENT_M0M1M2 ||
                 mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP ||
                 mom_type == GKYL_F_MOMENT_HAMILTONIAN ?
               0 :
               -1;
    case GKYL_BALANCE_GK_ENERGY:
      // The energy is the third component of the Hamiltonian moments.
      return mom_type == GKYL_F_MOMENT_HAMILTONIAN ? 2 : -1;
    default:
      assert(false); // Balance diagnostic not available.
      return -1;
  }
}

// Index of the first moment of a list containing the quantity balanced by a
// diagnostic, -1 if none does.
static int
gk_species_balance_find_mom(
  enum gkyl_gyrokinetic_balance_diag diag_type, int num_moms,
  const enum gkyl_distribution_moments *moms
)
{
  for (int m = 0; m < num_moms; m++) {
    if (gk_species_balance_mom_comp(diag_type, moms[m]) >= 0) {
      return m;
    }
  }
  return -1;
}

// Whether this rank counts the integrated bflux through the boundary at
// (dir, edge). The bflux diagnostics sum over the ranks of the boundary plane
// (comm_plane), so a single rank of the plane must contribute to the total.
static bool
gk_species_balance_owns_boundary(gkyl_gyrokinetic_app *app, int dir, enum gkyl_edge_loc edge)
{
  bool owns = edge == GKYL_LOWER_EDGE ? app->local.lower[dir] == app->global.lower[dir] :
                                        app->local.upper[dir] == app->global.upper[dir];
  for (int d = 0; d < app->cdim; d++) {
    if (d != dir) {
      owns = owns && app->local.lower[d] == app->global.lower[d];
    }
  }
  return owns;
}

// Last values appended to an integrated diagnostic, which must be computed at time tm.
static void
gk_species_balance_getlast(gkyl_dynvec diag, double tm, double *avals)
{
  bool has_last = gkyl_dynvec_getlast(diag, avals);
  assert(has_last && gkyl_dynvec_getlast_tm(diag) == tm);
}

// Balance of a quantity of the species, from the integrated diagnostics of f,
// df/dt, the source and the boundary fluxes.
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

  // Balanced quantity, from the integrated moments of f.
  int fm = bal->int_mom[GK_SPECIES_BALANCE_F][d];
  int f_comp = bal->int_comp[GK_SPECIES_BALANCE_F][d];
  double f_moms[gks->integ_moms[fm].num_mom];
  gk_species_balance_getlast(gks->integ_diag[fm], tm, f_moms);
  double mom = f_moms[f_comp];
  bal->mom[d] = mom;

  double dt = gks->fdot_dt;
  if (dt > 0.0) {
    // Rates are only defined once a step has been taken.
    double fdot_moms[gks->integ_moms[fm].num_mom];
    gk_species_balance_getlast(gks->fdot_integ_diag[fm], tm, fdot_moms);
    double mom_dot = fdot_moms[f_comp];

    // Source, including the BGK source terms.
    double src = 0.0;
    if (bal->has_source_bgk) {
      int bm = bal->int_mom[GK_SPECIES_BALANCE_SRC_BGK][d];
      double src_bgk_moms[gks->bgk_src.integ_moms[bm].num_mom];
      gk_species_balance_getlast(gks->bgk_src.integ_diag[bm], tm, src_bgk_moms);
      src += src_bgk_moms[bal->int_comp[GK_SPECIES_BALANCE_SRC_BGK][d]];
    }
    if (bal->has_source) {
      int sm = bal->int_mom[GK_SPECIES_BALANCE_SRC][d];
      double src_moms[gks->src.integ_moms[sm].num_mom];
      gk_species_balance_getlast(gks->src.integ_diag[sm], tm, src_moms);
      src += src_moms[bal->int_comp[GK_SPECIES_BALANCE_SRC][d]];
    }

    // Total flux through the boundaries.
    struct gk_boundary_fluxes *bflux = &gks->bflux;
    double bflux_tot_local = 0.0, bflux_tot = 0.0;
    for (int b = 0; b < bflux->num_boundaries; ++b) {
      if (gk_species_balance_owns_boundary(
            app, bflux->boundaries_dir[b], bflux->boundaries_edge[b]
          )) {
        int bm = bal->int_mom[GK_SPECIES_BALANCE_BFLUX][d];
        double bflux_moms[bflux->moms_op[bflux->diag_int_mom_idx[bm]].num_mom];
        gk_species_balance_getlast(bflux->intmom[b * bflux->num_diag_int_mom + bm], tm, bflux_moms);
        bflux_tot_local += bflux_moms[bal->int_comp[GK_SPECIES_BALANCE_BFLUX][d]];
      }
    }
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &bflux_tot_local, &bflux_tot);

    // Quantity gained/lost in the last step, also relative to mom.
    double mom_err = src - bflux_tot - mom_dot;
    bal_vals[0] = mom_err;
    bal_vals[1] = mom_err * dt / mom;
    bal_vals[2] = mom_dot;
    bal_vals[3] = src;
    bal_vals[4] = bflux_tot;
  }

  gkyl_dynvec_append(bal->diag[d], tm, bal_vals);
}

void
gk_species_balance_add_int_moms(
  struct gk_species_balance *bal, enum gk_species_balance_term term, int *num_moms,
  enum gkyl_distribution_moments *moms
)
{
  if (!bal->enabled) {
    return;
  }

  // Append the moment with every balanced quantity if one is missing.
  bool has_all = true;
  for (int d = 0; d < bal->num_diags; d++) {
    has_all = has_all && gk_species_balance_find_mom(bal->diag_types[d], *num_moms, moms) >= 0;
  }
  if (!has_all) {
    assert(*num_moms < 12);
    moms[(*num_moms)++] = bal->mom_type;
  }

  // Record where each balanced quantity is.
  for (int d = 0; d < bal->num_diags; d++) {
    int m = gk_species_balance_find_mom(bal->diag_types[d], *num_moms, moms);
    bal->int_mom[term][d] = m;
    bal->int_comp[term][d] = gk_species_balance_mom_comp(bal->diag_types[d], moms[m]);
  }
}

void
gk_species_balance_pre_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *bal,
  enum gkyl_species_bflux_type *bflux_type, struct gkyl_phase_diagnostics_inp *add_bflux_moms_inp
)
{
  // Collect the balances requested by the species (written) and by the app
  // (only computed, for the sum over species).
  *bal = (struct gk_species_balance){};
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

  // Moment with every balanced quantity, added to the integrated diagnostics missing one.
  bal->mom_type = GKYL_F_MOMENT_M0;
  for (int d = 0; d < bal->num_diags; d++) {
    if (bal->diag_types[d] == GKYL_BALANCE_GK_ENERGY) {
      bal->mom_type = GKYL_F_MOMENT_HAMILTONIAN;
    }
  }

  bal->has_source = info->source.source_id != 0;
  bal->has_source_bgk = info->source_bgk.source_bgk_id != 0;

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
    // Turn on the integrated bflux diagnostics, with the balanced quantities.
    // The bflux appends the additional moments to those requested.
    const struct gkyl_phase_diagnostics_inp *bflux_diags = &info->boundary_flux_diagnostics;
    int num_moms = 0;
    enum gkyl_distribution_moments moms[12];
    for (int m = 0; m < bflux_diags->num_integrated_diag_moments; m++) {
      moms[num_moms++] = bflux_diags->integrated_diag_moments[m];
    }
    for (int m = 0; m < add_bflux_moms_inp->num_integrated_diag_moments; m++) {
      moms[num_moms++] = add_bflux_moms_inp->integrated_diag_moments[m];
    }
    int num_moms_prev = num_moms;
    gk_species_balance_add_int_moms(bal, GK_SPECIES_BALANCE_BFLUX, &num_moms, moms);
    for (int m = num_moms_prev; m < num_moms; m++) {
      add_bflux_moms_inp
        ->integrated_diag_moments[add_bflux_moms_inp->num_integrated_diag_moments++] = moms[m];
    }
    *bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS;
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
  for (int d = 0; d < bal->num_diags; d++) {
    bal->diag[d] = gkyl_dynvec_new(GKYL_DOUBLE, bal->num_comp);
    bal->is_first_write[d] = true;
  }
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
    gkyl_dynvec_release(bal->diag[d]);
  }
}
