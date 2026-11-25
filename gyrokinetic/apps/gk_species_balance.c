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
    case GKYL_GK_BALANCE_PARTICLE:
      // The density is the first component of these moments.
      return mom_type == GKYL_F_MOMENT_M0 || mom_type == GKYL_F_MOMENT_M0M1M2 ||
                 mom_type == GKYL_F_MOMENT_M0M1M2PARM2PERP ||
                 mom_type == GKYL_F_MOMENT_HAMILTONIAN ?
               0 :
               -1;
    case GKYL_GK_BALANCE_ENERGY:
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

// Component comp of the last value appended to an integrated diagnostic, which
// must have been computed at time tm.
static double
gk_species_balance_getlast(gkyl_dynvec diag, double tm, int comp)
{
  double vals[gkyl_dynvec_ncomp(diag)];
  bool has_last = gkyl_dynvec_getlast(diag, vals);
  assert(has_last && gkyl_dynvec_getlast_tm(diag) == tm);
  return vals[comp];
}

// Component comp of the last value of the integrated moment mom_type.
static double
gk_species_balance_getlast_integ_mom(
  const struct gk_integ_moms *im, enum gkyl_distribution_moments mom_type, int comp, double tm
)
{
  int m = gk_integ_moms_find(im, mom_type);
  assert(m >= 0);
  return gk_species_balance_getlast(im->diag[m], tm, comp);
}

// Total flux of the quantity balanced by diagnostic d through the boundaries.
static double
gk_species_balance_bflux_tot(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *gksbal, int d,
  double tm
)
{
  struct gk_boundary_fluxes *bflux = &gks->bflux;
  int bm = -1;
  for (int m = 0; m < bflux->num_calc_integ_moms; m++) {
    enum gkyl_distribution_moments type = bflux->moms_op[bflux->calc_integ_mom_idx[m]].type;
    if (type == gksbal->integ_mom[GK_SPECIES_BALANCE_TERM_BFLUX][d]) {
      bm = m;
      break;
    }
  }
  assert(bm >= 0);

  double bflux_tot_local = 0.0, bflux_tot = 0.0;
  for (int b = 0; b < bflux->num_boundaries; ++b) {
    if (gk_species_balance_owns_boundary(app, bflux->boundaries_dir[b], bflux->boundaries_edge[b])) {
      bflux_tot_local += gk_species_balance_getlast(
        bflux->intmom[b * bflux->num_calc_integ_moms + bm], tm,
        gksbal->integ_comp[GK_SPECIES_BALANCE_TERM_BFLUX][d]
      );
    }
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &bflux_tot_local, &bflux_tot);
  return bflux_tot;
}

// Balance of a quantity of the species, from the integrated diagnostics of f,
// df/dt, the sources and the boundary fluxes.
static void
gk_species_balance_calc_diag(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *gksbal, int d,
  double tm, double dt
)
{
  double bal_vals[gksbal->num_comp];
  for (int k = 0; k < gksbal->num_comp; k++) {
    bal_vals[k] = 0.0;
  }

  // Balanced quantity, from the integrated moments of f.
  enum gkyl_distribution_moments f_mom = gksbal->integ_mom[GK_SPECIES_BALANCE_TERM_F][d];
  int f_comp = gksbal->integ_comp[GK_SPECIES_BALANCE_TERM_F][d];
  double mom = gk_species_balance_getlast_integ_mom(&gks->integ_moms, f_mom, f_comp, tm);
  gksbal->mom[d] = mom;

  if (dt > 0.0) {
    // Rates are only defined once a step has been taken.
    double mom_dot = gk_species_balance_getlast_integ_mom(&gks->fdot_integ_moms, f_mom, f_comp, tm);

    // Sources.
    double src = gksbal->has_term[GK_SPECIES_BALANCE_TERM_SRC] ?
                   gk_species_balance_getlast_integ_mom(
                     &gks->src.integ_moms, gksbal->integ_mom[GK_SPECIES_BALANCE_TERM_SRC][d],
                     gksbal->integ_comp[GK_SPECIES_BALANCE_TERM_SRC][d], tm
                   ) :
                   0.0;
    double src_bgk = gksbal->has_term[GK_SPECIES_BALANCE_TERM_SRC_BGK] ?
                       gk_species_balance_getlast_integ_mom(
                         &gks->bgk_src.integ_moms,
                         gksbal->integ_mom[GK_SPECIES_BALANCE_TERM_SRC_BGK][d],
                         gksbal->integ_comp[GK_SPECIES_BALANCE_TERM_SRC_BGK][d], tm
                       ) :
                       0.0;

    // Total flux through the boundaries.
    double bflux_tot = gksbal->has_term[GK_SPECIES_BALANCE_TERM_BFLUX] ?
                         gk_species_balance_bflux_tot(app, gks, gksbal, d, tm) :
                         0.0;

    // Quantity gained/lost in the last step, also relative to mom.
    double mom_err = src + src_bgk - bflux_tot - mom_dot;
    bal_vals[GK_SPECIES_BALANCE_COMP_ERR] = mom_err;
    bal_vals[GK_SPECIES_BALANCE_COMP_ERR_NORM] = mom_err * dt / mom;
    bal_vals[GK_SPECIES_BALANCE_COMP_MOM_DOT] = mom_dot;
    bal_vals[GK_SPECIES_BALANCE_COMP_SRC] = src;
    bal_vals[GK_SPECIES_BALANCE_COMP_SRC_BGK] = src_bgk;
    bal_vals[GK_SPECIES_BALANCE_COMP_BFLUX] = bflux_tot;
  }

  gkyl_dynvec_append(gksbal->diag[d], tm, bal_vals);
}

void
gk_species_balance_add_integ_moms(
  struct gk_species_balance *gksbal, enum gk_species_balance_term term, int *num_moms,
  enum gkyl_distribution_moments *moms
)
{
  if (gksbal->num_diags == 0) {
    return;
  }

  // Append the moment with every balanced quantity if one is missing.
  bool has_all = true;
  for (int d = 0; d < gksbal->num_diags; d++) {
    has_all = has_all && gk_species_balance_find_mom(gksbal->diag_types[d], *num_moms, moms) >= 0;
  }
  if (!has_all) {
    assert(*num_moms < GKYL_F_MOMENT_COUNT);
    moms[(*num_moms)++] = gksbal->mom_type;
  }

  // Record which moment, and which component of it, has each balanced quantity.
  for (int d = 0; d < gksbal->num_diags; d++) {
    int m = gk_species_balance_find_mom(gksbal->diag_types[d], *num_moms, moms);
    gksbal->integ_mom[term][d] = moms[m];
    gksbal->integ_comp[term][d] = gk_species_balance_mom_comp(gksbal->diag_types[d], moms[m]);
  }
}

void
gk_species_balance_add_fdot_moms(
  const struct gk_species_balance *gksbal, int *num_moms, enum gkyl_distribution_moments *moms
)
{
  for (int d = 0; d < gksbal->num_diags; d++) {
    enum gkyl_distribution_moments f_mom = gksbal->integ_mom[GK_SPECIES_BALANCE_TERM_F][d];
    bool found = false;
    for (int m = 0; m < *num_moms; m++) {
      found = found || moms[m] == f_mom;
    }
    if (!found) {
      assert(*num_moms < GKYL_F_MOMENT_COUNT);
      moms[(*num_moms)++] = f_mom;
    }
  }
}

void
gk_species_balance_pre_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *gksbal,
  enum gkyl_species_bflux_type *bflux_type, struct gkyl_phase_diagnostics_inp *add_bflux_moms_inp
)
{
  // Collect the balances requested by the species and by the app (for the sum over species).
  *gksbal = (struct gk_species_balance){};
  struct gkyl_gyrokinetic_species *info = &gks->info;
  if (info->is_static) {
    return;
  }
  for (int d = 0; d < info->num_balance_diags; d++) {
    gksbal->diag_types[gksbal->num_diags++] = info->balance_diags[d];
  }
  for (int d = 0; d < app->gkbal.num_diags; d++) {
    bool found = false;
    for (int k = 0; k < gksbal->num_diags; k++) {
      found = found || gksbal->diag_types[k] == app->gkbal.diag_types[d];
    }
    if (!found) {
      assert(gksbal->num_diags < GKYL_GK_BALANCE_COUNT);
      gksbal->diag_types[gksbal->num_diags++] = app->gkbal.diag_types[d];
    }
  }
  if (gksbal->num_diags == 0) {
    return;
  }

  // Moment with every balanced quantity, added to the integrated diagnostics missing one.
  gksbal->mom_type = GKYL_F_MOMENT_M0;
  for (int d = 0; d < gksbal->num_diags; d++) {
    if (gksbal->diag_types[d] == GKYL_GK_BALANCE_ENERGY) {
      gksbal->mom_type = GKYL_F_MOMENT_HAMILTONIAN;
    }
  }

  // Terms of the balance.
  gksbal->has_term[GK_SPECIES_BALANCE_TERM_F] = true;
  gksbal->has_term[GK_SPECIES_BALANCE_TERM_SRC] = info->source.source_id != 0;
  gksbal->has_term[GK_SPECIES_BALANCE_TERM_SRC_BGK] = info->source_bgk.source_bgk_id != 0;
  // Boundaries with a boundary flux (same as in gk_species_bflux_init).
  for (int d = 0; d < app->cdim; ++d) {
    if (gks->bc_is_np[d] && (gks->lower_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX ||
                             gks->upper_bc[d].type != GKYL_BC_GK_SPECIES_ZERO_FLUX)) {
      gksbal->has_term[GK_SPECIES_BALANCE_TERM_BFLUX] = true;
    }
  }

  if (gksbal->has_term[GK_SPECIES_BALANCE_TERM_BFLUX]) {
    // Turn on the integrated bflux diagnostics, with the balanced quantities.
    // The bflux appends the additional moments to those requested.
    const struct gkyl_phase_diagnostics_inp *bflux_diags = &info->boundary_flux_diagnostics;
    int num_moms = 0;
    enum gkyl_distribution_moments moms[GKYL_F_MOMENT_COUNT];
    for (int m = 0; m < bflux_diags->num_integrated_diag_moments; m++) {
      moms[num_moms++] = bflux_diags->integrated_diag_moments[m];
    }
    for (int m = 0; m < add_bflux_moms_inp->num_integrated_diag_moments; m++) {
      moms[num_moms++] = add_bflux_moms_inp->integrated_diag_moments[m];
    }
    int num_moms_prev = num_moms;
    gk_species_balance_add_integ_moms(gksbal, GK_SPECIES_BALANCE_TERM_BFLUX, &num_moms, moms);
    for (int m = num_moms_prev; m < num_moms; m++) {
      add_bflux_moms_inp
        ->integrated_diag_moments[add_bflux_moms_inp->num_integrated_diag_moments++] = moms[m];
    }
    *bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS;
  }
}

void
gk_species_balance_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *gksbal
)
{
  if (gksbal->num_diags == 0) {
    return;
  }

  gksbal->num_comp = GK_SPECIES_BALANCE_COMP_COUNT;
  for (int d = 0; d < gksbal->num_diags; d++) {
    gksbal->diag[d] = gkyl_dynvec_new(GKYL_DOUBLE, gksbal->num_comp);
    gksbal->is_first_write[d] = true;
  }
}

void
gk_species_balance_calc(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *gksbal, double tm,
  double dt
)
{
  for (int d = 0; d < gksbal->num_diags; d++) {
    gk_species_balance_calc_diag(app, gks, gksbal, d, tm, dt);
  }
}

void
gk_species_balance_write(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gk_species_balance *gksbal
)
{
  struct timespec wst = gkyl_wall_clock();
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  for (int d = 0; d < gksbal->num_diags; d++) {
    if (rank == 0) {
      // File named after the balance, e.g. <name>-<species>_balance_particle.gkyl.
      const char *fmt = "%s-%s_balance_%s.gkyl";
      const char *diag_name = gkyl_gyrokinetic_balance_diag_strs[gksbal->diag_types[d]];
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, diag_name);
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, diag_name);

      if (gksbal->is_first_write[d]) {
        struct gkyl_msgpack_map_elem io_meta_bal[] = {
          {
            .key = "Description",
            .elem_type = GKYL_MP_STRING,
            .cval = "Species balance: mom_err = src_mom + src_bgk_mom - bflux_tot_mom - mom_dot, "
                    "mom_err_norm = mom_err*dt/mom.",
          },
          {
            .key = "Columns",
            .elem_type = GKYL_MP_STRING,
            .cval = "mom_err, mom_err_norm, mom_dot, src_mom, src_bgk_mom, bflux_tot_mom",
          },
        };
        int io_meta_len[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 2};
        const struct gkyl_msgpack_map_elem *io_meta[] = {
          gks->io_meta_basic, app->gk_geom->io_meta_basic, io_meta_bal
        };
        struct gkyl_msgpack_data *mt =
          gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(gksbal->diag[d], fileNm, mt);
        gkyl_msgpack_data_release(mt);
      } else {
        gkyl_dynvec_awrite(gksbal->diag[d], fileNm);
      }
      app->stat.n_diag_io += 1;
    }
    gksbal->is_first_write[d] = false;
    gkyl_dynvec_clear(gksbal->diag[d]);
  }
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_balance_release(
  const struct gkyl_gyrokinetic_app *app, const struct gk_species_balance *gksbal
)
{
  for (int d = 0; d < gksbal->num_diags; d++) {
    gkyl_dynvec_release(gksbal->diag[d]);
  }
}

void
gk_app_balance_init(
  struct gkyl_gyrokinetic_app *app, struct gk_app_balance *gkbal,
  const struct gkyl_phase_diagnostics_inp *diags
)
{
  *gkbal = (struct gk_app_balance){};
  gkbal->num_diags = diags->num_balance_diags;
  for (int d = 0; d < gkbal->num_diags; d++) {
    gkbal->diag_types[d] = diags->balance_diags[d];
    // Species components, and field_dot, bg_dot and apar_dot for the energy.
    gkbal->num_comp[d] = gkbal->diag_types[d] == GKYL_GK_BALANCE_ENERGY ?
                           GK_APP_BALANCE_COMP_COUNT :
                           GK_SPECIES_BALANCE_COMP_COUNT;
    gkbal->diag[d] = gkyl_dynvec_new(GKYL_DOUBLE, gkbal->num_comp[d]);
    gkbal->is_first_write[d] = true;
  }
}

void
gk_app_balance_calc(
  struct gkyl_gyrokinetic_app *app, struct gk_app_balance *gkbal, double tm, double dt
)
{
  // Sum the balance of each dynamic species; the energy balance also includes
  // the field terms. Must be called after the species balances are computed.
  for (int d = 0; d < gkbal->num_diags; d++) {
    double bal_vals[gkbal->num_comp[d]];
    for (int k = 0; k < gkbal->num_comp[d]; k++) {
      bal_vals[k] = 0.0;
    }
    double mom_tot = 0.0;
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species_balance *gksbal = &app->species[i].gksbal;
      if (gksbal->num_diags == 0) {
        continue;
      }
      int sd = -1;
      for (int k = 0; k < gksbal->num_diags; k++) {
        sd = gksbal->diag_types[k] == gkbal->diag_types[d] ? k : sd;
      }
      assert(sd >= 0); // The species computes every balance of the app.
      double sbal_vals[gksbal->num_comp];
      gkyl_dynvec_getlast(gksbal->diag[sd], sbal_vals);
      // Sum every component but the normalized error, recomputed below.
      for (int k = 0; k < gksbal->num_comp; k++) {
        if (k != GK_SPECIES_BALANCE_COMP_ERR_NORM) {
          bal_vals[k] += sbal_vals[k];
        }
      }
      mom_tot += gksbal->mom[sd];
    }
    if (dt > 0.0) {
      if (gkbal->diag_types[d] == GKYL_GK_BALANCE_ENERGY && app->field->update_field) {
        // Field terms of the energy balance.
        double energy, field_dot, apar_energy, apar_dot, bg_energy, bg_dot;
        gk_field_calc_energy_balance(
          app, app->field, tm, &energy, &field_dot, &apar_energy, &apar_dot, &bg_energy, &bg_dot
        );
        // mom_err = src + src_bgk - bflux - (mom_dot - field_dot - apar_dot + bg_dot).
        bal_vals[GK_SPECIES_BALANCE_COMP_ERR] += field_dot + apar_dot - bg_dot;
        bal_vals[GK_APP_BALANCE_COMP_FIELD_DOT] = field_dot;
        bal_vals[GK_APP_BALANCE_COMP_BG_DOT] = bg_dot;
        bal_vals[GK_APP_BALANCE_COMP_APAR_DOT] = apar_dot;
        mom_tot += bg_energy - energy - apar_energy;
      }
      bal_vals[GK_SPECIES_BALANCE_COMP_ERR_NORM] =
        bal_vals[GK_SPECIES_BALANCE_COMP_ERR] * dt / mom_tot;
    }
    gkyl_dynvec_append(gkbal->diag[d], tm, bal_vals);
  }
}

void
gk_app_balance_write(struct gkyl_gyrokinetic_app *app, struct gk_app_balance *gkbal)
{
  struct timespec wst = gkyl_wall_clock();
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  for (int d = 0; d < gkbal->num_diags; d++) {
    if (rank == 0) {
      // File named after the balance, e.g. <name>-balance_particle.gkyl.
      const char *fmt = "%s-balance_%s.gkyl";
      const char *diag_name = gkyl_gyrokinetic_balance_diag_strs[gkbal->diag_types[d]];
      int sz = gkyl_calc_strlen(fmt, app->name, diag_name);
      char fileNm[sz + 1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, diag_name);

      if (gkbal->is_first_write[d]) {
        bool is_energy = gkbal->diag_types[d] == GKYL_GK_BALANCE_ENERGY;
        struct gkyl_msgpack_map_elem io_meta_bal[] = {
          {
            .key = "Description",
            .elem_type = GKYL_MP_STRING,
            .cval = is_energy ?
                      "Energy balance summed over species: mom_err = src_mom + src_bgk_mom - "
                      "bflux_tot_mom - (mom_dot - field_dot - apar_dot + bg_dot), with field_dot "
                      "= dE_field/dt, bg_dot = d/dt int rho_bg phi and apar_dot = "
                      "-d/dt int |grad_perp A_par|^2/(2 mu0)." :
                      "Balance summed over species: mom_err = src_mom + src_bgk_mom - "
                      "bflux_tot_mom - mom_dot.",
          },
          {
            .key = "Columns",
            .elem_type = GKYL_MP_STRING,
            .cval = is_energy ?
                      "mom_err, mom_err_norm, mom_dot, src_mom, src_bgk_mom, bflux_tot_mom, "
                      "field_dot, bg_dot, apar_dot" :
                      "mom_err, mom_err_norm, mom_dot, src_mom, src_bgk_mom, bflux_tot_mom",
          },
        };
        int io_meta_len[] = {app->gk_geom->io_meta_basic_len, 2};
        const struct gkyl_msgpack_map_elem *io_meta[] = {app->gk_geom->io_meta_basic, io_meta_bal};
        struct gkyl_msgpack_data *mt =
          gkyl_msgpack_create_union(sizeof(io_meta_len) / sizeof(int), io_meta_len, io_meta);

        gkyl_dynvec_write_wmeta(gkbal->diag[d], fileNm, mt);
        gkyl_msgpack_data_release(mt);
      } else {
        gkyl_dynvec_awrite(gkbal->diag[d], fileNm);
      }
      app->stat.n_diag_io += 1;
    }
    gkbal->is_first_write[d] = false;
    gkyl_dynvec_clear(gkbal->diag[d]);
  }
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_app_balance_release(const struct gk_app_balance *gkbal)
{
  for (int d = 0; d < gkbal->num_diags; d++) {
    gkyl_dynvec_release(gkbal->diag[d]);
  }
}
