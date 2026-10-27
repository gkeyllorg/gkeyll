#include <gkyl_gyrokinetic_priv.h>

// Each entry supplies a local rate and its definition. Storage and component
// metadata follow this registry; transit_min uses a minimum MPI reduction.
struct gk_frequency_descriptor {
  const char *name, *definition;
  double (*calc)(gkyl_gyrokinetic_app *, struct gk_species *);
};

static double
gk_frequency_omegaH(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  return gk_species_omegaH_frequency(app, gks, gks->f);
}

// Thermal rates use cell-average primitive moments and stationary geometry.
static double
gk_frequency_transit(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  return gks->frequency_rates[0];
}

static double
gk_frequency_transit_max(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  return gks->transit_frequency * gks->collisionless.scale_fac *
         gk_fdot_multiplier_get_time_dilation_scale_const(app, &gks->fdot_mult);
}

static double
gk_frequency_transit_min(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  return gks->frequency_rates[1];
}

static double
gk_frequency_drift(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  return gks->frequency_rates[2];
}

static double
gk_frequency_collision(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  return gks->frequency_rates[3];
}

static const struct gk_frequency_descriptor frequency_registry[] = {
  {"omega_H",
   "Local omega_H stability estimate, including equation scaling; zero for Boltzmann/adiabatic fields. No CFL factor.",
   gk_frequency_omegaH},
  {"transit",
   "max_x sqrt(T_parallel/m)*abs(<C/(JB)>)/L_z; cell-average primitive temperature, computational block length L_z, equation scaling, no 2*pi or CFL factor.",
   gk_frequency_transit},
  {"transit_max",
   "Velocity-domain upper transit estimate: max|v_parallel|*max_x abs(<C/(JB)>)/L_z; equation scaling, no 2*pi or CFL factor.",
   gk_frequency_transit_max},
  {"transit_min",
   "Minimum velocity-domain transit estimate: min|v_parallel|*min_x abs(<C/(JB)>)/L_z; zero if the continuous velocity interval contains zero. Equation scaling, no 2*pi or CFL factor.",
   gk_frequency_transit_min},
  {"omega_d",
   "Thermal magnetic drift estimate max_x sum_i (abs(v_curvature^i)+abs(v_gradB^i))/L_i over simulated perpendicular directions; v_d=(T_parallel/m)*curl(b)/Omega+(T_perp/qB)*b cross grad(log B). Cell-center geometry and cell-average temperatures, no mean-flow contribution, equation scaling, no 2*pi or CFL factor. Zero in 1x.",
   gk_frequency_drift},
  {"nu",
   "max_x cell-average sum of self and cross LBO/BGK collision frequencies, refreshed at diagnostic time; includes configured nu_frac, excludes collisionless and time-dilation scaling. Zero without collisions.",
   gk_frequency_collision}
};

enum { num_frequencies = sizeof frequency_registry / sizeof frequency_registry[0] };

void
gk_species_frequencies_init(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  gks->frequencies_max = gkyl_dynvec_new(GKYL_DOUBLE, num_frequencies);
  gks->is_first_frequencies_write_call = true;
  gk_species_moment_init(
    app, gks, &gks->frequency_moms,
    gks->grid_vel.ndim == 2 ? GKYL_F_MOMENT_BIMAXWELLIAN : GKYL_F_MOMENT_MAXWELLIAN, false
  );
  gks->frequency_geometry = mkarr(false, 3, app->local_ext.volume);
  gks->frequency_nu = mkarr(false, app->basis.num_basis, app->local_ext.volume);

  // Geometry and velocity mapping are stationary. Cache the drift/transit factors
  // on the host for both CPU and GPU runs.
  int nb = app->basis.num_basis;
  struct gkyl_array *parfac = mkarr(app->use_gpu, nb, app->local_ext.volume);
  gkyl_dg_mul_op_range(
    &app->basis, 0, parfac, 0, app->gk_geom->geo_int.cmag, 0, app->gk_geom->geo_int.jacobtot_inv,
    &app->local
  );

  // Reuse CPU geometry; only GPU runs need temporary host copies.
  struct gkyl_array *parfac_host = parfac;
  struct gkyl_array *jacobtot_inv_host = app->gk_geom->geo_int.jacobtot_inv;
  struct gkyl_array *bmag_host = app->gk_geom->geo_int.bmag;
  struct gkyl_array *b_i_host = app->gk_geom->geo_int.b_i;
  struct gkyl_array *dualcurlbhatoverB_host = app->gk_geom->geo_int.dualcurlbhatoverB;
  if (app->use_gpu) {
    parfac_host = mkarr(false, nb, app->local_ext.volume);
    jacobtot_inv_host = mkarr(false, nb, app->local_ext.volume);
    bmag_host = mkarr(false, nb, app->local_ext.volume);
    b_i_host = mkarr(false, 3 * nb, app->local_ext.volume);
    dualcurlbhatoverB_host = mkarr(false, 3 * nb, app->local_ext.volume);
    gkyl_array_copy(parfac_host, parfac);
    gkyl_array_copy(jacobtot_inv_host, app->gk_geom->geo_int.jacobtot_inv);
    gkyl_array_copy(bmag_host, app->gk_geom->geo_int.bmag);
    gkyl_array_copy(b_i_host, app->gk_geom->geo_int.b_i);
    gkyl_array_copy(dualcurlbhatoverB_host, app->gk_geom->geo_int.dualcurlbhatoverB);
  }
  double max_parfac = 0.0, min_parfac = DBL_MAX;
  double norm = pow(sqrt(2.0), app->cdim);
  double center[3] = {0.0};

  // Iterate over the local range to compute and cache frequency-related geometry.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);
    const double *coeff = gkyl_array_cfetch(parfac_host, loc);
    double *geom = gkyl_array_fetch(gks->frequency_geometry, loc);
    geom[0] = fabs(coeff[0]) / norm;
    max_parfac = fmax(max_parfac, geom[0]);
    min_parfac = fmin(min_parfac, geom[0]);
    const double *bc = gkyl_array_cfetch(bmag_host, loc);
    const double *bic = gkyl_array_cfetch(b_i_host, loc);
    const double *cc = gkyl_array_cfetch(dualcurlbhatoverB_host, loc);
    const double *ij = gkyl_array_cfetch(jacobtot_inv_host, loc);
    double b = app->basis.eval_expand(center, bc), grad[3] = {0.0};
    for (int d = 0; d < app->cdim; ++d) {
      int gd = d == app->cdim - 1 ? 2 : d;
      grad[gd] = app->basis.eval_grad_expand(d, center, bc) * 2.0 / app->grid.dx[d];
    }
    // Cache the perpendicular contravariant drift coefficients divided by block lengths.
    // The two coefficients are kept separately for anisotropic temperatures.
    geom[1] = geom[2] = 0.0;
    for (int d = 0; d < app->cdim - 1; ++d) {
      int j = (d + 1) % 3, k = (d + 2) % 3;
      double cross = app->basis.eval_expand(center, bic + j * nb) * grad[k] -
                     app->basis.eval_expand(center, bic + k * nb) * grad[j];
      double length = app->grid.upper[d] - app->grid.lower[d];
      // Sum magnitudes to provide a conservative characteristic estimate.
      geom[1] += fabs(app->basis.eval_expand(center, cc + d * nb)) / length;
      geom[2] += fabs(cross * app->basis.eval_expand(center, ij) / b) / length;
    }
  }
  int vdim = gks->grid_vel.ndim, zdir = app->cdim - 1;
  double vmax = fmax(fabs(gks->vel_map->vbounds[0]), fabs(gks->vel_map->vbounds[vdim]));
  gks->transit_frequency = vmax * max_parfac / (app->grid.upper[zdir] - app->grid.lower[zdir]);
  double vmin = gks->vel_map->vbounds[0] <= 0.0 && gks->vel_map->vbounds[vdim] >= 0.0 ?
                  0.0 :
                  fmin(fabs(gks->vel_map->vbounds[0]), fabs(gks->vel_map->vbounds[vdim]));
  gks->transit_frequency_min = vmin * min_parfac / (app->grid.upper[zdir] - app->grid.lower[zdir]);
  if (app->use_gpu) {
    gkyl_array_release(parfac_host);
    gkyl_array_release(jacobtot_inv_host);
    gkyl_array_release(bmag_host);
    gkyl_array_release(b_i_host);
    gkyl_array_release(dualcurlbhatoverB_host);
  }
  gkyl_array_release(parfac);

  // Flat, zero-based component keys are readable by existing msgpack readers.
  struct gkyl_msgpack_map_elem meta[3 + 2 * num_frequencies];
  char name_keys[num_frequencies][32], definition_keys[num_frequencies][32];
  meta[0] = gmpe_cval(
    "Description",
    "Characteristic frequencies reduced across MPI ranks within this block; maxima except transit_min (minimum)."
  );
  meta[1] = gmpe_cval("units", "1/s");
  meta[2] = gmpe_ival("numComponents", num_frequencies);
  for (int i = 0; i < num_frequencies; ++i) {
    snprintf(name_keys[i], sizeof name_keys[i], "component%d", i);
    snprintf(definition_keys[i], sizeof definition_keys[i], "component%dDefinition", i);
    meta[3 + 2 * i] = (struct gkyl_msgpack_map_elem){
      .key = name_keys[i],
      .elem_type = GKYL_MP_STRING,
      .cval = frequency_registry[i].name,
    };
    meta[4 + 2 * i] = (struct gkyl_msgpack_map_elem){
      .key = definition_keys[i],
      .elem_type = GKYL_MP_STRING,
      .cval = frequency_registry[i].definition,
    };
  }
  int lengths[] = {gks->io_meta_basic_len, app->gk_geom->io_meta_basic_len, 3 + 2 * num_frequencies};
  const struct gkyl_msgpack_map_elem *lists[] = {
    gks->io_meta_basic, app->gk_geom->io_meta_basic, meta
  };
  gks->frequencies_meta = gkyl_msgpack_create_union(3, lengths, lists);
}

void
gk_species_frequencies_calc(gkyl_gyrokinetic_app *app, struct gk_species *gks, double tm)
{
  gk_species_moment_calc(app, &gks->frequency_moms, gks->local, app->local, gks->f);
  gkyl_array_copy(gks->frequency_moms.marr_host, gks->frequency_moms.marr);
  // Cross-species rates require all partners' moments from the same state.
  if (gks->lbo.collision_id == GKYL_LBO_COLLISIONS ||
      gks->bgk.collision_id == GKYL_BGK_COLLISIONS) {
    for (int i = 0; i < app->num_species; ++i) {
      struct gk_species *sp = &app->species[i];
      gk_species_lbo_moms(app, sp, &sp->lbo, sp->f);
      gk_species_bgk_moms(app, sp, &sp->bgk, sp->f);
    }
    gk_species_lbo_cross_moms(app, gks, &gks->lbo, gks->f);
    gk_species_bgk_cross_moms(app, gks, &gks->bgk, gks->f);
  }
  gkyl_array_clear(gks->frequency_nu, 0.0);
  if (gks->lbo.collision_id == GKYL_LBO_COLLISIONS) {
    gkyl_array_copy(gks->frequency_nu, gks->lbo.nu_sum);
  }
  if (gks->bgk.collision_id == GKYL_BGK_COLLISIONS) {
    gkyl_array_copy(gks->frequency_nu, gks->bgk.nu_sum);
  }
  double scale = gks->collisionless.scale_fac *
                 gk_fdot_multiplier_get_time_dilation_scale_const(app, &gks->fdot_mult);
  double norm = pow(sqrt(2.0), app->cdim);
  int nb = app->basis.num_basis, zdir = app->cdim - 1;
  double length = app->grid.upper[zdir] - app->grid.lower[zdir];
  memset(gks->frequency_rates, 0, sizeof gks->frequency_rates);
  gks->frequency_rates[1] = scale * gks->transit_frequency_min;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);
    const double *mom = gkyl_array_cfetch(gks->frequency_moms.marr_host, loc);
    const double *geom = gkyl_array_cfetch(gks->frequency_geometry, loc);
    const double *nu = gkyl_array_cfetch(gks->frequency_nu, loc);
    double tpar = fmax(0.0, mom[2 * nb] / norm),
           tperp = gks->grid_vel.ndim == 2 ? fmax(0.0, mom[3 * nb] / norm) : 0.0;
    gks->frequency_rates[0] = fmax(gks->frequency_rates[0], scale * sqrt(tpar) * geom[0] / length);
    double drift = gks->info.charge == 0.0 ? 0.0 :
                                             scale * gks->info.mass / fabs(gks->info.charge) *
                                               (tpar * geom[1] + tperp * geom[2]);
    gks->frequency_rates[2] = fmax(gks->frequency_rates[2], drift);
    gks->frequency_rates[3] = fmax(gks->frequency_rates[3], nu[0] / norm);
  }
  double local[num_frequencies], global[num_frequencies];
  for (int i = 0; i < num_frequencies; ++i) {
    local[i] = frequency_registry[i].calc(app, gks);
  }
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, num_frequencies, local, global);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &local[3], &global[3]);
  gkyl_dynvec_append(gks->frequencies_max, tm, global);
}

void
gk_species_frequencies_write(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0) {
    const char *fmt = "%s-%s_frequencies_max.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name);
    char filename[sz + 1];
    snprintf(filename, sizeof filename, fmt, app->name, gks->info.name);
    // A restarted run may enable this diagnostic for the first time.
    if (gks->is_first_frequencies_write_call || !gkyl_check_file_exists(filename)) {
      gkyl_dynvec_write_wmeta(gks->frequencies_max, filename, gks->frequencies_meta);
      gks->is_first_frequencies_write_call = false;
    } else {
      gkyl_dynvec_awrite(gks->frequencies_max, filename);
    }
  }
  gkyl_dynvec_clear(gks->frequencies_max);
  app->stat.n_diag_io += 1;
}

void
gk_species_frequencies_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  gk_species_moment_release(app, &gks->frequency_moms);
  gkyl_array_release(gks->frequency_geometry);
  gkyl_array_release(gks->frequency_nu);
  gkyl_dynvec_release(gks->frequencies_max);
  gkyl_msgpack_data_release(gks->frequencies_meta);
}
