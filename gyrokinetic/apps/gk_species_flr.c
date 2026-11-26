#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_species_priv.h>
#include <gkyl_gyrokinetic_priv.h>

#include <assert.h>
#include <math.h>

// Finite Larmor radius (FLR) model of a species. The model is chosen in the field input, and
// provides the gyroaverage of the potential (which fills gyro_phi and the potential of the
// Hamiltonian, hamil_phi), the gyroaverage of the density (which fills m0_gyroavg), the
// gyroaverage of the particle flow of the electromagnetic field equations (m1_gyroavg), and the
// gyroaverage of A_parallel and its time derivative (gyro_apar, gyro_apardot):
// - none: nothing is gyroaveraged.
// - Pade of each mu cell: the potential, A_parallel and its time derivative are gyroaveraged in
//   each mu cell with the Pade approximant (1 - div_perp((rho_j^2/4) grad_perp))^{-1}, rho_j^2 =
//   2*m*mu_j/(q^2*B) being the squared gyroradius of the mu cell j, and the density and the flow
//   of each mu cell with the same (symmetric) operator. The field equations are those without
//   FLR effects (Ishizawa et al., Phys. Plasmas 26, 082301 (2019)).

static void
gk_species_flr_hamil_mu_independent(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // The Hamiltonian uses the gyroaveraged potential.
  gks->hamil_phi = gkyl_array_acquire(gks->gyro_phi);
  gks->hamil_num_mu = 1;
}

/** No FLR effects. **/

static void
gk_species_flr_none_gyroaverage(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
)
{
}

static void
gk_species_flr_none_gyroaverage_density(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_dens_in, struct gkyl_array *jac_dens_gyroavg
)
{
}

static void
gk_species_flr_none_gyroaverage_moment(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_mom_in, struct gkyl_array *jac_mom_gyroavg
)
{
}

static void
gk_species_flr_none_gyroaverage_em(gkyl_gyrokinetic_app *app, struct gk_species *species)
{
}

static void
gk_species_flr_none_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
}

static void
gk_species_flr_none_new(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // The gyroaveraged potential and density are the potential and the density.
  gks->gyro_phi = gkyl_array_acquire(app->field->phi_smooth);
  gks->m0_gyroavg = gkyl_array_acquire(gks->m0.marr);
  gk_species_flr_hamil_mu_independent(app, gks);

  gks->gyroaverage = gk_species_flr_none_gyroaverage;
  gks->gyroaverage_density = gk_species_flr_none_gyroaverage_density;
  gks->gyroaverage_moment = gk_species_flr_none_gyroaverage_moment;
  gks->flr.release_func = gk_species_flr_none_release;
}

/** Pade approximant of the gyroaverage with the gyroradius of each mu cell. **/

static void
gk_species_flr_pade_mu_smooth_coeff(
  gkyl_gyrokinetic_app *app, struct gkyl_array *coeff, struct gkyl_array *tmp
)
{
  // Smooth each component of an FLR operator coefficient along z, so the
  // z-nodal deflated solver samples it consistently with the smoothed density.
  int ncomp = coeff->ncomp / app->basis.num_basis;
  for (int c = 0; c < ncomp; c++) {
    gkyl_array_set_offset(tmp, 1.0, coeff, c * app->basis.num_basis);
    app->field->fem_projection_par_rho_func(app, app->field, tmp, tmp);
    gkyl_array_set_offset(coeff, 1.0, tmp, c * app->basis.num_basis);
  }
}

static void
gk_species_flr_pade_mu_grid(
  gkyl_gyrokinetic_app *app, const struct gk_species *gks, double *mu_cells, double *mu_avg_weight
)
{
  // Magnetic moment each mu cell is gyroaveraged at, and its weight in averages over mu: the
  // average magnetic moment in the cell of a Maxwellian in mu, exp(-mu*B/T), and the fraction
  // of that Maxwellian in the cell. T/B follows from the reference
  // gyroradius of the species, rho^2 = T*m/(q*B)^2.
  int num_mu = gkyl_range_shape(&gks->local_vel, 1);
  const struct gkyl_velocity_map *gvm = gks->vel_map;
  int vmap_num_basis = gvm->vmap_basis_ho.num_basis;
  double gyroradius_bmag = gks->info.flr.bmag ? gks->info.flr.bmag : app->bmag_ref;
  double temp_over_bmag =
    gks->flr_rhoSq_ref * pow(gks->info.charge, 2.0) * gyroradius_bmag / gks->info.mass;
  double weight_tot = 0.0;
  for (int j = 0; j < num_mu; j++) {
    int vidx[] = {gks->local_vel.lower[0], gks->local_vel.lower[1] + j};
    const double *vmap_d = gkyl_array_cfetch(gvm->vmap_ho, gkyl_range_idx(&gvm->local_vel, vidx));
    const double *mu_map_d = &vmap_d[vmap_num_basis];
    double cell_lo = -1.0, cell_up = 1.0;
    double mu_lo = fmax(0.0, gvm->vmap_basis_ho.eval_expand(&cell_lo, mu_map_d));
    double mu_up = gvm->vmap_basis_ho.eval_expand(&cell_up, mu_map_d);
    double cell_len = mu_up - mu_lo;
    // Fraction of the Maxwellian at the lower edge that is left at the upper edge.
    double decay = exp(-cell_len / temp_over_bmag);
    mu_cells[j] = mu_lo + temp_over_bmag - cell_len * decay / (1.0 - decay);
    mu_avg_weight[j] = exp(-mu_lo / temp_over_bmag) * (1.0 - decay);
    weight_tot += mu_avg_weight[j];
  }
  for (int j = 0; j < num_mu; j++) {
    mu_avg_weight[j] /= weight_tot;
  }
}

static void
gk_species_flr_pade_mu_gyroaverage(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
)
{
  // Gyroaverage the potential with the gyroradius of each mu cell, solving
  // J*phi_j - div_perp((rho_j^2/4) J g^ij grad_perp phi_j) = J*phi with phi as the boundary
  // value. The Hamiltonian uses this potential, constant in each mu cell. Also compute its
  // average over a Maxwellian in mu (used wherever a potential independent of mu is needed).
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_pade_mu *pm = &species->flr.pade_mu;
  int num_basis = app->basis.num_basis;
  gkyl_dg_mul_op_range(
    &app->basis, 0, pm->buff, 0, app->gk_geom->geo_int.jacobgeo, 0, field_in, &app->local
  );
  gkyl_array_clear_range(field_gyroavg, 0.0, &app->local);
  for (int j = 0; j < pm->num_mu; j++) {
    gkyl_deflated_fem_poisson_advance(pm->ops[j], pm->buff, field_in, pm->buff_gyroavg);
    gkyl_array_set_offset(species->hamil_phi, 1.0, pm->buff_gyroavg, j * num_basis);
    gkyl_array_accumulate_range(field_gyroavg, pm->mu_avg_weight[j], pm->buff_gyroavg, &app->local);
  }
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_pade_mu_gyroaverage_density(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_dens_in, struct gkyl_array *jac_dens_gyroavg
)
{
  // Gyroaverage the density (times J) in each mu cell with the operator that gyroaverages
  // the potential of that mu cell, and add them up. The operators are symmetric, so the
  // integral of f*q*<phi> over phase space equals the integral of phi times the gyroaveraged
  // charge density. The difference between the gyroaveraged and gyrocenter densities of the
  // initial condition is removed, see gk_species_flr_pade_mu_gyroaverage_density_init.
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_pade_mu *pm = &species->flr.pade_mu;
  gkyl_array_clear_range(jac_dens_gyroavg, 0.0, &app->local);
  for (int j = 0; j < pm->num_mu; j++) {
    // Density (times J) of the mu cell: the M0 moment over the phase space cells of that mu cell.
    gk_species_moment_calc(app, &pm->m0_mu, pm->mu_ranges[j], app->local, fin);
    // Smooth the density along z to get a consistent field for the deflated solve.
    app->field->fem_projection_par_rho_func(app, app->field, pm->m0_mu.marr, pm->buff_dens);
    // Boundary value n = (J*n)/J.
    gkyl_dg_mul_op_range(
      &app->basis, 0, pm->buff, 0, app->gk_geom->geo_int.jacobgeo_inv, 0, pm->buff_dens, &app->local
    );
    gkyl_deflated_fem_poisson_advance(pm->ops[j], pm->buff_dens, pm->buff, pm->buff_gyroavg);
    gkyl_dg_mul_op_range(
      &app->basis, 0, pm->buff_gyroavg, 0, app->gk_geom->geo_int.jacobgeo, 0, pm->buff_gyroavg,
      &app->local
    );
    gkyl_array_accumulate_range(jac_dens_gyroavg, 1.0, pm->buff_gyroavg, &app->local);
  }
  gkyl_array_accumulate_range(jac_dens_gyroavg, -1.0, pm->dens_gyroavg_init, &app->local);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_pade_mu_gyroaverage_density_init(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_dens_in, struct gkyl_array *jac_dens_gyroavg
)
{
  // First gyroaverage of the density, done with the initial condition. The gyroaveraged
  // density of an equilibrium with a density or temperature profile differs from its
  // gyrocenter density by O((rho/L)^2), and so do the charge densities of species with
  // different gyroradii. The long wavelength polarization turns that into an O(1) potential
  // (e*phi/T of the order of the relative profile variation). Hence store the difference
  // between the gyroaveraged and the gyrocenter densities of the initial condition, and
  // remove it from the gyroaveraged density from now on: the initial charge density then is
  // that of the gyrocenters, and FLR effects act on departures from the initial condition.
  struct gk_species_flr_pade_mu *pm = &species->flr.pade_mu;
  gkyl_array_clear(pm->dens_gyroavg_init, 0.0);
  gk_species_flr_pade_mu_gyroaverage_density(app, species, fin, jac_dens_in, jac_dens_gyroavg);
  gkyl_array_set_range(pm->dens_gyroavg_init, 1.0, jac_dens_gyroavg, &app->local);
  gkyl_array_accumulate_range(pm->dens_gyroavg_init, -1.0, jac_dens_in, &app->local);
  // Without the difference, the gyroaveraged density of the initial condition is its density.
  gkyl_array_set_range(jac_dens_gyroavg, 1.0, jac_dens_in, &app->local);

  species->gyroaverage_density = gk_species_flr_pade_mu_gyroaverage_density;
}

static void
gk_species_flr_pade_mu_gyroaverage_moment(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_mom_in, struct gkyl_array *jac_mom_gyroavg
)
{
  // Gyroaverage the particle flow (times J) of fin in each mu cell (the M1 moment over the
  // phase space cells of that mu cell) with the operator that gyroaverages A_parallel in that
  // mu cell, and add them up. The flow is not smoothed along z: the solves are done cell by
  // cell along z, as for A_parallel.
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_pade_mu *pm = &species->flr.pade_mu;
  gkyl_array_clear_range(jac_mom_gyroavg, 0.0, &app->local);
  for (int j = 0; j < pm->num_mu; j++) {
    gk_species_moment_calc(app, &pm->m1_mu, pm->mu_ranges[j], app->local, fin);
    // Boundary value: flow = (J*flow)/J.
    gkyl_dg_mul_op_range(
      &app->basis, 0, pm->buff, 0, app->gk_geom->geo_int.jacobgeo_inv, 0, pm->m1_mu.marr,
      &app->local
    );
    gkyl_deflated_fem_poisson_advance_dg(
      pm->ops[j], pm->m1_mu.marr, pm->buff, pm->buff_gyroavg, pm->buff_dg
    );
    gkyl_dg_mul_op_range(
      &app->basis, 0, pm->buff_gyroavg, 0, app->gk_geom->geo_int.jacobgeo, 0, pm->buff_gyroavg,
      &app->local
    );
    gkyl_array_accumulate_range(jac_mom_gyroavg, 1.0, pm->buff_gyroavg, &app->local);
  }
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_pade_mu_gyroaverage_em(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
)
{
  // Gyroaverage A_parallel or its time derivative with the gyroradius of each mu cell. These
  // fields are discontinuous along z (the time derivative of A_parallel comes from Ohm's law),
  // so the solves are done cell by cell along z.
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_pade_mu *pm = &species->flr.pade_mu;
  int num_basis = app->basis.num_basis;
  gkyl_dg_mul_op_range(
    &app->basis, 0, pm->buff, 0, app->gk_geom->geo_int.jacobgeo, 0, field_in, &app->local
  );
  for (int j = 0; j < pm->num_mu; j++) {
    gkyl_deflated_fem_poisson_advance_dg(
      pm->ops[j], pm->buff, field_in, pm->buff_gyroavg, pm->buff_dg
    );
    gkyl_array_set_offset(field_gyroavg, 1.0, pm->buff_gyroavg, j * num_basis);
  }
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_pade_mu_gyroaverage_apar(gkyl_gyrokinetic_app *app, struct gk_species *species)
{
  gk_species_flr_pade_mu_gyroaverage_em(
    app, species, species->flr.pade_mu.apar_src, species->gyro_apar
  );
}

static void
gk_species_flr_pade_mu_gyroaverage_apardot(gkyl_gyrokinetic_app *app, struct gk_species *species)
{
  gk_species_flr_pade_mu_gyroaverage_em(app, species, app->field->apardot, species->gyro_apardot);
}

static void
gk_species_flr_pade_mu_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  const struct gk_species_flr_pade_mu *pm = &gks->flr.pade_mu;
  for (int j = 0; j < pm->num_mu; j++) {
    gkyl_deflated_fem_poisson_release(pm->ops[j]);
  }
  gkyl_free(pm->ops);
  gkyl_free(pm->mu_avg_weight);
  gkyl_array_release(pm->buff);
  gkyl_array_release(pm->buff_dens);
  gkyl_array_release(pm->buff_gyroavg);
  gkyl_array_release(pm->buff_dg);
  gkyl_free(pm->mu_ranges);
  gk_species_moment_release(app, &pm->m0_mu);
  gkyl_array_release(pm->dens_gyroavg_init);
  if (pm->use_em) {
    gkyl_array_release(pm->apar_src);
    gk_species_moment_release(app, &pm->m1_mu);
  }
}

static void
gk_species_flr_pade_mu_new(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Only available in 2x and 3x with 2 velocity dimensions.
  assert(app->cdim > 1 && gks->local_vel.ndim == 2);

  struct gk_species_flr_pade_mu *pm = &gks->flr.pade_mu;
  // The objects for electromagnetic fields are created in gk_species_flr_em_new.
  pm->use_em = false;
  int num_basis = app->basis.num_basis;
  long size = app->local_ext.volume;
  double charge = gks->info.charge, mass = gks->info.mass;

  int num_mu = gkyl_range_shape(&gks->local_vel, 1);
  pm->num_mu = num_mu;
  double *mu_cells = gkyl_malloc(sizeof(double[num_mu]));
  pm->mu_avg_weight = gkyl_malloc(sizeof(double[num_mu]));
  gk_species_flr_pade_mu_grid(app, gks, mu_cells, pm->mu_avg_weight);

  pm->buff = mkarr(app->use_gpu, num_basis, size);
  pm->buff_dens = mkarr(app->use_gpu, num_basis, size);
  pm->buff_gyroavg = mkarr(app->use_gpu, num_basis, size);
  pm->buff_dg = mkarr(app->use_gpu, num_basis, size);

  // Laplacian weight of the operator of a mu cell: rho_j^2/4 = mu_j*m/(2*q^2*B) times
  // J*g^ij. The gyroradius varies with B inside the divergence, which keeps the operator
  // symmetric and the gyroaverage of a constant equal to that constant. Store J*g^ij/B,
  // smoothed along z, and scale it for each mu cell.
  int num_metric = 2 * (app->cdim - 1) - 1;
  const struct gkyl_array *metric[] = {
    app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj
  };
  struct gkyl_array *metric_over_bmag = mkarr(app->use_gpu, num_metric * num_basis, size);
  struct gkyl_array *rhoSqD4 = mkarr(app->use_gpu, num_metric * num_basis, size);
  for (int i = 0; i < num_metric; i++) {
    // J*g^ij/B = (J*g^ij) * (1/(J*B)) * J.
    gkyl_dg_mul_op_range(
      &app->basis, 0, pm->buff, 0, app->gk_geom->geo_int.jacobtot_inv, 0, metric[i], &app->local_ext
    );
    gkyl_dg_mul_op_range(
      &app->basis, i, metric_over_bmag, 0, app->gk_geom->geo_int.jacobgeo, 0, pm->buff,
      &app->local_ext
    );
  }
  // Factor multiplying the field in the FLR operator.
  struct gkyl_array *kSq = mkarr(app->use_gpu, num_basis, size);
  gkyl_array_set(kSq, -1.0, app->gk_geom->geo_int.jacobgeo);
  gk_species_flr_pade_mu_smooth_coeff(app, metric_over_bmag, pm->buff_dens);
  gk_species_flr_pade_mu_smooth_coeff(app, kSq, pm->buff_dens);

  // Gyroaverage BCs: the input field is its own boundary value at non-periodic boundaries.
  struct gkyl_poisson_bc flr_bc = app->field->poisson_bcs;
  for (int d = 0; d < app->cdim - 1; d++) {
    if (flr_bc.lo_type[d] != GKYL_POISSON_PERIODIC) {
      flr_bc.lo_type[d] = GKYL_POISSON_DIRICHLET_VARYING;
    }
    if (flr_bc.up_type[d] != GKYL_POISSON_PERIODIC) {
      flr_bc.up_type[d] = GKYL_POISSON_DIRICHLET_VARYING;
    }
  }
  // Deflated Poisson solves are performed on range assuming decomposition is *only* in z.
  pm->ops = gkyl_malloc(sizeof(struct gkyl_deflated_fem_poisson *[num_mu]));
  for (int j = 0; j < num_mu; j++) {
    gkyl_array_set(rhoSqD4, mu_cells[j] * mass / (2.0 * charge * charge), metric_over_bmag);
    pm->ops[j] = gkyl_deflated_fem_poisson_new(
      app->grid, app->basis_on_dev, app->basis, app->local, app->local, rhoSqD4, kSq, flr_bc, NULL,
      app->use_gpu
    );
  }
  gkyl_array_release(metric_over_bmag);
  gkyl_array_release(rhoSqD4);
  gkyl_array_release(kSq);
  gkyl_free(mu_cells);

  // Phase space range of each mu cell, and the moment that computes the density (times J) of
  // the distribution in one of them.
  int pdim = gks->local.ndim;
  pm->mu_ranges = gkyl_malloc(sizeof(struct gkyl_range[num_mu]));
  for (int j = 0; j < num_mu; j++) {
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    for (int d = 0; d < pdim; d++) {
      lower[d] = gks->local.lower[d];
      upper[d] = gks->local.upper[d];
    }
    lower[pdim - 1] = upper[pdim - 1] = gks->local.lower[pdim - 1] + j;
    gkyl_sub_range_init(&pm->mu_ranges[j], &gks->local, lower, upper);
  }
  gk_species_moment_init(app, gks, &pm->m0_mu, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_WEIGHT_NONE, false);
  pm->dens_gyroavg_init = mkarr(app->use_gpu, num_basis, size);
  gks->flr.dens_offset_from_ic = true;

  // Gyroaveraged potential and density, and the potential of the Hamiltonian in each mu cell.
  gks->gyro_phi = mkarr(app->use_gpu, num_basis, size);
  gks->m0_gyroavg = mkarr(app->use_gpu, num_basis, size);
  gks->hamil_num_mu = num_mu;
  gks->hamil_phi = mkarr(app->use_gpu, num_mu * num_basis, size);

  gks->gyroaverage = gk_species_flr_pade_mu_gyroaverage;
  gks->gyroaverage_density = gk_species_flr_pade_mu_gyroaverage_density_init;
  // Set in gk_species_flr_em_new for electromagnetic fields, the only users of the flow.
  gks->gyroaverage_moment = gk_species_flr_none_gyroaverage_moment;
  gks->flr.release_func = gk_species_flr_pade_mu_release;
}

/** Selection of the model. **/

void
gk_species_flr_new(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // The species is gyroaveraged if the field enables FLR effects and the species provides a
  // reference gyroradius, rho_s = sqrt(Tperp*m_s)/(|q_s|*bmag) unless given directly.
  gks->use_flr = app->field->info.flr.type != GKYL_GK_FLR_NONE &&
                 (gks->info.flr.gyroradius > 0.0 || gks->info.flr.Tperp > 0.0);
  gks->flr_rhoSq_ref = 0.0;
  gks->flr.dens_offset_from_ic = false;
  if (gks->use_flr) {
    double gyroradius_bmag = gks->info.flr.bmag ? gks->info.flr.bmag : app->bmag_ref;
    gks->flr_rhoSq_ref =
      gks->info.flr.gyroradius > 0.0 ?
        pow(gks->info.flr.gyroradius, 2.0) :
        gks->info.flr.Tperp * gks->info.mass / pow(gks->info.charge * gyroradius_bmag, 2.0);
  }

  switch (gks->use_flr ? app->field->info.flr.type : GKYL_GK_FLR_NONE) {
    case GKYL_GK_FLR_PADE_MU:
      gk_species_flr_pade_mu_new(app, gks);
      break;
    default:
      gk_species_flr_none_new(app, gks);
      break;
  }
  // A_parallel and its time derivative are gyroaveraged only with electromagnetic fields (see
  // gk_species_flr_em_new).
  gks->gyroaverage_apar = gk_species_flr_none_gyroaverage_em;
  gks->gyroaverage_apardot = gk_species_flr_none_gyroaverage_em;
}

void
gk_species_flr_em_new(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  int num_basis = app->basis.num_basis;
  long size = app->local_ext.volume;
  // Gyroaveraged particle flow (for the current density of the electromagnetic field equations).
  gks->m1_gyroavg = gks->use_flr ? mkarr(app->use_gpu, num_basis, size) :
                                   gkyl_array_acquire(gks->m1.marr);

  // A_parallel used by the collisionless terms. Only EM uses A_parallel at the current RK stage
  // (apar_curr only exists if EM).
  struct gkyl_array *apar = gks->info.collisionless.type == GKYL_GK_COLLISIONLESS_EM ?
                              app->field->apar_curr :
                              app->field->apar;
  bool pade_mu_em = gks->use_flr && app->field->info.flr.type == GKYL_GK_FLR_PADE_MU &&
                    app->field->info.gkfield_id == GKYL_GK_FIELD_EM;
  if (pade_mu_em) {
    // Pade of each mu cell: A_parallel and its time derivative in each mu cell, as the potential
    // of the Hamiltonian, and the particle flow of each mu cell gyroaveraged with the same
    // operator.
    struct gk_species_flr_pade_mu *pm = &gks->flr.pade_mu;
    pm->use_em = true;
    pm->apar_src = gkyl_array_acquire(apar);
    gks->gyro_apar = mkarr(app->use_gpu, pm->num_mu * num_basis, size);
    gks->gyro_apardot = mkarr(app->use_gpu, pm->num_mu * num_basis, size);
    gkyl_array_clear(gks->gyro_apar, 0.0);
    gkyl_array_clear(gks->gyro_apardot, 0.0);
    gk_species_moment_init(app, gks, &pm->m1_mu, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_WEIGHT_NONE, false);

    gks->gyroaverage_moment = gk_species_flr_pade_mu_gyroaverage_moment;
    gks->gyroaverage_apar = gk_species_flr_pade_mu_gyroaverage_apar;
    gks->gyroaverage_apardot = gk_species_flr_pade_mu_gyroaverage_apardot;
    return;
  }
  // A_parallel and its time derivative are not gyroaveraged.
  gks->gyro_apar = gkyl_array_acquire(apar);
  gks->gyro_apardot = gkyl_array_acquire(app->field->apardot);
}

void
gk_species_flr_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  gkyl_array_release(gks->gyro_phi);
  gkyl_array_release(gks->m0_gyroavg);
  gkyl_array_release(gks->hamil_phi);
  gks->flr.release_func(app, gks);
}
