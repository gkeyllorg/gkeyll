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
// - Pade: the potential, the density and the flow are gyroaveraged with
//   (1 - (rho^2/2) nabla_perp^2)^{-1}, rho being the reference gyroradius of the species;
//   A_parallel is not.
// - gyromatrix: the potential, A_parallel and its time derivative are gyroaveraged over the
//   gyro-ring of each mu cell, and the density and the flow in each mu cell with the adjoint
//   operator.
// - Pade of each mu cell: as the gyromatrix model, with the Pade approximant
//   (1 - div_perp((rho_j^2/4) grad_perp))^{-1} in place of the gyroaverage over the gyro-ring,
//   rho_j^2 = 2*m*mu_j/(q^2*B) being the squared gyroradius of the mu cell j. The field
//   equations are those without FLR effects (Ishizawa et al., Phys. Plasmas 26, 082301 (2019)).

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

/** Pade approximant of the gyroaverage. **/

static void
gk_species_flr_pade_gyroaverage(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
)
{
  // Solve the gyroaveraging problem J*(1 - (rho_s^2/2) nabla_perp^2) phi_g = J*phi_in
  // with phi_in as the boundary value, e.g. to get the potential used in the Hamiltonian.
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_pade *pade = &species->flr.pade;
  gkyl_dg_mul_op_range(
    &app->basis, 0, pade->buff, 0, app->gk_geom->geo_int.jacobgeo, 0, field_in, &app->local
  );
  gkyl_deflated_fem_poisson_advance(pade->op, pade->buff, field_in, field_gyroavg);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_pade_gyroaverage_moment(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_mom_in, struct gkyl_array *jac_mom_gyroavg
)
{
  // Same as gk_species_flr_pade_gyroaverage but for a moment that already
  // carries the Jacobian (e.g. the density J*n), returning J*n_g.
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_pade *pade = &species->flr.pade;
  // Smooth the moment along z to get a consistent field for the deflated solve.
  app->field->fem_projection_par_rho_func(app, app->field, jac_mom_in, pade->buff_dens);
  // Boundary value n = (J*n)/J.
  gkyl_dg_mul_op_range(
    &app->basis, 0, pade->buff, 0, app->gk_geom->geo_int.jacobgeo_inv, 0, pade->buff_dens,
    &app->local
  );
  gkyl_deflated_fem_poisson_advance(pade->op, pade->buff_dens, pade->buff, jac_mom_gyroavg);
  gkyl_dg_mul_op_range(
    &app->basis, 0, jac_mom_gyroavg, 0, app->gk_geom->geo_int.jacobgeo, 0, jac_mom_gyroavg,
    &app->local
  );
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_pade_gyroaverage_density(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_dens_in, struct gkyl_array *jac_dens_gyroavg
)
{
  // The Pade gyroaverage does not depend on mu: gyroaverage the density moment.
  gk_species_flr_pade_gyroaverage_moment(app, species, fin, jac_dens_in, jac_dens_gyroavg);
}

static void
gk_species_flr_pade_smooth_coeff(
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
gk_species_flr_pade_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  const struct gk_species_flr_pade *pade = &gks->flr.pade;
  gkyl_array_release(pade->rhoSqD2);
  gkyl_array_release(pade->kSq);
  gkyl_array_release(pade->buff);
  gkyl_array_release(pade->buff_dens);
  gkyl_deflated_fem_poisson_release(pade->op);
}

static void
gk_species_flr_pade_new(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  assert(app->cdim > 1);
  struct gk_species_flr_pade *pade = &gks->flr.pade;
  int num_basis = app->basis.num_basis;

  gks->gyro_phi = mkarr(app->use_gpu, num_basis, app->local_ext.volume);
  gks->m0_gyroavg = mkarr(app->use_gpu, num_basis, app->local_ext.volume);
  gk_species_flr_hamil_mu_independent(app, gks);
  pade->buff = mkarr(app->use_gpu, num_basis, app->local_ext.volume);
  pade->buff_dens = mkarr(app->use_gpu, num_basis, app->local_ext.volume);

  double flr_weight = gks->flr_rhoSq_ref;
  // Initialize the weight in the Laplacian operator.
  pade->rhoSqD2 = mkarr(app->use_gpu, (2 * (app->cdim - 1) - 1) * num_basis, app->local_ext.volume);
  gkyl_array_set_offset(pade->rhoSqD2, 0.5 * flr_weight, app->gk_geom->geo_int.gxxj, 0 * num_basis);
  if (app->cdim > 2) {
    gkyl_array_set_offset(
      pade->rhoSqD2, 0.5 * flr_weight, app->gk_geom->geo_int.gxyj, 1 * num_basis
    );
    gkyl_array_set_offset(
      pade->rhoSqD2, 0.5 * flr_weight, app->gk_geom->geo_int.gyyj, 2 * num_basis
    );
  }
  // Initialize the factor multiplying the field in the FLR operator.
  pade->kSq = mkarr(app->use_gpu, num_basis, app->local_ext.volume);
  gkyl_array_set(pade->kSq, -1.0, app->gk_geom->geo_int.jacobgeo);
  gk_species_flr_pade_smooth_coeff(app, pade->rhoSqD2, pade->buff_dens);
  gk_species_flr_pade_smooth_coeff(app, pade->kSq, pade->buff_dens);

  // Gyroaverage BCs: the input field is its own boundary value at
  // non-periodic boundaries.
  struct gkyl_poisson_bc flr_bc = app->field->poisson_bcs;
  for (int d = 0; d < app->cdim - 1; d++) {
    if (flr_bc.lo_type[d] != GKYL_POISSON_PERIODIC) {
      flr_bc.lo_type[d] = GKYL_POISSON_DIRICHLET_VARYING;
    }
    if (flr_bc.up_type[d] != GKYL_POISSON_PERIODIC) {
      flr_bc.up_type[d] = GKYL_POISSON_DIRICHLET_VARYING;
    }
  }
  // Deflated Poisson solve is performed on range assuming decomposition is *only* in z.
  pade->op = gkyl_deflated_fem_poisson_new(
    app->grid, app->basis_on_dev, app->basis, app->local, app->local, pade->rhoSqD2, pade->kSq,
    flr_bc, NULL, app->use_gpu
  );

  gks->gyroaverage = gk_species_flr_pade_gyroaverage;
  gks->gyroaverage_density = gk_species_flr_pade_gyroaverage_density;
  gks->gyroaverage_moment = gk_species_flr_pade_gyroaverage_moment;
  gks->flr.release_func = gk_species_flr_pade_release;
}

/** Gyromatrix: gyroaverage over the gyro-ring of each mu cell. **/

static void
gk_species_flr_gyromatrix_gyroaverage(
  gkyl_gyrokinetic_app *app, struct gk_species *species, struct gkyl_array *field_in,
  struct gkyl_array *field_gyroavg
)
{
  // Gyroaverage the potential over the gyro-ring of each mu cell. The Hamiltonian uses this
  // potential, constant in each mu cell. Also compute its average over a Maxwellian in mu
  // (used wherever a potential independent of mu is needed).
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_gyromatrix *gm = &species->flr.gyromatrix;
  gkyl_deflated_fem_gyroaverage_advance(gm->gyroavg, field_in, species->hamil_phi);
  gkyl_array_clear_range(field_gyroavg, 0.0, &app->local);
  for (int j = 0; j < gm->num_mu; j++) {
    gkyl_array_accumulate_offset_range(
      field_gyroavg, gm->mu_avg_weight[j], species->hamil_phi, j * app->basis.num_basis, &app->local
    );
  }
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_gyromatrix_gyroaverage_density(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_dens_in, struct gkyl_array *jac_dens_gyroavg
)
{
  // Gyroaverage the density (times J) with the adjoint of the gyroaverage of
  // the potential, applied to the density in each mu cell, so that the
  // integral of f*q*<phi> over phase space equals the integral of phi times
  // the gyroaveraged charge density. The difference between the gyroaveraged
  // and gyrocenter densities of the initial condition is removed, see
  // gk_species_flr_gyromatrix_gyroaverage_density_init.
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_gyromatrix *gm = &species->flr.gyromatrix;
  gkyl_gk_dens_mu_cells_advance(gm->dens_mu_op, &species->local, &app->local, fin, gm->dens_cells);
  gkyl_deflated_fem_gyroaverage_advance_adjoint(gm->gyroavg, gm->dens_cells, jac_dens_gyroavg);
  gkyl_array_accumulate_range(jac_dens_gyroavg, -1.0, gm->dens_gyroavg_init, &app->local);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_gyromatrix_gyroaverage_density_init(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_dens_in, struct gkyl_array *jac_dens_gyroavg
)
{
  // First gyroaverage of the density, done with the initial condition. The
  // gyroradius depends on position through the magnetic field, so even a
  // uniform gyrocenter density has a gyroaveraged density that differs from it
  // by O((rho/L_B)^2). This is of higher order than the model retains, yet the
  // long wavelength polarization turns it into an O(1) error in the potential.
  // Hence store the difference between the gyroaveraged and the gyrocenter
  // densities of the initial condition, and remove it from the gyroaveraged
  // density from now on; the initial charge density then is that of the
  // gyrocenters, and FLR effects act on departures from the initial condition.
  struct gk_species_flr_gyromatrix *gm = &species->flr.gyromatrix;
  gkyl_gk_dens_mu_cells_advance(gm->dens_mu_op, &species->local, &app->local, fin, gm->dens_cells);
  gkyl_deflated_fem_gyroaverage_advance_adjoint(gm->gyroavg, gm->dens_cells, gm->dens_gyroavg_init);
  gkyl_array_accumulate_range(gm->dens_gyroavg_init, -1.0, jac_dens_in, &app->local);

  species->gyroaverage_density = gk_species_flr_gyromatrix_gyroaverage_density;
  species->gyroaverage_density(app, species, fin, jac_dens_in, jac_dens_gyroavg);
}

static void
gk_species_flr_gyromatrix_gyroaverage_moment(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *jac_mom_in, struct gkyl_array *jac_mom_gyroavg
)
{
  // Gyroaverage the particle flow (times J) of fin with the adjoint of the gyroaverage of
  // A_parallel, applied to the flow in each mu cell, so that the integral of f*q*vpar*<A_par>
  // over phase space equals the integral of A_par times the gyroaveraged current density. The
  // flow in each mu cell is the M1 moment over the phase space cells of that mu cell. Unlike the
  // density no initial offset is removed: the flow of the initial condition, and the offset it
  // would have, are small, and the time derivative of the flow (Ohm's law) has none.
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_gyromatrix *gm = &species->flr.gyromatrix;
  int pdim = species->local.ndim, num_basis = app->basis.num_basis;
  for (int j = 0; j < gm->num_mu; j++) {
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    for (int d = 0; d < pdim; d++) {
      lower[d] = species->local.lower[d];
      upper[d] = species->local.upper[d];
    }
    lower[pdim - 1] = upper[pdim - 1] = species->local.lower[pdim - 1] + j;
    struct gkyl_range mu_range;
    gkyl_sub_range_init(&mu_range, &species->local, lower, upper);
    gk_species_moment_calc(app, &gm->m1_mu, mu_range, app->local, fin);
    gkyl_array_set_offset(gm->flow_cells, 1.0, gm->m1_mu.marr, j * num_basis);
  }
  gkyl_deflated_fem_gyroaverage_advance_adjoint(gm->gyroavg, gm->flow_cells, jac_mom_gyroavg);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_gyromatrix_gyroaverage_apar(gkyl_gyrokinetic_app *app, struct gk_species *species)
{
  // Gyroaverage A_parallel over the gyro-ring of each mu cell (cell by cell along z, so that the
  // adjoint is that of the flow, and a field discontinuous along z is not smoothed).
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_gyromatrix *gm = &species->flr.gyromatrix;
  gkyl_deflated_fem_gyroaverage_advance_dg(gm->gyroavg, gm->apar_src, species->gyro_apar);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_flr_gyromatrix_gyroaverage_apardot(
  gkyl_gyrokinetic_app *app, struct gk_species *species
)
{
  // Gyroaverage the time derivative of A_parallel from Ohm's law, which is discontinuous along
  // z, cell by cell along z (see gk_species_flr_gyromatrix_gyroaverage_apar).
  struct timespec wst = gkyl_wall_clock();
  struct gk_species_flr_gyromatrix *gm = &species->flr.gyromatrix;
  gkyl_deflated_fem_gyroaverage_advance_dg(gm->gyroavg, app->field->apardot, species->gyro_apardot);
  app->stat.species_gyroavg_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_flr_gyromatrix_accumulate_polarization(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, const struct gkyl_array *phi,
  struct gkyl_array *out
)
{
  struct gk_species_flr_gyromatrix *gm = &gks->flr.gyromatrix;

  // Weighted gyroaverages W_j G_j phi (the density in each mu cell is recomputed before it is
  // next used, so its array serves as a buffer here).
  gkyl_deflated_fem_gyroaverage_advance(gm->gyroavg, phi, gm->pol_gyroavg_buff);
  struct gkyl_array *weighted = gm->dens_cells;
  for (int j = 0; j < gm->num_mu; j++) {
    gkyl_dg_mul_op_range(
      &app->basis, j, weighted, j, gm->pol_wgt, j, gm->pol_gyroavg_buff, &app->local
    );
  }

  // out += pol_dens*phi - sum_j G_j^dagger W_j G_j phi.
  gkyl_deflated_fem_gyroaverage_advance_adjoint(gm->gyroavg, weighted, gm->pol_buff);
  gkyl_array_accumulate_range(out, -1.0, gm->pol_buff, &app->local);
  gkyl_dg_mul_op_range(&app->basis, 0, gm->pol_buff, 0, gm->pol_dens, 0, phi, &app->local);
  gkyl_array_accumulate_range(out, 1.0, gm->pol_buff, &app->local);
}

static void
gk_species_flr_gyromatrix_polarization_new(
  gkyl_gyrokinetic_app *app, struct gk_species *gks,
  const struct gkyl_deflated_fem_gyroaverage_inp *gyroavg_inp, const double *mu_edges
)
{
  // Polarization of the species, P phi = pol_dens*phi - sum_j G_j^dagger W_j G_j phi, where G_j
  // is the gyroaverage in the mu cell j, W_j is proportional to the fraction of a reference
  // Maxwellian F_0 in that mu cell, and pol_dens = sum_j G_j^dagger W_j so that P vanishes for
  // a constant potential. The thermal gyroradius of F_0 is the reference gyroradius rho of the
  // species everywhere, and the W_j add up to q^2 n/T = eps/rho^2 (times J), with
  // eps = n m/B^2 the polarization weight of the species, times the ratio of the thermal
  // magnetic moment of F_0 to the average over F_0 of the magnetic moments the gyroaverages
  // are done at. At long wavelength P then is the polarization, -div(eps grad_perp), of the
  // field equation without FLR effects, whatever the mu grid. This holds up to a non-periodic
  // boundary where the potential vanishes: the gyroaverage extends the potential beyond it as
  // an odd function, and the product of two potentials, which pol_dens weighs, is even, so the
  // adjoint in pol_dens is that of the gyroaverage that reflects the field about the boundary.
  struct gk_species_flr_gyromatrix *gm = &gks->flr.gyromatrix;
  int cdim = app->cdim, num_basis = app->basis.num_basis, num_mu = gm->num_mu;
  long size = app->local_ext.volume;

  double polarization_bmag =
    app->field->info.polarization_bmag ? app->field->info.polarization_bmag : app->bmag_ref;
  double pol_weight = gks->info.polarization_density * gks->info.mass / pow(polarization_bmag, 2.0);
  double pol_fac = pol_weight / gks->flr_rhoSq_ref;
  // Thermal magnetic moment of F_0 over the magnetic field: rho^2 = 2*m*mu/(q^2*B) averages to
  // 2*rho_ref^2 when mu averages to rho_ref^2*q^2*B/m.
  double mu_th_over_bmag = gks->flr_rhoSq_ref * pow(gks->info.charge, 2.0) / gks->info.mass;

  // The weights are computed on the host, at the nodes of the configuration space basis.
  struct gkyl_array *bmag_ho = mkarr(false, num_basis, size);
  struct gkyl_array *jacobgeo_ho = mkarr(false, num_basis, size);
  struct gkyl_array *pol_wgt_ho = mkarr(false, num_mu * num_basis, size);
  gkyl_array_copy(bmag_ho, app->gk_geom->geo_corn.bmag);
  gkyl_array_copy(jacobgeo_ho, app->gk_geom->geo_int.jacobgeo);

  double node_coords[GKYL_MAX_CDIM * num_basis];
  app->basis.node_list(node_coords);
  double *wgt_nod = gkyl_malloc(sizeof(double[num_mu * num_basis]));

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&app->local, iter.idx);
    const double *bmag_d = gkyl_array_cfetch(bmag_ho, linidx);
    const double *jacobgeo_d = gkyl_array_cfetch(jacobgeo_ho, linidx);

    for (int a = 0; a < num_basis; a++) {
      double bmag_a = app->basis.eval_expand(&node_coords[a * cdim], bmag_d);
      double jacobgeo_a = app->basis.eval_expand(&node_coords[a * cdim], jacobgeo_d);
      double mu_th = mu_th_over_bmag * bmag_a;

      // Integral of F_0 = exp(-mu/mu_th) over each mu cell, and of F_0 times the magnetic
      // moment of the gyroaverage of each mu cell.
      double weight_tot = 0.0, weight_mu_tot = 0.0;
      for (int j = 0; j < num_mu; j++) {
        double weight = exp(-mu_edges[j] / mu_th) - exp(-mu_edges[j + 1] / mu_th);
        wgt_nod[j * num_basis + a] = weight;
        weight_tot += weight;
        weight_mu_tot += weight * gyroavg_inp->mu[j];
      }
      double norm = pol_fac * jacobgeo_a * mu_th / weight_mu_tot;
      for (int j = 0; j < num_mu; j++) {
        wgt_nod[j * num_basis + a] *= norm;
      }
    }

    double *pol_wgt_d = gkyl_array_fetch(pol_wgt_ho, linidx);
    for (int j = 0; j < num_mu; j++) {
      app->basis.nodal_to_modal(&wgt_nod[j * num_basis], &pol_wgt_d[j * num_basis]);
    }
  }

  gm->pol_wgt = mkarr(app->use_gpu, num_mu * num_basis, size);
  gm->pol_dens = mkarr(app->use_gpu, num_basis, size);
  gm->pol_gyroavg_buff = mkarr(app->use_gpu, num_mu * num_basis, size);
  gm->pol_buff = mkarr(app->use_gpu, num_basis, size);
  gkyl_array_copy(gm->pol_wgt, pol_wgt_ho);

  // pol_dens = sum_j G_j^dagger W_j with the reflecting gyroaverage, only needed here (so it
  // is created on the host).
  struct gkyl_deflated_fem_gyroaverage_inp refl_inp = *gyroavg_inp;
  refl_inp.reflect_bnd = true;
  refl_inp.use_gpu = false;
  struct gkyl_deflated_fem_gyroaverage *gyroavg_refl =
    gkyl_deflated_fem_gyroaverage_inew(&refl_inp);
  struct gkyl_array *pol_dens_ho = mkarr(false, num_basis, size);
  gkyl_deflated_fem_gyroaverage_advance_adjoint(gyroavg_refl, pol_wgt_ho, pol_dens_ho);
  gkyl_array_copy(gm->pol_dens, pol_dens_ho);
  gkyl_array_release(pol_dens_ho);
  gkyl_deflated_fem_gyroaverage_release(gyroavg_refl);

  gkyl_free(wgt_nod);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(jacobgeo_ho);
  gkyl_array_release(pol_wgt_ho);
}

static void
gk_species_flr_gyromatrix_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  const struct gk_species_flr_gyromatrix *gm = &gks->flr.gyromatrix;
  gkyl_free(gm->mu_avg_weight);
  gkyl_deflated_fem_gyroaverage_release(gm->gyroavg);
  gkyl_gk_dens_mu_cells_release(gm->dens_mu_op);
  gkyl_array_release(gm->dens_cells);
  gkyl_array_release(gm->dens_gyroavg_init);
  gkyl_array_release(gm->pol_wgt);
  gkyl_array_release(gm->pol_dens);
  gkyl_array_release(gm->pol_gyroavg_buff);
  gkyl_array_release(gm->pol_buff);
  if (gm->use_em) {
    gkyl_array_release(gm->apar_src);
    gkyl_array_release(gm->flow_cells);
    gk_species_moment_release(app, &gm->m1_mu);
  }
}

static void
gk_species_flr_gyromatrix_new(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Only available in 2x and 3x with 2 velocity dimensions.
  assert(app->cdim > 1 && gks->local_vel.ndim == 2);

  struct gk_species_flr_gyromatrix *gm = &gks->flr.gyromatrix;
  // The objects for electromagnetic fields are created in gk_species_flr_em_new.
  gm->use_em = false;
  int cdim = app->cdim, num_basis = app->basis.num_basis;
  double charge = gks->info.charge, mass = gks->info.mass;

  // Magnetic moment at the boundaries of the mu cells.
  int num_mu = gkyl_range_shape(&gks->local_vel, 1);
  gm->num_mu = num_mu;
  double *mu_edges = gkyl_malloc(sizeof(double[num_mu + 1]));
  const struct gkyl_velocity_map *gvm = gks->vel_map;
  int vmap_num_basis = gvm->vmap_basis_ho.num_basis;
  for (int j = 0; j < num_mu; j++) {
    int vidx[] = {gks->local_vel.lower[0], gks->local_vel.lower[1] + j};
    const double *vmap_d = gkyl_array_cfetch(gvm->vmap_ho, gkyl_range_idx(&gvm->local_vel, vidx));
    const double *mu_map_d = &vmap_d[vmap_num_basis];
    double cell_lo = -1.0, cell_up = 1.0;
    mu_edges[j] = fmax(0.0, gvm->vmap_basis_ho.eval_expand(&cell_lo, mu_map_d));
    mu_edges[j + 1] = gvm->vmap_basis_ho.eval_expand(&cell_up, mu_map_d);
  }

  // The particles of a mu cell are all gyroaveraged over the same gyro-ring, that of the
  // average magnetic moment in the cell of a Maxwellian in mu, exp(-mu*B/T), which also gives
  // the weight of the mu cell in averages over mu. T/B follows from the reference gyroradius of
  // the species, rho^2 = T*m/(q*B)^2. With this choice the gyroaverage of such a Maxwellian is
  // exact to second order in the gyroradius, whatever the mu grid.
  double gyroradius_bmag = gks->info.flr.bmag ? gks->info.flr.bmag : app->bmag_ref;
  double temp_over_bmag = gks->flr_rhoSq_ref * pow(charge, 2.0) * gyroradius_bmag / mass;
  double *mu_cells = gkyl_malloc(sizeof(double[num_mu]));
  gm->mu_avg_weight = gkyl_malloc(sizeof(double[num_mu]));
  double weight_tot = 0.0;
  for (int j = 0; j < num_mu; j++) {
    double cell_len = mu_edges[j + 1] - mu_edges[j];
    // Fraction of the Maxwellian at the lower edge that is left at the upper edge.
    double decay = exp(-cell_len / temp_over_bmag);
    mu_cells[j] = mu_edges[j] + temp_over_bmag - cell_len * decay / (1.0 - decay);
    gm->mu_avg_weight[j] = exp(-mu_edges[j] / temp_over_bmag) * (1.0 - decay);
    weight_tot += gm->mu_avg_weight[j];
  }
  for (int j = 0; j < num_mu; j++) {
    gm->mu_avg_weight[j] /= weight_tot;
  }

  // Contravariant metric coefficients, smoothed along z so that the
  // gyroaverage on each plane of constant z does not depend on which side of
  // the plane the geometry is evaluated from.
  struct gkyl_array *metric[3];
  int metric_comp[] = {0, 1, 3}; // Location of g^xx, g^xy, g^yy in the metric.
  for (int i = 0; i < 3; i++) {
    metric[i] = mkarr(app->use_gpu, num_basis, app->local_ext.volume);
    gkyl_array_set_offset(metric[i], 1.0, app->gk_geom->geo_int.gij, metric_comp[i] * num_basis);
    app->field->fem_projection_par_rho_func(app, app->field, metric[i], metric[i]);
  }

  struct gkyl_deflated_fem_gyroaverage_inp gyroavg_inp = {
    .grid = &app->grid,
    .basis = &app->basis,
    .local = &app->local,
    .is_periodic =
      {app->field->poisson_bcs.lo_type[0] == GKYL_POISSON_PERIODIC,
       cdim > 2 && app->field->poisson_bcs.lo_type[1] == GKYL_POISSON_PERIODIC},
    .gxx = metric[0],
    .gxy = metric[1],
    .gyy = metric[2],
    .bmag = app->gk_geom->geo_corn.bmag,
    .mass = mass,
    .charge = charge,
    .num_mu = num_mu,
    .mu = mu_cells,
    .use_gpu = app->use_gpu,
  };
  gm->gyroavg = gkyl_deflated_fem_gyroaverage_inew(&gyroavg_inp);
  gk_species_flr_gyromatrix_polarization_new(app, gks, &gyroavg_inp, mu_edges);
  for (int i = 0; i < 3; i++) {
    gkyl_array_release(metric[i]);
  }
  gkyl_free(mu_edges);
  gkyl_free(mu_cells);

  gm->dens_mu_op =
    gkyl_gk_dens_mu_cells_new(&gks->grid, &app->basis, &gks->basis, mass, app->use_gpu);
  gm->dens_cells = mkarr(app->use_gpu, num_mu * num_basis, app->local_ext.volume);
  gm->dens_gyroavg_init = mkarr(app->use_gpu, num_basis, app->local_ext.volume);

  // Gyroaveraged potential and density, and the potential of the Hamiltonian in each mu cell.
  gks->gyro_phi = mkarr(app->use_gpu, num_basis, app->local_ext.volume);
  gks->m0_gyroavg = mkarr(app->use_gpu, num_basis, app->local_ext.volume);
  gks->hamil_num_mu = num_mu;
  gks->hamil_phi = mkarr(app->use_gpu, num_mu * num_basis, app->local_ext.volume);

  gks->gyroaverage = gk_species_flr_gyromatrix_gyroaverage;
  gks->gyroaverage_density = gk_species_flr_gyromatrix_gyroaverage_density_init;
  // Set in gk_species_flr_em_new for electromagnetic fields, the only users of the flow.
  gks->gyroaverage_moment = gk_species_flr_none_gyroaverage_moment;
  gks->flr.release_func = gk_species_flr_gyromatrix_release;
}

/** Pade approximant of the gyroaverage with the gyroradius of each mu cell. **/

static void
gk_species_flr_pade_mu_grid(
  gkyl_gyrokinetic_app *app, const struct gk_species *gks, double *mu_cells, double *mu_avg_weight
)
{
  // Magnetic moment each mu cell is gyroaveraged at, and its weight in averages over mu: the
  // average magnetic moment in the cell of a Maxwellian in mu, exp(-mu*B/T), and the fraction
  // of that Maxwellian in the cell (as in the gyromatrix model). T/B follows from the reference
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
    gkyl_array_accumulate_range(
      field_gyroavg, pm->mu_avg_weight[j], pm->buff_gyroavg, &app->local
    );
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
  int num_basis = app->basis.num_basis;
  gkyl_gk_dens_mu_cells_advance(pm->dens_mu_op, &species->local, &app->local, fin, pm->dens_cells);
  gkyl_array_clear_range(jac_dens_gyroavg, 0.0, &app->local);
  for (int j = 0; j < pm->num_mu; j++) {
    gkyl_array_set_offset(pm->buff_gyroavg, 1.0, pm->dens_cells, j * num_basis);
    // Smooth the density along z to get a consistent field for the deflated solve.
    app->field->fem_projection_par_rho_func(app, app->field, pm->buff_gyroavg, pm->buff_dens);
    // Boundary value n = (J*n)/J.
    gkyl_dg_mul_op_range(
      &app->basis, 0, pm->buff, 0, app->gk_geom->geo_int.jacobgeo_inv, 0, pm->buff_dens,
      &app->local
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
  int pdim = species->local.ndim;
  gkyl_array_clear_range(jac_mom_gyroavg, 0.0, &app->local);
  for (int j = 0; j < pm->num_mu; j++) {
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    for (int d = 0; d < pdim; d++) {
      lower[d] = species->local.lower[d];
      upper[d] = species->local.upper[d];
    }
    lower[pdim - 1] = upper[pdim - 1] = species->local.lower[pdim - 1] + j;
    struct gkyl_range mu_range;
    gkyl_sub_range_init(&mu_range, &species->local, lower, upper);
    gk_species_moment_calc(app, &pm->m1_mu, mu_range, app->local, fin);
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
  gkyl_gk_dens_mu_cells_release(pm->dens_mu_op);
  gkyl_array_release(pm->dens_cells);
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
  // smoothed along z as in the Pade model, and scale it for each mu cell.
  int num_metric = 2 * (app->cdim - 1) - 1;
  const struct gkyl_array *metric[] = {
    app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj
  };
  struct gkyl_array *metric_over_bmag = mkarr(app->use_gpu, num_metric * num_basis, size);
  struct gkyl_array *rhoSqD4 = mkarr(app->use_gpu, num_metric * num_basis, size);
  for (int i = 0; i < num_metric; i++) {
    // J*g^ij/B = (J*g^ij) * (1/(J*B)) * J.
    gkyl_dg_mul_op_range(
      &app->basis, 0, pm->buff, 0, app->gk_geom->geo_int.jacobtot_inv, 0, metric[i],
      &app->local_ext
    );
    gkyl_dg_mul_op_range(
      &app->basis, i, metric_over_bmag, 0, app->gk_geom->geo_int.jacobgeo, 0, pm->buff,
      &app->local_ext
    );
  }
  // Factor multiplying the field in the FLR operator.
  struct gkyl_array *kSq = mkarr(app->use_gpu, num_basis, size);
  gkyl_array_set(kSq, -1.0, app->gk_geom->geo_int.jacobgeo);
  gk_species_flr_pade_smooth_coeff(app, metric_over_bmag, pm->buff_dens);
  gk_species_flr_pade_smooth_coeff(app, kSq, pm->buff_dens);

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
      app->grid, app->basis_on_dev, app->basis, app->local, app->local, rhoSqD4, kSq, flr_bc,
      NULL, app->use_gpu
    );
  }
  gkyl_array_release(metric_over_bmag);
  gkyl_array_release(rhoSqD4);
  gkyl_array_release(kSq);
  gkyl_free(mu_cells);

  pm->dens_mu_op =
    gkyl_gk_dens_mu_cells_new(&gks->grid, &app->basis, &gks->basis, mass, app->use_gpu);
  pm->dens_cells = mkarr(app->use_gpu, num_mu * num_basis, size);
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
    case GKYL_GK_FLR_PADE:
      gk_species_flr_pade_new(app, gks);
      break;
    case GKYL_GK_FLR_GYROMATRIX:
      gk_species_flr_gyromatrix_new(app, gks);
      break;
    case GKYL_GK_FLR_PADE_MU:
      gk_species_flr_pade_mu_new(app, gks);
      break;
    default:
      gk_species_flr_none_new(app, gks);
      break;
  }
  // A_parallel and its time derivative are gyroaveraged only by the gyromatrix and Pade of each
  // mu cell models with electromagnetic fields (see gk_species_flr_em_new).
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
  bool gyromatrix_em = gks->use_flr && app->field->info.flr.type == GKYL_GK_FLR_GYROMATRIX &&
                       app->field->info.gkfield_id == GKYL_GK_FIELD_EM;
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
    gk_species_moment_init(
      app, gks, &pm->m1_mu, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_WEIGHT_NONE, false
    );

    gks->gyroaverage_moment = gk_species_flr_pade_mu_gyroaverage_moment;
    gks->gyroaverage_apar = gk_species_flr_pade_mu_gyroaverage_apar;
    gks->gyroaverage_apardot = gk_species_flr_pade_mu_gyroaverage_apardot;
    return;
  }
  if (!gyromatrix_em) {
    // A_parallel and its time derivative are not gyroaveraged.
    gks->gyro_apar = gkyl_array_acquire(apar);
    gks->gyro_apardot = gkyl_array_acquire(app->field->apardot);
    return;
  }

  // Gyromatrix: A_parallel and its time derivative in each mu cell, as the potential of the
  // Hamiltonian, and the particle flow gyroaveraged with the adjoint operator.
  struct gk_species_flr_gyromatrix *gm = &gks->flr.gyromatrix;
  gm->use_em = true;
  gm->apar_src = gkyl_array_acquire(apar);
  gks->gyro_apar = mkarr(app->use_gpu, gm->num_mu * num_basis, size);
  gks->gyro_apardot = mkarr(app->use_gpu, gm->num_mu * num_basis, size);
  gkyl_array_clear(gks->gyro_apar, 0.0);
  gkyl_array_clear(gks->gyro_apardot, 0.0);
  gm->flow_cells = mkarr(app->use_gpu, gm->num_mu * num_basis, size);
  gkyl_array_clear(gm->flow_cells, 0.0);
  gk_species_moment_init(
    app, gks, &gm->m1_mu, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_WEIGHT_NONE, false
  );

  gks->gyroaverage_moment = gk_species_flr_gyromatrix_gyroaverage_moment;
  gks->gyroaverage_apar = gk_species_flr_gyromatrix_gyroaverage_apar;
  gks->gyroaverage_apardot = gk_species_flr_gyromatrix_gyroaverage_apardot;
}

void
gk_species_flr_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  gkyl_array_release(gks->gyro_phi);
  gkyl_array_release(gks->m0_gyroavg);
  gkyl_array_release(gks->hamil_phi);
  gks->flr.release_func(app, gks);
}
