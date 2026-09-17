#include <assert.h>
#include <gkyl_vlasov_priv.h>

#include <gkyl_array_rio.h>

void
vm_species_fpo_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_fpo_collisions *fpo)
{
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;
  struct gkyl_basis surf_basis;

  fpo->write_diagnostics = s->info.collisions.write_diagnostics;

  // initialize surface basis for potentials on velocity space edges
  if (app->poly_order == 1) {
    // JMR 1/25/26: Polynomial order 1 not supported for the time being.
    // The p=1 kernels in vlasov/ker/fpo use pure p=1 Serendipity, not hybrid
    gkyl_exit("p=1 Serendipity not supported for the FPO! Use p=2 instead.");
  }
  else {
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, app->poly_order);
  }

  // allocate gamma and initialize it
  fpo->gamma = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *gamma_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);

  int num_quad = app->poly_order == 1 ? app->poly_order+2 : app->poly_order+1;
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    num_quad, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, gamma_host);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(fpo->gamma, gamma_host);
  gkyl_array_release(gamma_host);

  // initialize the potentials and solver for potentials
  fpo->h = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  fpo->g = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
 
  fpo->h_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->g_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->dhdv_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->dgdv_surf = mkarr(app->use_gpu, 2*vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->d2gdv2_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);

  fpo->pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(&s->grid,
    &app->local_ext, &s->local_ext, &app->confBasis, &app->basis, app->poly_order+1, app->use_gpu);

  // allocate moments needed for FPO update
  vm_species_moment_init(app, s, &fpo->lte_moms, GKYL_F_MOMENT_LTE, false);
  vm_species_moment_init(app, s, &fpo->moms, GKYL_F_MOMENT_M0M1M2, false);

  // need number of surface quad points for drag coeff sign information arrays
  int num_surf_quad_nodes;
  if (app->poly_order == 1)
    num_surf_quad_nodes = pow(app->poly_order+1, cdim)*pow(app->poly_order+2, vdim-1);
  else
    num_surf_quad_nodes = pow(app->poly_order+1, pdim-1);

  // initialize drag and diffusion coefficient arrays
  fpo->drag_coeff = mkarr(app->use_gpu, vdim*app->basis.num_basis, s->local_ext.volume);
  fpo->drag_coeff_surf = mkarr(app->use_gpu, vdim*surf_basis.num_basis, s->local_ext.volume);
  fpo->sgn_drag_coeff_surf = mkarr(app->use_gpu, vdim*num_surf_quad_nodes, s->local_ext.volume);
  fpo->const_sgn_drag_coeff_surf = mkarr(app->use_gpu, vdim, s->local_ext.volume);

  fpo->diff_coeff = mkarr(app->use_gpu, vdim*vdim*app->basis.num_basis, s->local_ext.volume);
  fpo->diff_coeff_surf = mkarr(app->use_gpu, 2*vdim*vdim*surf_basis.num_basis, s->local_ext.volume);

  // utility struct for computing drag and diffusion coefficients
  fpo->coeff_recovery = gkyl_fpo_vlasov_coeff_recovery_new(&s->grid, &app->basis, &s->local_ext, fpo->offsets, app->use_gpu);

  // host-side arrays for writing out H, G, a, D
  fpo->h_host = fpo->h;
  fpo->g_host = fpo->g;
  fpo->drag_coeff_host = fpo->drag_coeff;
  fpo->diff_coeff_host = fpo->diff_coeff;
  if (app->use_gpu) {
    fpo->h_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    fpo->g_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    fpo->drag_coeff_host = mkarr(false, vdim*app->basis.num_basis, s->local_ext.volume);
    fpo->diff_coeff_host = mkarr(false, vdim*vdim*app->basis.num_basis, s->local_ext.volume);
  }

  // velocity space boundary corrections for momentum and energy
  fpo->fpo_moms = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  fpo->boundary_corrections = mkarr(app->use_gpu, 2*(vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
  fpo->drag_diff_coeff_corrs = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);

  // updater to compute moments for conservation corrections
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&app->confBasis,
    &app->basis, &s->local, app->use_gpu);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = { .a = fpo->drag_coeff, .D = fpo->diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  fpo->fpo_mom_calc = gkyl_mom_calc_new(&s->grid, fpo_mom_type, app->use_gpu);

  // updater to compute boundary corrections for conservation corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = s->info.lower[d];
    v_bounds[d + vdim] = s->info.upper[d];
  }
  fpo->bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&s->grid,
    &app->confBasis, &app->basis, &s->local, v_bounds, fpo->diff_coeff, app->use_gpu);

  // updater object for drag and diffusion coefficient correction
  fpo->coeff_correct_calc = gkyl_fpo_coeff_correct_new(&s->grid, &app->confBasis, &app->local, app->use_gpu);

  // initialize FPO updater
  fpo->coll_slvr = gkyl_dg_updater_fpo_vlasov_new(&s->grid, &app->basis, &s->local, app->use_gpu);
}

// computes drag coefficient and diffusion tensor
void
vm_species_fpo_drag_diff_coeffs(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_fpo_collisions *fpo, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  // calculate needed moments
  vm_species_moment_calc(&fpo->lte_moms, s->local, app->local, fin);
  vm_species_moment_calc(&fpo->moms, s->local, app->local, fin);

  // calculate maxwellian potentials
  gkyl_proj_maxwellian_pots_on_basis_advance(fpo->pot_slvr, &s->local, &app->local,
    fpo->lte_moms.marr, fpo->h, fpo->g, fpo->h_surf, fpo->g_surf,
    fpo->dhdv_surf, fpo->dgdv_surf, fpo->d2gdv2_surf);

  // calculate drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(fpo->coeff_recovery, &s->grid, app->basis,
    &s->local, &app->local, fpo->gamma,
    fpo->h, fpo->dhdv_surf, fpo->drag_coeff, fpo->drag_coeff_surf, app->use_gpu);

  gkyl_calc_fpo_diff_coeff_recovery(fpo->coeff_recovery, &s->grid, app->basis,
    &s->local, &app->local, fpo->gamma,
    fpo->g, fpo->g_surf, fpo->dgdv_surf, fpo->d2gdv2_surf,
    fpo->diff_coeff, fpo->diff_coeff_surf, app->use_gpu);

  // Calculate corrections for momentum and energy conservation
  if (app->use_gpu) {
    // calculate volume corrections
    gkyl_mom_calc_advance_cu(fpo->fpo_mom_calc, &s->local, &app->local, fin, fpo->fpo_moms);
  }
  else {
    // calculate volume corrections
    gkyl_mom_calc_advance(fpo->fpo_mom_calc, &s->local, &app->local, fin, fpo->fpo_moms);
  }
  // calculate boundary corrections
  gkyl_mom_calc_bcorr_advance(fpo->bcorr_calc,
    &s->local, &app->local, fin, fpo->boundary_corrections);

  // solve linear system for corrections and accumulate onto drag/diff coefficients
  gkyl_fpo_coeff_correct_advance(fpo->coeff_correct_calc,
    &app->local, &s->local, fpo->fpo_moms, fpo->boundary_corrections,
    fpo->moms.marr, fpo->drag_diff_coeff_corrs,
    fpo->drag_coeff, fpo->drag_coeff_surf,
    fpo->diff_coeff, fpo->diff_coeff_surf, app->use_gpu);

  // compute sign information at interfaces with corrected drag coefficient
  gkyl_calc_fpo_sgn_drag_coeff(fpo->coeff_recovery, app->basis, &s->local,
    fpo->drag_coeff_surf, fpo->sgn_drag_coeff_surf,
    fpo->const_sgn_drag_coeff_surf, app->use_gpu);

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

// updates the collision terms in the rhs
void
vm_species_fpo_rhs(gkyl_vlasov_app *app, const struct vm_species *s,
  struct vm_fpo_collisions *fpo, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();

  wst = gkyl_wall_clock();
 
  // accumulate update due to collisions onto rhs
  gkyl_dg_updater_fpo_vlasov_advance(fpo->coll_slvr,
    fpo->offsets, &s->local,
    fpo->drag_coeff, fpo->drag_coeff_surf,
    fpo->sgn_drag_coeff_surf, fpo->const_sgn_drag_coeff_surf,
    fpo->diff_coeff, fpo->diff_coeff_surf, fin, s->cflrate, rhs);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_species_fpo_release(const struct gkyl_vlasov_app *app, const struct vm_fpo_collisions *fpo)
{
  gkyl_array_release(fpo->gamma);
  gkyl_array_release(fpo->h);
  gkyl_array_release(fpo->g);
  if (app->use_gpu) {
    gkyl_array_release(fpo->h_host);
    gkyl_array_release(fpo->g_host);
    gkyl_array_release(fpo->drag_coeff_host);
    gkyl_array_release(fpo->diff_coeff_host);
  }

  gkyl_array_release(fpo->h_surf);
  gkyl_array_release(fpo->g_surf);
  gkyl_array_release(fpo->dhdv_surf);
  gkyl_array_release(fpo->dgdv_surf);
  gkyl_array_release(fpo->d2gdv2_surf);

  gkyl_proj_maxwellian_pots_on_basis_release(fpo->pot_slvr);

  vm_species_moment_release(app, &fpo->lte_moms);
  vm_species_moment_release(app, &fpo->moms);

  gkyl_array_release(fpo->drag_coeff);
  gkyl_array_release(fpo->diff_coeff);
  gkyl_array_release(fpo->drag_coeff_surf);
  gkyl_array_release(fpo->diff_coeff_surf);
  gkyl_array_release(fpo->sgn_drag_coeff_surf);
  gkyl_array_release(fpo->const_sgn_drag_coeff_surf);

  gkyl_array_release(fpo->fpo_moms);
  gkyl_array_release(fpo->boundary_corrections);
  gkyl_array_release(fpo->drag_diff_coeff_corrs);

  gkyl_mom_calc_release(fpo->fpo_mom_calc);
  gkyl_mom_calc_bcorr_release(fpo->bcorr_calc);

  gkyl_fpo_vlasov_coeff_correct_release(fpo->coeff_correct_calc);
  gkyl_fpo_vlasov_coeff_recovery_release(fpo->coeff_recovery);

  gkyl_dg_updater_fpo_vlasov_release(fpo->coll_slvr);
}

