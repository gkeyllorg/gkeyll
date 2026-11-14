#include <assert.h>
#include <gkyl_vlasov_priv.h>

static void
vmfpo_moms_disabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_fpo_collisions *fpo,
  const struct gkyl_array *fin
)
{
  // Empty method.
}

// Compute the drag coefficient and diffusion tensor, including the corrections
// that make the discrete operator conserve momentum and energy.
static void
vmfpo_moms_enabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_fpo_collisions *fpo,
  const struct gkyl_array *fin
)
{
  struct timespec wst = gkyl_wall_clock();

  // Compute the LTE moments (n, V_drift, T/m) of the Maxwellian the Rosenbluth
  // potentials are computed from, and the M0, M1, M2 moments.
  vm_species_lte_moms(app, vms, &vms->lte, fin);
  vm_species_moment_calc(&fpo->moms, vms->local, app->local, fin);

  // Project the Rosenbluth potentials of a Maxwellian, and their derivatives at the
  // edges of velocity space.
  gkyl_proj_maxwellian_pots_on_basis_advance(
    fpo->pot_slvr, &vms->local, &app->local, vms->lte.moms.marr, fpo->h, fpo->g, fpo->h_surf,
    fpo->g_surf, fpo->dhdv_surf, fpo->dgdv_surf, fpo->d2gdv2_surf
  );

  // Compute the drag coefficient and diffusion tensor with recovery.
  gkyl_calc_fpo_drag_coeff_recovery(
    fpo->coeff_recovery, &vms->grid, vms->basis, &vms->local, &app->local, fpo->gamma, fpo->h,
    fpo->dhdv_surf, fpo->drag_coeff, fpo->drag_coeff_surf, app->use_gpu
  );
  gkyl_calc_fpo_diff_coeff_recovery(
    fpo->coeff_recovery, &vms->grid, vms->basis, &vms->local, &app->local, fpo->gamma, fpo->g,
    fpo->g_surf, fpo->dgdv_surf, fpo->d2gdv2_surf, fpo->diff_coeff, fpo->diff_coeff_surf,
    app->use_gpu
  );

  // Volume and boundary moments in the momentum and energy conservation constraints.
  if (app->use_gpu) {
    gkyl_mom_calc_advance_cu(fpo->fpo_mom_calc, &vms->local, &app->local, fin, fpo->fpo_moms);
  } else {
    gkyl_mom_calc_advance(fpo->fpo_mom_calc, &vms->local, &app->local, fin, fpo->fpo_moms);
  }
  gkyl_mom_calc_bcorr_advance(
    fpo->bcorr_calc, &vms->local, &app->local, fin, fpo->boundary_corrections
  );

  // Solve the linear system for the corrections and accumulate them onto the drag
  // coefficient and diffusion tensor.
  gkyl_fpo_coeff_correct_advance(
    fpo->coeff_correct_calc, &app->local, &vms->local, fpo->fpo_moms, fpo->boundary_corrections,
    fpo->moms.marr, fpo->drag_diff_coeff_corrs, fpo->drag_coeff, fpo->drag_coeff_surf,
    fpo->diff_coeff, fpo->diff_coeff_surf, app->use_gpu
  );

  // Sign of the corrected drag coefficient at cell interfaces (for upwinding).
  gkyl_calc_fpo_sgn_drag_coeff(
    fpo->coeff_recovery, vms->basis, &vms->local, fpo->drag_coeff_surf, fpo->sgn_drag_coeff_surf,
    fpo->const_sgn_drag_coeff_surf, app->use_gpu
  );

  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);
}

static void
vmfpo_rhs_disabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_fpo_collisions *fpo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  // Empty method.
}

static void
vmfpo_rhs_enabled(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_fpo_collisions *fpo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();

  // Accumulate update due to collisions onto rhs.
  gkyl_dg_updater_fpo_vlasov_advance(
    fpo->coll_slvr, fpo->offsets, &vms->local, fpo->drag_coeff, fpo->drag_coeff_surf,
    fpo->sgn_drag_coeff_surf, fpo->const_sgn_drag_coeff_surf, fpo->diff_coeff, fpo->diff_coeff_surf,
    fin, vms->cflrate, rhs
  );

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

static void
vmfpo_write_mom_disabled(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame)
{
  // Empty method.
}

static void
vmfpo_write_mom_enabled(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame)
{
  struct vm_fpo_collisions *fpo = &vms->fpo;

  // Recompute the potentials and coefficients from the current distribution function.
  fpo->moms_func(app, vms, fpo, vms->f);

  struct timespec wtm = gkyl_wall_clock();
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = vms->basis.id,
  });

  // Copy data from device to host before writing it out.
  if (app->use_gpu) {
    gkyl_array_copy(fpo->h_host, fpo->h);
    gkyl_array_copy(fpo->g_host, fpo->g);
    gkyl_array_copy(fpo->drag_coeff_host, fpo->drag_coeff);
    gkyl_array_copy(fpo->diff_coeff_host, fpo->diff_coeff);
  }

  // Write out the Rosenbluth potentials, drag coefficient and diffusion tensor.
  const char *fmts[] = {
    "%s-%s_fpo_H_%d.gkyl", "%s-%s_fpo_G_%d.gkyl", "%s-%s_fpo_drag_coeff_%d.gkyl",
    "%s-%s_fpo_diff_coeff_%d.gkyl"
  };
  const struct gkyl_array *arrs[] = {
    fpo->h_host, fpo->g_host, fpo->drag_coeff_host, fpo->diff_coeff_host
  };
  for (int i = 0; i < 4; ++i) {
    int sz = gkyl_calc_strlen(fmts[i], app->name, vms->name, frame);
    char fileNm[sz + 1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmts[i], app->name, vms->name, frame);
    gkyl_comm_array_write(vms->comm, &vms->grid, &vms->local, mt, arrs[i], fileNm);
  }
  app->stat.n_diag_io += 4;

  vlasov_array_meta_release(mt);
  app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
}

void
vm_species_fpo_init(
  struct gkyl_vlasov_app *app, struct vm_species *vms, struct vm_fpo_collisions *fpo
)
{
  fpo->collision_id = vms->info.collisions.collision_id;
  fpo->write_coll_diagnostics = vms->info.collisions.write_coll_diagnostics;

  // Empty methods.
  fpo->moms_func = vmfpo_moms_disabled;
  fpo->rhs_func = vmfpo_rhs_disabled;
  fpo->write_mom_func = vmfpo_write_mom_disabled;

  if (fpo->collision_id == GKYL_FPO_COLLISIONS) {
    int cdim = app->cdim, vdim = app->vdim;
    int pdim = cdim + vdim;
    int num_conf_basis = app->basis.num_basis;
    int num_phase_basis = vms->basis.num_basis;

    // The FPO is the collision operator of a non-relativistic species in three
    // velocity dimensions, and its kernels assume a uniform velocity grid (no
    // velocity-space Jacobian).
    if (vdim != 3 || vms->model_id != GKYL_MODEL_DEFAULT) {
      gkyl_exit("vm_species_fpo: the FPO needs a non-relativistic species with vdim=3.");
    }
    assert(vms->vel_map->is_identity);
    // Only self collisions with the Rosenbluth potentials of a Maxwellian.
    if (vms->info.collisions.num_cross_collisions > 0) {
      gkyl_exit("vm_species_fpo: the FPO does not support cross-species collisions.");
    }
    if (!vms->info.collisions.self_nu) {
      gkyl_exit("vm_species_fpo: the FPO needs self_nu (the prefactor Gamma of the operator).");
    }

    // Basis for the expansions on velocity-space surfaces.
    struct gkyl_basis surf_basis;
    if (vms->basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY && app->poly_order == 2) {
      gkyl_cart_modal_serendip(&surf_basis, pdim - 1, app->poly_order);
    } else {
      gkyl_exit("vm_species_fpo: the FPO only supports the p=2 Serendipity basis.");
    }
    int num_surf_basis = surf_basis.num_basis;
    // Number of surface quadrature nodes, for the sign of the drag coefficient.
    int num_surf_quad_nodes = pow(app->poly_order + 1, pdim - 1);

    // Project the prefactor of the operator, Gamma = (q^2/(eps0 m))^2 log(Lambda)/(4 pi),
    // which the user provides through self_nu.
    double nu_frac = vms->info.collisions.nu_frac ? vms->info.collisions.nu_frac : 1.0;
    fpo->gamma = mkarr(app->use_gpu, num_conf_basis, app->local_ext.volume);
    struct gkyl_array *gamma_ho = mkarr(false, num_conf_basis, app->local_ext.volume);
    vm_species_proj_conf_func(
      app, vms, vms->info.collisions.self_nu, vms->info.collisions.self_nu_ctx, gamma_ho
    );
    gkyl_array_copy(fpo->gamma, gamma_ho);
    gkyl_array_release(gamma_ho);
    gkyl_array_scale(fpo->gamma, nu_frac);

    // Rosenbluth potentials and their projection updater.
    fpo->h = mkarr(app->use_gpu, num_phase_basis, vms->local_ext.volume);
    fpo->g = mkarr(app->use_gpu, num_phase_basis, vms->local_ext.volume);

    // Potentials and their derivatives at the edges of velocity space, the boundary
    // conditions of the recovery (dG/dv at each edge in the two transverse directions).
    fpo->h_surf = mkarr(app->use_gpu, vdim * num_surf_basis, vms->local_ext.volume);
    fpo->g_surf = mkarr(app->use_gpu, vdim * num_surf_basis, vms->local_ext.volume);
    fpo->dhdv_surf = mkarr(app->use_gpu, vdim * num_surf_basis, vms->local_ext.volume);
    fpo->dgdv_surf = mkarr(app->use_gpu, 2 * vdim * num_surf_basis, vms->local_ext.volume);
    fpo->d2gdv2_surf = mkarr(app->use_gpu, vdim * num_surf_basis, vms->local_ext.volume);

    fpo->pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(
      &vms->grid, &app->local_ext, &vms->local_ext, &app->basis, &vms->basis, app->poly_order + 1,
      app->use_gpu
    );

    // Moment calculator to get M0, M1, M2 for the conservation corrections.
    vm_species_moment_init(app, vms, &fpo->moms, GKYL_F_MOMENT_M0M1M2, false);

    // Drag coefficient and diffusion tensor, their surface expansions at the lower
    // cell boundary, and the sign of the drag coefficient there.
    fpo->drag_coeff = mkarr(app->use_gpu, vdim * num_phase_basis, vms->local_ext.volume);
    fpo->drag_coeff_surf = mkarr(app->use_gpu, vdim * num_surf_basis, vms->local_ext.volume);
    fpo->sgn_drag_coeff_surf =
      mkarr(app->use_gpu, vdim * num_surf_quad_nodes, vms->local_ext.volume);
    fpo->const_sgn_drag_coeff_surf = mkarr(app->use_gpu, vdim, vms->local_ext.volume);

    fpo->diff_coeff = mkarr(app->use_gpu, vdim * vdim * num_phase_basis, vms->local_ext.volume);
    fpo->diff_coeff_surf =
      mkarr(app->use_gpu, 2 * vdim * vdim * num_surf_basis, vms->local_ext.volume);

    // Updater for the drag coefficient and diffusion tensor. It also sets the stencil
    // offsets it shares with the collision solver.
    fpo->coeff_recovery = gkyl_fpo_vlasov_coeff_recovery_new(
      &vms->grid, &vms->basis, &vms->local_ext, fpo->offsets, app->use_gpu
    );

    // Host-side copies for I/O.
    if (fpo->write_coll_diagnostics) {
      if (app->use_gpu) {
        fpo->h_host = mkarr(false, fpo->h->ncomp, fpo->h->size);
        fpo->g_host = mkarr(false, fpo->g->ncomp, fpo->g->size);
        fpo->drag_coeff_host = mkarr(false, fpo->drag_coeff->ncomp, fpo->drag_coeff->size);
        fpo->diff_coeff_host = mkarr(false, fpo->diff_coeff->ncomp, fpo->diff_coeff->size);
      } else {
        fpo->h_host = fpo->h;
        fpo->g_host = fpo->g;
        fpo->drag_coeff_host = fpo->drag_coeff;
        fpo->diff_coeff_host = fpo->diff_coeff;
      }
    }

    // Volume and boundary moments in the momentum and energy conservation constraints,
    // and the resulting corrections to the drag coefficient and diffusion tensor.
    fpo->fpo_moms = mkarr(app->use_gpu, (vdim + 1) * num_conf_basis, app->local_ext.volume);
    fpo->boundary_corrections =
      mkarr(app->use_gpu, 2 * (vdim + 1) * num_conf_basis, app->local_ext.volume);
    fpo->drag_diff_coeff_corrs =
      mkarr(app->use_gpu, (vdim + 1) * num_conf_basis, app->local_ext.volume);

    struct gkyl_mom_type *fpo_mom_type =
      gkyl_mom_fpo_vlasov_new(&app->basis, &vms->basis, &vms->local, app->use_gpu);
    struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
      .a = fpo->drag_coeff,
      .D = fpo->diff_coeff,
    };
    gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
    fpo->fpo_mom_calc = gkyl_mom_calc_new(&vms->grid, fpo_mom_type, app->use_gpu);
    // The calculator holds its own reference to the moment type.
    gkyl_mom_type_release(fpo_mom_type);

    double v_bounds[2 * GKYL_MAX_DIM] = {0.0};
    for (int d = 0; d < vdim; ++d) {
      v_bounds[d] = vms->info.lower[d];
      v_bounds[d + vdim] = vms->info.upper[d];
    }
    fpo->bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(
      &vms->grid, &app->basis, &vms->basis, &vms->local, v_bounds, fpo->diff_coeff, app->use_gpu
    );

    fpo->coeff_correct_calc =
      gkyl_fpo_coeff_correct_new(&vms->grid, &app->basis, &app->local, app->use_gpu);

    // FPO updater.
    fpo->coll_slvr =
      gkyl_dg_updater_fpo_vlasov_new(&vms->grid, &vms->basis, &vms->local, app->use_gpu);

    // Methods chosen at runtime.
    fpo->moms_func = vmfpo_moms_enabled;
    fpo->rhs_func = vmfpo_rhs_enabled;
    if (fpo->write_coll_diagnostics) {
      fpo->write_mom_func = vmfpo_write_mom_enabled;
    }
  }
}

void
vm_species_fpo_moms(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_fpo_collisions *fpo,
  const struct gkyl_array *fin
)
{
  fpo->moms_func(app, vms, fpo, fin);
}

void
vm_species_fpo_rhs(
  gkyl_vlasov_app *app, const struct vm_species *vms, struct vm_fpo_collisions *fpo,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  fpo->rhs_func(app, vms, fpo, fin, rhs);
}

void
vm_species_fpo_write_mom(gkyl_vlasov_app *app, struct vm_species *vms, double tm, int frame)
{
  vms->fpo.write_mom_func(app, vms, tm, frame);
}

void
vm_species_fpo_release(const struct gkyl_vlasov_app *app, const struct vm_fpo_collisions *fpo)
{
  if (fpo->collision_id == GKYL_FPO_COLLISIONS) {
    gkyl_dg_updater_fpo_vlasov_release(fpo->coll_slvr);

    gkyl_fpo_vlasov_coeff_correct_release(fpo->coeff_correct_calc);
    gkyl_mom_calc_bcorr_release(fpo->bcorr_calc);
    gkyl_mom_calc_release(fpo->fpo_mom_calc);

    gkyl_array_release(fpo->drag_diff_coeff_corrs);
    gkyl_array_release(fpo->boundary_corrections);
    gkyl_array_release(fpo->fpo_moms);

    if (fpo->write_coll_diagnostics) {
      if (app->use_gpu) {
        gkyl_array_release(fpo->h_host);
        gkyl_array_release(fpo->g_host);
        gkyl_array_release(fpo->drag_coeff_host);
        gkyl_array_release(fpo->diff_coeff_host);
      }
    }

    gkyl_fpo_vlasov_coeff_recovery_release(fpo->coeff_recovery);

    gkyl_array_release(fpo->diff_coeff_surf);
    gkyl_array_release(fpo->diff_coeff);
    gkyl_array_release(fpo->const_sgn_drag_coeff_surf);
    gkyl_array_release(fpo->sgn_drag_coeff_surf);
    gkyl_array_release(fpo->drag_coeff_surf);
    gkyl_array_release(fpo->drag_coeff);

    vm_species_moment_release(app, &fpo->moms);

    gkyl_proj_maxwellian_pots_on_basis_release(fpo->pot_slvr);

    gkyl_array_release(fpo->d2gdv2_surf);
    gkyl_array_release(fpo->dgdv_surf);
    gkyl_array_release(fpo->dhdv_surf);
    gkyl_array_release(fpo->g_surf);
    gkyl_array_release(fpo->h_surf);

    gkyl_array_release(fpo->g);
    gkyl_array_release(fpo->h);
    gkyl_array_release(fpo->gamma);
  }
}
