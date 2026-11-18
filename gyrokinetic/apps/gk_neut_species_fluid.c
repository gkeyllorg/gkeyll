#include <gkyl_gk_neut_species_priv.h>
#include <gk_neut_species_fluid_cu.h>

// True if this rank owns the physical boundary at (dir, edge). Communication
// has already populated ghosts shared with neighboring ranks.
static bool
gk_neut_species_fluid_owns_edge(
  const struct gk_neut_species *species, int dir, enum gkyl_edge_loc edge
)
{
  return edge == GKYL_LOWER_EDGE ? species->local.lower[dir] == species->global.lower[dir] :
                                   species->local.upper[dir] == species->global.upper[dir];
}

static void
gk_neut_species_fluid_diffusion_dirichlet_ghost(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, int dir,
  enum gkyl_edge_loc edge, const struct gkyl_range *ghost, double number_density,
  struct gkyl_array *density, struct gkyl_array *tensor
)
{
  if (!gk_neut_species_fluid_owns_edge(species, dir, edge)) {
    return;
  }

  const int skin_idx_dir = edge == GKYL_LOWER_EDGE ? species->local.lower[dir] :
                                                     species->local.upper[dir];
#ifdef GKYL_HAVE_CUDA
  if (app->use_gpu) {
    gk_neut_fluid_diffusion_dirichlet_ghost_cu(
      app->cdim, dir, skin_idx_dir, species->info.mass, number_density, &species->local_ext, ghost,
      density, tensor
    );
    return;
  }
#endif
  const int nb = app->basis.num_basis;
  const int nK = app->cdim * app->cdim;
  const double rho_b_modal = pow(sqrt(2.0), app->cdim) * species->info.mass * number_density;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost);
  while (gkyl_range_iter_next(&iter)) {
    int skin_idx[GKYL_MAX_CDIM];
    for (int d = 0; d < app->cdim; ++d) {
      skin_idx[d] = iter.idx[d];
    }
    skin_idx[dir] = skin_idx_dir;

    const double *rho_skin =
      gkyl_array_cfetch(density, gkyl_range_idx(&species->local_ext, skin_idx));
    double *rho_ghost = gkyl_array_fetch(density, gkyl_range_idx(&species->local_ext, iter.idx));

    // Mirror the skin polynomial and choose its exterior value so the
    // arithmetic trace at the physical face is rho_b:
    //   rho_g = 2 rho_b - R(rho_s).
    app->basis.flip_odd_sign(dir, rho_skin, rho_ghost);
    for (int k = 0; k < nb; ++k) {
      rho_ghost[k] = -rho_ghost[k];
    }
    rho_ghost[0] += 2.0 * rho_b_modal;

    const double *K_skin = gkyl_array_cfetch(tensor, gkyl_range_idx(&species->local_ext, skin_idx));
    double *K_ghost = gkyl_array_fetch(tensor, gkyl_range_idx(&species->local_ext, iter.idx));
    for (int k = 0; k < nK; ++k) {
      app->basis.flip_odd_sign(dir, &K_skin[k * nb], &K_ghost[k * nb]);
    }
  }
}

// Impose zero diffusive flux through the physical face by mirroring the skin
// cell into the ghost cell:
//   rho_g = R(rho_s),  K_g = -R(K_s),
// where R reflects across the face. With these ghosts the interior surface
// kernel reduces exactly to the one-sided zero-flux boundary kernel, including
// the off-diagonal (mapped-geometry) tensor terms; see the unit test
// conf_diffusion_zero_flux_boundary_kernels. Unlike the hyper_dg zero-flux
// flag, this acts only on the given ghost range (e.g. the SOL part of a
// parallel boundary whose confined part is periodic).
static void
gk_neut_species_fluid_diffusion_zero_flux_ghost(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, int dir,
  enum gkyl_edge_loc edge, const struct gkyl_range *ghost, struct gkyl_array *density,
  struct gkyl_array *tensor
)
{
  if (!gk_neut_species_fluid_owns_edge(species, dir, edge)) {
    return;
  }

  const int skin_idx_dir = edge == GKYL_LOWER_EDGE ? species->local.lower[dir] :
                                                     species->local.upper[dir];
#ifdef GKYL_HAVE_CUDA
  if (app->use_gpu) {
    gk_neut_fluid_diffusion_zero_flux_ghost_cu(
      app->cdim, dir, skin_idx_dir, &species->local_ext, ghost, density, tensor
    );
    return;
  }
#endif
  const int nb = app->basis.num_basis;
  const int nK = app->cdim * app->cdim;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost);
  while (gkyl_range_iter_next(&iter)) {
    int skin_idx[GKYL_MAX_CDIM];
    for (int d = 0; d < app->cdim; ++d) {
      skin_idx[d] = iter.idx[d];
    }
    skin_idx[dir] = skin_idx_dir;
    long skin_loc = gkyl_range_idx(&species->local_ext, skin_idx);
    long ghost_loc = gkyl_range_idx(&species->local_ext, iter.idx);

    app->basis.flip_odd_sign(
      dir, gkyl_array_cfetch(density, skin_loc), gkyl_array_fetch(density, ghost_loc)
    );
    const double *K_skin = gkyl_array_cfetch(tensor, skin_loc);
    double *K_ghost = gkyl_array_fetch(tensor, ghost_loc);
    for (int k = 0; k < nK; ++k) {
      app->basis.flip_even_sign(dir, &K_skin[k * nb], &K_ghost[k * nb]);
    }
  }
}

// Recycled neutral mass flux at a target, stored in the ghost cells in the
// same representation as the ion boundary-flux moments:
//   F_n = m_n R sum_s g_s,
// where g_s is the J-weighted M0 rate of impacting species s deposited in the
// ghost cell (its ghost-cell integral is the particle loss rate through the
// face). Within an explicit RK stage g_s is the rate of the current stage
// (bflux.current_moms). In the operator-split implicit stage it must be the
// RK-weighted rate over the whole step (bflux.f), i.e. what the species lost
// during the step; current_moms then points to an RK work array that the
// stage combination has since rescaled.
static void
gk_neut_species_fluid_recycling_flux(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, int dir,
  enum gkyl_edge_loc edge, const struct gkyl_range *ghost, bool step_averaged
)
{
  const int e = edge == GKYL_LOWER_EDGE ? 0 : 1;
  const struct gkyl_gyrokinetic_bc *bc = edge == GKYL_LOWER_EDGE ? &species->lower_bc[dir] :
                                                                   &species->upper_bc[dir];
  gkyl_array_clear_range(species->diffusion_recycling_flux, 0.0, ghost);
  for (int i = 0; i < species->diffusion_recycling_num_species[e]; ++i) {
    struct gk_species *ion = &app->species[species->diffusion_recycling_ion_idx[e][i]];
    struct gkyl_array **ion_flux_moms = step_averaged ? ion->bflux.f : ion->bflux.current_moms;
    assert(ion_flux_moms);
    gkyl_array_clear_range(species->diffusion_recycling_flux_tmp, 0.0, ghost);
    gk_species_bflux_get_flux_mom(
      &ion->bflux, dir, edge, GKYL_F_MOMENT_M0, ion_flux_moms,
      species->diffusion_recycling_flux_tmp, ghost
    );
    gkyl_array_accumulate_range(
      species->diffusion_recycling_flux, species->info.mass * bc->emission.recycling_frac,
      species->diffusion_recycling_flux_tmp, ghost
    );
  }
}

// Add the recycled flux entering through the physical face to the density
// component (the first num_basis coefficients) of rhs in the skin cells. The
// ghost-cell projection of a face flux and the skin-cell projection of the same
// face flux are mirror images, so the skin increment is R(F_n); it has the same
// cell integral, so the neutrals gain exactly R times the ion loss rate.
static void
gk_neut_species_fluid_recycling_source(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, int dir,
  enum gkyl_edge_loc edge, const struct gkyl_range *ghost, struct gkyl_array *rhs
)
{
  if (!gk_neut_species_fluid_owns_edge(species, dir, edge)) {
    return;
  }

  const int skin_idx_dir = edge == GKYL_LOWER_EDGE ? species->local.lower[dir] :
                                                     species->local.upper[dir];
#ifdef GKYL_HAVE_CUDA
  if (app->use_gpu) {
    gk_neut_fluid_recycling_source_cu(
      app->cdim, dir, skin_idx_dir, &species->local_ext, ghost, species->diffusion_recycling_flux,
      rhs
    );
    return;
  }
#endif
  const int nb = app->basis.num_basis;
  double src[8] = {0.0};

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost);
  while (gkyl_range_iter_next(&iter)) {
    int skin_idx[GKYL_MAX_CDIM];
    for (int d = 0; d < app->cdim; ++d) {
      skin_idx[d] = iter.idx[d];
    }
    skin_idx[dir] = skin_idx_dir;

    app->basis.flip_odd_sign(
      dir,
      gkyl_array_cfetch(
        species->diffusion_recycling_flux, gkyl_range_idx(&species->local_ext, iter.idx)
      ),
      src
    );
    double *rhs_skin = gkyl_array_fetch(rhs, gkyl_range_idx(&species->local_ext, skin_idx));
    for (int k = 0; k < nb; ++k) {
      rhs_skin[k] += src[k];
    }
  }
}

// Ghost ranges of the parallel boundaries. In an LCFS-spanning domain only the
// SOL part of each parallel boundary is a material target.
static void
gk_neut_species_fluid_par_ghost_ranges(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  const struct gkyl_range *ghost[2]
)
{
  const int par_dir = app->cdim - 1;
  ghost[0] =
    app->gk_geom->has_LCFS ? &app->local_lower_ghost_par_sol : &species->local_lower_ghost[par_dir];
  ghost[1] =
    app->gk_geom->has_LCFS ? &app->local_upper_ghost_par_sol : &species->local_upper_ghost[par_dir];
}

// Fill the parallel ghosts of the diffusion density and tensor. Recycling
// targets get a zero diffusive flux (the recycled flux is added separately as
// a boundary source); fixed-density targets get a Dirichlet trace, which is
// zero for the homogeneous part of the implicit operator.
static void
gk_neut_species_fluid_diffusion_par_ghosts(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gkyl_array *density, bool homogeneous_bc
)
{
  const int par_dir = app->cdim - 1;
  const struct gkyl_range *ghost[2];
  gk_neut_species_fluid_par_ghost_ranges(app, species, ghost);
  for (int e = 0; e < 2; ++e) {
    const enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    const struct gkyl_gyrokinetic_bc *bc = e == 0 ? &species->lower_bc[par_dir] :
                                                    &species->upper_bc[par_dir];
    const enum gkyl_gyrokinetic_neut_fluid_diffusion_bc_type diff_bc =
      e == 0 ? species->info.diffusion.lower_bc_type : species->info.diffusion.upper_bc_type;
    if (bc->type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_fluid_diffusion_zero_flux_ghost(
        app, species, par_dir, edge, ghost[e], density, species->diffusion_tensor
      );
    } else if (diff_bc == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET) {
      const double n_b =
        e == 0 ? species->info.diffusion.lower_bc_density : species->info.diffusion.upper_bc_density;
      gk_neut_species_fluid_diffusion_dirichlet_ghost(
        app, species, par_dir, edge, ghost[e], homogeneous_bc ? 0.0 : n_b, density,
        species->diffusion_tensor
      );
    }
  }
}

// Add the recycled neutral flux at the parallel targets to rhs. step_averaged
// selects the step-averaged impacting flux (implicit stage) instead of the
// current RK-stage flux (explicit RHS).
static void
gk_neut_species_fluid_diffusion_add_recycling(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *rhs,
  bool step_averaged
)
{
  const int par_dir = app->cdim - 1;
  const struct gkyl_range *ghost[2];
  gk_neut_species_fluid_par_ghost_ranges(app, species, ghost);
  for (int e = 0; e < 2; ++e) {
    const enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    const struct gkyl_gyrokinetic_bc *bc = e == 0 ? &species->lower_bc[par_dir] :
                                                    &species->upper_bc[par_dir];
    if (bc->type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_fluid_recycling_flux(app, species, par_dir, edge, ghost[e], step_averaged);
      gk_neut_species_fluid_recycling_source(app, species, par_dir, edge, ghost[e], rhs);
    }
  }
}

static void
gk_neut_species_fluid_apply_bc_dynamic(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *f
)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_comm_array_per_sync(
    species->comm, &species->local, &species->local_ext, app->num_periodic_dir, app->periodic_dirs,
    f
  );
  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// Squared neutral thermal speed T_n/m_n in the numerator of
// D = vt^2/(ne*<sigma v>_cx), from the neutral LTE moments
// (n, u_e1, u_e2, u_bhat, T/m) computed for the current stage.
static const struct gkyl_array *
fluid_diffusion_vt_sq(gkyl_gyrokinetic_app *app, struct gk_neut_species *species)
{
  gkyl_array_set_offset(
    species->diffusion_vt_sq, 1.0, species->lte.moms.marr, 4 * app->basis.num_basis
  );
  return species->diffusion_vt_sq;
}

static double
gk_neut_species_fluid_rhs_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(species->cflrate, 0.0);
  if (species->has_diffusion) {
    gkyl_array_clear(species->diffusion_cflrate, 0.0);
  }
  gkyl_array_clear(rhs, 0.0);

  // Collisionless terms.
  struct timespec wst = gkyl_wall_clock();
  if (species->has_diffusion && !species->implicit_diffusion) {
    if (species->use_reaction_rate_diffusion) {
      struct gk_react *react = &species->react_neut;
      int iz_idx = species->ionization_react_idx;
      int cx_idx = species->diffusion_cx_react_idx;
      const struct gkyl_array *vt_sq_D = fluid_diffusion_vt_sq(app, species);
      struct gk_species *elc = &app->species[react->elc_idx[iz_idx]];

      // Evaluate <sigma v>_cx for the diffusion closure independently of the
      // neutral-density guard used by the reaction source.
      gkyl_dg_cx_coll_rate(
        react->cx[cx_idx], app->species[react->ion_idx[cx_idx]].lte.moms.marr,
        species->lte.moms.marr, react->upar_ion[cx_idx], species->diffusion_coeff
      );

      // denominator = ne*<sigma v>_cx.
      gkyl_dg_mul_op_range(
        &app->basis, 0, species->diffusion_moment_ratio, 0, elc->lte.moms.marr, 0,
        species->diffusion_coeff, &species->local
      );

      // Construct D at the p1 interpolation nodes. Flooring the collision
      // frequency at every node (rather than only its cell average) prevents
      // a positive average with a negative/near-zero within-cell value from
      // producing an unbounded diffusion coefficient. A multilinear p1 field
      // that is positive at every vertex is positive throughout the cell.
      const double nu_min = species->info.diffusion.min_collision_frequency;
      if (nu_min > 0.0) {
#ifdef GKYL_HAVE_CUDA
        if (app->use_gpu) {
          gk_neut_fluid_diffusion_floor_cu(
            app->cdim, nu_min, &species->local, vt_sq_D,
            species->diffusion_moment_ratio, species->diffusion_coeff
          );
        } else
#endif
        {
          double nodes[8 * GKYL_MAX_CDIM] = {0.0};
          double D_nodal[8] = {0.0};
          app->basis.node_list(nodes);
          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, &species->local);
          while (gkyl_range_iter_next(&iter)) {
            long loc = gkyl_range_idx(&species->local, iter.idx);
            const double *vt_sq = gkyl_array_cfetch(vt_sq_D, loc);
            const double *nu = gkyl_array_cfetch(species->diffusion_moment_ratio, loc);
            double *D = gkyl_array_fetch(species->diffusion_coeff, loc);
            for (int q = 0; q < app->basis.num_basis; ++q) {
              double vt_sq_q = app->basis.eval_expand(&nodes[q * app->cdim], vt_sq);
              double nu_q = app->basis.eval_expand(&nodes[q * app->cdim], nu);
              D_nodal[q] = fmax(vt_sq_q, 0.0) / fmax(nu_q, nu_min);
            }
            app->basis.nodal_to_modal(D_nodal, D);
          }
        }
      } else {
        // With no floor requested, retain the original weak quotient.
        gkyl_dg_div_op_range(
          species->diffusion_div_mem, &app->basis, 0, species->diffusion_coeff, 0,
          vt_sq_D, 0, species->diffusion_moment_ratio, &species->local
        );
      }

      // The variable-coefficient surface terms read neighboring cells.
      // Populate coefficient ghosts independently of the state BC update.
      gkyl_comm_array_per_sync(
        species->comm, &species->local, &species->local_ext, species->num_periodic_dir,
        species->periodic_dirs, species->diffusion_coeff
      );
      gkyl_comm_array_sync(
        species->comm, &species->local, &species->local_ext, species->diffusion_coeff
      );
    }

    // The evolved density is N=J*rho. Form the physical density rho=N/J and
    // the mapped tensor K^{ij}=J*D*g_neut^{ij} used by
    //   d_t N = d_i(K^{ij} d_j rho).
    // Use weak division by J instead of multiplication by its separately
    // projected inverse. This preserves a spatially constant physical density
    // exactly when the stored state is N=J*rho.
    gkyl_dg_div_op_range(
      species->diffusion_div_mem, &app->basis, 0, species->diffusion_density, 0, fin, 0,
      app->gk_geom->geo_int.jacobgeo, &species->local
    );
    gkyl_dg_mul_op_range(
      &app->basis, 0, species->diffusion_geom_factor, 0, species->diffusion_coeff, 0,
      app->gk_geom->geo_int.jacobgeo, &species->local
    );

    const int metric_idx_1x[1][1] = {{5}};
    const int metric_idx_2x[2][2] = {{0, 2}, {2, 5}};
    const int metric_idx_3x[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
    for (int i = 0; i < app->cdim; ++i) {
      for (int j = 0; j < app->cdim; ++j) {
        int gij_idx = app->cdim == 1 ? metric_idx_1x[i][j] :
                      app->cdim == 2 ? metric_idx_2x[i][j] :
                                       metric_idx_3x[i][j];
        gkyl_dg_mul_op_range(
          &app->basis, i * app->cdim + j, species->diffusion_tensor, 0,
          species->diffusion_geom_factor, gij_idx, app->gk_geom->geo_int.gij_neut, &species->local
        );
      }
    }
    gkyl_comm_array_per_sync(
      species->comm, &species->local, &species->local_ext, species->num_periodic_dir,
      species->periodic_dirs, species->diffusion_density
    );
    gkyl_comm_array_sync(
      species->comm, &species->local, &species->local_ext, species->diffusion_density
    );
    gkyl_comm_array_per_sync(
      species->comm, &species->local, &species->local_ext, species->num_periodic_dir,
      species->periodic_dirs, species->diffusion_tensor
    );
    gkyl_comm_array_sync(
      species->comm, &species->local, &species->local_ext, species->diffusion_tensor
    );

    // Fixed-density or recycling open-field-line boundaries.
    gk_neut_species_fluid_diffusion_par_ghosts(app, species, species->diffusion_density, false);

    gkyl_dg_updater_conf_diffusion_advance(
      species->diffusion_slvr, &species->local, species->diffusion_density,
      species->diffusion_cflrate, rhs
    );
    gkyl_array_accumulate_range(species->cflrate, 1.0, species->diffusion_cflrate, &species->local);
    gk_neut_species_fluid_diffusion_add_recycling(app, species, rhs, false);

    // Density is the only independently diffused quantity. Propagate its
    // change to the dependent conserved moments while holding their specific
    // values (moment/rho) fixed.
    for (int m = 1; m < species->num_moments; ++m) {
      gkyl_dg_div_op_range(
        species->diffusion_div_mem, &app->basis, m, species->diffusion_moment_ratio, m, fin, 0, fin,
        &species->local
      );
      gkyl_dg_mul_op_range(
        &app->basis, m, rhs, m, species->diffusion_moment_ratio, 0, rhs, &species->local
      );
    }
  }
  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);

  // Ionization loss and other enabled reactions with charged species.
  gk_neut_species_react_rhs(app, species, &species->react_neut, fin, rhs);

  // Add the explicit ionization frequency ne*<sigma v>_iz to the CFL rate.
  if (species->ionization_react_idx >= 0) {
    struct gk_react *react = &species->react_neut;
    for (int i = 0; i < react->num_react; ++i) {
      if (react->react_id[i] != GKYL_REACT_IZ) {
        continue;
      }
      struct gk_species *elc = &app->species[react->elc_idx[i]];
      gkyl_dg_mul_op_range(
        &app->basis, 0, species->ionization_rate, 0, elc->lte.moms.marr, 0, react->coeff_react[i],
        &species->local
      );
      gkyl_dg_calc_average_range(
        &app->basis, 0, species->ionization_cflrate, 0, species->ionization_rate, species->local
      );
      gkyl_array_accumulate_range(
        species->cflrate, 1.0, species->ionization_cflrate, &species->local
      );
    }
  }

  // Compute volume-integrated reactions in sca.
  gk_neut_species_scaling_rhs(app, species, &species->sca, fin, rhs);

  app->stat.n_neut_species_omega_cfl += 1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }

  omega_cfl = omega_cfl_ho[0];

  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl / omega_cfl;
}

static double
fluid_diffusion_dot(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, const struct gkyl_array *a,
  const struct gkyl_array *b
)
{
#ifdef GKYL_HAVE_CUDA
  if (app->use_gpu) {
    double local_sum =
      gk_neut_fluid_diffusion_dot_cu(&species->local, a, b, species->diffusion_implicit_dot);
    double global_sum = 0.0;
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &local_sum, &global_sum);
    return global_sum;
  }
#endif
  double local_sum = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &species->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&species->local, iter.idx);
    const double *ap = gkyl_array_cfetch(a, loc);
    const double *bp = gkyl_array_cfetch(b, loc);
    for (int k = 0; k < app->basis.num_basis; ++k) {
      local_sum += ap[k] * bp[k];
    }
  }
  double global_sum = 0.0;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, &local_sum, &global_sum);
  return global_sum;
}

static void
fluid_diffusion_vec_set(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *species, struct gkyl_array *out,
  double a, const struct gkyl_array *x, double b, const struct gkyl_array *y
)
{
#ifdef GKYL_HAVE_CUDA
  if (app->use_gpu) {
    gk_neut_fluid_diffusion_vec_set_cu(&species->local, out, a, x, b, y);
    return;
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &species->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&species->local, iter.idx);
    double *op = gkyl_array_fetch(out, loc);
    const double *xp = gkyl_array_cfetch(x, loc);
    const double *yp = y ? gkyl_array_cfetch(y, loc) : 0;
    for (int k = 0; k < app->basis.num_basis; ++k) {
      op[k] = a * xp[k] + (yp ? b * yp[k] : 0.0);
    }
  }
}

static void
fluid_diffusion_prepare_implicit(gkyl_gyrokinetic_app *app, struct gk_neut_species *species)
{
  if (species->use_reaction_rate_diffusion) {
    struct gk_react *react = &species->react_neut;
    int iz_idx = species->ionization_react_idx;
    int cx_idx = species->diffusion_cx_react_idx;
    const struct gkyl_array *vt_sq_D = fluid_diffusion_vt_sq(app, species);
    struct gk_species *elc = &app->species[react->elc_idx[iz_idx]];
    gkyl_dg_cx_coll_rate(
      react->cx[cx_idx], app->species[react->ion_idx[cx_idx]].lte.moms.marr, species->lte.moms.marr,
      react->upar_ion[cx_idx], species->diffusion_coeff
    );
    gkyl_dg_mul_op_range(
      &app->basis, 0, species->diffusion_moment_ratio, 0, elc->lte.moms.marr, 0,
      species->diffusion_coeff, &species->local
    );

    const double nu_min = species->info.diffusion.min_collision_frequency;
    if (nu_min > 0.0) {
#ifdef GKYL_HAVE_CUDA
      if (app->use_gpu) {
        gk_neut_fluid_diffusion_floor_cu(
          app->cdim, nu_min, &species->local, vt_sq_D,
          species->diffusion_moment_ratio, species->diffusion_coeff
        );
      } else
#endif
      {
        double nodes[8 * GKYL_MAX_CDIM] = {0.0}, D_nodal[8] = {0.0};
        app->basis.node_list(nodes);
        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, &species->local);
        while (gkyl_range_iter_next(&iter)) {
          long loc = gkyl_range_idx(&species->local, iter.idx);
          const double *vt_sq = gkyl_array_cfetch(vt_sq_D, loc);
          const double *nu = gkyl_array_cfetch(species->diffusion_moment_ratio, loc);
          double *D = gkyl_array_fetch(species->diffusion_coeff, loc);
          for (int q = 0; q < app->basis.num_basis; ++q) {
            double vtq = app->basis.eval_expand(&nodes[q * app->cdim], vt_sq);
            double nuq = app->basis.eval_expand(&nodes[q * app->cdim], nu);
            D_nodal[q] = fmax(vtq, 0.0) / fmax(nuq, nu_min);
          }
          app->basis.nodal_to_modal(D_nodal, D);
        }
      }
    } else {
      gkyl_dg_div_op_range(
        species->diffusion_div_mem, &app->basis, 0, species->diffusion_coeff, 0,
        vt_sq_D, 0, species->diffusion_moment_ratio, &species->local
      );
    }
  }

  gkyl_dg_mul_op_range(
    &app->basis, 0, species->diffusion_geom_factor, 0, species->diffusion_coeff, 0,
    app->gk_geom->geo_int.jacobgeo, &species->local
  );
  const int metric_idx_1x[1][1] = {{5}};
  const int metric_idx_2x[2][2] = {{0, 2}, {2, 5}};
  const int metric_idx_3x[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
  for (int i = 0; i < app->cdim; ++i) {
    for (int j = 0; j < app->cdim; ++j) {
      int gij_idx = app->cdim == 1 ? metric_idx_1x[i][j] :
                    app->cdim == 2 ? metric_idx_2x[i][j] :
                                     metric_idx_3x[i][j];
      gkyl_dg_mul_op_range(
        &app->basis, i * app->cdim + j, species->diffusion_tensor, 0,
        species->diffusion_geom_factor, gij_idx, app->gk_geom->geo_int.gij_neut, &species->local
      );
    }
  }
  gkyl_comm_array_per_sync(
    species->comm, &species->local, &species->local_ext, species->num_periodic_dir,
    species->periodic_dirs, species->diffusion_tensor
  );
  gkyl_comm_array_sync(
    species->comm, &species->local, &species->local_ext, species->diffusion_tensor
  );
}

static void
fluid_diffusion_apply_L(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *species, struct gkyl_array *rho,
  struct gkyl_array *out, bool homogeneous_bc
)
{
  gkyl_comm_array_per_sync(
    species->comm, &species->local, &species->local_ext, species->num_periodic_dir,
    species->periodic_dirs, rho
  );
  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, rho);

  gk_neut_species_fluid_diffusion_par_ghosts(app, species, rho, homogeneous_bc);

  gkyl_array_clear(out, 0.0);
  gkyl_array_clear(species->diffusion_cflrate, 0.0);
  gkyl_dg_updater_conf_diffusion_advance(
    species->diffusion_slvr, &species->local, rho, species->diffusion_cflrate, out
  );
  // The recycled flux is independent of the neutral state, so it belongs to
  // the inhomogeneous part c of the affine operator L(rho)=L_h(rho)+c.
  if (!homogeneous_bc) {
    gk_neut_species_fluid_diffusion_add_recycling(app, species, out, true);
  }
}

static void
fluid_diffusion_apply_A(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *species, double dt, struct gkyl_array *x,
  struct gkyl_array *out
)
{
  gkyl_dg_mul_op_range(
    &app->basis, 0, out, 0, x, 0, app->gk_geom->geo_int.jacobgeo, &species->local
  );
  fluid_diffusion_apply_L(app, species, x, species->diffusion_implicit_work, true);
  gkyl_array_accumulate_range(out, -dt, species->diffusion_implicit_work, &species->local);
}

// Limit the p1 slopes of physical density without changing its cell average.
// A multilinear p1 field reaches its extrema at the cell vertices, so scaling
// every non-constant modal coefficient by one factor makes rho nonnegative
// throughout the cell. A negative cell average cannot be repaired
// conservatively and is set to zero.
static void
fluid_diffusion_enforce_positivity(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *species, struct gkyl_array *rho
)
{
  if (app->basis.poly_order != 1) {
    return;
  }
#ifdef GKYL_HAVE_CUDA
  if (app->use_gpu) {
    gk_neut_fluid_diffusion_positivity_cu(app->cdim, &species->local, rho);
    return;
  }
#endif

  double nodes[8 * GKYL_MAX_CDIM] = {0.0};
  double center[GKYL_MAX_CDIM] = {0.0};
  app->basis.node_list(nodes);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &species->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&species->local, iter.idx);
    double *rhop = gkyl_array_fetch(rho, loc);
    double rho_avg = app->basis.eval_expand(center, rhop);
    if (rho_avg <= 0.0) {
      for (int k = 0; k < app->basis.num_basis; ++k) {
        rhop[k] = 0.0;
      }
      continue;
    }

    double rho_min = DBL_MAX;
    for (int q = 0; q < app->basis.num_basis; ++q) {
      rho_min = fmin(rho_min, app->basis.eval_expand(&nodes[q * app->cdim], rhop));
    }
    if (rho_min < 0.0) {
      double theta = fmin(1.0, rho_avg / (rho_avg - rho_min));
      // Move infinitesimally inside the admissible set to avoid a negative
      // vertex caused solely by roundoff.
      theta *= 1.0 - 1.0e-14;
      for (int k = 1; k < app->basis.num_basis; ++k) {
        rhop[k] *= theta;
      }
    }
  }
}

static double
gk_neut_species_fluid_rhs_implicit_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);
  if (!species->implicit_diffusion) {
    gkyl_array_set(rhs, 1.0, fin);
    return DBL_MAX;
  }

  fluid_diffusion_prepare_implicit(app, species);
  struct gkyl_array *b = species->diffusion_implicit_rhs;
  struct gkyl_array *x = species->diffusion_implicit_x;
  struct gkyl_array *r = species->diffusion_implicit_r;
  struct gkyl_array *rh = species->diffusion_implicit_rhat;
  struct gkyl_array *p = species->diffusion_implicit_p;
  struct gkyl_array *v = species->diffusion_implicit_v;
  struct gkyl_array *s = species->diffusion_implicit_s;
  struct gkyl_array *t = species->diffusion_implicit_t;

  // Initial guess is the pre-diffusion physical mass density.
  gkyl_dg_div_op_range(
    species->diffusion_div_mem, &app->basis, 0, x, 0, fin, 0, app->gk_geom->geo_int.jacobgeo,
    &species->local
  );

  // The inhomogeneous Dirichlet trace and the recycled flux make L affine. If
  // L(rho)=L_h(rho)+c, backward Euler is
  // (J-dt*L_h)rho_new=N_old+dt*c.
  gkyl_array_clear(s, 0.0);
  fluid_diffusion_apply_L(app, species, s, v, false);
#ifdef GKYL_HAVE_CUDA
  if (app->use_gpu) {
    gk_neut_fluid_diffusion_vec_set_cu(&species->local, b, 1.0, fin, dt, v);
  } else
#endif
  {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &species->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&species->local, iter.idx);
      double *bp = gkyl_array_fetch(b, loc);
      const double *fp = gkyl_array_cfetch(fin, loc);
      const double *cp = gkyl_array_cfetch(v, loc);
      for (int k = 0; k < app->basis.num_basis; ++k) {
        bp[k] = fp[k] + dt * cp[k];
      }
    }
  }

  fluid_diffusion_apply_A(app, species, dt, x, t);
  fluid_diffusion_vec_set(app, species, r, 1.0, b, -1.0, t);
  fluid_diffusion_vec_set(app, species, rh, 1.0, r, 0.0, 0);
  gkyl_array_clear(p, 0.0);
  gkyl_array_clear(v, 0.0);

  const double norm_b = sqrt(fmax(fluid_diffusion_dot(app, species, b, b), DBL_MIN));
  const double tol =
    species->info.diffusion.implicit_tol > 0.0 ? species->info.diffusion.implicit_tol : 1.0e-10;
  const int max_iter =
    species->info.diffusion.implicit_max_iter > 0 ? species->info.diffusion.implicit_max_iter : 200;
  double rho_old = 1.0, alpha = 1.0, omega = 1.0;
  double residual = sqrt(fluid_diffusion_dot(app, species, r, r)) / norm_b;
  int niter = 0;
  for (; niter < max_iter && residual > tol; ++niter) {
    double rho_new = fluid_diffusion_dot(app, species, rh, r);
    if (fabs(rho_new) < DBL_MIN) {
      break;
    }
    double beta = (rho_new / rho_old) * (alpha / omega);
    // p = r + beta*(p-omega*v).
#ifdef GKYL_HAVE_CUDA
    if (app->use_gpu) {
      gk_neut_fluid_diffusion_p_update_cu(&species->local, p, r, v, beta, omega);
    } else
#endif
    {
      struct gkyl_range_iter vit;
      gkyl_range_iter_init(&vit, &species->local);
      while (gkyl_range_iter_next(&vit)) {
        long loc = gkyl_range_idx(&species->local, vit.idx);
        double *pp = gkyl_array_fetch(p, loc);
        const double *rp = gkyl_array_cfetch(r, loc);
        const double *vp = gkyl_array_cfetch(v, loc);
        for (int k = 0; k < app->basis.num_basis; ++k) {
          pp[k] = rp[k] + beta * (pp[k] - omega * vp[k]);
        }
      }
    }
    fluid_diffusion_apply_A(app, species, dt, p, v);
    double rhv = fluid_diffusion_dot(app, species, rh, v);
    if (fabs(rhv) < DBL_MIN) {
      break;
    }
    alpha = rho_new / rhv;
    fluid_diffusion_vec_set(app, species, s, 1.0, r, -alpha, v);
    double snorm = sqrt(fluid_diffusion_dot(app, species, s, s)) / norm_b;
    if (snorm <= tol) {
      gkyl_array_accumulate_range(x, alpha, p, &species->local);
      residual = snorm;
      ++niter;
      break;
    }
    fluid_diffusion_apply_A(app, species, dt, s, t);
    double tt = fluid_diffusion_dot(app, species, t, t);
    if (tt < DBL_MIN) {
      break;
    }
    omega = fluid_diffusion_dot(app, species, t, s) / tt;
    gkyl_array_accumulate_range(x, alpha, p, &species->local);
    gkyl_array_accumulate_range(x, omega, s, &species->local);
    fluid_diffusion_vec_set(app, species, r, 1.0, s, -omega, t);
    residual = sqrt(fluid_diffusion_dot(app, species, r, r)) / norm_b;
    if (fabs(omega) < DBL_MIN) {
      break;
    }
    rho_old = rho_new;
  }
  species->diffusion_implicit_last_iter = niter;
  species->diffusion_implicit_last_residual = residual;
  if (residual > tol) {
    gkyl_gyrokinetic_app_cout(
      app, stderr,
      "WARNING: implicit neutral diffusion for %s did not converge: iter=%d residual=%.6e\n",
      species->info.name, niter, residual
    );
  }

  fluid_diffusion_enforce_positivity(app, species, x);

  // Convert physical rho back to N=J*rho and retain the specific values of
  // dependent moments while updating their density factor.
  gkyl_dg_mul_op_range(
    &app->basis, 0, rhs, 0, x, 0, app->gk_geom->geo_int.jacobgeo, &species->local
  );
  for (int m = 1; m < species->num_moments; ++m) {
    gkyl_dg_div_op_range(
      species->diffusion_div_mem, &app->basis, m, species->diffusion_moment_ratio, m, fin, 0, fin,
      &species->local
    );
    gkyl_dg_mul_op_range(
      &app->basis, m, rhs, m, species->diffusion_moment_ratio, 0, rhs, &species->local
    );
  }
  return DBL_MAX;
}

static void
gk_neut_species_fluid_release_dynamic(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *ns
)
{
  // Release memory allocated for dynamic neutrals.
  gkyl_array_release(ns->cflrate);

  if (ns->has_diffusion) {
    gkyl_array_release(ns->diffusion_coeff);
    if (ns->diffusion_coeff_host) {
      gkyl_array_release(ns->diffusion_coeff_host);
    }
    gkyl_array_release(ns->diffusion_tensor);
    gkyl_array_release(ns->diffusion_geom_factor);
    gkyl_array_release(ns->diffusion_density);
    gkyl_array_release(ns->diffusion_cflrate);
    if (ns->diffusion_cflrate_host) {
      gkyl_array_release(ns->diffusion_cflrate_host);
    }
    if (ns->implicit_diffusion) {
      gkyl_array_release(ns->diffusion_implicit_rhs);
      gkyl_array_release(ns->diffusion_implicit_x);
      gkyl_array_release(ns->diffusion_implicit_r);
      gkyl_array_release(ns->diffusion_implicit_rhat);
      gkyl_array_release(ns->diffusion_implicit_p);
      gkyl_array_release(ns->diffusion_implicit_v);
      gkyl_array_release(ns->diffusion_implicit_s);
      gkyl_array_release(ns->diffusion_implicit_t);
      gkyl_array_release(ns->diffusion_implicit_work);
      if (app->use_gpu) {
        gkyl_cu_free(ns->diffusion_implicit_dot);
      }
    }
    if (ns->diffusion_recycling_flux) {
      gkyl_array_release(ns->diffusion_recycling_flux);
    }
    if (ns->diffusion_recycling_flux_tmp) {
      gkyl_array_release(ns->diffusion_recycling_flux_tmp);
    }
    gkyl_array_release(ns->diffusion_moment_ratio);
    if (ns->diffusion_vt_sq) {
      gkyl_array_release(ns->diffusion_vt_sq);
    }
    gkyl_dg_bin_op_mem_release(ns->diffusion_div_mem);
    gkyl_dg_updater_conf_diffusion_release(ns->diffusion_slvr);
  }
  if (ns->ionization_react_idx >= 0) {
    gkyl_array_release(ns->ionization_rate);
    gkyl_array_release(ns->ionization_cflrate);
    if (ns->ionization_cflrate_host) {
      gkyl_array_release(ns->ionization_cflrate_host);
    }
  }

  if (app->use_gpu) {
    gkyl_cu_free(ns->omega_cfl);
  } else {
    gkyl_free(ns->omega_cfl);
  }

  // Release integrated mom data.
  gk_neut_species_moment_release(app, &ns->integ_moms);

  // Release integrated mom diag data.
  gkyl_dynvec_release(ns->integ_diag);

  if (app->use_gpu) {
    gkyl_cu_free(ns->red_integ_diag);
    gkyl_cu_free(ns->red_integ_diag_global);
  } else {
    gkyl_free(ns->red_integ_diag);
    gkyl_free(ns->red_integ_diag_global);
  }
}

static void
gk_neut_species_fluid_release(const gkyl_gyrokinetic_app *app, const struct gk_neut_species *ns)
{
  // Release resources for fluid neutral species.
  gkyl_msgpack_map_elem_release(ns->io_meta_basic_len, ns->io_meta_basic);
  gkyl_msgpack_map_elem_release(ns->io_meta_conf_len, ns->io_meta_conf);

  gkyl_array_release(ns->f);
  gkyl_array_release(ns->f1);
  gkyl_array_release(ns->fnew);
  gkyl_array_release(ns->f_host);

  if (ns->info.init_from_file.type == 0) {
    gk_neut_species_projection_release(app, &ns->proj_init);
  }
  gkyl_comm_release(ns->comm);

  for (int i = 0; i < ns->info.num_diag_moments; ++i) {
    gk_neut_species_moment_release(app, &ns->moms[i]);
  }
  gkyl_free(ns->moms);

  gk_neut_species_bgk_release(app, &ns->bgk);

  gk_neut_species_react_release(app, &ns->react_neut);

  // Free boundary flux memory.
  gk_neut_species_bflux_release(app, ns, &ns->bflux);

  gk_neut_species_lte_release(app, &ns->lte);

  // Free memory for the object that scales the species according to a balance
  // between recycling and reactions.
  gk_neut_species_scaling_release(app, &ns->sca);

  ns->release_is_static_func(app, ns);
}

static void
gk_neut_species_fluid_init_dynamic(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns
)
{
  int cdim = app->cdim;

  // Allocate additional moment arrays for time stepping.
  ns->f1 = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);
  ns->fnew = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);

  // Allocate cflrate (scalar array).
  ns->cflrate = mkarr(app->use_gpu, 1, ns->local_ext.volume);

  ns->omega_cfl = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));

  // Initialize diffusion of neutral density (component 0) only. The same
  // scalar coefficient is used in every configuration-space direction.
  ns->use_reaction_rate_diffusion = ns->info.diffusion.use_reaction_rates;
  ns->has_diffusion = ns->info.diffusion.D > 0.0 || ns->use_reaction_rate_diffusion;
  ns->implicit_diffusion = ns->has_diffusion && ns->info.diffusion.is_implicit;
  ns->ionization_react_idx = -1;
  ns->diffusion_cx_react_idx = -1;
  ns->diffusion_coeff = 0;
  ns->diffusion_coeff_host = 0;
  ns->diffusion_tensor = 0;
  ns->diffusion_geom_factor = 0;
  ns->diffusion_density = 0;
  ns->diffusion_cflrate = 0;
  ns->diffusion_cflrate_host = 0;
  ns->diffusion_recycling_flux = 0;
  ns->diffusion_recycling_flux_tmp = 0;
  ns->diffusion_moment_ratio = 0;
  ns->diffusion_vt_sq = 0;
  ns->diffusion_div_mem = 0;
  ns->diffusion_slvr = 0;
  ns->diffusion_recycling_num_species[0] = 0;
  ns->diffusion_recycling_num_species[1] = 0;
  ns->ionization_rate = 0;
  ns->ionization_cflrate = 0;
  ns->ionization_cflrate_host = 0;

  for (int i = 0; i < ns->info.react_neut.num_react; ++i) {
    enum gkyl_react_id react_id = ns->info.react_neut.react_type[i].react_id;
    if (react_id == GKYL_REACT_IZ && ns->ionization_react_idx < 0) {
      ns->ionization_react_idx = i;
    } else if (react_id == GKYL_REACT_CX && ns->diffusion_cx_react_idx < 0) {
      ns->diffusion_cx_react_idx = i;
    }
  }
  if (ns->ionization_react_idx >= 0) {
    ns->ionization_rate = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    ns->ionization_cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  }
  if (ns->has_diffusion) {
    // Configuration-space restriction of the single-integration-by-parts plus
    // recovery-DG anomalous-diffusion scheme. Generated CPU kernels currently
    // cover 1x-3x serendipity p1.
    assert(cdim > 0 && cdim <= 3);
    assert(app->basis.poly_order == 1);
    ns->diffusion_coeff = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    ns->diffusion_tensor =
      mkarr(app->use_gpu, cdim * cdim * app->basis.num_basis, app->local_ext.volume);
    ns->diffusion_geom_factor = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    ns->diffusion_density = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    ns->diffusion_cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
    if (ns->info.diffusion.write_diagnostics) {
      ns->diffusion_coeff_host = app->use_gpu ?
                                   mkarr(false, app->basis.num_basis, app->local_ext.volume) :
                                   gkyl_array_acquire(ns->diffusion_coeff);
      ns->diffusion_cflrate_host = app->use_gpu ? mkarr(false, 1, app->local_ext.volume) :
                                                  gkyl_array_acquire(ns->diffusion_cflrate);
      if (ns->ionization_cflrate) {
        ns->ionization_cflrate_host = app->use_gpu ? mkarr(false, 1, app->local_ext.volume) :
                                                     gkyl_array_acquire(ns->ionization_cflrate);
      }
    }
    if (ns->implicit_diffusion) {
      ns->diffusion_implicit_rhs = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_x = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_r = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_rhat =
        mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_p = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_v = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_s = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_t = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_work =
        mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_implicit_dot = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : 0;
      ns->diffusion_implicit_last_iter = 0;
      ns->diffusion_implicit_last_residual = 0.0;
    }

    const int par_dir = cdim - 1;
    const bool recycle_lo = ns->lower_bc[par_dir].type == GKYL_BC_GK_SPECIES_RECYCLE;
    const bool recycle_up = ns->upper_bc[par_dir].type == GKYL_BC_GK_SPECIES_RECYCLE;
    if (recycle_lo || recycle_up) {
      ns->diffusion_recycling_flux = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      ns->diffusion_recycling_flux_tmp =
        mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
      gkyl_array_clear(ns->diffusion_recycling_flux, 0.0);
      gkyl_array_clear(ns->diffusion_recycling_flux_tmp, 0.0);
    }
    for (int e = 0; e < 2; ++e) {
      const struct gkyl_gyrokinetic_bc *bc = e == 0 ? &ns->lower_bc[par_dir] :
                                                      &ns->upper_bc[par_dir];
      if (bc->type != GKYL_BC_GK_SPECIES_RECYCLE) {
        continue;
      }
      assert(bc->emission.num_species > 0);
      ns->diffusion_recycling_num_species[e] = bc->emission.num_species;
      for (int i = 0; i < bc->emission.num_species; ++i) {
        ns->diffusion_recycling_ion_idx[e][i] =
          gk_find_species_idx(app, bc->emission.in_species[i]);
        assert(ns->diffusion_recycling_ion_idx[e][i] >= 0);
      }
    }
    gkyl_array_clear(ns->diffusion_coeff, 0.0);
    gkyl_array_clear(ns->diffusion_tensor, 0.0);
    gkyl_array_clear(ns->diffusion_geom_factor, 0.0);
    gkyl_array_clear(ns->diffusion_density, 0.0);
    gkyl_array_clear(ns->diffusion_cflrate, 0.0);
    if (!ns->use_reaction_rate_diffusion) {
      // Modal constant coefficient in cdim dimensions.
      gkyl_array_shiftc(ns->diffusion_coeff, pow(sqrt(2.0), cdim) * ns->info.diffusion.D, 0);
    }

    if (ns->use_reaction_rate_diffusion) {
      assert(ns->ionization_react_idx >= 0);
      assert(ns->diffusion_cx_react_idx >= 0);
      ns->diffusion_vt_sq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    }

    ns->diffusion_moment_ratio = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);
    ns->diffusion_div_mem =
      app->use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis) :
                     gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);

    int diffusion_zero_flux_flags[2 * GKYL_MAX_CDIM] = {0};
    for (int d = 0; d < cdim; ++d) {
      diffusion_zero_flux_flags[d] = diffusion_zero_flux_flags[d + cdim] = ns->bc_is_np[d];
    }
    if (recycle_lo || ns->info.diffusion.lower_bc_type == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET) {
      diffusion_zero_flux_flags[par_dir] = 0;
    }
    if (recycle_up || ns->info.diffusion.upper_bc_type == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET) {
      diffusion_zero_flux_flags[par_dir + cdim] = 0;
    }
    ns->diffusion_slvr = gkyl_dg_updater_conf_diffusion_new(
      &app->grid, &app->basis, &app->local_ext, ns->diffusion_tensor,
      app->gk_geom->geo_int.jacobgeo_inv, diffusion_zero_flux_flags, app->use_gpu
    );
  }

  // Allocate data for integrated moments.
  gk_neut_species_moment_init(app, ns, &ns->integ_moms, GKYL_F_MOMENT_M0M1M2, true);

  // Allocate data for integrated diagnostics.
  if (app->use_gpu) {
    ns->red_integ_diag = gkyl_cu_malloc(sizeof(double[ns->integ_moms.num_mom]));
    ns->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[ns->integ_moms.num_mom]));
  } else {
    ns->red_integ_diag = gkyl_malloc(sizeof(double[ns->integ_moms.num_mom]));
    ns->red_integ_diag_global = gkyl_malloc(sizeof(double[ns->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  ns->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, ns->integ_moms.num_mom);
  ns->is_first_integ_write_call = true;

  // Set function pointers
  ns->rhs_func = gk_neut_species_fluid_rhs_dynamic;
  ns->rhs_implicit_func = gk_neut_species_fluid_rhs_implicit_dynamic;
  ns->bc_func = gk_neut_species_fluid_apply_bc_dynamic;
  ns->release_func = gk_neut_species_fluid_release;
  ns->release_is_static_func = gk_neut_species_fluid_release_dynamic;
  ns->step_f_func = gk_neut_species_step_f_dynamic;
  ns->combine_func = gk_neut_species_combine_dynamic;
  ns->copy_func = gk_neut_species_copy_range_dynamic;
  ns->write_func = gk_neut_species_write_dynamic;
  ns->write_mom_func =
    gk_neut_species_write_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->calc_integrated_mom_func =
    gk_neut_species_calc_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->write_integrated_mom_func =
    gk_neut_species_write_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

static void
gk_neut_species_fluid_init_static(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s
)
{
  // Set pointers for RK methods.
  s->f1 = gkyl_array_acquire(s->f);
  s->fnew = gkyl_array_acquire(s->f);

  // Set function pointers
  s->rhs_func = gk_neut_species_rhs_static;
  s->rhs_implicit_func = gk_neut_species_rhs_implicit_static;
  s->bc_func = gk_neut_species_apply_bc_static;
  s->release_func = gk_neut_species_fluid_release;
  s->release_is_static_func = gk_neut_species_release_static;
  s->step_f_func = gk_neut_species_step_f_static;
  s->combine_func = gk_neut_species_combine_static;
  s->copy_func = gk_neut_species_copy_range_static;
  s->write_func = gk_neut_species_write_init_only;
  s->write_mom_func = gk_neut_species_write_mom_init_only;
  s->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_static;
  s->write_integrated_mom_func = gk_neut_species_write_integrated_mom_static;
  s->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

void
gk_neut_species_fluid_init(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns
)
{
  ns->is_fluid = true; // Fluid neutrals.
  assert(
    ns->info.vdim == 0
  ); // Ensure user provided vdim=0 in input file, or didn't provide it at all.

  ns->model_id = GKYL_MODEL_DEFAULT;
  ns->field_id = GKYL_FIELD_NULL;

  int cdim = app->cdim;

  // Number of moments depends on eqn_type.
  ns->num_moments = 5; // rho, rho*ux, rho*uy, rho*uy, totE

  // Use the same basis as conf-space.
  ns->basis = app->basis;
  ns->basis_on_dev = app->basis_on_dev;

  // Use the same grid as conf-space.
  ns->grid = app->grid;
  ns->global_ext = app->global_ext;
  ns->global = app->global;

  // Use the same communicator as conf-space.
  ns->comm = gkyl_comm_acquire(app->comm);

  // Use the same local range as conf-space.
  ns->local_ext = app->local_ext;
  ns->local = app->local;

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // add the parallel direction in case TS BCs are needed.
  ns->num_periodic_dir = app->num_periodic_dir;
  for (int d = 0; d < ns->num_periodic_dir; ++d) {
    ns->periodic_dirs[d] = app->periodic_dirs[d];
  }

  // In an LCFS-spanning domain, populate parallel ghosts periodically first.
  // The diffusion RHS later overwrites only the SOL subset with recycled or
  // fixed Dirichlet data, leaving the confined subset periodic.
  if (app->gk_geom->has_LCFS) {
    ns->periodic_dirs[ns->num_periodic_dir] = app->cdim - 1;
    ns->num_periodic_dir += 1;
  }

  for (int d = 0; d < app->cdim; ++d) {
    ns->bc_is_np[d] = true;
  }
  for (int d = 0; d < ns->num_periodic_dir; ++d) {
    ns->bc_is_np[ns->periodic_dirs[d]] = false;
  }

  // Store the BCs from the input file.
  for (int d = 0; d < app->cdim; ++d) {
    struct gkyl_gyrokinetic_bc *bc_lo =
      gk_fetch_bc_with_dir_edge(ns->info.bcs, 2 * app->cdim, d, GKYL_LOWER_EDGE);
    if (bc_lo != 0) {
      ns->lower_bc[d] = *bc_lo;
    } else {
      ns->lower_bc[d].type = GKYL_BC_GK_SKIP;
    }

    struct gkyl_gyrokinetic_bc *bc_up =
      gk_fetch_bc_with_dir_edge(ns->info.bcs, 2 * app->cdim, d, GKYL_UPPER_EDGE);
    if (bc_up != 0) {
      ns->upper_bc[d] = *bc_up;
    } else {
      ns->upper_bc[d].type = GKYL_BC_GK_SKIP;
    }
  }

  // Species properties metadata.
  struct gkyl_msgpack_map_elem io_meta_sprop[] = {
    {.key = "mass", .elem_type = GKYL_MP_DOUBLE, .dval = ns->info.mass},
    {.key = "charge", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "gas_gamma", .elem_type = GKYL_MP_DOUBLE, .dval = ns->info.gas_gamma},
    {.key = "vdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = ns->info.vdim},
  };

  // Metadata for integrated quantities.
  const struct gkyl_msgpack_map_elem *io_meta_basic_union[] = {app->io_meta_basic, io_meta_sprop};
  int io_meta_basic_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0])
  };
  ns->io_meta_basic = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_basic_union) / sizeof(io_meta_basic_union[0]), io_meta_basic_union_len,
    io_meta_basic_union, &ns->io_meta_basic_len
  );

  // Metadata for conf-space quantities.
  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    {.key = "value_form", .elem_type = GKYL_MP_STRING, .cval = "modal"},
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = app->basis.id},
    {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0},
  };
  const struct gkyl_msgpack_map_elem *io_meta_conf_union[] = {
    app->io_meta_basic, io_meta_sprop, io_meta_conf
  };
  int io_meta_conf_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0]),
    sizeof(io_meta_conf) / sizeof(io_meta_conf[0])
  };
  ns->io_meta_conf = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_conf_union) / sizeof(io_meta_conf_union[0]), io_meta_conf_union_len,
    io_meta_conf_union, &ns->io_meta_conf_len
  );

  // Metadata for phase-space quantities.
  ns->io_meta_phase = ns->io_meta_conf;
  ns->io_meta_phase_len = ns->io_meta_conf_len;

  // Allocate distribution function array for initialization and I/O.
  ns->f = mkarr(app->use_gpu, ns->num_moments * ns->basis.num_basis, ns->local_ext.volume);
  ns->f_host = app->use_gpu ? mkarr(false, ns->f->ncomp, ns->f->size) : gkyl_array_acquire(ns->f);

  // Create skin/ghost ranges fir applying BCs. Only used for dynamic neutrals but included here to avoid
  // code duplication since the "ghost" array is needed.
  int ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    ghost[d] = 1;
  }

  // Create skin/ghost ranges.
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &ns->local_lower_skin[dir], &ns->local_lower_ghost[dir], dir, GKYL_LOWER_EDGE, &ns->local_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &ns->local_upper_skin[dir], &ns->local_upper_ghost[dir], dir, GKYL_UPPER_EDGE, &ns->local_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &ns->global_lower_skin[dir], &ns->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &ns->global_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &ns->global_upper_skin[dir], &ns->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &ns->local_ext, ghost
    );
  }

  // Initialize projection routine for initial conditions.
  if (ns->info.init_from_file.type == 0) {
    gk_neut_species_projection_init(app, ns, ns->info.projection, &ns->proj_init);
  }

  // Allocate objects for computing diagnostic moments.
  int ndm = ns->info.num_diag_moments;
  ns->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m = 0; m < ndm; ++m) {
    gk_neut_species_moment_init(app, ns, &ns->moms[m], ns->info.diag_moments[m], false);
  }

  // Initialize boundary fluxes.
  ns->bflux = (struct gk_boundary_fluxes){};
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp){};
  // Set the operation type for the bflux app.
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  gk_neut_species_bflux_init(app, ns, &ns->bflux, bflux_type, add_bflux_moms_inp);

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  ns->lte = (struct gk_lte){};
  struct correct_all_moms_inp corr_inp = {
    .correct_all_moms = false,
    .max_iter = 0,
    .iter_eps = 10,
    .use_last_converged = false,
  };
  gk_neut_species_lte_init(app, ns, &ns->lte, corr_inp);

  // Initialize the object that scales the species according to a balance
  // between recycling and reactions.
  ns->sca = (struct gk_scaling){};
  gk_neut_species_scaling_init(app, ns, &ns->sca);

  // Initialize BGK collisions with null type (not applicable to fluids).
  ns->bgk = (struct gk_bgk_collisions){};
  ns->info.collisions.collision_id = 0;
  gk_neut_species_bgk_init(app, ns, &ns->bgk);

  // Initialize reactions with charged species. Fluid neutrals use the
  // ionization reaction as a density loss; CX supplies the diffusion rate.
  ns->react_neut = (struct gk_react){};
  gk_neut_species_react_init(app, ns, ns->info.react_neut, &ns->react_neut);

  ns->src = (struct gk_source){};
  if (!ns->info.is_static) {
    gk_neut_species_fluid_init_dynamic(gk, app, ns);
  } else {
    gk_neut_species_fluid_init_static(gk, app, ns);
  }
}
