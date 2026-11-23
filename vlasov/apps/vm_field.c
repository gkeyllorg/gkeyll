#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_priv.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static bool
vm_field_is_fixed_func_bc(const struct vm_field *field, int d)
{
  return field->lower_bc[d] == GKYL_FIELD_FIXED_FUNC || field->upper_bc[d] == GKYL_FIELD_FIXED_FUNC;
}

static struct gkyl_array *
vm_field_get_fixed_func_bc_buffer(
  const struct vm_field *field, const struct gkyl_array *f, int d, enum gkyl_edge_loc edge
)
{
  bool use_no_J_buffer = field->field_id == GKYL_FIELD_GR_D_B && f == field->em_no_J;

  if (edge == GKYL_LOWER_EDGE) {
    return use_no_J_buffer ? field->bc_buffer_lo_fixed_no_J[d] : field->bc_buffer_lo_fixed[d];
  }

  return use_no_J_buffer ? field->bc_buffer_up_fixed_no_J[d] : field->bc_buffer_up_fixed[d];
}

void
vm_field_buffer_fixed_func_bc(gkyl_vlasov_app *app, struct vm_field *field)
{
  bool has_fixed_func_bc = false;
  for (int d = 0; d < app->cdim; ++d) {
    has_fixed_func_bc = has_fixed_func_bc || vm_field_is_fixed_func_bc(field, d);
  }

  if (!has_fixed_func_bc) {
    return;
  }

  if (field->field_id == GKYL_FIELD_GR_D_B) {
    gkyl_dg_gr_maxwell_divide_Jc(
      &app->basis, &app->local_ext, app->vm_geom->det_h, field->em, field->em_no_J, app->use_gpu
    );
  }

  for (int d = 0; d < app->cdim; ++d) {
    if (field->lower_bc[d] == GKYL_FIELD_FIXED_FUNC) {
      gkyl_array_copy_to_buffer(field->bc_buffer_lo_fixed[d]->data, field->em, &app->lower_ghost[d]);
      if (field->field_id == GKYL_FIELD_GR_D_B) {
        gkyl_array_copy_to_buffer(
          field->bc_buffer_lo_fixed_no_J[d]->data, field->em_no_J, &app->lower_ghost[d]
        );
      }
    }

    if (field->upper_bc[d] == GKYL_FIELD_FIXED_FUNC) {
      gkyl_array_copy_to_buffer(field->bc_buffer_up_fixed[d]->data, field->em, &app->upper_ghost[d]);
      if (field->field_id == GKYL_FIELD_GR_D_B) {
        gkyl_array_copy_to_buffer(
          field->bc_buffer_up_fixed_no_J[d]->data, field->em_no_J, &app->upper_ghost[d]
        );
      }
    }
  }
}

// initialize field object
struct vm_field *
vm_field_new(struct gkyl_vm *vm, struct gkyl_vlasov_app *app)
{
  struct vm_field *f = gkyl_malloc(sizeof(struct vm_field));

  f->info = vm->field;
  f->field_id = f->info.field_id;
  for (int d = 0; d < 3; ++d) {
    f->bc_buffer_lo_fixed[d] = f->bc_buffer_up_fixed[d] = 0;
    f->bc_buffer_lo_fixed_no_J[d] = f->bc_buffer_up_fixed_no_J[d] = 0;
  }

  // A PKPM species takes the magnetic unit vector b = B/|B|, its gradient and d_t B/|B| of the
  // stage from the field, and reads the physical fields of the stage in em_no_J.
  f->has_bvar = false;
  for (int i = 0; i < vm->num_species; ++i) {
    f->has_bvar = f->has_bvar || vm->species[i].type == GKYL_SPECIES_PKPM;
  }
  if (f->has_bvar && f->field_id != GKYL_FIELD_E_B) {
    gkyl_exit("vm_field: a PKPM species needs the standard Maxwell field (E_B).");
  }

  // Standard E_B Maxwell on a non-identity position map evolves the J-weighted
  // fields J*E, J*B (so the mapped curl kernels and the J-weighted current are
  // consistent) and keeps the physical E, B in em_no_J for the Lorentz force and
  // I/O; with PKPM species it does so on every map (J = 1 on the identity). The GR
  // path uses the same em_no_J array (un-weighted by det_h instead).
  f->weight_by_pos_jacob = (f->field_id == GKYL_FIELD_E_B) &&
                           (!app->pos_map->is_identity || f->has_bvar);
  if (f->field_id == GKYL_FIELD_GR_D_B) {
    // Geometry of GR-Maxwell, owned by the app.
    f->geom = app->vm_geom;
    f->em_no_J = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);
    f->em_no_J_host = app->use_gpu ? mkarr(false, f->em_no_J->ncomp, f->em_no_J->size) :
                                     gkyl_array_acquire(f->em_no_J);
  } else if (f->weight_by_pos_jacob) {
    f->em_no_J = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);
    f->em_no_J_host = app->use_gpu ? mkarr(false, f->em_no_J->ncomp, f->em_no_J->size) :
                                     gkyl_array_acquire(f->em_no_J);
  }

  // allocate EM arrays
  f->em = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);
  f->em1 = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);
  f->emnew = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);

  // Host EM array for  I/O.
  f->em_host = app->use_gpu ? mkarr(false, f->em->ncomp, f->em->size) : gkyl_array_acquire(f->em);

  f->em_energy = mkarr(app->use_gpu, 6, app->local_ext.volume);
  if (app->use_gpu) {
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[6]));
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 6);
  f->is_first_energy_write_call = true;

  // Resistive layer: E is damped at the rate sigma(x), sampled at the (p+1)^cdim Gauss-Legendre
  // points of each cell, explicitly on the stage right-hand side (vm_field_complete_update) or
  // exactly at each node of the source solve of the split (vm_fluid_em_coupling). B is not damped.
  // sigma multiplies the evolved (J-weighted on a mapped mesh) field, so it is sampled at physical
  // coordinates.
  int num_quad_tot = 1;
  for (int d = 0; d < app->cdim; ++d) {
    num_quad_tot *= app->poly_order + 1;
  }
  f->sigma = mkarr(app->use_gpu, num_quad_tot, app->local_ext.volume);
  gkyl_array_clear(f->sigma, 0.0);
  f->has_sigma = false;
  if (f->info.sigma) {
    f->has_sigma = true;
    f->sigma_slvr = gkyl_dg_calc_resistive_layer_inew(&(struct gkyl_dg_calc_resistive_layer_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .range = &app->local_ext,
      .eval = f->info.sigma,
      .ctx = f->info.sigma_ctx,
      .c2p_func = vm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
      .sigma = f->sigma,
      .use_gpu = app->use_gpu,
    });
  }
  // The largest rate of the layer on this rank, for the time-step of its explicit increment.
  f->sigma_max = 0.0;
  if (f->has_sigma) {
    double sigma_red_ho[num_quad_tot];
    if (app->use_gpu) {
      double *sigma_red = gkyl_cu_malloc(sizeof(double[num_quad_tot]));
      gkyl_array_reduce_range(sigma_red, f->sigma, GKYL_MAX, &app->local);
      gkyl_cu_memcpy(sigma_red_ho, sigma_red, sizeof(double[num_quad_tot]), GKYL_CU_MEMCPY_D2H);
      gkyl_cu_free(sigma_red);
    } else {
      gkyl_array_reduce_range(sigma_red_ho, f->sigma, GKYL_MAX, &app->local);
    }
    for (int q = 0; q < num_quad_tot; ++q) {
      f->sigma_max = fmax(f->sigma_max, sigma_red_ho[q]);
    }
  }

  // Initialize external EM fields (always used by implicit fluid sources, so always initialize)
  f->ext_em = mkarr(app->use_gpu, 6 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->ext_em, 0.0);
  f->has_ext_em = false;
  f->ext_em_evolve = false;
  // setup external electromagnetic field
  if (f->info.ext_em) {
    f->has_ext_em = true;
    if (f->info.ext_em_evolve) {
      f->ext_em_evolve = f->info.ext_em_evolve;
    }

    f->ext_em_host = app->use_gpu ? mkarr(false, f->ext_em->ncomp, f->ext_em->size) :
                                    gkyl_array_acquire(f->ext_em);
    // The external field enters the Lorentz force as a physical field, sampled at physical
    // coordinates.
    f->ext_em_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = app->basis.poly_order + 1,
      .num_ret_vals = 6,
      .eval = f->info.ext_em,
      .ctx = f->info.ext_em_ctx,
      .c2p_func = vm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // No external potentials in Vlasov-Maxwell.
  f->has_ext_pot = f->ext_pot_evolve = false;

  // Initialize applied currents (always used by implicit fluid sources, so always initialize)
  f->app_current = mkarr(app->use_gpu, 3 * app->basis.num_basis, app->local_ext.volume);
  gkyl_array_clear(f->app_current, 0.0);
  // app_current drives the J-weighted field equation and carries J; the implicit fluid-EM
  // coupling works with the physical current.
  f->app_current_no_J = 0;
  if (f->weight_by_pos_jacob) {
    f->app_current_no_J = mkarr(app->use_gpu, 3 * app->basis.num_basis, app->local_ext.volume);
    gkyl_array_clear(f->app_current_no_J, 0.0);
  }
  f->has_app_current = false;
  f->app_current_evolve = false;
  // setup external currents
  if (f->info.app_current) {
    if (f->field_id == GKYL_FIELD_GR_D_B) {
      gkyl_exit("vm_field: applied currents are not supported with GR-Maxwell fields.");
    }
    f->has_app_current = true;
    if (f->info.app_current_evolve) {
      f->app_current_evolve = f->info.app_current_evolve;
    }

    // Separate host array (also on CPUs): it holds the physical applied current for I/O, while
    // app_current carries the conf-Jacobian weight on a mapped mesh.
    f->app_current_host = mkarr(false, f->app_current->ncomp, f->app_current->size);
    // Sampled at physical coordinates; vm_field_calc_app_current weights it by the conf
    // Jacobian when the field evolves J*E.
    f->app_current_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &app->grid,
      .basis = &app->basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = app->basis.poly_order + 1,
      .num_ret_vals = 3,
      .eval = f->info.app_current,
      .ctx = f->info.app_current_ctx,
      .c2p_func = vm_pos_c2p,
      .c2p_func_ctx = app->pos_map,
    });
  }

  // allocate cflrate (scalar array)
  f->cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  if (app->use_gpu) {
    f->omega_cfl = gkyl_cu_malloc(sizeof(double));
  } else {
    f->omega_cfl = gkyl_malloc(sizeof(double));
  }

  // equation object
  double c = 1.0;
  if (f->field_id == GKYL_FIELD_GR_D_B) {
    // Check that the speed of light is in natural units for GR.
    assert(fabs(f->info.epsilon0 * f->info.mu0 - 1.0) < 1e-12);
  } else {
    c = 1 / sqrt(f->info.epsilon0 * f->info.mu0);
  }
  double ef = f->info.elcErrorSpeedFactor, mf = f->info.mgnErrorSpeedFactor;
  double K_phi = f->info.K_phi, K_psi = f->info.K_psi;

  struct gkyl_dg_eqn *eqn;

  // Modal surface expansion of the configuration-space flux (GR-Maxwell).
  if (f->field_id == GKYL_FIELD_GR_D_B) {
    // Number of surface nodes of the tensor basis.
    f->num_surf_conf_nodes = pow(app->poly_order + 1, app->cdim - 1);

    f->conf_flux_surf =
      mkarr(app->use_gpu, app->cdim * 8 * f->num_surf_conf_nodes, app->local_ext.volume);
    struct gkyl_dg_gr_maxwell_conf_flux_surf_inp inp_conf_flux = {
      .conf_basis = &app->basis,
      .conf_grid = &app->grid,
      .pos_map = app->pos_map,
      .field_id = f->field_id,
      .theta_pole_lo = app->vm_geom->theta_pole_lo,
      .theta_pole_up = app->vm_geom->theta_pole_up,
      .chi = ef,
      .gamma = mf,
      .use_gpu = app->use_gpu,
    };
    f->calc_conf_flux = gkyl_dg_gr_maxwell_conf_flux_surf_inew(&inp_conf_flux);
  }

  // Input structure for building the dg eqn object
  struct gkyl_dg_maxwell_inp inp_dg_maxwell = {
    .cbasis = &app->basis,
    .crange = &app->local,
    .jacob_pos =
      app->pos_map->jacob_pos, // position-map Jacobian for the mapped curl (identity => no effect)
    .conf_flux_surf = (f->field_id == GKYL_FIELD_GR_D_B) ? f->conf_flux_surf : 0,
    .lapse = (f->field_id == GKYL_FIELD_GR_D_B) ? app->vm_geom->lapse : 0,
    .shift = (f->field_id == GKYL_FIELD_GR_D_B) ? app->vm_geom->shift : 0,
    .h_ij = (f->field_id == GKYL_FIELD_GR_D_B) ? app->vm_geom->h_ij : 0,
    .h_ij_inv = (f->field_id == GKYL_FIELD_GR_D_B) ? app->vm_geom->h_ij_inv : 0,
    .det_h = (f->field_id == GKYL_FIELD_GR_D_B) ? app->vm_geom->det_h : 0,
    .lightSpeed = c,
    .field_id = f->field_id,
    .elcErrorSpeedFactor = ef,
    .mgnErrorSpeedFactor = mf,
    .use_gpu = app->use_gpu,
  };
  eqn = gkyl_dg_maxwell_inew(&inp_dg_maxwell);

  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_flags[2 * GKYL_MAX_DIM] = {false};

  // Absorbing edges: the solver updates the cells of the edge with the boundary surface kernels of
  // the equation (the outgoing waves only), on the ranks that own the edge. The BC updaters below
  // fill the ghost cells of the edge for the operations that read them (limiter).
  for (int d = 0; d < app->cdim; ++d) {
    const enum gkyl_field_bc_type *bc = d == 0 ? f->info.bcx : (d == 1 ? f->info.bcy : f->info.bcz);
    bool is_periodic = false;
    for (int i = 0; i < app->num_periodic_dir; ++i) {
      is_periodic = is_periodic || app->periodic_dirs[i] == d;
    }
    bool absorb_lo = !is_periodic && bc[0] == GKYL_FIELD_ABSORB;
    bool absorb_up = !is_periodic && bc[1] == GKYL_FIELD_ABSORB;
    if ((absorb_lo || absorb_up) && f->field_id != GKYL_FIELD_E_B) {
      gkyl_exit(
        "vm_field: absorbing boundaries are only supported with the standard Maxwell field (E_B)."
      );
    }
    zero_flux_flags[d] = absorb_lo && app->local.lower[d] == app->global.lower[d];
    zero_flux_flags[d + app->cdim] = absorb_up && app->local.upper[d] == app->global.upper[d];
  }

  // Maxwell solver
  f->slvr = gkyl_hyper_dg_new(
    &app->grid, &app->basis, eqn, app->cdim, up_dirs, zero_flux_flags, 1, app->use_gpu
  );

  // Check if limiter_fac is specified for adjusting how much diffusion is applied through slope limiter
  // If not specified, set to 0.0 and updater sets default behavior (1/sqrt(3); see gkyl_dg_calc_em_vars.h)
  double limiter_fac = f->info.limiter_fac == 0 ? 0.0 : f->info.limiter_fac;
  f->limit_em = f->info.limit_em == 0 ? false : true;
  if (f->limit_em && f->has_bvar) {
    // The slope limiter does not keep the discrete divergence of B that the Maxwell solver
    // preserves, and b = B/|B| of a PKPM species follows the limited field.
    gkyl_exit(
      "vm_field: the limiter of the electromagnetic fields (limit_em) is not supported with PKPM "
      "species."
    );
  }

  struct gkyl_wv_eqn *maxwell = gkyl_wv_maxwell_new(c, ef, mf, app->use_gpu);
  // Create updaters for limiting EM fields
  f->calc_em_vars = gkyl_dg_calc_em_vars_new(
    &app->grid, &app->basis, &app->local_ext, maxwell, app->geom, limiter_fac, 0, app->use_gpu
  );
  gkyl_wv_eqn_release(maxwell);

  // determine which directions are not periodic
  int num_periodic_dir = app->num_periodic_dir, is_np[3] = {1, 1, 1};
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np[app->periodic_dirs[d]] = 0;
  }

  for (int dir = 0; dir < app->cdim; ++dir) {
    f->lower_bc[dir] = f->upper_bc[dir] = GKYL_FIELD_COPY;
    if (is_np[dir]) {
      const enum gkyl_field_bc_type *bc;
      if (dir == 0) {
        bc = f->info.bcx;
      } else if (dir == 1) {
        bc = f->info.bcy;
      } else {
        bc = f->info.bcz;
      }

      f->lower_bc[dir] = bc[0];
      f->upper_bc[dir] = bc[1];
    }
  }

  f->use_ghost_current = false;
  if (f->info.use_ghost_current) {
    if (app->cdim != 1 || app->num_periodic_dir != 1) {
      gkyl_exit("vm_field: ghost currents need one periodic configuration-space dimension.");
    }
    f->use_ghost_current = true;
    f->ghost_current = mkarr(app->use_gpu, 1, app->local_ext.volume);
    if (app->use_gpu) {
      f->red_ghost_current = gkyl_cu_malloc(sizeof(double[1]));
    }
  }

  // The magnetic unit vector b = B/|B| (p=1) by nodal normalization, its weak gradient, and the
  // physical time derivative of the fields of the stage with d_t B/|B| (zero for a static field).
  if (f->has_bvar) {
    int num_basis = app->basis.num_basis;
    long num_cells = app->local_ext.volume;
    f->tot_em = mkarr(app->use_gpu, 8 * num_basis, num_cells);
    f->bvar = mkarr(app->use_gpu, 3 * num_basis, num_cells);
    f->grad_b = mkarr(app->use_gpu, app->cdim * 3 * num_basis, num_cells);
    gkyl_array_clear(f->grad_b, 0.0);
    f->dt_em_no_J = mkarr(app->use_gpu, 8 * num_basis, num_cells);
    gkyl_array_clear(f->dt_em_no_J, 0.0);
    f->dtB_over_magB = mkarr(app->use_gpu, 3 * num_basis, num_cells);
    gkyl_array_clear(f->dtB_over_magB, 0.0);
    f->calc_bvar = gkyl_dg_calc_pkpm_bvar_inew(&(struct gkyl_dg_calc_pkpm_bvar_inp){
      .conf_grid = &app->grid,
      .conf_basis = &app->basis,
      .pos_map = app->pos_map,
      .use_gpu = app->use_gpu,
    });
  }

  f->use_geom_sources = false;
  f->geom_source = 0;
  f->calc_geom_source = 0;

  // Geometric sources are only active for GR_D_B fields. For that model,
  // use the input flag, which defaults to true in Lua.
  bool use_geom_sources = f->field_id == GKYL_FIELD_GR_D_B && f->info.use_geom_sources;
  if (use_geom_sources) {
    f->use_geom_sources = true;
    f->geom_source = mkarr(app->use_gpu, 8 * app->basis.num_basis, app->local_ext.volume);
    struct gkyl_dg_gr_maxwell_geom_source_inp inp_geom_source = {
      .conf_basis = &app->basis,
      .conf_grid = &app->grid,
      .field_id = f->field_id,
      .chi = ef,
      .gamma = mf,
      .K_phi = K_phi,
      .K_psi = K_psi,
      .use_gpu = app->use_gpu,
    };
    f->calc_geom_source = gkyl_dg_gr_maxwell_geom_source_inew(&inp_geom_source);
  }

  // allocate buffer for applying BCs
  long buff_sz = 0;
  // compute buffer size needed
  for (int dir = 0; dir < app->cdim; ++dir) {
    long vol = GKYL_MAX2(app->lower_skin[dir].volume, app->upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  f->bc_buffer = mkarr(app->use_gpu, 8 * app->basis.num_basis, buff_sz);

  for (int d = 0; d < app->cdim; ++d) {
    // Lower BC updater. Copy BCs by default.
    enum gkyl_bc_basic_type bctype = GKYL_BC_COPY;
    if (f->lower_bc[d] == GKYL_FIELD_COPY) {
      bctype = GKYL_BC_COPY;
    } else if (f->lower_bc[d] == GKYL_FIELD_PEC_WALL) {
      bctype = GKYL_BC_MAXWELL_PEC;
    } else if (f->lower_bc[d] == GKYL_FIELD_SYM_WALL) {
      bctype = GKYL_BC_MAXWELL_SYM;
    } else if (f->lower_bc[d] == GKYL_FIELD_RESERVOIR) {
      bctype = GKYL_BC_MAXWELL_RESERVOIR;
    } else if (f->lower_bc[d] == GKYL_FIELD_FIXED_FUNC) {
      bctype = GKYL_BC_FIXED_FUNC;
    } else if (f->lower_bc[d] == GKYL_FIELD_THETA_POLE) {
      bctype = GKYL_BC_MAXWELL_THETA_POLE;
    }

    f->bc_lo[d] = gkyl_bc_basic_new(
      d, GKYL_LOWER_EDGE, bctype, app->basis_on_dev, &app->lower_skin[d], &app->lower_ghost[d],
      f->em->ncomp, app->cdim, app->use_gpu
    );

    // Upper BC updater. Copy BCs by default.
    bctype = GKYL_BC_COPY;
    if (f->upper_bc[d] == GKYL_FIELD_COPY) {
      bctype = GKYL_BC_COPY;
    } else if (f->upper_bc[d] == GKYL_FIELD_PEC_WALL) {
      bctype = GKYL_BC_MAXWELL_PEC;
    } else if (f->upper_bc[d] == GKYL_FIELD_SYM_WALL) {
      bctype = GKYL_BC_MAXWELL_SYM;
    } else if (f->upper_bc[d] == GKYL_FIELD_RESERVOIR) {
      bctype = GKYL_BC_MAXWELL_RESERVOIR;
    } else if (f->upper_bc[d] == GKYL_FIELD_FIXED_FUNC) {
      bctype = GKYL_BC_FIXED_FUNC;
    } else if (f->upper_bc[d] == GKYL_FIELD_THETA_POLE) {
      bctype = GKYL_BC_MAXWELL_THETA_POLE;
    }
    f->bc_up[d] = gkyl_bc_basic_new(
      d, GKYL_UPPER_EDGE, bctype, app->basis_on_dev, &app->upper_skin[d], &app->upper_ghost[d],
      f->em->ncomp, app->cdim, app->use_gpu
    );

    if (f->lower_bc[d] == GKYL_FIELD_FIXED_FUNC) {
      f->bc_buffer_lo_fixed[d] = mkarr(app->use_gpu, f->em->ncomp, app->lower_ghost[d].volume);
      if (f->field_id == GKYL_FIELD_GR_D_B) {
        f->bc_buffer_lo_fixed_no_J[d] =
          mkarr(app->use_gpu, f->em->ncomp, app->lower_ghost[d].volume);
      }
    }
    if (f->upper_bc[d] == GKYL_FIELD_FIXED_FUNC) {
      f->bc_buffer_up_fixed[d] = mkarr(app->use_gpu, f->em->ncomp, app->upper_ghost[d].volume);
      if (f->field_id == GKYL_FIELD_GR_D_B) {
        f->bc_buffer_up_fixed_no_J[d] =
          mkarr(app->use_gpu, f->em->ncomp, app->upper_ghost[d].volume);
      }
    }
  }

  gkyl_dg_eqn_release(eqn);

  // Set the type-specific dispatch methods (Vlasov-Maxwell).
  f->update_func = f->has_bvar ? vm_field_update_bvar : vm_field_update;
  f->combine_func = vm_field_combine;
  f->copy_range_func = vm_field_copy_range;
  f->apply_ic_func = vm_field_apply_ic;
  f->apply_bc_func = vm_field_apply_bc;
  f->limiter_func = vm_field_limiter;
  f->complete_update_func = vm_field_complete_update;
  f->calc_ext_em_func = vm_field_calc_ext_em;
  f->calc_app_current_func = vm_field_calc_app_current;
  f->calc_ext_pot_func = vm_field_calc_ext_pot;
  f->calc_energy_func = vm_field_calc_energy;
  f->write_func = vm_field_write;
  f->write_energy_func = vm_field_write_energy;
  f->from_file_func = vm_field_from_file;
  f->release_func = vm_field_release;

  return f;
}

void
vm_field_apply_ic(
  gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *fin[], double t0
)
{
  (void)fin; // unused: the Vlasov-Maxwell IC comes from the field init function.
  if (!app->has_field) {
    return;
  }

  int poly_order = app->poly_order;
  // Project the (physical-space) EM field IC on the physical coordinates of the
  // (possibly mapped) conf mesh via the position-map c2p (identity => no effect).
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->basis,
    .qtype = GKYL_GAUSS_QUAD,
    .num_quad = poly_order + 1,
    .num_ret_vals = 8,
    .eval = field->info.init,
    .ctx = field->info.ctx,
    .c2p_func = vm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });

  // run updater; need to project onto extended range for ease of handling
  // subsequent operations over extended range such as magnetic field unit vector computation
  // This is needed to fill the corner cells as the corner cells may not be filled by
  // boundary conditions and we cannot divide by 0 anywhere or the weak divisions will fail
  gkyl_proj_on_basis_advance(proj, t0, &app->local_ext, field->em_host);
  gkyl_proj_on_basis_release(proj);

  if (field->field_id == GKYL_FIELD_GR_D_B) {
    // The input function is specified in primitive variables. Preserve that in em_no_J
    // and rescale to the conservative J-weighted fields used for evolution.
    gkyl_array_copy(field->em_no_J, field->em_host);
    gkyl_dg_gr_maxwell_rescale_Jc(
      &app->basis, &app->local_ext, app->vm_geom->det_h, field->em_no_J, field->em, app->use_gpu
    );
    // On a mapped conf mesh the evolved field also carries the (cell-constant)
    // position-map Jacobian: J_pos*J_c*(D,B).
    if (!app->pos_map->is_identity) {
      gkyl_vlasov_position_map_rescale_jacobpos_conf(
        app->pos_map, &app->local_ext, field->em, field->em
      );
    }
  } else if (field->weight_by_pos_jacob) {
    // Physical E, B were projected into em_host/em_no_J; store the J-weighted
    // fields evolved by the mapped Maxwell solver (em = J*E, J*B).
    gkyl_array_copy(field->em_no_J, field->em_host);
    gkyl_vlasov_position_map_rescale_jacobpos_conf(
      app->pos_map, &app->local_ext, field->em_no_J, field->em
    );
  } else if (app->use_gpu) {
    gkyl_array_copy(field->em, field->em_host);
  }
  vm_field_buffer_fixed_func_bc(app, field);
  // Apply limiter at t=0 to insure slopes are well-behaved at beginning of simulation
  vm_field_limiter(app, field, field->em);

  // pre-compute external EM field and applied current if present
  // pre-computation necessary in case external EM field or applied current
  // are time-independent and not computed in the time-stepping loop
  vm_field_calc_ext_em(app, field, t0);
  vm_field_calc_app_current(app, field, t0);

  if (field->has_bvar) {
    // b of the initial field (which includes the external field, so after it), over the extended
    // range: the boundary conditions of the species read it in the ghost cells. The boundary
    // conditions may not fill the corner cells, which hold the projection of the initial field.
    vm_field_apply_bc(app, field, field->em);
    vm_field_calc_bvar(app, field, field->em);
  }
}

void
vm_field_calc_ext_pot(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  // No external potentials in Vlasov-Maxwell.
}

void
vm_field_calc_ext_em(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  if (field->has_ext_em) {
    gkyl_proj_on_basis_advance(field->ext_em_proj, tm, &app->local_ext, field->ext_em_host);
    if (app->use_gpu) {
      // Note: ext_em_host is same as ext_em when not on GPUs.
      gkyl_array_copy(field->ext_em, field->ext_em_host);
    }
  }
}

void
vm_field_calc_app_current(gkyl_vlasov_app *app, struct vm_field *field, double tm)
{
  if (field->has_app_current) {
    gkyl_proj_on_basis_advance(
      field->app_current_proj, tm, &app->local_ext, field->app_current_host
    );
    gkyl_array_copy(field->app_current, field->app_current_host);
    if (field->app_current_no_J) {
      gkyl_array_copy(field->app_current_no_J, field->app_current);
    }
    // On a mapped mesh the field evolves J*E, so the current it is driven by is J*J_app (like the
    // species current, a moment of J_x*J_v*f); app_current_host keeps the physical current.
    if (field->weight_by_pos_jacob) {
      gkyl_vlasov_position_map_rescale_jacobpos_conf(
        app->pos_map, &app->local_ext, field->app_current, field->app_current
      );
    }
  }
}

void
vm_field_accumulate_current(
  gkyl_vlasov_app *app, const struct gkyl_array *fin[], const struct gkyl_array *fluidin[],
  struct gkyl_array *emout
)
{
  // Accumulate the current of each species onto the E rows (kinetic: -q/eps0 M1; Euler fluids
  // coupled explicitly: -q/(m eps0) rho u).
  int num_species = app->num_species;
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_accumulate_field_coupling(app, &app->species[i], fin[i], fluidin[i], emout);
  }

  // Accumulate the applied current onto the E rows, unless the fluid-EM coupling is implicit: its
  // source step applies it (vm_fluid_em_coupling.c).
  if (app->field->has_app_current && !app->has_fluid_em_coupling) {
    gkyl_array_accumulate_range(
      emout, -1.0 / app->field->info.epsilon0, app->field->app_current, &app->local
    );
  }
}

void
vm_field_accumulate_geom_sources(
  gkyl_vlasov_app *app, const struct gkyl_array *emin, const struct vm_geom *vm_geom,
  struct gkyl_array *emout
)
{
  // Accumulate the geometric source terms onto the fields .
  // Accumulate *only* if there is geometry sources to accumulate.
  if (app->field->use_geom_sources) {
    gkyl_array_clear(app->field->geom_source, 0.0);
    gkyl_dg_gr_maxwell_geom_source_advance(
      app->field->calc_geom_source, &app->local, vm_geom->geom_factor_con, emin,
      app->field->geom_source
    );
    gkyl_array_accumulate_range(emout, 1.0, app->field->geom_source, &app->local);
  }
}

void
vm_field_limiter(gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *em)
{
  if (field->limit_em) {
    // Limit the slopes of the solution
    gkyl_dg_calc_em_vars_limiter(field->calc_em_vars, &app->local, em);

    // Apply boundary conditions after limiting solution
    vm_field_apply_bc(app, field, em);
  }
}

// Combine the RK stages of the field state (out = c1*arr1 + c2*arr2) over the extended range.
void
vm_field_combine(
  gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *out, double c1,
  const struct gkyl_array *arr1, double c2, const struct gkyl_array *arr2
)
{
  array_combine(out, c1, arr1, c2, arr2, &app->local_ext);
}

// Copy the field state (out = inp) over the extended range.
void
vm_field_copy_range(
  gkyl_vlasov_app *app, struct vm_field *field, struct gkyl_array *out, const struct gkyl_array *inp
)
{
  gkyl_array_copy_range(out, inp, &app->local_ext);
}

// Vlasov-Maxwell field update: compute the RHS of Maxwell's equations.
double
vm_field_update(
  gkyl_vlasov_app *app, int stage, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  return vm_field_rhs(app, app->field, emin, emout);
}

// Transport step of the staggered operator split: the right-hand side of the stage with zero E
// rows. The transport step does not advance E (the source step does), so B drifts with the E of
// the source step; the terms of the B rows that depend on B (divergence cleaning of the magnetic
// field, the flux of an absorbing edge) follow the stages.
static double
vm_field_update_split(
  gkyl_vlasov_app *app, int stage, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  double dt_out = vm_field_rhs(app, app->field, emin, emout);
  gkyl_array_set_offset(emout, 0.0, app->field->app_current, 0);
  return dt_out;
}

void
vm_field_calc_bvar(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em)
{
  struct timespec tm = gkyl_wall_clock();

  // Physical fields of the stage (the evolved em carries the position-map Jacobian).
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local_ext, em->ncomp, em, field->em_no_J
  );
  gkyl_array_clear(field->tot_em, 0.0);
  gkyl_array_set(field->tot_em, 1.0, field->em_no_J);
  if (field->has_ext_em) {
    gkyl_array_accumulate_range(field->tot_em, 1.0, field->ext_em, &app->local_ext);
  }
  // Assumes magnetic field boundary conditions applied so the unit vector is
  // defined everywhere in the extended domain.
  gkyl_dg_calc_pkpm_bvar_advance(field->calc_bvar, &app->local_ext, field->tot_em, field->bvar);
  // Weak gradient of b for the mirror terms (div b), from b in the ghost cells.
  gkyl_dg_calc_pkpm_bvar_grad(field->calc_bvar, &app->local, field->bvar, field->grad_b);

  app->stat.field_em_vars_tm += gkyl_time_diff_now_sec(tm);
}

// Physical time derivative of the fields from a right-hand side of the evolved ones (which carry
// the position-map Jacobian), and d_t B/|B| at the nodes of b (tot_em holds the physical fields
// the right-hand side was computed from). A static field keeps both zero.
static void
vm_field_calc_dtB(gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *rhs)
{
  if (field->info.is_static) {
    return;
  }
  gkyl_vlasov_position_map_divide_jacobpos_conf(
    app->pos_map, &app->local, rhs->ncomp, rhs, field->dt_em_no_J
  );
  gkyl_dg_calc_pkpm_bvar_dtB_over_magB(
    field->calc_bvar, &app->local, field->tot_em, field->dt_em_no_J, field->dtB_over_magB
  );
}

void
vm_field_calc_diag_vars(gkyl_vlasov_app *app, struct vm_field *field)
{
  vm_field_calc_bvar(app, field, field->em);
  if (!field->info.is_static) {
    // emnew is free scratch between time steps.
    gkyl_array_clear(field->cflrate, 0.0);
    gkyl_array_clear(field->emnew, 0.0);
    gkyl_hyper_dg_advance(field->slvr, &app->local, field->em, field->cflrate, field->emnew);
    vm_field_calc_dtB(app, field, field->emnew);
  }
}

// Field update of a stage with PKPM species: b of the stage fields first (the species read it),
// then the right-hand side of Maxwell's equations, which does not depend on the species and whose
// magnetic part is the d_t B of the first-harmonic force terms.
double
vm_field_update_bvar(
  gkyl_vlasov_app *app, int stage, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  vm_field_calc_bvar(app, app->field, emin);
  double dt_out = vm_field_rhs(app, app->field, emin, emout);
  vm_field_calc_dtB(app, app->field, emout);
  return dt_out;
}

// The same in the transport step of the staggered operator split (zero E rows).
static double
vm_field_update_split_bvar(
  gkyl_vlasov_app *app, int stage, double tcurr, const struct gkyl_array *fin[],
  const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  vm_field_calc_bvar(app, app->field, emin);
  double dt_out = vm_field_update_split(app, stage, tcurr, fin, emin, emout);
  vm_field_calc_dtB(app, app->field, emout);
  return dt_out;
}

// Completion of a stage of the transport step of the staggered operator split: the currents of
// all species (and the applied current) enter the source step, so the stage only finalizes
// emout = emin + dt*RHS.
static void
vm_field_complete_update_split(
  gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  gkyl_array_accumulate(gkyl_array_scale(emout, dt), 1.0, emin);
}

double
vm_field_ghost_current_mean(
  gkyl_vlasov_app *app, struct vm_field *field, double fac, const struct gkyl_array *arr
)
{
  double avals_ghost_current[1], avals_ghost_current_global[1];
  // First set the scalar ghost current array to fac / nx times the leading coefficient, where nx
  // is the number of x cells.
  gkyl_array_set_range(field->ghost_current, fac / app->grid.cells[0], arr, &app->local);
  // Sum over the whole domain for the global average.
  if (app->use_gpu) {
    gkyl_array_reduce_range(field->red_ghost_current, field->ghost_current, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(
      avals_ghost_current, field->red_ghost_current, sizeof(double[1]), GKYL_CU_MEMCPY_D2H
    );
  } else {
    gkyl_array_reduce_range(avals_ghost_current, field->ghost_current, GKYL_SUM, &app->local);
  }
  gkyl_comm_allreduce_host(
    app->comm, GKYL_DOUBLE, GKYL_SUM, 1, avals_ghost_current, avals_ghost_current_global
  );
  return avals_ghost_current_global[0];
}

void
vm_field_accumulate_ghost_current(
  gkyl_vlasov_app *app, struct vm_field *field, double fac, const struct gkyl_array *current,
  struct gkyl_array *target
)
{
  if (field->use_ghost_current) {
    double mean = vm_field_ghost_current_mean(app, field, fac, current);
    // Set the scalar ghost current array to the global average current and accumulate to the
    // electric field, weighted by the position-map Jacobian like the field it drives (a uniform
    // current is then no source of E on a position map).
    gkyl_array_clear(field->ghost_current, mean);
    if (field->weight_by_pos_jacob) {
      gkyl_vlasov_position_map_rescale_jacobpos_conf(
        app->pos_map, &app->local, field->ghost_current, field->ghost_current
      );
    }
    gkyl_array_accumulate_range(target, 1.0, field->ghost_current, &app->local);
  }
}

void
vm_field_set_source_split(gkyl_vlasov_app *app, struct vm_field *field)
{
  // The split advances E in the source step only (with the resistive layer): the correction
  // potential of div(E) would be driven by an E the transport step does not advance.
  if (field->info.elcErrorSpeedFactor != 0.0) {
    gkyl_exit(
      "vm_field: divergence cleaning of the electric field is not supported with the implicit "
      "fluid-EM coupling."
    );
  }
  field->update_func = field->has_bvar ? vm_field_update_split_bvar : vm_field_update_split;
  field->complete_update_func = vm_field_complete_update_split;
}

// Vlasov-Maxwell completion of the field update: accumulate the species current
// onto the field RHS (unless the field is static), then finalize the explicit
// step emout = emin + dt*RHS.
void
vm_field_complete_update(
  gkyl_vlasov_app *app, double dt, const struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *emout
)
{
  struct timespec wst = gkyl_wall_clock();

  // (can't accumulate current when field is static)
  if (!app->field->info.is_static) {
    // Accumulate the species and applied currents onto the E rows, and the damping of E by the
    // resistive layer.
    vm_field_accumulate_current(app, fin, fluidin, emout);
    if (app->field->has_sigma && app->field->field_id != GKYL_FIELD_GR_D_B) {
      gkyl_dg_calc_resistive_layer_advance(app->field->sigma_slvr, &app->local, emin, emout);
    }
    app->stat.current_tm += gkyl_time_diff_now_sec(wst);
  }

  // Accumulate geometric source terms onto the field RHS.
  if (app->field->use_geom_sources) {
    vm_field_accumulate_geom_sources(app, emin, app->vm_geom, emout);
  }

  // complete update of field (even when field is static, it is
  // safest to do this accumulate as it ensure emout = emin)
  gkyl_array_accumulate(gkyl_array_scale(emout, dt), 1.0, emin);
}

// Compute the RHS for field update, returning maximum stable time-step.
double
vm_field_rhs(
  gkyl_vlasov_app *app, struct vm_field *field, const struct gkyl_array *em, struct gkyl_array *rhs
)
{
  struct timespec wst = gkyl_wall_clock();

  double dt_out = DBL_MAX;

  gkyl_array_clear(field->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  if (field->field_id == GKYL_FIELD_GR_D_B) {
    // Divide out configuration-space Jacobian.
    gkyl_dg_gr_maxwell_divide_Jc(
      &app->basis, &app->local, app->vm_geom->det_h, em, field->em_no_J, app->use_gpu
    );

    // Apply BCs after dividing out J so ghost cells are populated for
    // conf_flux_surf, which references the ghost cells for the flux (except on
    // theta-pole faces, whose flux is built from the cell next to the pole). The
    // input em already has its BCs applied by the caller.
    vm_field_apply_bc(app, field, field->em_no_J);

    // Compute the surface expansion of the GR-Maxwell flux at configuration-space surfaces.
    gkyl_dg_gr_maxwell_conf_flux_surf_advance(
      field->calc_conf_flux, &app->local, &app->local_ext, field->geom->lapse, field->geom->shift,
      field->geom->h_ij, field->geom->h_ij_inv, field->geom->det_h, em, field->em_no_J,
      field->cflrate, field->conf_flux_surf
    );
  }

  if (!field->info.is_static) {
    gkyl_hyper_dg_advance(field->slvr, &app->local, em, field->cflrate, rhs);

    gkyl_array_reduce_range(field->omega_cfl, field->cflrate, GKYL_MAX, &app->local);

    app->stat.n_field_omega_cfl += 1;
    struct timespec tm = gkyl_wall_clock();

    double omega_cfl_ho[1];
    if (app->use_gpu) {
      gkyl_cu_memcpy(omega_cfl_ho, field->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
    } else {
      omega_cfl_ho[0] = field->omega_cfl[0];
    }
    dt_out = app->cfl / omega_cfl_ho[0];

    app->stat.field_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  }

  app->stat.field_rhs_tm += gkyl_time_diff_now_sec(wst);

  return dt_out;
}

// Determine which directions are periodic and which directions are not periodic,
// and then apply boundary conditions for EM fields
void
vm_field_apply_bc(gkyl_vlasov_app *app, const struct vm_field *field, struct gkyl_array *f)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    app->comm, &app->local, &app->local_ext, num_periodic_dir, app->periodic_dirs, f
  );

  int is_np_bc[3] = {1, 1, 1}; // flags to indicate if direction is not periodic
  for (int d = 0; d < num_periodic_dir; ++d) {
    is_np_bc[app->periodic_dirs[d]] = 0;
  }

  for (int d = 0; d < cdim; ++d) {
    if (is_np_bc[d]) {
      switch (field->lower_bc[d]) {
        case GKYL_FIELD_COPY:
        case GKYL_FIELD_ABSORB:
        case GKYL_FIELD_PEC_WALL:
        case GKYL_FIELD_SYM_WALL:
        case GKYL_FIELD_THETA_POLE:
        case GKYL_FIELD_RESERVOIR:
          gkyl_bc_basic_advance(field->bc_lo[d], field->bc_buffer, f);
          break;

        case GKYL_FIELD_FIXED_FUNC:
          gkyl_bc_basic_advance(
            field->bc_lo[d], vm_field_get_fixed_func_bc_buffer(field, f, d, GKYL_LOWER_EDGE), f
          );
          break;

        default:
          break;
      }

      switch (field->upper_bc[d]) {
        case GKYL_FIELD_COPY:
        case GKYL_FIELD_ABSORB:
        case GKYL_FIELD_PEC_WALL:
        case GKYL_FIELD_SYM_WALL:
        case GKYL_FIELD_THETA_POLE:
        case GKYL_FIELD_RESERVOIR:
          gkyl_bc_basic_advance(field->bc_up[d], field->bc_buffer, f);
          break;

        case GKYL_FIELD_FIXED_FUNC:
          gkyl_bc_basic_advance(
            field->bc_up[d], vm_field_get_fixed_func_bc_buffer(field, f, d, GKYL_UPPER_EDGE), f
          );
          break;

        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, f);

  app->stat.field_bc_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_field_write(gkyl_vlasov_app *app, double tm, int frame, const struct gkyl_array *fin[])
{
  struct timespec wst = gkyl_wall_clock();

  // The field frame also stores the steps of the staggered operator split (zero without it).
  struct gkyl_msgpack_data *mt = vlasov_array_meta_new((struct vlasov_output_meta){
    .frame = frame,
    .stime = tm,
    .poly_order = app->poly_order,
    .basis_type = app->basis.id,
    .dt_prev = app->op_split_dt_prev,
    .dt_suggested = app->op_split_dt_suggested,
  });

  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz + 1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  const struct gkyl_array *field_to_write = app->field->em_host;
  if (app->field->field_id == GKYL_FIELD_GR_D_B) {
    // For GR Maxwell, write primitive D/B fields by default.
    gkyl_dg_gr_maxwell_divide_Jc(
      &app->basis, &app->local, app->vm_geom->det_h, app->field->em, app->field->em_no_J,
      app->use_gpu
    );

    if (app->use_gpu) {
      gkyl_array_copy(app->field->em_no_J_host, app->field->em_no_J);
    }
    // On a mapped mesh em_no_J still carries the position-map Jacobian.
    if (!app->pos_map->is_identity) {
      gkyl_vlasov_position_map_divide_jacobpos_conf(
        app->pos_map, &app->local, app->field->em_no_J_host->ncomp, app->field->em_no_J_host,
        app->field->em_no_J_host
      );
    }
    field_to_write = app->field->em_no_J_host;
  } else if (app->field->weight_by_pos_jacob) {
    // Write the physical E, B (divide the stored J*E, J*B by the conf Jacobian).
    if (app->use_gpu) {
      gkyl_array_copy(app->field->em_host, app->field->em);
    }
    gkyl_vlasov_position_map_divide_jacobpos_conf(
      app->pos_map, &app->local, app->field->em_host->ncomp, app->field->em_host,
      app->field->em_no_J_host
    );
    field_to_write = app->field->em_no_J_host;
  } else if (app->use_gpu) {
    // Copy data from device to host before writing it out.
    gkyl_array_copy(app->field->em_host, app->field->em);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, field_to_write, fileNm);

  if (app->field->has_ext_em) {
    // Only write out external fields at t=0 or if they are time-dependent.
    if (frame == 0 || app->field->ext_em_evolve) {
      const char *fmt_ext_em = "%s-field_ext_em_%d.gkyl";
      int sz_ext_em = gkyl_calc_strlen(fmt_ext_em, app->name, frame);
      char fileNm_ext_em[sz_ext_em + 1]; // Ensures no buffer overflow.
      snprintf(fileNm_ext_em, sizeof fileNm_ext_em, fmt_ext_em, app->name, frame);

      // External EM field computed with project on basis, so just use host copy.
      vm_field_calc_ext_em(app, app->field, tm);

      gkyl_comm_array_write(
        app->comm, &app->grid, &app->local, mt, app->field->ext_em_host, fileNm_ext_em
      );
    }
  }
  if (app->field->has_app_current) {
    // Only write the applied current at t=0 or if it is time-dependent.
    if (frame == 0 || app->field->app_current_evolve) {
      const char *fmt_app_current = "%s-field_app_current_%d.gkyl";
      int sz_app_current = gkyl_calc_strlen(fmt_app_current, app->name, frame);
      char fileNm_app_current[sz_app_current + 1]; // Ensures no buffer overflow.
      snprintf(fileNm_app_current, sizeof fileNm_app_current, fmt_app_current, app->name, frame);

      // Applied current computed with project on basis, so just use host copy.
      vm_field_calc_app_current(app, app->field, tm);

      gkyl_comm_array_write(
        app->comm, &app->grid, &app->local, mt, app->field->app_current_host, fileNm_app_current
      );
    }
  }

  vlasov_array_meta_release(mt);

  app->stat.field_io_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_field_io += 1;
}

void
vm_field_calc_energy(
  gkyl_vlasov_app *app, double tm, struct vm_field *field, const struct gkyl_array *fin[]
)
{
  struct timespec wst = gkyl_wall_clock();

  for (int i = 0; i < 6; ++i) {
    gkyl_dg_calc_l2_range(&app->basis, i, field->em_energy, i, field->em, app->local);
  }
  gkyl_array_scale_range(field->em_energy, app->grid.cellVolume, &app->local);
  // The field stores J*E, J*B, so the L2 carries J^2 while the (computational)
  // cellVolume under-counts the physical volume by J; the net is one extra factor
  // of J. Divide it out so the energy is the physical int J E^2, int J B^2
  // (= int E^2, B^2 dx_phys). Identity map => J=1 => unchanged.
  if (field->weight_by_pos_jacob ||
      (field->field_id == GKYL_FIELD_GR_D_B && !app->pos_map->is_identity)) {
    gkyl_vlasov_position_map_divide_jacobpos_conf(
      app->pos_map, &app->local, field->em_energy->ncomp, field->em_energy, field->em_energy
    );
  }

  double energy[6] = {0.0};
  if (app->use_gpu) {
    gkyl_array_reduce_range(field->em_energy_red, field->em_energy, GKYL_SUM, &app->local);
    gkyl_cu_memcpy(energy, field->em_energy_red, sizeof(double[6]), GKYL_CU_MEMCPY_D2H);
  } else {
    gkyl_array_reduce_range(energy, field->em_energy, GKYL_SUM, &app->local);
  }

  double energy_global[6] = {0.0};
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, energy, energy_global);

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

void
vm_field_write_energy(gkyl_vlasov_app *app)
{
  struct timespec wst = gkyl_wall_clock();

  // Write out integrated field energy.
  const char *fmt = "%s-field-energy.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name);
  char fileNm[sz + 1]; // Ensures no buffer overflow.
  snprintf(fileNm, sizeof fileNm, fmt, app->name);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    if (app->field->is_first_energy_write_call) {
      // Write to a new file (this ensure previous output is removed).
      gkyl_dynvec_write(app->field->integ_energy, fileNm);
      app->field->is_first_energy_write_call = false;
    } else {
      // Append to existing file.
      gkyl_dynvec_awrite(app->field->integ_energy, fileNm);
    }
  }
  gkyl_dynvec_clear(app->field->integ_energy);

  app->stat.n_field_diag_io += 1;
  app->stat.field_diag_io_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_app_restart_status
vm_field_from_file(gkyl_vlasov_app *app, struct vm_field *field, const char *fname)
{
  struct gkyl_app_restart_status rstat = vlasov_header_from_file(app, fname);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Fixed-function field BCs are frozen from the initial conditions, so seed
    // those buffers before the local solution is overwritten by restart data.
    vm_field_apply_ic(app, field, 0, 0.0);

    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, field->em_host, fname);
    if (app->use_gpu) {
      gkyl_array_copy(field->em, field->em_host);
    }
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      // For GR, rescale the primitive fields to the evolved quantities by
      // multiplying by Jc.
      if (field->field_id == GKYL_FIELD_GR_D_B) {
        gkyl_array_copy(field->em_no_J, field->em_host);
        gkyl_dg_gr_maxwell_rescale_Jc(
          &app->basis, &app->local_ext, app->vm_geom->det_h, field->em_no_J, field->em, app->use_gpu
        );
        // On a mapped conf mesh the evolved field also carries the position-map
        // Jacobian; the file holds the physical fields.
        if (!app->pos_map->is_identity) {
          gkyl_vlasov_position_map_rescale_jacobpos_conf(
            app->pos_map, &app->local, field->em, field->em
          );
        }
      } else if (field->weight_by_pos_jacob) {
        // Restart files hold the physical E, B; rescale to the evolved J*E, J*B on
        // the interior (BCs re-fill the ghost cells below).
        gkyl_array_copy(field->em_no_J, field->em_host);
        gkyl_vlasov_position_map_rescale_jacobpos_conf(
          app->pos_map, &app->local, field->em_no_J, field->em
        );
      }

      vm_field_apply_bc(app, field, field->em);
    }
  }

  // Recompute the external EM field and applied current: they are not part of
  // the restart file, and the time-stepping loop only recomputes them when they
  // are time dependent.
  vm_field_calc_ext_em(app, field, rstat.stime);
  vm_field_calc_app_current(app, field, rstat.stime);

  if (field->has_bvar && rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // b of the restart field (which includes the external field, so after it), over the
    // extended range (the reflecting walls of the species read it in the ghost cells)
    vm_field_calc_bvar(app, field, field->em);
  }

  return rstat;
}

// release resources for field
void
vm_field_release(const gkyl_vlasov_app *app, struct vm_field *f)
{
  gkyl_array_release(f->em);
  gkyl_array_release(f->em1);
  gkyl_array_release(f->emnew);
  if (f->weight_by_pos_jacob) {
    gkyl_array_release(f->em_no_J);
    gkyl_array_release(f->em_no_J_host);
  }
  if (f->field_id == GKYL_FIELD_GR_D_B) {
    gkyl_array_release(f->em_no_J);
    gkyl_array_release(f->em_no_J_host);
    gkyl_dg_gr_maxwell_conf_flux_surf_release(f->calc_conf_flux);
    gkyl_array_release(f->conf_flux_surf);
    if (f->use_geom_sources) {
      gkyl_dg_gr_maxwell_geom_source_release(f->calc_geom_source);
      gkyl_array_release(f->geom_source);
    }
  }
  gkyl_array_release(f->em_host);

  gkyl_array_release(f->bc_buffer);
  gkyl_array_release(f->cflrate);
  gkyl_array_release(f->em_energy);
  gkyl_dynvec_release(f->integ_energy);

  gkyl_array_release(f->sigma);
  if (f->has_sigma) {
    gkyl_dg_calc_resistive_layer_release(f->sigma_slvr);
  }

  gkyl_array_release(f->ext_em);
  if (f->has_ext_em) {
    gkyl_array_release(f->ext_em_host);
    gkyl_proj_on_basis_release(f->ext_em_proj);
  }
  gkyl_array_release(f->app_current);
  if (f->has_app_current) {
    gkyl_array_release(f->app_current_host);
    gkyl_proj_on_basis_release(f->app_current_proj);
  }
  if (f->app_current_no_J) {
    gkyl_array_release(f->app_current_no_J);
  }
  if (f->has_bvar) {
    gkyl_array_release(f->tot_em);
    gkyl_array_release(f->bvar);
    gkyl_array_release(f->grad_b);
    gkyl_array_release(f->dt_em_no_J);
    gkyl_array_release(f->dtB_over_magB);
    gkyl_dg_calc_pkpm_bvar_release(f->calc_bvar);
  }

  gkyl_hyper_dg_release(f->slvr);

  gkyl_dg_calc_em_vars_release(f->calc_em_vars);

  if (app->use_gpu) {
    gkyl_cu_free(f->omega_cfl);
    gkyl_cu_free(f->em_energy_red);
  } else {
    gkyl_free(f->omega_cfl);
  }

  if (f->use_ghost_current) {
    gkyl_array_release(f->ghost_current);
    if (app->use_gpu) {
      gkyl_cu_free(f->red_ghost_current);
    }
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d = 0; d < app->cdim; ++d) {
    if (f->bc_buffer_lo_fixed[d]) {
      gkyl_array_release(f->bc_buffer_lo_fixed[d]);
    }
    if (f->bc_buffer_up_fixed[d]) {
      gkyl_array_release(f->bc_buffer_up_fixed[d]);
    }
    if (f->bc_buffer_lo_fixed_no_J[d]) {
      gkyl_array_release(f->bc_buffer_lo_fixed_no_J[d]);
    }
    if (f->bc_buffer_up_fixed_no_J[d]) {
      gkyl_array_release(f->bc_buffer_up_fixed_no_J[d]);
    }
    gkyl_bc_basic_release(f->bc_lo[d]);
    gkyl_bc_basic_release(f->bc_up[d]);
  }

  gkyl_free(f);
}
