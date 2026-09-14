#include <acutest.h>

#include <float.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_loss_cone_mask_gyrokinetic.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_velocity_map.h>

struct loss_cone_mask_test_ctx {
  int cdim;
  double eV;
  double R_m;
  double B_m;
  double z_m;
  double mass, charge;
  double n0, T0, B0;
  double phi_fac;
  double z_max, vpar_max, mu_max;
  int Nz, Nvpar, Nmu;
};

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                 : gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
bmag_func_1x(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[0];
  struct loss_cone_mask_test_ctx *params = ctx;

  fout[0] = params->B_m
    * (1.0 - ((params->R_m - 1.0) / params->R_m) * pow(cos(z), 2.0));
}

static void
phi_func_1x_zero(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

static void
phi_func_1x_nonzero(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xc[0];
  struct loss_cone_mask_test_ctx *params = ctx;

  // fout[0] = params->phi_fac * params->T0 / params->eV * (1.0 - pow(z / params->z_max, 2.0));
  fout[0] = params->phi_fac * params->T0 / params->eV * (1.0 - pow(cos(4*z), 2.0));
}

static double
analytic_bmag(const struct loss_cone_mask_test_ctx *ctx, double z)
{
  return ctx->B_m * (1.0 - ((ctx->R_m - 1.0) / ctx->R_m) * pow(cos(z), 2.0));
}

static double
analytic_phi(const struct loss_cone_mask_test_ctx *ctx, double z,
  bool use_nonzero_phi)
{
  return use_nonzero_phi
    ? ctx->phi_fac * ctx->T0 / ctx->eV * (1.0 - pow(cos(4.0 * z), 2.0))
    : 0.0;
}

static double
analytic_effective_potential(const struct loss_cone_mask_test_ctx *ctx,
  double z, double mu, bool use_nonzero_phi)
{
  return mu * analytic_bmag(ctx, z)
    + ctx->charge * analytic_phi(ctx, z, use_nonzero_phi);
}

/**
 * Build the expected mask directly from the analytic B(z), phi(z), Cartesian
 * velocity-cell endpoints, and H=mv_parallel^2/2+mu*B+q*phi. The P1 fields
 * are nodally projected, so their cell-corner values equal these analytic
 * functions. This oracle deliberately does not evaluate the production DG
 * arrays or call any loss-cone helper.
 */
static void
build_analytic_reference_mask_1x2v(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *phase_range, const struct loss_cone_mask_test_ctx *ctx,
  bool use_nonzero_phi, struct gkyl_array *mask_ref)
{
  struct gkyl_range corner_range;
  gkyl_range_init_from_shape(&corner_range, 3, (int[3]) { 2, 2, 2 });

  struct gkyl_range_iter phase_iter;
  gkyl_range_iter_init(&phase_iter, phase_range);
  while (gkyl_range_iter_next(&phase_iter)) {
    bool cell_trapped = true;
    struct gkyl_range_iter corner_iter;
    gkyl_range_iter_init(&corner_iter, &corner_range);
    while (cell_trapped && gkyl_range_iter_next(&corner_iter)) {
      int z_vertex = phase_iter.idx[0] - 1 + corner_iter.idx[0];
      double z = grid->lower[0] + z_vertex * grid->dx[0];
      double vpar = grid->lower[1]
        + (phase_iter.idx[1] - 1 + corner_iter.idx[1]) * grid->dx[1];
      double mu = grid->lower[2]
        + (phase_iter.idx[2] - 1 + corner_iter.idx[2]) * grid->dx[2];
      double u_curr = analytic_effective_potential(ctx, z, mu, use_nonzero_phi);
      double h_curr = 0.5 * ctx->mass * vpar * vpar + u_curr;

      double barrier_left = -DBL_MAX, barrier_right = -DBL_MAX;
      for (int vertex = 0; vertex <= ctx->Nz; ++vertex) {
        double z_scan = grid->lower[0] + vertex * grid->dx[0];
        double u_scan = analytic_effective_potential(ctx, z_scan, mu,
          use_nonzero_phi);
        if (vertex <= z_vertex) {
          barrier_left = GKYL_MAX2(barrier_left, u_scan);
        }
        if (vertex >= z_vertex) {
          barrier_right = GKYL_MAX2(barrier_right, u_scan);
        }
      }

      cell_trapped = h_curr < GKYL_MIN2(barrier_left, barrier_right);
    }

    long linidx_phase = gkyl_range_idx(phase_range, phase_iter.idx);
    double *m = gkyl_array_fetch(mask_ref, linidx_phase);
    m[0] = cell_trapped ? 1.0 : 0.0;
  }

}

static void
run_case_1x2v(int poly_order, bool use_gpu, bool use_nonzero_phi)
{
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mass_proton = GKYL_PROTON_MASS;

  struct loss_cone_mask_test_ctx ctx = {
    .cdim = 1,
    .eV = eV,
    .R_m = 8.0,
    .B_m = 4.0,
    .z_m = GKYL_PI / 2.0,
    .mass = 2.014 * mass_proton,
    .charge = eV,
    .n0 = 1e18,
    .T0 = 100 * eV,
    .phi_fac = 5.0,
    .z_max = GKYL_PI - 0.5,
    .Nz = 64,
    .Nvpar = 16,
    .Nmu = 16,
  };
  ctx.B0 = ctx.B_m / 2.0;
  ctx.vpar_max = 6.0 * sqrt(ctx.T0 / ctx.mass);
  ctx.mu_max = 0.5 * ctx.mass * pow(ctx.vpar_max, 2) / ctx.B0;

  double lower[] = { -ctx.z_max, -ctx.vpar_max, 0.0 };
  double upper[] = { ctx.z_max, ctx.vpar_max, ctx.mu_max };
  int cells[] = { ctx.Nz, ctx.Nvpar, ctx.Nmu };
  const int ndim = sizeof(cells) / sizeof(cells[0]);
  const int cdim = ctx.cdim;
  const int vdim = ndim - cdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d = 0; d < cdim; ++d) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }

  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d = 0; d < vdim; ++d) {
    lower_vel[d] = lower[cdim + d];
    upper_vel[d] = upper[cdim + d];
    cells_vel[d] = cells[cdim + d];
  }

  struct gkyl_rect_grid grid, grid_conf, grid_vel;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_ext_conf;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_ext_vel;
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[GKYL_MAX_DIM] = { 0 };
  for (int d = 0; d < cdim; ++d) {
    ghost[d] = ghost_conf[d];
  }
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);
  struct gkyl_velocity_map *gvm_ho = use_gpu
    ? gkyl_velocity_map_new(c2p_in, grid, grid_vel, local, local_ext,
        local_vel, local_ext_vel, false)
    : 0;

  struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *bmag_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);

  // Nodal P1 projection makes the DG corner values exactly equal to the
  // analytic B(z) values used by the independent reference below.
  gkyl_eval_on_nodes *evbmag = gkyl_eval_on_nodes_new(&grid_conf, &basis_conf,
    1, bmag_func_1x, &ctx);
  gkyl_eval_on_nodes_advance(evbmag, 0.0, &local_conf, bmag_ho);
  gkyl_eval_on_nodes_release(evbmag);
  gkyl_array_copy(bmag, bmag_ho);

  evalf_t phi_func = use_nonzero_phi ? phi_func_1x_nonzero : phi_func_1x_zero;
  gkyl_eval_on_nodes *evphi = gkyl_eval_on_nodes_new(&grid_conf, &basis_conf, 1, phi_func, &ctx);
  gkyl_eval_on_nodes_advance(evphi, 0.0, &local_conf, phi_ho);
  gkyl_eval_on_nodes_release(evphi);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj = {
    .conf_basis = &basis_conf,
    .vel_map = gvm,
    .use_gpu = use_gpu,
    .mass = ctx.mass,
    .charge = ctx.charge,
  };
  struct gkyl_loss_cone_mask_gyrokinetic *proj_mask =
    gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj);

  struct gkyl_array *mask = mkarr(use_gpu, 1, local_ext.volume);
  struct gkyl_array *mask_ho = mkarr(false, 1, local_ext.volume);
  struct gkyl_array *mask_cpu = use_gpu ? mkarr(false, 1, local_ext.volume) : 0;
  struct gkyl_array *mask_ref = mkarr(false, 1, local_ext.volume);

  gkyl_loss_cone_mask_gyrokinetic_advance(proj_mask, &local, &local_conf,
    bmag, phi, 0, 0, mask);
  gkyl_array_copy(mask_ho, mask);

  gkyl_array_clear(mask_ref, 0.0);
  build_analytic_reference_mask_1x2v(&grid, &local, &ctx, use_nonzero_phi,
    mask_ref);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    struct gkyl_loss_cone_mask_gyrokinetic_inp inp_proj_cpu = {
      .conf_basis = &basis_conf,
      .vel_map = gvm_ho,
      .use_gpu = false,
      .mass = ctx.mass,
      .charge = ctx.charge,
    };
    struct gkyl_loss_cone_mask_gyrokinetic *proj_mask_cpu =
      gkyl_loss_cone_mask_gyrokinetic_inew(&inp_proj_cpu);

    gkyl_loss_cone_mask_gyrokinetic_advance(proj_mask_cpu, &local, &local_conf,
      bmag_ho, phi_ho, 0, 0, mask_cpu);

    int gpu_cpu_mismatches = 0;
    struct gkyl_range_iter iter_cmp;
    gkyl_range_iter_init(&iter_cmp, &local);
    while (gkyl_range_iter_next(&iter_cmp)) {
      long linidx = gkyl_range_idx(&local, iter_cmp.idx);
      const double *mg = gkyl_array_cfetch(mask_ho, linidx);
      const double *mc = gkyl_array_cfetch(mask_cpu, linidx);

      bool same = fabs(mg[0] - mc[0]) < 1e-12;
      TEST_CHECK(same);
      if (!same && gpu_cpu_mismatches < 8) {
        printf("gpu/cpu mismatch idx=(%d,%d,%d): gpu=%g cpu=%g\n",
          iter_cmp.idx[0], iter_cmp.idx[1], iter_cmp.idx[2], mg[0], mc[0]);
        gpu_cpu_mismatches++;
      }
    }

    gkyl_loss_cone_mask_gyrokinetic_release(proj_mask_cpu);
  }
#endif

  int mismatches = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *m = gkyl_array_cfetch(mask_ho, linidx);
    const double *r = gkyl_array_cfetch(mask_ref, linidx);

    bool ok = fabs(m[0] - r[0]) < 1e-12;
    TEST_CHECK(ok);
    if (!ok && mismatches < 8) {
      printf("mask mismatch idx=(%d,%d,%d): got=%g ref=%g\n",
        iter.idx[0], iter.idx[1], iter.idx[2], m[0], r[0]);
      mismatches++;
    }
  }

  gkyl_array_release(bmag);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(mask);
  gkyl_array_release(mask_ho);
  if (mask_cpu) {
    gkyl_array_release(mask_cpu);
  }
  gkyl_array_release(mask_ref);

  gkyl_loss_cone_mask_gyrokinetic_release(proj_mask);
  gkyl_velocity_map_release(gvm);
  if (gvm_ho) {
    gkyl_velocity_map_release(gvm_ho);
  }
}

struct constant_boundary_case {
  const char *name;
  double charge;
  double phi_plasma;
  double phi_wall_lo;
  double phi_wall_up;
  bool use_wall_lo;
  bool use_wall_up;
  enum gkyl_gk_trapped_passing_orbit_type lower_orbit;
  enum gkyl_gk_trapped_passing_orbit_type upper_orbit;
};

static void
set_constant_1x_p1(struct gkyl_array *field, const struct gkyl_range *range,
  double value)
{
  gkyl_array_clear(field, 0.0);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(range, iter.idx);
    double *field_d = gkyl_array_fetch(field, linidx);
    field_d[0] = sqrt(2.0) * value;
  }
}

static double
constant_boundary_energy(double charge, double phi_plasma, double phi_wall,
  bool use_wall, enum gkyl_gk_trapped_passing_orbit_type orbit)
{
  if (orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL) {
    return DBL_MAX;
  }
  return orbit == GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH && use_wall
    ? GKYL_MAX2(0.0, charge * (phi_wall - phi_plasma)) : 0.0;
}

static void
run_constant_boundary_cases(bool use_gpu)
{
  // With constant B and phi, mu*B cancels between H and each wall barrier.
  // These cases therefore have an exact, independent criterion in each
  // v_parallel cell: every velocity endpoint must satisfy K < Delta U at
  // both boundaries.
  const struct constant_boundary_case cases[] = {
    {
      .name = "grounded_wall_electron",
      .charge = -1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      .name = "grounded_wall_ion_no_barrier",
      .charge = 1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      .name = "grounded_wall_ion_reversed_phi",
      .charge = 1.0, .phi_plasma = -4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      .name = "asymmetric_biased_wall_electron",
      .charge = -1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 3.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      .name = "asymmetric_biased_wall_ion",
      .charge = 1.0, .phi_plasma = 0.0,
      .phi_wall_lo = 4.0, .phi_wall_up = 1.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      // This is the asymmetric electron case shifted everywhere by +7 V.
      .name = "gauge_shifted_asymmetric_electron",
      .charge = -1.0, .phi_plasma = 11.0,
      .phi_wall_lo = 7.0, .phi_wall_up = 10.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      // Delta U=2 and a v_parallel endpoint is exactly |v|=2. Equality
      // reaches the wall and must be classified as passing/absorbed.
      .name = "sheath_cutoff_equality",
      .charge = -1.0, .phi_plasma = 2.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      .name = "passing_orbits",
      .charge = -1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_PASSING,
      .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_PASSING,
    },
    {
      .name = "passing_and_wall_trapped_orbits",
      .charge = -1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_PASSING,
      .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL,
    },
    {
      .name = "wall_trapped_and_sheath_trapped_orbits",
      .charge = -1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL,
      .use_wall_up = true, .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
    },
    {
      .name = "sheath_trapped_and_wall_trapped_orbits",
      .charge = -1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .use_wall_lo = true, .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_SHEATH,
      .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL,
    },
    {
      .name = "wall_trapped_orbits",
      .charge = -1.0, .phi_plasma = 4.0,
      .phi_wall_lo = 0.0, .phi_wall_up = 0.0,
      .lower_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL,
      .upper_orbit = GKYL_GK_TRAP_PASS_ORBIT_TRAPPED_WALL,
    },
  };

  double lower[] = { -1.0, -4.0, 0.0 };
  double upper[] = { 1.0, 4.0, 1.0 };
  int cells[] = { 2, 8, 1 };
  double lower_conf[] = { lower[0] }, upper_conf[] = { upper[0] };
  double lower_vel[] = { lower[1], lower[2] };
  double upper_vel[] = { upper[1], upper[2] };
  int cells_conf[] = { cells[0] }, cells_vel[] = { cells[1], cells[2] };

  struct gkyl_rect_grid grid, grid_conf, grid_vel;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  gkyl_rect_grid_init(&grid_conf, 1, lower_conf, upper_conf, cells_conf);
  gkyl_rect_grid_init(&grid_vel, 2, lower_vel, upper_vel, cells_vel);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, 1, 1);

  int ghost_conf[] = { 1 };
  int ghost_vel[] = { 0, 0 };
  int ghost[] = { 1, 0, 0 };
  struct gkyl_range local_conf, local_ext_conf, local_vel, local_ext_vel;
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_ext_vel, use_gpu);

  struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_wall_lo = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_wall_up = mkarr(use_gpu, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *mask = mkarr(use_gpu, 1, local_ext.volume);

  struct gkyl_array *bmag_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_wall_lo_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *phi_wall_up_ho = mkarr(false, basis_conf.num_basis, local_ext_conf.volume);
  struct gkyl_array *mask_ho = mkarr(false, 1, local_ext.volume);

  set_constant_1x_p1(bmag_ho, &local_ext_conf, 2.0);
  gkyl_array_copy(bmag, bmag_ho);

  for (int c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
    const struct constant_boundary_case *test = &cases[c];

    set_constant_1x_p1(phi_ho, &local_ext_conf, test->phi_plasma);
    set_constant_1x_p1(phi_wall_lo_ho, &local_ext_conf, test->phi_wall_lo);
    set_constant_1x_p1(phi_wall_up_ho, &local_ext_conf, test->phi_wall_up);
    gkyl_array_copy(phi, phi_ho);
    gkyl_array_copy(phi_wall_lo, phi_wall_lo_ho);
    gkyl_array_copy(phi_wall_up, phi_wall_up_ho);

    struct gkyl_loss_cone_mask_gyrokinetic *up =
      gkyl_loss_cone_mask_gyrokinetic_inew(&(struct gkyl_loss_cone_mask_gyrokinetic_inp) {
        .conf_basis = &basis_conf,
        .vel_map = gvm,
        .use_gpu = use_gpu,
        .mass = 1.0,
        .charge = test->charge,
        .lower_orbit = test->lower_orbit,
        .upper_orbit = test->upper_orbit,
      });

    const struct gkyl_array *wall_lo = test->use_wall_lo ? phi_wall_lo : 0;
    const struct gkyl_array *wall_up = test->use_wall_up ? phi_wall_up : 0;
    gkyl_loss_cone_mask_gyrokinetic_advance(up, &local, &local_conf,
      bmag, phi, wall_lo, wall_up, mask);
    gkyl_array_copy(mask_ho, mask);

    double barrier_lo = constant_boundary_energy(test->charge, test->phi_plasma,
      test->phi_wall_lo, test->use_wall_lo, test->lower_orbit);
    double barrier_up = constant_boundary_energy(test->charge, test->phi_plasma,
      test->phi_wall_up, test->use_wall_up, test->upper_orbit);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      int ivpar = iter.idx[1];
      double vpar_lo = grid.lower[1] + (ivpar - 1) * grid.dx[1];
      double vpar_up = vpar_lo + grid.dx[1];
      double kinetic_max = 0.5 * GKYL_MAX2(vpar_lo * vpar_lo, vpar_up * vpar_up);
      double expected = kinetic_max < barrier_lo && kinetic_max < barrier_up ? 1.0 : 0.0;

      long linidx = gkyl_range_idx(&local, iter.idx);
      const double *mask_d = gkyl_array_cfetch(mask_ho, linidx);
      double actual = mask_d[0];
      bool ok = fabs(actual - expected) < 1e-12;
      TEST_CHECK(ok);
      if (!ok) {
        TEST_MSG("%s idx=(%d,%d,%d): got=%g expected=%g Kmax=%g barriers=(%g,%g)",
          test->name, iter.idx[0], iter.idx[1], iter.idx[2], actual, expected,
          kinetic_max, barrier_lo, barrier_up);
      }
    }

    gkyl_loss_cone_mask_gyrokinetic_release(up);
  }

  gkyl_array_release(bmag);
  gkyl_array_release(phi);
  gkyl_array_release(phi_wall_lo);
  gkyl_array_release(phi_wall_up);
  gkyl_array_release(mask);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phi_wall_lo_ho);
  gkyl_array_release(phi_wall_up_ho);
  gkyl_array_release(mask_ho);
  gkyl_velocity_map_release(gvm);
}

void
test_1x2v_p1_gk_ho(void)
{
  // Magnetic-mirror barriers only: compare every host mask cell with the
  // analytic Hamiltonian reference built above.
  run_case_1x2v(1, false, false);
}

void
test_1x2v_p1_nonzero_phi_gk_ho(void)
{
  // Add a spatially varying electrostatic potential to check that q*phi is
  // included in both directional barriers.
  run_case_1x2v(1, false, true);
}

void
test_1x2v_p1_constant_boundaries_gk_ho(void)
{
  // Check exact wall-potential thresholds, open/closed trajectories, gauge
  // invariance, and equality at the absorbing cutoff.
  run_constant_boundary_cases(false);
}

#ifdef GKYL_HAVE_CUDA
// The device cases use the same independent expectations as the host cases
// and also compare the device mask cell-by-cell with a host updater run.
void
test_1x2v_p1_gk_dev(void)
{
  run_case_1x2v(1, true, false);
}

void
test_1x2v_p1_nonzero_phi_gk_dev(void)
{
  run_case_1x2v(1, true, true);
}

void
test_1x2v_p1_constant_boundaries_gk_dev(void)
{
  run_constant_boundary_cases(true);
}
#endif

TEST_LIST = {
  { "test_1x2v_p1_gk_ho", test_1x2v_p1_gk_ho },
  { "test_1x2v_p1_nonzero_phi_gk_ho", test_1x2v_p1_nonzero_phi_gk_ho },
  { "test_1x2v_p1_constant_boundaries_gk_ho", test_1x2v_p1_constant_boundaries_gk_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x2v_p1_gk_dev", test_1x2v_p1_gk_dev },
  { "test_1x2v_p1_nonzero_phi_gk_dev", test_1x2v_p1_nonzero_phi_gk_dev },
  { "test_1x2v_p1_constant_boundaries_gk_dev", test_1x2v_p1_constant_boundaries_gk_dev },
#endif
  { NULL, NULL },
};
