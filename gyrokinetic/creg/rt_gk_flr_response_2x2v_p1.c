// Single-mode validation of the Pade FLR chain (see flr_effects.tex).
//
// Straight slab, uniform B, Dirichlet walls in x (the fully periodic
// perpendicular domain is a known unresolved issue in fem_poisson_perp),
// periodic in z, and a single ion density mode
//   n_i = n0*(1 + A*sin(kx*x)),  n_e = n0,  kx = pi/Lx,
// which is an eigenfunction of every operator in the chain. The initial
// field solve is linear and diagonal in kx, so with b = (kx*rho_i)^2 the
// potential with FLR effects on and off satisfies
//   phi_flr/phi_noflr = (1 + b)/(1 + b/2),
// from the density gyroaverage (Gamma_1 = 1/(1+b/2)) and the FLR
// inversion of phi (A = 1 + b). No time stepping is needed; this test
// applies ICs, solves the field, and checks amplitudes at frame 0.
//
// We sweep kx*rho_i (and thus b)
// from small to large to check two limits of the Pade chain:
//   - b -> 0: rho_i -> 0 limit, ratio -> 1 (FLR effects vanish).
//   - b >> 1: ratio -> 2, the large-argument saturation of the Pade form.
//
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_const.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <rt_arg_parse.h>

struct flr_ctx {
  int cdim, vdim;
  double me, qe, mi, qi;
  double n0, Te0, Ti0, B0;
  double rho_i; // Ion gyroradius, sets the FLR scale.
  double kx; // Perturbation wavenumber.
  double pert_amp; // Relative amplitude of the density perturbation.
  double Lx, Lz;
  double vpar_max_elc, mu_max_elc, vpar_max_ion, mu_max_ion;
  int cells[GKYL_MAX_DIM];
};

static struct flr_ctx
create_ctx(double kx_rho)
{
  double eV = GKYL_ELEMENTARY_CHARGE;
  double me = GKYL_ELECTRON_MASS, mi = GKYL_PROTON_MASS;
  double qe = -eV, qi = eV;

  double n0 = 1.0e19;
  double Te0 = 100.0*eV, Ti0 = 100.0*eV;
  double B0 = 1.0;

  double vte = sqrt(Te0/me), vti = sqrt(Ti0/mi);
  double rho_i = sqrt(Ti0*mi)/(qi*B0);

  double kx = kx_rho/rho_i;
  double Lx = M_PI/kx;
  double Lz = 1.0;

  struct flr_ctx ctx = {
    .cdim = 2, .vdim = 2,
    .me = me, .qe = qe, .mi = mi, .qi = qi,
    .n0 = n0, .Te0 = Te0, .Ti0 = Ti0, .B0 = B0,
    .rho_i = rho_i,
    .kx = kx,
    .pert_amp = 1.0e-2,
    .Lx = Lx, .Lz = Lz,
    .vpar_max_elc = 5.0*vte,
    .mu_max_elc = me*pow(5.0*vte,2)/(2.0*B0),
    .vpar_max_ion = 5.0*vti,
    .mu_max_ion = mi*pow(5.0*vti,2)/(2.0*B0),
    .cells = {32, 4, 12, 8},
  };
  return ctx;
}

static void
density_elc(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = app->n0;
}

static void
density_ion(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  double x = xn[0];
  fout[0] = app->n0*(1.0 + app->pert_amp*sin(app->kx*x));
}

static void
temp_elc(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = app->Te0;
}

static void
temp_ion(double t, const double *xn, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = app->Ti0;
}

static void
zero_func(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 0.0;
}

static void
mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

static void
bfield_func(double t, const double *xc, double *fout, void *ctx)
{
  struct flr_ctx *app = ctx;
  fout[0] = 0.0; fout[1] = 0.0; fout[2] = app->B0;
}

// Run the initial field solve with or without FLR effects and write frame 0.
static void
run_case(struct flr_ctx *ctx, struct gkyl_app_args *app_args, struct gkyl_comm *comm,
  bool use_flr, const char *name)
{
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx->qe, .mass = ctx->me,
    .vdim = ctx->vdim,
    .lower = { -ctx->vpar_max_elc, 0.0},
    .upper = {  ctx->vpar_max_elc, ctx->mu_max_elc},
    .cells = { ctx->cells[2], ctx->cells[3] },
    .polarization_density = ctx->n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = ctx,
      .ctx_upar = ctx,
      .ctx_temp = ctx,
      .density = density_elc,
      .upar = zero_func,
      .temp = temp_elc,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
    },
    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_M0},
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx->qi, .mass = ctx->mi,
    .vdim = ctx->vdim,
    .lower = { -ctx->vpar_max_ion, 0.0},
    .upper = {  ctx->vpar_max_ion, ctx->mu_max_ion},
    .cells = { ctx->cells[2], ctx->cells[3] },
    .polarization_density = ctx->n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = ctx,
      .ctx_upar = ctx,
      .ctx_temp = ctx,
      .density = density_ion,
      .upar = zero_func,
      .temp = temp_ion,
    },

    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
    },
    .num_diag_moments = 1,
    .diag_moments = {GKYL_F_MOMENT_M0},
  };

  if (use_flr) {
    ion.flr.type = GKYL_GK_FLR_PADE_CONST;
    ion.flr.Tperp = ctx->Ti0;
  }

  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ES,
    .polarization_bmag = ctx->B0,
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    },
  };

  if (use_flr) {
    field.flr.type = GKYL_GK_FLR_PADE_CONST;
    field.flr.gyroradius = ctx->rho_i;
  }

  struct gkyl_gk gk = {
    .cdim = ctx->cdim,
    .lower = { 0.0, -ctx->Lz/2.0 },
    .upper = { ctx->Lx, ctx->Lz/2.0 },
    .cells = { ctx->cells[0], ctx->cells[1] },
    .poly_order = 1,
    .basis_type = app_args->basis_type,

    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MAPC2P,
      .world = {0.0},
      .mapc2p = mapc2p,
      .c2p_ctx = ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = ctx,
    },

    .num_periodic_dir = 1,
    .periodic_dirs = {1},

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .comm = comm,
      .cuts = {app_args->cuts[0], app_args->cuts[1]},
      .use_gpu = app_args->use_gpu,
    },
  };
  strcpy(gk.name, name);

  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);
  gkyl_gyrokinetic_app_apply_ic(app, 0.0);
  gkyl_gyrokinetic_app_write_field(app, 0.0, 0);
  gkyl_gyrokinetic_app_release(app);
}

// L2 norm of a DG field, sqrt(int f^2 dx).
static double
calc_l2(struct gkyl_rect_grid grid, struct gkyl_range range, struct gkyl_range range_ext,
  struct gkyl_basis basis, struct gkyl_array *f)
{
  struct gkyl_array *l2 = gkyl_array_new(GKYL_DOUBLE, 1, range_ext.volume);
  gkyl_dg_calc_l2_range(&basis, 0, l2, 0, f, range);
  gkyl_array_scale_range(l2, grid.cellVolume, &range);
  double l2red[1];
  gkyl_array_reduce_range(l2red, l2, GKYL_SUM, &range);
  gkyl_array_release(l2);
  return sqrt(l2red[0]);
}

static int
run_and_check_case(int case_idx, double kx_rho, struct gkyl_app_args *app_args, struct gkyl_comm *comm)
{
  struct flr_ctx ctx = create_ctx(kx_rho);

  char name_off[64], name_on[64];
  snprintf(name_off, sizeof(name_off), "rt_gk_flr_response_2x2v_p1_case%d_flroff", case_idx);
  snprintf(name_on, sizeof(name_on), "rt_gk_flr_response_2x2v_p1_case%d_flron", case_idx);

  run_case(&ctx, app_args, comm, false, name_off);
  run_case(&ctx, app_args, comm, true, name_on);

  // Read back both potentials.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ctx.cdim, (double[]) { 0.0, -ctx.Lz/2.0 },
    (double[]) { ctx.Lx, ctx.Lz/2.0 }, (int[]) { ctx.cells[0], ctx.cells[1] });
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ctx.cdim, 1);
  int nghost[GKYL_MAX_CDIM] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, nghost, &local_ext, &local);

  struct gkyl_array *phi_off = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_array *phi_on = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  struct gkyl_rect_grid grid_read;
  char fname_off[96], fname_on[96];
  snprintf(fname_off, sizeof(fname_off), "%s-field_0.gkyl", name_off);
  snprintf(fname_on, sizeof(fname_on), "%s-field_0.gkyl", name_on);
  gkyl_grid_sub_array_read(&grid_read, &local, phi_off, fname_off);
  gkyl_grid_sub_array_read(&grid_read, &local, phi_on, fname_on);

  // Measured amplitudes and expected factors.
  double l2_off = calc_l2(grid, local, local_ext, basis, phi_off);
  double l2_on = calc_l2(grid, local, local_ext, basis, phi_on);

  double b = pow(ctx.kx*ctx.rho_i, 2.0);
  double ratio_expected = (1.0 + b)/(1.0 + b/2.0);
  double ratio = l2_on/l2_off;

  // Analytic amplitude of the FLR-off solve: phi = qi*A*n0/(eps*kx^2)*cos(kx*x).
  double eps_pol = ctx.n0*(ctx.mi + ctx.me)/pow(ctx.B0, 2.0);
  double phi_amp = ctx.qi*ctx.pert_amp*ctx.n0/(eps_pol*pow(ctx.kx, 2.0));
  double l2_off_expected = phi_amp*sqrt(ctx.Lx*ctx.Lz/2.0);

  // Mode shape: phi_on must be the same mode scaled by the expected ratio.
  struct gkyl_array *diff = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  gkyl_array_set(diff, 1.0, phi_on);
  gkyl_array_accumulate(diff, -ratio_expected, phi_off);
  double l2_shape = calc_l2(grid, local, local_ext, basis, diff);

  double err_ratio = fabs(ratio/ratio_expected - 1.0);
  double err_abs = fabs(l2_off/l2_off_expected - 1.0);
  double err_shape = l2_shape/l2_on;

  printf("\nFLR single-mode response test (kx*rho_i = %g, b = %g):\n", ctx.kx*ctx.rho_i, b);
  printf("  |phi| FLR-off measured/expected: %e / %e (err = %.3e)\n", l2_off, l2_off_expected, err_abs);
  printf("  |phi| FLR-on:                    %e\n", l2_on);
  printf("  amplitude ratio measured/expected: %f / %f (err = %.3e)\n", ratio, ratio_expected, err_ratio);
  printf("  mode-shape error: %.3e\n", err_shape);

  int num_fail = 0;
  num_fail += err_ratio > 0.02;
  num_fail += err_abs > 0.05;
  num_fail += err_shape > 0.02;
  printf("%s\n", num_fail == 0? "PASSED" : "FAILED");

  gkyl_array_release(phi_off);
  gkyl_array_release(phi_on);
  gkyl_array_release(diff);

  return num_fail;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Sweep kx*rho_i from the rho_i->0 limit (ratio -> 1) to the large-b Pade saturation (ratio -> 2).
  double kx_rho_list[] = { 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0 };
  int num_cases = sizeof(kx_rho_list)/sizeof(kx_rho_list[0]);

  int num_fail = 0;
  for (int i = 0; i < num_cases; i++)
    num_fail += run_and_check_case(i, kx_rho_list[i], &app_args, comm);

  printf("\n%d/%d cases failed.\n", num_fail, num_cases);

  gkyl_gyrokinetic_comms_release(comm);

  return num_fail;
}
